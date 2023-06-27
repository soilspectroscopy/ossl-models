
## Packages
library("tidyverse")
library("lubridate")
library("mlr3verse")
library("mlr3extralearners")
library("mlr3pipelines")
library("yardstick")
library("qs")

## Folders
mnt.dir <- "~/projects/temp/ossl_models/comparison/"

## Importing data
modeling.combinations <- read_csv(file.path(mnt.dir, "modeling_combinations.csv"))

## Modeling pipeline

i=1
for(i in 1:nrow(modeling.combinations)) {

  imodel_type <- modeling.combinations[[i,"model_type"]]
  ibase_predictions <- modeling.combinations[[i,"base_predictions"]]
  isoil_property <- modeling.combinations[[i,"soil_property"]]
  ispectra <- modeling.combinations[[i,"spectra"]]
  igeocovariates <- modeling.combinations[[i,"geocovariates"]]
  ipca_compression <- modeling.combinations[[i,"pca_compression"]]

  lgr::get_logger("mlr3")$set_threshold("warn")
  future::plan("multisession")

  if(imodel_type == "cubist") {

    # Define learner
    lrn_cubist <- lrn("regr.cubist",
                      committees = to_tune(5, 20),
                      neighbors = 0, unbiased = FALSE, seed = 1993) # Tuning not working for this HP

    # Define autotuner
    at_cubist = auto_tuner(tuner = tnr("random_search"),
                           learner = lrn_cubist,
                           resampling = rsmp("cv", folds = 5),
                           measure = msr("regr.rmse"),
                           terminator = trm("evals", n_evals = 10),
                           store_models = TRUE)

    if(ispectra == "mir") {

      data <- qread(file.path(mnt.dir, "pca_scores_mir_kssl_subset_v1.2.qs")) %>%
        mutate(!!isoil_property := log1p(!!as.name(isoil_property)))

      if(ipca_compression == "cumvar99perc"){

        pca.model <- qread(file.path(dirname(mnt.dir), "mpca_mir_mlr..eml_kssl_v1.2.qs"))

        # Defining the number of PCs by cumulative percent variance (99%)
        n.comps <- which(cumsum(pca.model$sdev/sum(pca.model$sdev)) > 0.99)[1]

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      } else if(ipca_compression == "n120") {

        n.comps <- 120

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      }


    } else if(ispectra == "visnir") {

      data <- qread(file.path(mnt.dir, "pca_scores_visnir_kssl_subset_v1.2.qs")) %>%
        mutate(!!isoil_property := log1p(!!as.name(isoil_property)))

      if(ipca_compression == "cumvar99perc"){

        pca.model <- qread(file.path(dirname(mnt.dir), "mpca_visnir_mlr..eml_kssl_v1.2.qs"))

        # Defining the number of PCs by cumulative percent variance (99%)
        n.comps <- which(cumsum(pca.model$sdev/sum(pca.model$sdev)) > 0.99)[1]

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      } else if(ipca_compression == "n120") {

        n.comps <- 120

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      }

    }

    # Export regression matrix = task
    qsave(sel.data, paste0(mnt.dir,
                           "autotune/",
                           paste("task",
                                 imodel_type,
                                 ibase_predictions,
                                 isoil_property,
                                 ispectra,
                                 igeocovariates,
                                 ipca_compression,
                                 sep = "_"),
                           ".qs"))

    # Create regression task
    task <- as_task_regr(sel.data, id = "train", target = isoil_property, type = "regression")

    # Defining id column
    task$set_col_roles("id.layer_uuid_txt", roles = "name")

    # For block cv. If 'id.tile' not in the data.frame, default to random CV
    if("id.tile" %in% colnames(sel.data)) {
      task$set_col_roles("id.tile", roles = "group")
    }

    # Fit
    at_cubist$train(task)

    # Summary
    # at_cubist$tuning_result
    # at_cubist$tuning_instance

    # Saving model to disk
    qsave(at_cubist, paste0(mnt.dir,
                            "autotune/",
                            paste("model",
                                  imodel_type,
                                  ibase_predictions,
                                  isoil_property,
                                  ispectra,
                                  igeocovariates,
                                  ipca_compression,
                                  sep = "_"),
                            ".qs"))

  } else if(imodel_type == "ensemble") {

    # Learners
    learner_glmnet = lrn("regr.glmnet", predict_type = "response")

    learner_ranger = lrn("regr.ranger", predict_type = "response",
                         replace = TRUE, num.threads = 1)

    learner_xgboost = lrn("regr.xgboost", predict_type = "response",
                          booster = "gbtree", nthread = 1,
                          subsample = 0.67)

    learner_cubist = lrn("regr.cubist", predict_type = "response",
                         neighbors = 0, unbiased = FALSE, seed = 1993)

    # Base learners
    base_learners = list(learner_glmnet, learner_ranger, learner_xgboost, learner_cubist)

    # Super learner: linear model of base learners
    learner_lm = lrn("regr.lm", predict_type = "response")

    if(ibase_predictions == "insample") {

      meta_learner = pipeline_stacking(base_learners, learner_lm,
                                       method = "insample",
                                       use_features = FALSE)

    } else if(ibase_predictions == "cv5") {

      meta_learner = pipeline_stacking(base_learners, learner_lm,
                                       method = "cv", folds = 5,
                                       use_features = FALSE)

    }

    # Setting as a learner
    learner_ensemble = as_learner(meta_learner)
    learner_ensemble$id = "ensemble"
    learner_ensemble$predict_type = "response"

    # Hyperparameters space, all crossed, i.e. not tuned separately
    search_space_ensemble = ps(
      regr.glmnet.alpha = p_dbl(0, 1),
      regr.glmnet.lambda = p_dbl(0.001, 0.1),
      regr.ranger.num.trees = p_int(20, 100),
      regr.ranger.min.node.size = p_int(5, 20),
      regr.xgboost.nrounds = p_int(20, 100),
      regr.xgboost.eta = p_dbl(0.3, 0.5),
      regr.xgboost.max_depth = p_int(5, 20),
      regr.cubist.committees = p_int(5, 10)
    )

    if(ispectra == "mir") {

      data <- qread(file.path(mnt.dir, "pca_scores_mir_kssl_subset_v1.2.qs")) %>%
        mutate(!!isoil_property := log1p(!!as.name(isoil_property)))

      if(ipca_compression == "cumvar99perc"){

        pca.model <- qread(file.path(dirname(mnt.dir), "mpca_mir_mlr..eml_kssl_v1.2.qs"))

        # Defining the number of PCs by cumulative percent variance (99%)
        n.comps <- which(cumsum(pca.model$sdev/sum(pca.model$sdev)) > 0.99)[1]

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      } else if(ipca_compression == "n120") {

        n.comps <- 120

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      }


    } else if(ispectra == "visnir") {

      data <- qread(file.path(mnt.dir, "pca_scores_visnir_kssl_subset_v1.2.qs")) %>%
        mutate(!!isoil_property := log1p(!!as.name(isoil_property)))

      if(ipca_compression == "cumvar99perc"){

        pca.model <- qread(file.path(dirname(mnt.dir), "mpca_visnir_mlr..eml_kssl_v1.2.qs"))

        # Defining the number of PCs by cumulative percent variance (99%)
        n.comps <- which(cumsum(pca.model$sdev/sum(pca.model$sdev)) > 0.99)[1]

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      } else if(ipca_compression == "n120") {

        n.comps <- 120

        selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

        if(igeocovariates == "ll") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps), # Compressed spectra
                   starts_with("clm")) # Geocovariates

        } else if(igeocovariates == "na") {

          sel.data <- data %>%
            select(id.layer_uuid_txt, # ID column
                   id.tile, # ID for block CV
                   all_of(isoil_property), # Target
                   all_of(selected.comps)) # Only compressed spectra

        }

      }

    }

    # Exporting regression matrix = task
    qsave(sel.data, paste0(mnt.dir,
                           "autotune/",
                           paste("task",
                                 imodel_type,
                                 ibase_predictions,
                                 isoil_property,
                                 ispectra,
                                 igeocovariates,
                                 ipca_compression,
                                 sep = "_"),
                           ".qs"))
    # Create regression task
    task <- as_task_regr(sel.data, id = "train", target = isoil_property, type = "regression")

    # Defining id column
    task$set_col_roles("id.layer_uuid_txt", roles = "name")

    # For block CV. If 'id.tile' not in the data.frame, default to random CV
    if("id.tile" %in% colnames(sel.data)) {
      task$set_col_roles("id.tile", roles = "group")
    }

    # Inner resampling for HPO
    inner_resampling = rsmp("cv", folds = 5)

    # Auto tuner
    at = auto_tuner(tuner = tnr("random_search"),
                    learner = learner_ensemble,
                    resampling = inner_resampling,
                    measure = msr("regr.rmse"),
                    search_space = search_space_ensemble,
                    terminator = trm("evals", n_evals = 20),
                    store_models = TRUE)

    # Fit
    at$train(task)

    # Summary
    # at$tuning_result
    # at$tuning_instance
    # summary(at$model$learner$model$regr.lm$model)

    # Saving model to disk
    qsave(at, paste0(mnt.dir,
                     "autotune/",
                     paste("model",
                           imodel_type,
                           ibase_predictions,
                           isoil_property,
                           ispectra,
                           igeocovariates,
                           ipca_compression,
                           sep = "_"),
                     ".qs"))

  }

}
