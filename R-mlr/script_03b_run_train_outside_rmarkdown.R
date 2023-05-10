
# Packages
packages <- c("tidyverse", "qs", "parallel", "lubridate",
              "mlr3verse", "mlr3extralearners", "mlr3pipelines",
              "yardstick")

invisible(lapply(packages, library, character.only = TRUE))

# Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

# Modeling combinations
modeling.combinations <- read_csv("../out/modeling_combinations_v1.2.csv")
count.table <- read_csv("../out/tab_dataset_count.csv")

# Defining available data from ossl-import
count.table <- count.table %>%
  filter(dataset %in% c("KSSL.SSL", "OSSL")) %>%
  rename(subset = dataset) %>%
  mutate(subset = recode(subset,
                         "KSSL.SSL" = "kssl",
                         "OSSL" = "ossl")) %>%
  pivot_longer(-all_of(c("soil_property", "subset")),
               names_to = "spectra_type", values_to = "count") %>%
  mutate(spectra_type = recode(spectra_type,
                               "n_mir" = "mir",
                               "n_visnir" = "visnir",
                               "n_neospectra" = "nir.neospectra"))

# Defining models with at least 500 observations
modeling.combinations <- left_join(modeling.combinations,
                                   count.table,
                                   by = c("soil_property", "spectra_type", "subset"))

modeling.combinations <- modeling.combinations %>%
  filter(count > 500) %>%
  filter(!(soil_property == "efferv_usda.a479_class"))

# Filtering already fitted models
modeling.combinations <- modeling.combinations %>%
  mutate(fitted = file.exists(paste0(dir, export_name, "/model_", model_name, ".qs")))

modeling.combinations <- modeling.combinations  %>%
  filter(!fitted)

## Parallelization is done inside the the autotuner
lgr::get_logger("mlr3")$set_threshold("warn")

i=1
for(i in 1:length(modeling.combinations)) {

  isoil_property = modeling.combinations[[i,"soil_property"]]
  imodel_name = modeling.combinations[[i,"model_name"]]
  iexport_name = modeling.combinations[[i,"export_name"]]
  ispectra_type = modeling.combinations[[i,"spectra_type"]]
  isubset = modeling.combinations[[i,"subset"]]
  igeo = modeling.combinations[[i,"geo"]]

  # Learners
  learner_glmnet = lrn("regr.glmnet", predict_type = "response")

  learner_ranger = lrn("regr.ranger", predict_type = "response",
                       replace = TRUE, num.threads = 1, verbose = FALSE)

  learner_xgboost = lrn("regr.xgboost", predict_type = "response",
                        booster = "gbtree", nthread = 1,
                        subsample = 0.67)

  learner_cubist = lrn("regr.cubist", predict_type = "response",
                       neighbors = 0, unbiased = FALSE, seed = 1993)

  # Base learners
  base_learners = list(learner_glmnet, learner_ranger, learner_xgboost, learner_cubist)

  # Meta learner: linear model of base learners
  learner_lm = lrn("regr.lm", predict_type = "response")

  meta_learner = pipeline_stacking(base_learners, learner_lm,
                                   # method = "cv", folds = 5,
                                   method = "insample",
                                   use_features = FALSE)

  # Setting ensemble as a learner
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

  # PCA scores
  n.comps <- 120
  selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

  data <- qread(paste0(dir, "pca.ossl/pca_scores_", ispectra_type, "_mlr3..eml_", isubset, "_v1.2.qs"))

  # Apply log transform to soil property
  if(grepl("log..", iexport_name)){

    data <- data %>%
      mutate(!!isoil_property := log1p(!!as.name(isoil_property)))

  }

  # Defining train data
  if(igeo == "ll") {

    sel.data <- data %>%
      select(id.layer_uuid_txt, # ID column
             id.tile, # ID for block CV
             all_of(isoil_property), # Target
             all_of(selected.comps), # Compressed spectra
             starts_with("clm")) %>% # Geocovariates
      filter(!is.na(!!as.name(isoil_property)))

  } else if(igeo == "na") {

    sel.data <- data %>%
      select(id.layer_uuid_txt, # ID column
             all_of(isoil_property), # Target
             all_of(selected.comps)) %>% # Only compressed spectra
      filter(!is.na(!!as.name(isoil_property)))

  }

  # Subset for speeding up HPO
  if(nrow(sel.data) >= 2000) {

    set.seed(1993)
    sel.data.hpo <- sel.data %>%
      sample_n(2000)

  } else {

    sel.data.hpo <- sel.data

  }

  # Exporting train data
  qsave(sel.data, paste0(dir,
                         iexport_name,
                         "/task_",
                         imodel_name,
                         ".qs"))

  # Create regression task
  task.hpo <- as_task_regr(sel.data.hpo, id = "hpo", target = isoil_property, type = "regression")
  task.train <- as_task_regr(sel.data, id = "train", target = isoil_property, type = "regression")

  # Defining id column
  task.hpo$set_col_roles("id.layer_uuid_txt", roles = "name")
  task.train$set_col_roles("id.layer_uuid_txt", roles = "name")

  # For block CV. If 'id.tile' not present in the data.frame, default to random CV
  if("id.tile" %in% colnames(sel.data)) {
    task.hpo$set_col_roles("id.tile", roles = "group")
    task.train$set_col_roles("id.tile", roles = "group")
  }

  # Inner resampling for HPO with 5-fold cv
  inner_resampling = rsmp("cv", folds = 5)

  # Auto tuner
  at = auto_tuner(tuner = tnr("random_search", batch_size = 3), # batch_size X 5 folds = 50 cores
                  learner = learner_ensemble,
                  resampling = inner_resampling,
                  measure = msr("regr.rmse"),
                  search_space = search_space_ensemble,
                  terminator = trm("evals", n_evals = 20),
                  store_models = FALSE)

  # Multicore processing
  future::plan("multisession")

  # Fit autotuner
  at$train(task.hpo)

  # # Overview
  # at$tuning_result
  # at$tuning_instance

  # Final model from best HPO
  final.model <- learner_ensemble
  final.model$param_set$values = at$tuning_result$learner_param_vals[[1]]
  final.model$train(task.train)

  future:::ClusterRegistry("stop")

  # # Overview
  # summary(final.model$model$regr.lm$model)

  # Saving trained final model to disk
  qsave(final.model, paste0(dir,
                            iexport_name,
                            "/model_",
                            imodel_name,
                            ".qs"))

}
