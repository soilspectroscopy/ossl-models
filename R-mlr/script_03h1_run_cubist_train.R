
# Packages
packages <- c("tidyverse", "qs", "parallel", "lubridate",
              "mlr3verse", "mlr3extralearners", "mlr3pipelines",
              "yardstick")

invisible(lapply(packages, library, character.only = TRUE))

# Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

# Modeling combinations
modeling.combinations <- read_csv("../out/modeling_combinations_v1.2.csv", show_col_types = FALSE)
count.table <- read_csv("../out/tab_dataset_count.csv", show_col_types = FALSE)

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

# # Filtering already fitted models
# modeling.combinations <- modeling.combinations %>%
#   mutate(fitted = file.exists(paste0(dir, export_name, "/model_", model_name, ".qs")))
#
# modeling.combinations <- modeling.combinations  %>%
#   filter(!fitted)

# # Filtering corrupted files
# redo.train <- read_csv("../R-mlr/sandbox/corrupted/redo_train.csv", show_col_types = F)
#
# modeling.combinations <- modeling.combinations %>%
#   filter(model_name %in% unique(redo.train$model_name)) %>%
#   filter(soil_property %in% unique(redo.train$soil_property))

# target.dirs <- paste0(dir, unique(modeling.combinations$export_name))
#
# invisible(sapply(target.dirs, function(x) {
#   if(!dir.exists(x)){dir.create(x)}
# }))

## Parallelization is done inside the the autotuner
lgr::get_logger("mlr3")$set_threshold("warn")

i=1
for(i in 1:nrow(modeling.combinations)) {

  isoil_property = modeling.combinations[[i,"soil_property"]]
  imodel_name = modeling.combinations[[i,"model_name"]]
  iexport_name = modeling.combinations[[i,"export_name"]]
  ispectra_type = modeling.combinations[[i,"spectra_type"]]
  isubset = modeling.combinations[[i,"subset"]]
  igeo = modeling.combinations[[i,"geo"]]

  # Learner
  learner_cubist = lrn("regr.cubist", predict_type = "response",
                       neighbors = 0, unbiased = FALSE, seed = 1993)

  # Hyperparameters space
  search_space_ensemble = ps(
    committees = p_int(1, 20)
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
  sel.data <- data %>%
    select(id.layer_uuid_txt, # ID column
           all_of(isoil_property), # Target
           all_of(selected.comps)) %>% # Only compressed spectra
    filter(!is.na(!!as.name(isoil_property)))

  # Subset for speeding up HPO
  if(nrow(sel.data) >= 2000) {

    set.seed(1993)
    sel.data.hpo <- sel.data %>%
      sample_n(2000)

  } else {

    sel.data.hpo <- sel.data

  }

  # Exporting train data

  export.data <- paste0(dir, iexport_name, "/task_", imodel_name, ".qs")
  export.data.alt <- gsub("\\.qs", "\\.rds", export.data)

  if(file.exists(export.data)){file.remove(export.data)}
  if(file.exists(export.data.alt)){file.remove(export.data.alt)}

  tryCatch(
    expr = {qsave(sel.data, export.data)},
    error = function(e){saveRDS(sel.data, export.data.alt)}
  )

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
  # Grid search
  # Resolution 5 means that within the range, 5 equidistant values will be tested
  at = auto_tuner(tuner = tnr("grid_search", resolution = 5, batch_size = 1),
                  learner = learner_cubist,
                  resampling = inner_resampling,
                  measure = msr("regr.rmse"),
                  search_space = search_space_ensemble,
                  terminator = trm("none"),
                  store_models = FALSE)

  # Multicore processing
  future::plan("multisession")

  # Fit autotuner
  at$train(task.hpo)

  # # Overview
  # at$tuning_result
  # at$tuning_instance

  # Final model with best HPO
  final.model <- learner_cubist
  final.model$param_set$values = at$tuning_result$learner_param_vals[[1]]
  final.model$train(task.train)

  future:::ClusterRegistry("stop")

  # # Overview
  # summary(final.model$model$regr.lm$model)

  # Saving trained final model to disk
  export.model <- paste0(dir, iexport_name, "/model_", imodel_name, ".qs")
  export.model.alt <- gsub("\\.qs", "\\.rds", export.model)

  if(file.exists(export.model)){file.remove(export.model)}
  if(file.exists(export.model.alt)){file.remove(export.model.alt)}

  tryCatch(
    expr = {qsave(final.model, export.model)},
    error = function(e){saveRDS(final.model, export.model.alt)}
  )

}
