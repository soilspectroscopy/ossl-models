
# Packages
packages <- c("tidyverse", "qs", "parallel", "lubridate",
              "mlr3verse", "mlr3extralearners", "mlr3pipelines",
              "yardstick")

invisible(lapply(packages, library, character.only = TRUE))

library("conflicted")
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(recipes::prep)

# Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

# Modeling combinations
modeling.combinations <- read_csv("../out/fitted_modeling_combinations_v1.2_cubist.csv",
                                  show_col_types = FALSE)

# # Filtering already fitted models
# modeling.combinations <- modeling.combinations %>%
#   mutate(fitted = file.exists(paste0(dir, export_name, "/model_", model_name, ".qs")))
#
# modeling.combinations <- modeling.combinations  %>%
#   filter(!fitted)

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

  # PCA scores
  n.comps <- 120
  selected.comps <- paste0("PC", seq(1, n.comps, by = 1))

  data <- qread(paste0(dir, "pca.ossl/pca_scores_", ispectra_type, "_mlr3..eml_", isubset, "_v1.2.qs"))

  # CV predictions
  predictions <- qread(paste0(dir, iexport_name, "/cvpred_", imodel_name, ".qs"))

  predictions <- predictions %>%
    mutate(residual = abs(truth-response)) %>%
    select(id.layer_uuid_txt, residual)

  data <- right_join(data, predictions, by = "id.layer_uuid_txt") %>%
    relocate(residual, .after = id.layer_uuid_txt)

  # Defining train data
  sel.data <- data %>%
    select(id.layer_uuid_txt, # ID column
           residual, # Target
           all_of(selected.comps)) %>% # Only compressed spectra
    filter(!is.na(residual))

  # Create regression task
  task.error <- as_task_regr(sel.data, id = "train", target = "residual", type = "regression")

  # Defining id column
  task.error$set_col_roles("id.layer_uuid_txt", roles = "name")

  # For block CV. If 'id.tile' not present in the data.frame, default to random CV
  if("id.tile" %in% colnames(sel.data)) {
    task.error$set_col_roles("id.tile", roles = "group")
  }

  # Importing the best HP from prediction model and setting to error model
  error.model <- qread(paste0(dir, iexport_name, "/model_", imodel_name, ".qs"))

  # Train
  error.model$train(task.error)

  # Saving trained error model to disk
  export.file <- paste0(dir, iexport_name, "/error_model_", imodel_name, ".qs")
  export.file.alt <- gsub("\\.qs", "\\.rds", export.file)

  if(file.exists(export.file)){file.remove(export.file)}
  if(file.exists(export.file.alt)){file.remove(export.file.alt)}

  tryCatch(
    expr = {qsave(error.model, export.file)},
    error = function(e){saveRDS(error.model, export.file.alt)}
  )

  # Getting residual predictions and saving to disk
  results <- predict(error.model, newdata = as.data.table(task.error))

  pred.export <- sel.data %>%
    select(-starts_with("PC")) %>%
    mutate(pred_residual = results) %>%
    mutate(alpha_scores = residual/pred_residual)

  export.error <- paste0(dir, iexport_name, "/error_pred_", imodel_name, ".qs")
  export.error.alt <- gsub("\\.qs", "\\.rds", export.error)

  if(file.exists(export.error)){file.remove(export.error)}
  if(file.exists(export.error.alt)){file.remove(export.error.alt)}

  tryCatch(
    expr = {qsave(pred.export, export.error)},
    error = function(e){saveRDS(pred.export, export.error.alt)}
  )

  cat(paste0("Exported error model at ", now(), "\n\n"))

}
