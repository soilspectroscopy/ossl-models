
# Packages
packages <- c("tidyverse", "tidymodels",
              "lubridate", "qs", "mdatools",
              "future", "furrr")

lapply(packages, library, character.only = TRUE)

library("conflicted")
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(recipes::prep)
conflicts_prefer(mdatools::pls)

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

# Filtering target models for PLSR experiment
# Testing only with mir and a few soil properties
modeling.combinations <- modeling.combinations  %>%
  filter(spectra_type == "mir") %>%
  filter(soil_property %in% c("oc_usda.c729_w.pct", "c.tot_usda.a622_w.pct",
                              "clay.tot_usda.a334_w.pct", "sand.tot_usda.c60_w.pct",
                              "silt.tot_usda.c62_w.pct", "bd_usda.a4_g.cm3",
                              "ph.h2o_usda.a268_index", "caco3_usda.a54_w.pct",
                              "k.ext_usda.a725_cmolc.kg"))

modeling.combinations <- modeling.combinations %>%
  mutate(model_name = gsub("mlr3..eml", "plsr", model_name))

modeling.combinations

# # Filtering already fitted models
# modeling.combinations <- modeling.combinations %>%
#   mutate(fitted = file.exists(paste0(dir, export_name, "/model_", model_name, ".qs")))
#
# modeling.combinations <- modeling.combinations  %>%
#   filter(!fitted)

## Running model tuning and fitting

i=1
for(i in 1:nrow(modeling.combinations)) {

  isoil_property = modeling.combinations[[i,"soil_property"]]
  imodel_name = modeling.combinations[[i,"model_name"]]
  iexport_name = modeling.combinations[[i,"export_name"]]
  ispectra_type = modeling.combinations[[i,"spectra_type"]]
  isubset = modeling.combinations[[i,"subset"]]
  igeo = modeling.combinations[[i,"geo"]]

  cat(paste0("Running iteration ", paste0(i, "/", nrow(modeling.combinations)), " at ", now(), "\n"))

  # Loading dataset inside loop for memory management
  mir.spectral.range <- paste0("scan_mir.", seq(600, 4000, by = 2), "_abs")
  visnir.spectral.range <- paste0("scan_visnir.", seq(400, 2500, by = 2), "_ref")
  nir.neospectra.spectral.range <- paste0("scan_nir.", seq(1350, 2550, by = 2), "_ref")

  column.ids <- c("id.layer_uuid_txt", "dataset.code_ascii_txt")

  if(ispectra_type == "mir"){

    data <- qread(paste0(db.dir, "ossl_all_L1_v1.2.qs")) %>%
      select(all_of(column.ids), all_of(isoil_property), all_of(mir.spectral.range)) %>%
      filter(!is.na(!!as.name(isoil_property))) %>%
      filter(!is.na(scan_mir.1500_abs))

  } else if(ispectra_type == "visnir") {

    data <- qread(paste0(db.dir, "ossl_all_L1_v1.2.qs")) %>%
      select(all_of(column.ids), all_of(isoil_property), all_of(visnir.spectral.range)) %>%
      filter(!is.na(!!as.name(isoil_property))) %>%
      filter(!is.na(scan_visnir.1500_ref))

  } else if(ispectra_type == "nir.neospectra") {

    data <- qread(paste0(db.dir, "ossl_all_L1_v1.2.qs")) %>%
      select(all_of(column.ids), all_of(isoil_property), all_of(nir.neospectra.spectral.range)) %>%
      filter(!is.na(!!as.name(isoil_property))) %>%
      filter(!is.na(scan_nir.1500_ref))

  }

  if(isubset == "kssl") {

    data <- data %>%
      filter(dataset.code_ascii_txt == "KSSL.SSL")

  }

  # Preprocessing

  preprocessed <- data %>%
    select(-all_of(c(column.ids, isoil_property))) %>%
    as.matrix() %>%
    # Standard Normal Variate
    prospectr::standardNormalVariate(X = .) %>%
    as_tibble() %>%
    # Rebinding sample id
    bind_cols({data %>%
        select(all_of(c(column.ids, isoil_property)))}, .)

  # CV predictions for error model

  cv.predictions <- qread(paste0(dir, iexport_name, "/cvpred_", imodel_name, ".qs"))

  # Metrics and best HPO

  performance.metrics <- cv.predictions %>%
    pivot_longer(starts_with("pred"), names_to = "components", values_to = "response") %>%
    rename(truth := !!isoil_property) %>%
    group_by(components) %>%
    summarise(n = n(),
              rmse = rmse_vec(truth = truth, estimate = response),
              bias = msd_vec(truth = truth, estimate = response),
              rsq = rsq_vec(truth = truth, estimate = response),
              ccc = ccc_vec(truth = truth, estimate = response, bias = T),
              rpiq = rpiq_vec(truth = truth, estimate = response),
              .groups = "drop") %>%
    arrange(rmse)

  best.component <- performance.metrics %>%
    filter(row_number() == 1) %>%
    pull(components)

  n.comps <- as.numeric(gsub("pred|comp", "", best.component))

  cv.predictions <- cv.predictions %>%
    select(all_of(c(column.ids, isoil_property, best.component))) %>%
    rename(predicted := !!best.component) %>%
    mutate(residual = abs(!!as.name(isoil_property)-predicted)) %>%
    select(-all_of(c(isoil_property)))

  preprocessed <- preprocessed %>%
    left_join(cv.predictions, by = c("id.layer_uuid_txt", "dataset.code_ascii_txt")) %>%
    relocate(residual, .after = all_of(column.ids)) %>%
    select(-all_of(c(isoil_property)))

  cat(paste0("Imported data at ", now(), "\n"))

  # Recipe model

  if(!grepl("log", iexport_name)) {

    recipe.model <- function(dataset){
      dataset %>%
        recipe() %>%
        update_role(everything()) %>%
        update_role(all_of(c(column.ids, "predicted")), new_role = "id") %>%
        update_role(all_of("residual"), new_role = "outcome") %>%
        prep(strings_as_factors = FALSE)
    }

  } else if(grepl("log", iexport_name)) {

    recipe.model <- function(dataset){
      dataset %>%
        recipe() %>%
        update_role(everything()) %>%
        update_role(all_of(c(column.ids, "predicted")), new_role = "id") %>%
        update_role(all_of("residual"), new_role = "outcome") %>%
        step_log(all_outcomes(), offset = 1, id = "log") %>%
        prep(strings_as_factors = FALSE)
    }

  }

  # Fitting error model and saving to disk

  training.outcome <- juice(recipe.model(preprocessed), composition = "matrix", all_outcomes())

  training.predictors <- juice(recipe.model(preprocessed), composition = "matrix", all_predictors())

  pls.model <- mdatools::pls(x = training.predictors, y = training.outcome, ncomp = n.comps,
                             scale = T, center = T, cv = NULL)

  export.model <- paste0(dir, iexport_name, "/error_model_", imodel_name, ".qs")
  export.model.alt <- gsub("\\.qs", "\\.rds", export.model)

  if(file.exists(export.model)){file.remove(export.model)}
  if(file.exists(export.model.alt)){file.remove(export.model.alt)}

  tryCatch(
    expr = {qsave(pls.model, export.model)},
    error = function(e){saveRDS(pls.model, export.model.alt)}
  )

  # Estimating conformity scores

  error.pred <- bind_cols(cv.predictions,
                          {tibble(pred_residual = pls.model$res$cal$y.pred[,n.comps,1])})

  error.pred <- error.pred %>%
    mutate(alpha_scores = residual/pred_residual)

  export.error <- paste0(dir, iexport_name, "/error_pred_", imodel_name, ".qs")
  export.error.alt <- gsub("\\.qs", "\\.rds", export.error)

  if(file.exists(export.error)){file.remove(export.error)}
  if(file.exists(export.error.alt)){file.remove(export.error.alt)}

  tryCatch(
    expr = {qsave(error.pred, export.error)},
    error = function(e){saveRDS(error.pred, export.error.alt)}
  )

  # Cleaning iteration and freeing memory

  keep.objects <- c("dir", "db.dir",
                    "modeling.combinations")

  remove.objects <- ls()[-grep(paste(keep.objects, collapse = "|"), ls())]
  rm(list = remove.objects)
  gc()

  cat(paste0("Exported final model at ", now(), "\n\n"))

}
