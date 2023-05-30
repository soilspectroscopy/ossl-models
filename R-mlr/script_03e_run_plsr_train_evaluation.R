
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

  preprocessed <- data %>%
    select(-all_of(c(column.ids, isoil_property))) %>%
    as.matrix() %>%
    # Standard Normal Variate
    prospectr::standardNormalVariate(X = .) %>%
    as_tibble() %>%
    # Rebinding sample id
    bind_cols({data %>%
        select(all_of(c(column.ids, isoil_property)))}, .)

  cat(paste0("Imported data at ", now(), "\n"))

  # Splitting into 10-folds for HPO
  set.seed(1993)
  modeling.folds <- preprocessed %>%
    filter(!is.na(!!as.name(isoil_property))) %>%
    vfold_cv(v = 10, repeats = 1) %>%
    unite(idfull, starts_with("id"), sep = "_")

  # Recipe model

  if(!grepl("log", iexport_name)) {

    recipe.model <- function(dataset){
      dataset %>%
        recipe() %>%
        update_role(everything()) %>%
        update_role(all_of(column.ids), new_role = "id") %>%
        update_role(all_of(isoil_property), new_role = "outcome") %>%
        prep(strings_as_factors = FALSE)
    }

  } else if(grepl("log", iexport_name)) {

    recipe.model <- function(dataset){
      dataset %>%
        recipe() %>%
        update_role(everything()) %>%
        update_role(all_of(column.ids), new_role = "id") %>%
        update_role(all_of(isoil_property), new_role = "outcome") %>%
        step_log(all_outcomes(), offset = 1, id = "log") %>%
        prep(strings_as_factors = FALSE)
    }

  }

  # Prediction function

  model.prediction.folds <- function(maxcomps = 30, split, id){

    # maxcomps = 30
    # split = modeling.folds[["splits"]][[1]]
    # id=1

    # Preparing pls matrices

    training.set <- analysis(split)

    training.outcome <- juice(recipe.model(training.set), composition = "matrix", all_outcomes())

    training.predictors <- juice(recipe.model(training.set), composition = "matrix", all_predictors())

    pls.model <- mdatools::pls(x = training.predictors, y = training.outcome, ncomp = maxcomps,
                               scale = T, center = T, cv = NULL)

    # Evaluation

    testing.set <- assessment(split)

    testing.outcome <- bake(recipe.model(training.set),
                            new_data = testing.set,
                            composition = "tibble", all_of(column.ids), all_outcomes())

    testing.predictors <- bake(recipe.model(training.set),
                               new_data = testing.set,
                               composition = "matrix", all_predictors())

    predictions <- predict(pls.model, x = testing.predictors)

    predictions$y.pred %>%
      as.data.frame() %>%
      as_tibble() %>%
      rename_with(~paste0("pred", seq(1, maxcomps, by=1), "comp"), everything()) %>%
      bind_cols(testing.outcome, .)

  }

  # future::plan("multisession", workers = 3, gc = FALSE)

  # cv.results <- future_map2_dfr(.x = modeling.folds$splits,
  #                               .y = modeling.folds$idfull,
  #                               ~model.prediction.folds(maxcomps = 30, split = .x, id = .y),
  #                               .options = furrr_options(seed = T))

  cv.results <- map2_dfr(.x = modeling.folds$splits,
                         .y = modeling.folds$idfull,
                         ~model.prediction.folds(maxcomps = 30, split = .x, id = .y))

  # future:::ClusterRegistry("stop")

  # Exporting results
  export.file <- paste0(dir,
                        iexport_name,
                        "/cvpred_",
                        imodel_name,
                        ".qs")

  export.file.alt <- gsub("\\.qs", "\\.csv", export.file)

  if(file.exists(export.file)){file.remove(export.file)}
  if(file.exists(export.file.alt)){file.remove(export.file.alt)}

  tryCatch(
    expr = {qsave(cv.results, export.file)},
    error = function(e){write_csv(cv.results, export.file.alt)}
    )

  cat(paste0("Exported CV results at ", now(), "\n"))

  # Metrics and best HPO

  performance.metrics <- cv.results %>%
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

  perfomance.annotation <- paste0("Lin's CCC = ", round(performance.metrics[[1,"ccc"]], 3),
                                  "\nRMSE = ", round(performance.metrics[[1,"rmse"]], 3))

  performance.metrics <- performance.metrics %>%
    mutate(soil_property = isoil_property,
           model_name = imodel_name,
           .before = 1)

  export.metrics <- paste0(dir,
                           iexport_name,
                           "/perfmetrics_",
                           imodel_name,
                           ".csv")

  if(file.exists(export.metrics)){file.remove(export.metrics)}

  write_csv(performance.metrics, export.metrics)

  best.hps <- performance.metrics %>%
    filter(row_number() == 1) %>%
    pull(components)

  # Plots

  if(grepl("log..", iexport_name)) {

    p.hex <- cv.results %>%
      select(all_of(isoil_property), all_of(best.hps)) %>%
      rename(truth := !!isoil_property, response := !!best.hps) %>%
      ggplot(aes(x = truth, y = response)) +
      geom_hex(bins = 30, alpha = 0.75) +
      geom_abline(intercept = 0, slope = 1) +
      geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2),
                label = perfomance.annotation) +
      scale_fill_viridis_c(trans = "log10") +
      labs(x = "log(observed)", y = "log(predicted)", fill = bquote(log[10](count)),
           title = iexport_name) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.key.size = unit(0.35, "cm"))

    r.max <- max(layer_scales(p.hex)$x$range$range)
    r.min <- min(layer_scales(p.hex)$x$range$range)

    s.max <-max(layer_scales(p.hex)$y$range$range)
    s.min <-min(layer_scales(p.hex)$y$range$range)

    t.max <-round(max(r.max,s.max),1)
    t.min <-round(min(r.min,s.min),1)

    p.hex <- p.hex + coord_equal(xlim=c(t.min,t.max),ylim=c(t.min,t.max))

    export.plot <- paste0(dir,
                         iexport_name,
                         "/valplot_",
                         imodel_name,
                         ".png")

    if(file.exists(export.plot)){file.remove(export.plot)}

    ggsave(export.plot,
           p.hex, dpi = 200, width = 5, height = 5, units = "in", scale = 1)

  } else {

    p.hex <- cv.results %>%
      select(all_of(isoil_property), all_of(best.hps)) %>%
      rename(truth := !!isoil_property, response := !!best.hps) %>%
      ggplot(aes(x = truth, y = response)) +
      geom_hex(bins = 30, alpha = 0.75) +
      geom_abline(intercept = 0, slope = 1) +
      geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2),
                label = perfomance.annotation) +
      scale_fill_viridis_c(trans = "log10") +
      labs(x = "observed", y = "predicted", fill = bquote(log[10](count)),
           title = iexport_name) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.key.size = unit(0.35, "cm"))

    r.max <- max(layer_scales(p.hex)$x$range$range)
    r.min <- min(layer_scales(p.hex)$x$range$range)

    s.max <-max(layer_scales(p.hex)$y$range$range)
    s.min <-min(layer_scales(p.hex)$y$range$range)

    t.max <-round(max(r.max,s.max),1)
    t.min <-round(min(r.min,s.min),1)

    p.hex <- p.hex + coord_equal(xlim=c(t.min,t.max),ylim=c(t.min,t.max))

    export.plot <- paste0(dir,
                          iexport_name,
                          "/valplot_",
                          imodel_name,
                          ".png")

    if(file.exists(export.plot)){file.remove(export.plot)}

    ggsave(export.plot,
           p.hex, dpi = 200, width = 5, height = 5, units = "in", scale = 1)

  }

  cat(paste0("Exported CV performance and plot at ", now(), "\n"))

  # Saving final model to disk

  training.outcome <- juice(recipe.model(preprocessed), composition = "matrix", all_outcomes())

  training.predictors <- juice(recipe.model(preprocessed), composition = "matrix", all_predictors())

  final.n.comps <- as.numeric(gsub("pred|comp", "", best.hps))

  pls.model <- mdatools::pls(x = training.predictors, y = training.outcome, ncomp = final.n.comps,
                             scale = T, center = T, cv = NULL)

  # # Reducing size to export
  # pls.model <- pls.model[c("Xmeans", "Ymeans", "ncomp", "center", "scale", "coefficients", "terms")]
  # class(pls.model) <- "mvr"

  export.model <- paste0(dir, iexport_name, "/model_", imodel_name, ".qs")
  export.model.alt <- gsub("\\.qs", "\\.rds", export.model)

  if(file.exists(export.model)){file.remove(export.model)}
  if(file.exists(export.model.alt)){file.remove(export.model.alt)}

  tryCatch(
    expr = {qsave(pls.model, export.model)},
    error = function(e){saveRDS(pls.model, export.model.alt)}
    )

  # Cleaning iteration and freeing memory

  keep.objects <- c("dir", "db.dir",
                    "modeling.combinations", "n.cores")

  remove.objects <- ls()[-grep(paste(keep.objects, collapse = "|"), ls())]
  rm(list = remove.objects)
  gc()

  cat(paste0("Exported final model at ", now(), "\n\n"))

}
