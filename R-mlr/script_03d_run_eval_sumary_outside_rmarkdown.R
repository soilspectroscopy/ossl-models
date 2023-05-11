
## Packages
library("tidyverse")
library("yardstick")
library("qs")

## Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

## Fitted models
fitted.modeling.combinations <- read_csv("./out/fitted_modeling_combinations_v1.2.csv",
                                         show_col_types = FALSE)

## Re-exporting plots due to inconsistencies

i=1
for(i in 1:nrow(fitted.modeling.combinations)) {

  # Parameters
  isoil_property = fitted.modeling.combinations[[i,"soil_property"]]
  imodel_name = fitted.modeling.combinations[[i,"model_name"]]
  iexport_name = fitted.modeling.combinations[[i,"export_name"]]
  ispectra_type = fitted.modeling.combinations[[i,"spectra_type"]]
  isubset = fitted.modeling.combinations[[i,"subset"]]
  igeo = fitted.modeling.combinations[[i,"geo"]]

  cat(paste0("Run ", i, "/", nrow(fitted.modeling.combinations), " at ", lubridate::now(), "\n"))

  # Predictions
  cv.results <- qread(paste0(dir,
                              iexport_name,
                              "/cvpred_",
                              imodel_name,
                              ".qs"))

  # Metrics
  performance.metrics <- cv.results %>%
    summarise(n = n(),
              rmse = rmse_vec(truth = truth, estimate = response),
              bias = msd_vec(truth = truth, estimate = response),
              rsq = rsq_vec(truth = truth, estimate = response),
              ccc = ccc_vec(truth = truth, estimate = response, bias = T),
              rpiq = rpiq_vec(truth = truth, estimate = response))

  perfomance.annotation <- paste0("Lin's CCC = ", round(performance.metrics[[1,"ccc"]], 3),
                                  "\nRMSE = ", round(performance.metrics[[1,"rmse"]], 3))

  # Plot
  if(grepl("log..", iexport_name)) {

    p.hex <- ggplot(cv.results, aes(x = truth, y = response)) +
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

    ggsave(paste0(dir,
                  iexport_name,
                  "/valplot_",
                  imodel_name,
                  ".png"),
           p.hex, dpi = 200, width = 5, height = 5, units = "in", scale = 1)

  } else {

    p.hex <- ggplot(cv.results, aes(x = truth, y = response)) +
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

    ggsave(paste0(dir,
                  iexport_name,
                  "/valplot_",
                  imodel_name,
                  ".png"),
           p.hex, dpi = 100, width = 5, height = 5, units = "in", scale = 1)

  }

}

## Full table with final performance

metrics.files <- list.files(dir, pattern = "perfmetrics_",
                            full.names = T, recursive = T)

performance <- metrics.files %>%
  purrr::map_dfr(.f = read_csv, show_col_types = FALSE)

write_csv(performance, "./out/fitted_models_performance_v1.2.csv")

# missing.files <- left_join(fitted.modeling.combinations,
#                            performance,
#                            by = c("soil_property",
#                                   "model_name"))
#
# missing.files <- missing.files %>%
#   filter(is.na(ccc))
#
# i=1
# for(i in 1:nrow(missing.files)) {
#
#   # Parameters
#   isoil_property = missing.files[[i,"soil_property"]]
#   imodel_name = missing.files[[i,"model_name"]]
#   iexport_name = missing.files[[i,"export_name"]]
#   ispectra_type = missing.files[[i,"spectra_type"]]
#   isubset = missing.files[[i,"subset"]]
#   igeo = missing.files[[i,"geo"]]
#
#   cat(paste0("Run ", i, "/", nrow(missing.files), " at ", lubridate::now(), "\n"))
#
#   # Predictions
#   cv.results <- qread(paste0(dir,
#                              iexport_name,
#                              "/cvpred_",
#                              imodel_name,
#                              ".qs"))
#
#   performance.metrics <- cv.results %>%
#     summarise(n = n(),
#               rmse = rmse_vec(truth = truth, estimate = response),
#               bias = msd_vec(truth = truth, estimate = response),
#               rsq = rsq_vec(truth = truth, estimate = response),
#               ccc = ccc_vec(truth = truth, estimate = response, bias = T),
#               rpiq = rpiq_vec(truth = truth, estimate = response))
#
#   performance.metrics <- performance.metrics %>%
#     mutate(soil_property = isoil_property,
#            model_name = imodel_name,
#            .before = 1)
#
#   export.file <- paste0(dir,
#                        iexport_name,
#                        "/perfmetrics_",
#                        imodel_name,
#                        ".csv")
#
#   if(file.exists(export.file)){file.remove(export.file)}
#
#   write_csv(performance.metrics, export.file)
#
# }


