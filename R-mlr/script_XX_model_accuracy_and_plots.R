# zecojls@gmail.com

options(timeout=6000)

library("tidyverse")
library("lubridate")
library("tidymodels")
library("ggpubr")

listed.models <- read_csv("out/list_ossl_models_v1.2.csv")

## Goodness-of-fit statistics

names(listed.models)

performance.metrics.list <- list()

i=1
for(i in 1:nrow(listed.models)) {

  isoil_property <- listed.models[[i,"variable"]]
  imodel_name <- listed.models[[i,"model_name"]]
  imodel_url <- listed.models[[i,"model_url"]]

  in.rds = readRDS(url(imodel_url), "rb")

  predictions <- tibble(observed = (in.rds$learner.model$super.model$learner.model$residuals+
                                      in.rds$learner.model$super.model$learner.model$fitted.values),
                        predicted = in.rds$learner.model$super.model$learner.model$fitted.values)

  performance.metrics <- predictions %>%
    summarise(n = n(),
              rmse = rmse_vec(truth = observed, estimate = predicted),
              bias = msd_vec(truth = observed, estimate = predicted),
              rsq = rsq_vec(truth = observed, estimate = predicted),
              ccc = ccc_vec(truth = observed, estimate = predicted, bias = T),
              rpiq = rpiq_vec(truth = observed, estimate = predicted)) %>%
    mutate(soil_property = isoil_property,
           model_name = imodel_name,
           .before = 1)

  performance.metrics.list[[i]] <- performance.metrics

  rm(in.rds)
  gc()

  cat(paste0("Run ", isoil_property, ", ", imodel_name, " at ", now(), "\n"))

}

performance.metrics <- Reduce(bind_rows, performance.metrics.list)
performance.metrics

# Export results
write_csv(performance.metrics, "out/ossl_models_v1.2_performance.csv")

# performance.metrics <- read_csv("out/ossl_models_v1.2_performance.csv")
# clipr::write_clip(performance.metrics)

## Accuracy plots
listed.models <- read_csv("out/list_ossl_models_v1.2.csv")
performance.metrics <- read_csv("out/ossl_models_v1.2_performance.csv")

performance.metrics <- performance.metrics %>%
  mutate_if(is.numeric, round, 3)

i=1
for(i in 1:nrow(performance.metrics)) {

  isoil_property <- listed.models[[i,"variable"]]
  imodel_name <- listed.models[[i,"model_name"]]
  imodel_url <- listed.models[[i,"model_url"]]

  iperformance.metrics <- performance.metrics %>%
    filter(soil_property == isoil_property) %>%
    filter(model_name == imodel_name) %>%
    select(all_of(c("n", "rmse", "bias", "rsq", "ccc", "rpiq"))) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    mutate(concat = paste0(metric, ": ", value)) %>%
    pull(concat) %>%
    paste(., collapse = ", ")

  in.rds = readRDS(url(imodel_url), "rb")

  predictions <- tibble(observed = (in.rds$learner.model$super.model$learner.model$residuals+
                                      in.rds$learner.model$super.model$learner.model$fitted.values),
                        predicted = in.rds$learner.model$super.model$learner.model$fitted.values)

  p.hex <- ggplot(predictions, aes(x = observed, y = predicted)) +
    geom_hex(bins = 30, alpha = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    scale_fill_viridis_c(trans = "log10") +
    labs(x = "Observed", y = "Predicted", fill = bquote(log[10](count)),
         title = isoil_property, subtitle = iperformance.metrics) +
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

  ggsave(paste0("out/validation_plots/", isoil_property, "..", imodel_name, ".png"),
         p.hex, width = 6, height = 6, units = "in", scale = 1)

  rm(in.rds)
  gc()

  cat(paste0("Run ", isoil_property, ", ", imodel_name, " at ", now(), "\n"))

}


## Global layers ----
# cog.lst = list.files("/data/WORLDCLIM", ".tif", full.names = TRUE)
# write.csv(data.frame(filename=basename(cog.lst)), "./out/global_layers1km.csv")

# ## Prepare sample data
# new.data = vroom::vroom("/mnt/diskstation/data/ossl/dataset/validation/JSset_KSSL.csv", n_max = 20)
# dim(new.data)
# names(new.data)[1:40]
# mir.raw = new.data[,35:ncol(new.data)]
# #write.csv(new.data[,35:ncol(new.data)], "/home/tomislav/Documents/git/ossl_modeling/api/server/R-code/models/sample_mir_data.csv")
# #write.csv(new.data[,c(1,4,5,32,33)], "/home/tomislav/Documents/git/ossl_modeling/api/server/R-code/models/sample_soilsite_data.csv")
#
# new.nirdata = vroom::vroom("/mnt/diskstation/data/ossl/dataset/KSSL/VNIR_Spectra_Library_spectra.csv", n_max = 20)
# dim(new.nirdata)
# names(new.nirdata)[1:40]
# visnir.raw = new.nirdata[,-2]
# write.csv(visnir.raw, "./api/server/R-code/models/sample_visnir_data.csv")
#
# ## Validation data
# mir.raw = read.csv("./api/server/R-code/models/sample_mir_data.csv")[,-1]
# new.soil = read.csv("./api/server/R-code/models/sample_soilsite_data.csv")
# lon = new.soil$longitude.std.decimal.degrees
# lat = new.soil$latitude.std.decimal.degrees
# hzn_depth = new.soil$lay_depth_to_top + (new.soil$lay_depth_to_bottom - new.soil$lay_depth_to_top)/2
# visnir.raw = read.csv("./api/server/R-code/models/sample_visnir_data.csv")[,-c(1:2)]
# str(visnir.raw[,1:5])
#
# ## load models:
# ossl.pca.mir = readRDS.gz("/mnt/landmark/ossl_mlr/pca.ossl/mpca_mir_kssl_v1.rds")
# ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/mir_mlr..eml_kssl_na_v1.rds")
# ## 1: predict using mir only:
# pred.oc = predict.ossl(t.var="log..oc_usda.calc_wpct", mir.raw=mir.raw, ossl.model=ossl.model, ossl.pca.mir=ossl.pca.mir, ylim=c(0,100))
# str(pred.oc, max.level = 1)
# ## 2: predict using mir with geo:
# ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/mir_mlr..eml_kssl_ll_v1.rds")
# pred.oc2 = predict.ossl(t.var="log..oc_usda.calc_wpct", mir.raw=mir.raw, ossl.model=ossl.model, ylim=c(0,100),
#              ossl.pca.mir=ossl.pca.mir, geo.type="ll", lon=lon, lat=lat, hzn_depth=hzn_depth)
# ## 3: predict using visnir:
# ossl.pca.visnir = readRDS.gz("/mnt/landmark/ossl_mlr/pca.ossl/mpca_visnir_kssl_v1.rds")
# ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/visnir_mlr..eml_ossl_na_v1.rds")
# pred.oc3 = predict.ossl(t.var="log..oc_usda.calc_wpct", visnir.raw=visnir.raw, ossl.model=ossl.model,
#                         ylim=c(0,100), spc.type = "visnir", ossl.pca.visnir = ossl.pca.visnir)
# ## 4: predict using mir (ossl subset) with geo:
# ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/mir_mlr..eml_ossl_ll_v1.rds")
# str(ossl.model$features)
# pred.oc3 = predict.ossl(t.var="log..oc_usda.calc_wpct", mir.raw=mir.raw, ossl.model=ossl.model, ylim=c(0,100),
#                         ossl.pca.mir=ossl.pca.mir, geo.type="ll", lon=lon, lat=lat, hzn_depth=hzn_depth)
#
# ## test overlay using COG on Wasabi
# for(j in 1:length(cog.lst)){
#   #x = terra::extract(terra::rast(paste0("/vsicurl/", url, cog.lst[j])), terra::vect(pnts[1,]))
#   x = lapply(1:length(cog.lst), function(j) terra::extract(terra::rast(paste0("/vsicurl/", url, cog.lst[j])), data.frame(lon, lat)[5,]))
#   print(x)
# }
#
# ## https://stackoverflow.com/questions/29958561/how-to-plot-points-on-hexbin-graph-in-r
# ## large file
# library(hexbin)
# library(grid)
# ossl.pcax.mir = readRDS.gz("/mnt/landmark/ossl_mlr/pca.ossl/pca_mir_kssl_v1.rds")
# reds = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[-1])
# hb <- hexbin(ossl.pcax.mir$x[,1:2], xbins=60)
# p <- plot(hb, colramp = reds, main='PCA MIR KSSL')
# pushHexport(p$plot.vp)
# grid.points(X1.pc$mir.PC1, X1.pc$mir.PC2, pch=17)
#
# hb2 <- hexbin(ossl.pcax.mir$x[,3:4], xbins=60)
# #op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
# p2 <- plot(hb2, colramp = reds, main='PCA MIR KSSL')
# pushHexport(p2$plot.vp)
# grid.points(X1.pc$mir.PC3, X1.pc$mir.PC4, pch=17)
# upViewport()
#
# ## move files
# #in.rds = list.files("/mnt/soilspec4gg/ossl/ossl_models", full.names = TRUE, recursive = TRUE)
# ## 314
# #try( file.copy(from=gsub("/mnt/soilspec4gg/ossl/ossl_models/", "/mnt/landmark/ossl_mlr/models/", gsub("v1", "v1.1", in.rds)), to=in.rds, overwrite = TRUE) )
