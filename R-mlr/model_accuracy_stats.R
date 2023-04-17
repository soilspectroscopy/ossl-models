# zecojls@gmail.com

library("tidyverse")
library("lubridate")
library("tidymodels")

## Goodness-of-fit statistics

# library("hexbin")
# library("plotKML")
# library("lattice")
# library("plyr")
# source("./R-mlr/SSL_functions.R")

t.vars = c("log..oc_usda.c729_w.pct",
           "log..n.tot_usda.a623_w.pct",
           "silt.tot_usda.c62_w.pct",
           "clay.tot_usda.a334_w.pct",
           "sand.tot_usda.a60_w.pct",
           "log..cec_usda.a723_cmolc.kg",
           "ph.h2o_usda.a268_index",
           "ph.cacl2_usda.a481_index",
           "log..al.ox_usda.a59_w.pct",
           "log..caco3_usda.a54_w.pct",
           "log..k.ext_usda.a725_cmolc.kg",
           "log..mg.ext_usda.a724_cmolc.kg",
           "log..ca.ext_usda.a722_cmolc.kg",
           "acidity_usda.a795_cmolc.kg",
           "bd_usda.a4_g.cm3")

mn.lst = c("mir_mlr..eml_kssl_ll_v1.2.rds",
           "mir_mlr..eml_kssl_na_v1.2.rds",
           "mir_mlr..eml_ossl_ll_v1.2.rds",
           "nir.neospectra_mlr..eml_ossl_ll_v1.2.rds")

model.combinations <- crossing(soil_property = t.vars,
                               model_type = mn.lst)

performance.metrics.list <- list()

i=1
for(i in 1:nrow(model.combinations)) {

  isoil_property <- model.combinations[[i,"soil_property"]]
  imodel_type <- model.combinations[[i,"model_type"]]

  base.url <- "http://s3.us-east-1.wasabisys.com/soilspectroscopy/ossl_models/"

  in.rds = readRDS(url(paste0(base.url, isoil_property, "/", imodel_type, ".rds"), "rb"))

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
           model_type = imodel_type,
           .before = 1)

  performance.metrics.list[[i]] <- performance.metrics

  rm(in.rds)
  gc()

  cat(paste0("Run ", isoil_property, ", ", imodel_type, " at ", now(), "\n"))

}

# Export results

## Accuracy plots

# for(i in 1:nrow(acc.mat)){
#   for(j in 1:length(mn.lst)){
#     in.rds = paste0("/mnt/landmark/ossl_mlr/models/", t.vars[i], "/", mn.lst[j], ".rds")
#     out.file = paste0("/mnt/landmark/ossl_mlr/models/", t.vars[i], "/ap.", mn.lst[j], ".rds.png")
#     if(file.exists(in.rds) & !file.exists(out.file)){
#       t.m = readRDS.gz(in.rds)
#       yh = t.m$learner.model$super.model$learner.model$fitted.values
#       meas = t.m$learner.model$super.model$learner.model$model[,t.vars[i]]
#       t.var.breaks = quantile(meas, c(0.001, 0.01, 0.999), na.rm=TRUE)
#       plot_hexbin(varn=t.vars[i], breaks=c(t.var.breaks[1], seq(t.var.breaks[2], t.var.breaks[3], length=25)), meas=ifelse(meas<0, 0, meas), pred=ifelse(yh<0, 0, yh), main=t.vars[i], out.file=out.file, log.plot=FALSE, colorcut=c(0,0.01,0.02,0.03,0.06,0.12,0.20,0.35,1.0))
#       #gc()
#     }
#   }
# }
#
# ## global layers ----
# cog.lst = list.files("/data/WORLDCLIM", ".tif", full.names = TRUE)
# write.csv(data.frame(filename=basename(cog.lst)), "./out/global_layers1km.csv")
#
# ## prepare sample data
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
