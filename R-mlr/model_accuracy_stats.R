## Produce summary statistics and plots
## tom.hengl@opengeohub.org

library(hexbin)
library(plotKML)
library(lattice)
library(plyr)
source("./R-mlr/SSL_functions.R")

t.vars = c("log..oc_usda.calc_wpct", "log..n.tot_usda.4h2_wpct", "silt.tot_usda.3a1_wpct",
           "clay.tot_usda.3a1_wpct", "sand.tot_usda.3a1_wpct", "log..ecec_usda.4b4_cmolkg",
           "ph.h2o_usda.4c1_index", "ph.cacl2_usda.4c1_index", "log..al.kcl_usda.4b3_cmolkg",
           "log..k.ext_usda.4b1_cmolkg", "log..caco3_usda.4e1_wpct", "log..mg.ext_usda.4b1_cmolkg",
           "log..ca.ext_usda.4b1_cmolkg", "log..gyp_usda.4e2_wpct",
           "log..cec.ext_usda.4b1_cmolkg", "bsat_usda.4b4_wpct", "bd.od_usda.3b2_gcm3")
mn.lst = c("mir_mlr..eml_kssl_na_v1.1", "mir_mlr..eml_kssl_ll_v1.1", "visnir_mlr..eml_ossl_na_v1.1",
           "mir_mlr..eml_ossl_ll_v1.1", "visnir.mir_mlr..eml_ossl_na_v1.1", "visnir.mir_mlr..eml_ossl_ll_v1.1")
p.lst = c("RMSE", "R.square", "N.tot", "N.outliers")
acc.mat = data.frame(matrix(nrow=length(t.vars), ncol=2+length(mn.lst)*4))
colnames(acc.mat) = c("variable", "std", sapply(mn.lst, function(i){paste0(i, "_", p.lst)}))
acc.mat$variable = t.vars
acc.mat$std = sapply(t.vars, function(i){ if(length(grep("log..", i))>0) { sd(log1p(rm.ossl[,gsub("log..","",i)]), na.rm = TRUE) } else { sd(rm.ossl[,i], na.rm = TRUE) } })
for(i in 1:nrow(acc.mat)){
  for(j in 1:length(mn.lst)){
    in.rds = paste0("/mnt/landmark/ossl_mlr/models/", t.vars[i], "/", mn.lst[j], ".rds")
    if(file.exists(in.rds)){
      if(is.na(acc.mat[i,paste0(mn.lst[j], "_RMSE")])){
        t.m = readRDS.gz(in.rds)
        x.s = summary(t.m$learner.model$super.model$learner.model)
        RMSE = signif(sqrt(sum(t.m$learner.model$super.model$learner.model$residuals^2) / t.m$learner.model$super.model$learner.model$df.residual), 3)
        acc.mat[i,paste0(mn.lst[j], "_RMSE")] = RMSE
        acc.mat[i,paste0(mn.lst[j], "_R.square")] = round(x.s$adj.r.squared, 3)
        acc.mat[i,paste0(mn.lst[j], "_N.tot")] = x.s$df[2]
        acc.mat[i,paste0(mn.lst[j], "_N.outliers")] = sum(abs(x.s$residuals) > 3*RMSE)
        gc()
      }
    }
  }
}
write.csv(acc.mat, "out/accuracy_matrix_ossl_models.csv")
acc.matL = data.frame(matrix(nrow=length(t.vars)*5, ncol=5))
names(acc.matL) = c("variable", "model", "R.square", "RMSE", "N.tot")
acc.matL$variable = as.vector(sapply(t.vars, function(i){rep(i, 5)}))
acc.matL$model = as.vector(rep(mn.lst[-5], length(t.vars)))
for(i in 1:nrow(acc.matL)){
  acc.matL[i,3] = acc.mat[which(acc.mat$variable==acc.matL$variable[i]), which(names(acc.mat)==paste0(acc.matL$model[i], "_R.square"))]
  acc.matL[i,4] = acc.mat[which(acc.mat$variable==acc.matL$variable[i]), which(names(acc.mat)==paste0(acc.matL$model[i], "_RMSE"))]
  acc.matL[i,5] = acc.mat[which(acc.mat$variable==acc.matL$variable[i]), which(names(acc.mat)==paste0(acc.matL$model[i], "_N.tot"))]
}
write.csv(acc.matL, "out/accuracy_list_long.csv")

## Accuracy plots ----
for(i in 1:nrow(acc.mat)){
  for(j in 1:length(mn.lst)){
    in.rds = paste0("/mnt/landmark/ossl_mlr/models/", t.vars[i], "/", mn.lst[j], ".rds")
    out.file = paste0("/mnt/landmark/ossl_mlr/models/", t.vars[i], "/ap.", mn.lst[j], ".rds.png")
    if(file.exists(in.rds) & !file.exists(out.file)){
      t.m = readRDS.gz(in.rds)
      yh = t.m$learner.model$super.model$learner.model$fitted.values
      meas = t.m$learner.model$super.model$learner.model$model[,t.vars[i]]
      t.var.breaks = quantile(meas, c(0.001, 0.01, 0.999), na.rm=TRUE)
      plot_hexbin(varn=t.vars[i], breaks=c(t.var.breaks[1], seq(t.var.breaks[2], t.var.breaks[3], length=25)), meas=ifelse(meas<0, 0, meas), pred=ifelse(yh<0, 0, yh), main=t.vars[i], out.file=out.file, log.plot=FALSE, colorcut=c(0,0.01,0.02,0.03,0.06,0.12,0.20,0.35,1.0))
      #gc()
    }
  }
}

## global layers ----
cog.lst = list.files("/data/WORLDCLIM", ".tif", full.names = TRUE)
write.csv(data.frame(filename=basename(cog.lst)), "./out/global_layers1km.csv")

## prepare sample data
new.data = vroom::vroom("/mnt/diskstation/data/ossl/dataset/validation/JSset_KSSL.csv", n_max = 20)
dim(new.data)
names(new.data)[1:40]
mir.raw = new.data[,35:ncol(new.data)]
#write.csv(new.data[,35:ncol(new.data)], "/home/tomislav/Documents/git/ossl_modeling/api/server/R-code/models/sample_mir_data.csv")
#write.csv(new.data[,c(1,4,5,32,33)], "/home/tomislav/Documents/git/ossl_modeling/api/server/R-code/models/sample_soilsite_data.csv")

new.nirdata = vroom::vroom("/mnt/diskstation/data/ossl/dataset/KSSL/VNIR_Spectra_Library_spectra.csv", n_max = 20)
dim(new.nirdata)
names(new.nirdata)[1:40]
visnir.raw = new.nirdata[,-2]
write.csv(visnir.raw, "./api/server/R-code/models/sample_visnir_data.csv")

## Validation data
mir.raw = read.csv("./api/server/R-code/models/sample_mir_data.csv")[,-1]
new.soil = read.csv("./api/server/R-code/models/sample_soilsite_data.csv")
lon = new.soil$longitude.std.decimal.degrees
lat = new.soil$latitude.std.decimal.degrees
hzn_depth = new.soil$lay_depth_to_top + (new.soil$lay_depth_to_bottom - new.soil$lay_depth_to_top)/2
visnir.raw = read.csv("./api/server/R-code/models/sample_visnir_data.csv")[,-c(1:2)]
str(visnir.raw[,1:5])

## load models:
ossl.pca.mir = readRDS.gz("/mnt/landmark/ossl_mlr/pca.ossl/mpca_mir_kssl_v1.rds")
ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/mir_mlr..eml_kssl_na_v1.rds")
## 1: predict using mir only:
pred.oc = predict.ossl(t.var="log..oc_usda.calc_wpct", mir.raw=mir.raw, ossl.model=ossl.model, ossl.pca.mir=ossl.pca.mir, ylim=c(0,100))
str(pred.oc, max.level = 1)
## 2: predict using mir with geo:
ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/mir_mlr..eml_kssl_ll_v1.rds")
pred.oc2 = predict.ossl(t.var="log..oc_usda.calc_wpct", mir.raw=mir.raw, ossl.model=ossl.model, ylim=c(0,100),
             ossl.pca.mir=ossl.pca.mir, geo.type="ll", lon=lon, lat=lat, hzn_depth=hzn_depth)
## 3: predict using visnir:
ossl.pca.visnir = readRDS.gz("/mnt/landmark/ossl_mlr/pca.ossl/mpca_visnir_kssl_v1.rds")
ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/visnir_mlr..eml_ossl_na_v1.rds")
pred.oc3 = predict.ossl(t.var="log..oc_usda.calc_wpct", visnir.raw=visnir.raw, ossl.model=ossl.model,
                        ylim=c(0,100), spc.type = "visnir", ossl.pca.visnir = ossl.pca.visnir)
## 4: predict using mir (ossl subset) with geo:
ossl.model = readRDS.gz("/mnt/landmark/ossl_mlr/models/log..oc_usda.calc_wpct/mir_mlr..eml_ossl_ll_v1.rds")
str(ossl.model$features)
pred.oc3 = predict.ossl(t.var="log..oc_usda.calc_wpct", mir.raw=mir.raw, ossl.model=ossl.model, ylim=c(0,100),
                        ossl.pca.mir=ossl.pca.mir, geo.type="ll", lon=lon, lat=lat, hzn_depth=hzn_depth)

## test overlay using COG on Wasabi
for(j in 1:length(cog.lst)){
  #x = terra::extract(terra::rast(paste0("/vsicurl/", url, cog.lst[j])), terra::vect(pnts[1,]))
  x = lapply(1:length(cog.lst), function(j) terra::extract(terra::rast(paste0("/vsicurl/", url, cog.lst[j])), data.frame(lon, lat)[5,]))
  print(x)
}

## https://stackoverflow.com/questions/29958561/how-to-plot-points-on-hexbin-graph-in-r
## large file
library(hexbin)
library(grid)
ossl.pcax.mir = readRDS.gz("/mnt/landmark/ossl_mlr/pca.ossl/pca_mir_kssl_v1.rds")
reds = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[-1])
hb <- hexbin(ossl.pcax.mir$x[,1:2], xbins=60)
p <- plot(hb, colramp = reds, main='PCA MIR KSSL')
pushHexport(p$plot.vp)
grid.points(X1.pc$mir.PC1, X1.pc$mir.PC2, pch=17)

hb2 <- hexbin(ossl.pcax.mir$x[,3:4], xbins=60)
#op <- par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
p2 <- plot(hb2, colramp = reds, main='PCA MIR KSSL')
pushHexport(p2$plot.vp)
grid.points(X1.pc$mir.PC3, X1.pc$mir.PC4, pch=17)
upViewport()

## move files
#in.rds = list.files("/mnt/soilspec4gg/ossl/ossl_models", full.names = TRUE, recursive = TRUE)
## 314
#try( file.copy(from=gsub("/mnt/soilspec4gg/ossl/ossl_models/", "/mnt/landmark/ossl_mlr/models/", gsub("v1", "v1.1", in.rds)), to=in.rds, overwrite = TRUE) )
