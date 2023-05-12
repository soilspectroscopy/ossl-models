
## Packages
library("tidyverse")
library("qs")

# https://github.com/pierreroudier/asdreader
library("asdreader")

# https://github.com/spectral-cockpit/opusreader2
# install.packages("opusreader2", repos = c(
#   spectralcockpit = 'https://spectral-cockpit.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))
library("opusreader2")


## Samples
opus.1 <- read_opus_single(dsn = "sample-data/235157XS01.0")
plot(as.numeric(opus.1$ab$data), type = 'l')

opus.2 <- read_opus_single(dsn = "sample-data/icr056141.0")
plot(as.numeric(opus.2$ab$data), type = 'l')

asd.1 <- get_spectra("sample-data/101453MD01.asd")
asd.1
plot(as.numeric(asd.1), type = 'l')

asd.2 <- get_spectra("sample-data/235157MD01.asd")
asd.2
dim(asd.2)
class(asd.2)
plot(as.numeric(asd.2), type = 'l')

## csv samples
mir <- read_csv("sample-data/sample_mir_data.csv", show_col_types = F)
visnir <- read_csv("sample-data/sample_visnir_data.csv", show_col_types = F)
nir <- read_csv("sample-data/sample_neospectra_data.csv", show_col_types = F)

## Predict
predict.ossl <- function(target = "clay.tot_usda.a334_w.pct",
                         spectra = NULL,
                         spectra.type = "mir",
                         subset.type = "ossl",
                         geo.type = "na",
                         ncomps = 120,
                         models.dir = "~/mnt-ossl/ossl_models/"){

  ## Reference table
  reference <- read_csv("")

  ## Check if input scans pass some minimum checks
  if(spectra.type == "mir"){
    if(!any(class(spectra.raw) == "data.frame")){
      stop("Input dataset '*.raw' not a correctly formated scan file. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
    if(nrow(mir.raw) > 1000 | ncol(mir.raw) < 1400){
      stop("Input dataset '*.raw' dimensions invalid. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
  } else if(spc.type == "visnir"){
    if(!any(class(visnir.raw)=="data.frame")){
      stop("Input dataset '*.raw' not a correctly formated scan file. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
    if(nrow(visnir.raw) > 1000 | ncol(visnir.raw) < 1000){
      stop("Input dataset '*.raw' dimensions invalid. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
  } else if(spc.type == "nir.neospectra"){
    if(!any(class(visnir.raw)=="data.frame")){
      stop("Input dataset '*.raw' not a correctly formated scan file. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
    if(nrow(visnir.raw) > 1000 | ncol(visnir.raw) < 1000){
      stop("Input dataset '*.raw' dimensions invalid. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
  }

  if(missing(ossl.model)){
    model.url = paste0("http://s3.us-east-1.wasabisys.com/soilspectroscopy/ossl_models/", t.var, "/", spc.type, "_mlr..eml_", subset.type, "_", geo.type, "_v1.rds")
    ossl.model = readRDS(url(model.rds, "rb"))
  }

  ## convert to PC scores
  if(spc.type == "mir" | spc.type == "visnir.mir"){
    wn = as.numeric(gsub("X", "", names(mir.raw)))
    spc = as.matrix(mir.raw)
    #colnames(spc) = paste(wn)
    spc = as.data.frame(prospectr::resample(spc, wn, seq(600, 4000, by=2), interpol = "spline"))
    spc = lapply(spc, function(j){ round(ifelse(j<0, NA, ifelse(j>3, 3, j))*1000) })
    spc = as.data.frame(do.call(cbind, spc))
    names(spc) = paste0("scan_mir.", seq(600, 4000, by=2), "_abs")
    ## convert to 1st der:
    spc = prospectr::gapDer(spc, m=1, w=5, s=1, delta.wav=2)
    class(ossl.pca.mir) = "prcomp"
    X1.pc = as.data.frame(predict(ossl.pca.mir, newdata=as.data.frame(spc)))[,1:n.spc]
    colnames(X1.pc) = paste0("mir.PC", 1:n.spc)
  } else {
    X1.pc = NA
  }
  if(spc.type == "visnir" | spc.type == "visnir.mir"){
    wn = as.numeric(gsub("X", "", names(visnir.raw)))
    spc = as.matrix(visnir.raw)
    #colnames(spc) = paste(wn)
    spc = as.data.frame(prospectr::resample(spc, wn, seq(350, 2500, by=2), interpol = "spline"))
    spc = lapply(spc, function(j){ round(ifelse(j<0, NA, ifelse(j>1, 1, j))*100, 1) })
    spc = as.data.frame(do.call(cbind, spc))
    names(spc) = paste0("scan_visnir.", seq(350, 2500, by=2), "_pcnt")
    ## convert to SNV
    spc = prospectr::standardNormalVariate(spc)
    class(ossl.pca.visnir) = "prcomp"
    X2.pc = as.data.frame(predict(ossl.pca.visnir, newdata=as.data.frame(spc)))[,1:n.spc]
    colnames(X2.pc) = paste0("visnir.PC", 1:n.spc)
  } else {
    X2.pc = NA
  }
  ## obtain GeoTIFF values
  if(geo.type=="ll"){
    pnts = SpatialPoints(data.frame(lon, lat), proj4string = CRS("EPSG:4326"))
    cog.lst = paste0(cog.dir, ossl.model$features[grep("clm_", ossl.model$features)], ".tif")
    ov.tmp = extract.cog(pnts, cog.lst)
  } else {
    ov.tmp = NA
  }
  ## Bind all covariates together
  X = do.call(cbind, list(X1.pc, X2.pc, ov.tmp, data.frame(hzn_depth=hzn_depth)))
  X = X[,which(unlist(lapply(X, function(x) !all(is.na(x)))))]
  X$dataset.code_ascii_c = factor(rep(dataset.code_ascii_c, nrow(X)), levels = c("NEON.SSL", "KSSL.SSL", "CAF.SSL", "AFSIS1.SSL", "AFSIS2.SSL", "ICRAF.ISRIC", "LUCAS.SSL"))
  X <- fastDummies::dummy_cols(X, select_columns = "dataset.code_ascii_c")
  ## predict
  ossl.model$features[which(!ossl.model$features %in% names(X))]
  pred = predict(ossl.model, newdata=X[,ossl.model$features])
  ## uncertainty
  if(sd==TRUE){
    out.c <- as.matrix(as.data.frame(mlr::getStackedBaseLearnerPredictions(ossl.model, newdata=X[,ossl.model$features])))
    cf = eml.cf(ossl.model)
    model.error <- sqrt(matrixStats::rowSds(out.c, na.rm=TRUE)^2 * cf)
  }
  ## Return result as a data.frame:
  out = data.frame(pred.mean=pred$data$response, pred.error=model.error)
  ## back-transform
  if(length(grep("log..", t.var))>0){
    out$tpred.mean = expm1(out$pred.mean)
    if(!is.null(ylim)) {
      out$tpred.mean = ifelse(out$tpred.mean< ylim[1], ylim[1], ifelse(out$tpred.mean > ylim[2], ylim[2], out$tpred.mean))
    }
    out$lower.1std = expm1(out$pred.mean - out$pred.error)
    out$upper.1std = expm1(out$pred.mean + out$pred.error)
  } else {
    out$lower.1std = out$pred.mean - out$pred.error
    out$upper.1std = out$pred.mean + out$pred.error
  }
  if(!is.null(ylim)) {
    out$lower.1std = ifelse(out$lower.1std < ylim[1], ylim[1], ifelse(out$lower.1std > ylim[2], ylim[2], out$lower.1std))
    out$upper.1std = ifelse(out$upper.1std < ylim[1], ylim[1], ifelse(out$upper.1std > ylim[2], ylim[2], out$upper.1std))
  }
  return(list(pred=out, x=X, model=ossl.model$learner.model$super.model$learner.model$model, cf=cf))
}
