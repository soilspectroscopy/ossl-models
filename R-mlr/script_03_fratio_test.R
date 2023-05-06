
## Packages
library("tidyverse")
library("readr")
library("qs")
library("microbenchmark")

source("R-mlr/SSL_functions.R")

neospectra.models.dir <- "~/projects/mnt-neospectra/ossl_models/"
neospectra.data.dir <- "~/projects/mnt-neospectra/preprocessed/"
ossl.models.dir <- "~/projects/temp/ossl_models/"

# ## Test sample
# neospectra.raw <- qread(paste0(neospectra.data.dir, "neospectra_nir_avgSample_raw.qs"))
#
# neospectra.raw.sample <- neospectra.raw %>%
#   sample_n(10)
#
# neospectra.raw.sample <- neospectra.raw.sample %>%
#   rename(sample_id = kssl_id) %>%
#   mutate(sample_id = row_number())
#
# neospectra.raw.sample
#
# write_csv(neospectra.raw.sample, "sample-data/sample_neospectra_data.csv")

mir.test <- read_csv("sample-data/sample_mir_data.csv")
# site.test <- read_csv("sample-data/sample_soilsite_data.csv")
# visnir.test <- read_csv("sample-data/sample_visnir_data.csv")
# neospectra.test <- read_csv("sample-data/sample_neospectra_data.csv")

mir.test

## OSSL models
model.pca <- readRDS.gz(paste0(ossl.models.dir, "mpca_mir_kssl_v1.2.rds"))
model.target.na <- readRDS.gz(paste0(ossl.models.dir, "mir_mlr..eml_kssl_na_v1.2.rds"))
model.target.ll <- readRDS.gz(paste0(ossl.models.dir, "mir_mlr..eml_kssl_ll_v1.2.rds"))

# Basic formatting
mir.raw <- as.data.frame(mir.test)
ids <- mir.raw[,1]
wn = as.numeric(gsub("X", "", names(mir.raw[,-1])))
spc = as.matrix(mir.raw[,-1])

# Wavenumber headers can be provided in increasing or decreasing order
# For interpolation, the correct sequence must be provided
if(diff(wn)[1] < 0) {
   new.wn <- seq(4000, 600, by=-2)
} else if(diff(wn)[1] > 0) {
  new.wn <- seq(600, 4000, by=2)
}

# Resample
spc = as.data.frame(prospectr::resample(spc, wn, new.wn, interpol = "spline"))
names(spc) = paste0("scan_mir.", new.wn, "_abs")
spc[1:5, 1:5]

# Preprocess with SNV
spc = prospectr::standardNormalVariate(spc)

# Compress by PCA
ossl.pca.mir <- model.pca
n.spc <- 120
class(ossl.pca.mir) = "prcomp"
X1.pc = as.data.frame(predict(ossl.pca.mir, newdata=as.data.frame(spc)))[,1:n.spc]
colnames(X1.pc) = paste0("mir.PC", 1:n.spc)
X1.pc[1:5,1:5]

# Other info
X2.pc <- NA
ov.tmp <- NA
hzn_depth <- NA

# Bind all predictors together
X = do.call(cbind, list(X1.pc, X2.pc, ov.tmp, hzn_depth))
head(names(X))
tail(names(X))
X[1:5,(ncol(X)-5):ncol(X)]
X = X[,which(unlist(lapply(X, function(x) !all(is.na(x)))))]

# Predict
ossl.model <- model.target.na
ossl.model$features[which(!ossl.model$features %in% names(X))]
pred = predict(ossl.model, newdata=X[,ossl.model$features])

# Uncertainty
out.c <- as.matrix(as.data.frame(mlr::getStackedBaseLearnerPredictions(ossl.model, newdata=X[,ossl.model$features])))
cf = eml.cf(ossl.model)
model.error <- sqrt(matrixStats::rowSds(out.c, na.rm=TRUE)^2 * cf)

# Result as a data.frame:
out = data.frame(id = ids, pred.mean=pred$data$response, pred.error=model.error)

## Back-transform
t.var <- "clay.tot_usda.a334_w.pct"
ylim=NULL
out <- out2
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

list(pred=out, x=X, cf=cf)
