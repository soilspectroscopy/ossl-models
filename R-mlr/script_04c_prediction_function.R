
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

mir <- read_csv("sample-data/sample_mir_data.csv", show_col_types = F)
visnir <- read_csv("sample-data/sample_visnir_data.csv", show_col_types = F)
nir <- read_csv("sample-data/sample_neospectra_data.csv", show_col_types = F)

## Predict
predict.ossl <- function(target = "clay.tot_usda.a334_w.pct",
                         # spectra.file = "sample-data/235157XS01.0",
                         # spectra.file = "sample-data/icr056141.0",
                         # spectra.file = "sample-data/101453MD01.asd",
                         # spectra.file = "sample-data/235157MD01.asd",
                         # spectra.file = "sample-data/sample_mir_data.csv",
                         # spectra.file = "sample-data/sample_visnir_data.csv",
                         spectra.file = "sample-data/sample_neospectra_data.csv",
                         # spectra.type = "mir",
                         # spectra.type = "visnir",
                         spectra.type = "nir.neospectra",
                         subset.type = "ossl",
                         geo.type = "na",
                         ncomps = 120,
                         models.dir = "~/mnt-ossl/ossl_models/"){

  ## Reference table
  reference <- read_csv("https://github.com/soilspectroscopy/ossl-models/raw/main/out/front_end_table.csv",
                        show_col_types = F)

  ## Formatting and preprocessing
  mir.spectral.range <- paste0("scan_mir.", seq(600, 4000, by = 2), "_abs")
  visnir.spectral.range <- paste0("scan_visnir.", seq(400, 2500, by = 2), "_ref")
  nir.neospectra.spectral.range <- paste0("scan_nir.", seq(1350, 2550, by = 2), "_ref")

  ## Check if input scans pass some minimum checks, load and preprocess
  if(spectra.type == "mir"){

    if(grepl(".csv", spectra.file)) {

      # Load csv
      data <- readr::read_csv(spectra.file, show_col_types = F)

      # Original range
      old.range <- as.numeric(names(data[,-1]))

      # Check if has a valid range with 5 units of tolerance
      if(min(old.range)-5 > 600 & max(old.range)+5 < 4000) {
        stop(paste0("Input dataset should have a valid spectral range.",
                    "\nCheck https://soilspectroscopy.github.io/ossl-manual/ for examples."))
      }

      # Defining if spectra is in increasing or decreasing order
      if(diff(old.range)[1] < 0) {
        new.range <- rev(seq(600, 4000, by = 2))
      } else {
        new.range <- seq(600, 4000, by = 2)
      }

      # Preparation
      data.prep <- data %>%
        dplyr::select(-1) %>%
        as.matrix() %>%
        # Resampling
        prospectr::resample(X = ., wav = old.range, new.wav = new.range, interpol = "spline") %>%
        # Standard Normal Variate
        prospectr::standardNormalVariate(X = .) %>%
        dplyr::as_tibble() %>%
        # Rebinding sample id
        dplyr::bind_cols({data %>%
            dplyr::select(1)}, .) %>%
        # Renaming to OSSL format
        dplyr::select(1, all_of(as.character(seq(600, 4000, by = 2)))) %>%
        dplyr::rename_with(~mir.spectral.range, as.character(seq(600, 4000, by = 2)))

    } else if(grepl(".0", spectra.file)){

      # Load opus file
      data <- opusreader2::read_opus_single(dsn = spectra.file)

      # Get absorbance from internal storage and format to wide table
      data <- as.data.frame(t(data$ab$data)) %>%
        tibble::rownames_to_column("wavenumber") %>%
        dplyr::mutate(wavenumber = round(as.numeric(wavenumber), 3)) %>%
        dplyr::rename(absorbance = 2) %>%
        dplyr::as_tibble() %>%
        tidyr::pivot_wider(id_cols = everything(), names_from = "wavenumber",
                    values_from = "absorbance") %>%
        dplyr::mutate(sample_id = 1, .before = 1)

      # Original range
      old.range <- as.numeric(names(data[,-1]))

      # Check if has a valid range with 5 units of tolerance
      if(min(old.range)-5 > 600 & max(old.range)+5 < 4000) {
        stop(paste0("Input dataset should have a valid spectral range.",
                    "\nCheck https://soilspectroscopy.github.io/ossl-manual/ for examples."))
      }

      # Defining if spectra is in increasing or decreasing order
      if(diff(old.range)[1] < 0) {
        new.range <- rev(seq(600, 4000, by = 2))
      } else {
        new.range <- seq(600, 4000, by = 2)
      }

      # Preparing
      data.prep <- data %>%
        dplyr::select(-1) %>%
        as.matrix() %>%
        # Resampling
        prospectr::resample(X = ., wav = old.range, new.wav = new.range, interpol = "spline") %>%
        # Standard Normal Variate
        prospectr::standardNormalVariate(X = .) %>%
        dplyr::as_tibble() %>%
        # Rebinding sample id
        dplyr::bind_cols({data %>%
                     select(1)}, .) %>%
        # Renaming to OSSL format
        dplyr::select(1, all_of(as.character(seq(600, 4000, by = 2)))) %>%
        dplyr::rename_with(~mir.spectral.range, as.character(seq(600, 4000, by = 2)))

    } else {

      stop(paste0("Input dataset not correctly imported (must be a csv, opus [.0], or asd file).",
                  "\nSee https://soilspectroscopy.github.io/ossl-manual/ for examples."))

    }

    # Check spectra
    # data.prep %>%
    #   pivot_longer(-1) %>%
    #   ggplot(aes(x = as.numeric(gsub("scan_mir.|_abs", "", name)), y = value)) +
    #   geom_line()

  } else if(spectra.type == "visnir"){

    if(grepl(".csv", spectra.file)) {

      data <- readr::read_csv(spectra.file, show_col_types = F)

      # Original range
      old.range <- as.numeric(names(data[,-1]))

      # Selecting spectra from 400 nm
      old.range <- old.range[which(old.range > 400)]

      # Check if has a valid range with 5 units of tolerance
      if(min(old.range)-5 > 400 & max(old.range)+5 < 2500) {
        stop(paste0("Input dataset should have a valid spectral range.",
                    "\nCheck https://soilspectroscopy.github.io/ossl-manual/ for examples."))
      }

      # Defining if spectra is in increasing or decreasing order
      if(diff(old.range)[1] < 0) {
        new.range <- rev(seq(400, 2500, by = 2))
      } else {
        new.range <- seq(400, 2500, by = 2)
      }

      # Preparing
      data.prep <- data %>%
        dplyr::select(all_of(as.character(old.range))) %>%
        as.matrix() %>%
        # Resampling
        prospectr::resample(X = ., wav = old.range, new.wav = new.range, interpol = "spline") %>%
        # Standard Normal Variate
        prospectr::standardNormalVariate(X = .) %>%
        dplyr::as_tibble() %>%
        # Rebinding sample id
        dplyr::bind_cols({data %>%
            dplyr::select(1)}, .) %>%
        # Renaming to OSSL format
        dplyr::select(1, all_of(as.character(seq(400, 2500, by = 2)))) %>%
        dplyr::rename_with(~visnir.spectral.range, as.character(seq(400, 2500, by = 2)))

    } else if(grepl(".asd", spectra.file)){

      data <- asdreader::get_spectra(spectra.file)

      data <- as.data.frame(t(data)) %>%
        tibble::rownames_to_column("wavelength") %>%
        dplyr::mutate(wavelength = round(as.numeric(wavelength), 3)) %>%
        dplyr::rename(reflectance = 2) %>%
        dplyr::as_tibble() %>%
        tidyr::pivot_wider(id_cols = everything(), names_from = "wavelength",
                    values_from = "reflectance") %>%
        dplyr::mutate(sample_id = 1, .before = 1)

      # Original range
      old.range <- as.numeric(names(data[,-1]))

      # Selecting spectra from 400 nm
      old.range <- old.range[which(old.range > 400)]

      # Check if has a valid range with 5 units of tolerance
      if(min(old.range)-5 > 400 & max(old.range)+5 < 2500) {
        stop(paste0("Input dataset should have a valid spectral range.",
                    "\nCheck https://soilspectroscopy.github.io/ossl-manual/ for examples."))
      }

      # Defining if spectra is in increasing or decreasing order
      if(diff(old.range)[1] < 0) {
        new.range <- rev(seq(400, 2500, by = 2))
      } else {
        new.range <- seq(400, 2500, by = 2)
      }

      # Preparing
      data.prep <- data %>%
        dplyr::select(all_of(as.character(old.range))) %>%
        as.matrix() %>%
        # Resampling
        prospectr::resample(X = ., wav = old.range, new.wav = new.range, interpol = "spline") %>%
        # Standard Normal Variate
        prospectr::standardNormalVariate(X = .) %>%
        dplyr::as_tibble() %>%
        # Rebinding sample id
        dplyr::bind_cols({data %>%
            dplyr::select(1)}, .) %>%
        # Renaming to OSSL format
        dplyr::select(1, all_of(as.character(seq(400, 2500, by = 2)))) %>%
        dplyr::rename_with(~visnir.spectral.range, as.character(seq(400, 2500, by = 2)))

    } else {

      stop(paste0("Input dataset not correctly imported (must be a csv, opus [.0], or asd file).",
                  "\nSee https://soilspectroscopy.github.io/ossl-manual/ for examples."))

    }

    # Check spectra
    # data.prep %>%
    #   pivot_longer(-1) %>%
    #   ggplot(aes(x = as.numeric(gsub("scan_visnir.|_ref", "", name)), y = value)) +
    #   geom_line()

  } else if(spectra.type == "nir.neospectra"){

    if(grepl(".csv", spectra.file)) {

      data <- readr::read_csv(spectra.file, show_col_types = F)

      # Original range
      old.range <- as.numeric(names(data[,-1]))

      # Check if has a valid range with 5 units of tolerance
      if(min(old.range)-5 > 1350 & max(old.range)+5 < 2550) {
        stop(paste0("Input dataset should have a valid spectral range.",
                    "\nCheck https://soilspectroscopy.github.io/ossl-manual/ for examples."))
      }

      # Defining if spectra is in increasing or decreasing order
      if(diff(old.range)[1] < 0) {
        new.range <- rev(seq(1350, 2550, by = 2))
      } else {
        new.range <- seq(1350, 2550, by = 2)
      }

      # Preparing
      data.prep <- data %>%
        dplyr::select(all_of(as.character(old.range))) %>%
        as.matrix() %>%
        # Resampling
        prospectr::resample(X = ., wav = old.range, new.wav = new.range, interpol = "spline") %>%
        # Standard Normal Variate
        prospectr::standardNormalVariate(X = .) %>%
        dplyr::as_tibble() %>%
        # Rebinding sample id
        dplyr::bind_cols({data %>%
            dplyr::select(1)}, .) %>%
        # Renaming to OSSL format
        dplyr::select(1, all_of(as.character(seq(1350, 2550, by = 2)))) %>%
        dplyr::rename_with(~nir.neospectra.spectral.range, as.character(seq(1350, 2550, by = 2)))

    } else {

      stop(paste0("Input dataset not correctly imported (must be a csv, opus [.0], or asd file).",
                  "\nSee https://soilspectroscopy.github.io/ossl-manual/ for examples."))

    }

    # Check spectra
    # data.prep %>%
    #   pivot_longer(-1) %>%
    #   ggplot(aes(x = as.numeric(gsub("scan_nir.|_ref", "", name)), y = value)) +
    #   geom_line()

  }

  ## Loading PCA model and compressing

  ## Loading prediction model
  if(missing(ossl.model)){
    model.url = paste0("http://s3.us-east-1.wasabisys.com/soilspectroscopy/ossl_models/", t.var, "/", spc.type, "_mlr..eml_", subset.type, "_", geo.type, "_v1.rds")
    ossl.model = readRDS(url(model.rds, "rb"))
  }

  ## Loading PCA model and compressing
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
