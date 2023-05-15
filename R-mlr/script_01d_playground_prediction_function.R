
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
library("matrixStats")

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

  ## Loading PCA model for compression
  spectra.type = "nir.neospectra"
  subset.type = "ossl"
  geo.type = "na"
  ncomps = 120
  models.dir = "~/mnt-ossl/ossl_models/"

  pca.model <- qs::qread(paste0(models.dir,
                            "pca.ossl/mpca_",
                            spectra.type,
                            "_mlr3..eml_",
                            subset.type,
                            "_v1.2.qs"))

  class(pca.model) <- "prcomp"

  data.scores <- dplyr::bind_cols(data.prep[,1], predict(pca.model, data.prep)) %>%
    dplyr::select(1, all_of(paste0("PC", seq(1, ncomps, 1))))

  ## Loading prediction model
  prediction.model <- qs::qread(paste0(models.dir,
                                   target,
                                   "/model_",
                                   spectra.type,
                                   "_mlr3..eml_",
                                   subset.type,
                                   "_",
                                   geo.type,
                                   "_v1.2.qs"))

  ## Preparing data
  task <- as.data.table(data.scores)

  ## Prediction
  data.prediction <- as.data.table(prediction.model$predict_newdata(task))

  data.prediction <- dplyr::bind_cols(task[,1], {
    data.prediction %>%
      dplyr::select(response) %>%
      rename(!!target := response)}) %>%
    as_tibble()

  ## Uncertainty (CI 95%) based on the standard error of prediction
  ## from the meta-learner linear model (intercept and 4 base learners, df = 4)
  prediction.model$predict_type = "se"

  crit <- qt(p = 0.05/2, df = 4, lower.tail = FALSE)

  confidence.predictions <- as.data.table(prediction.model$predict_newdata(task)) %>%
    dplyr::select(se) %>%
    dplyr::rename(std_error = se) %>%
    dplyr::mutate(critical_value = crit)

  out <- dplyr::bind_cols(data.prediction, confidence.predictions) %>%
    dplyr::mutate(lower_CI95 = !!as.name(target)-(std_error*critical_value),
                  upperCI_CI95 = !!as.name(target)+(std_error*critical_value)) %>%
    dplyr::select(-critical_value)

  ## Back-transforming
  if(grepl("log..", target)){
    out <- out %>%
      dplyr::mutate(!!target := expm1(!!as.name(target)),
                    std_error = expm1(std_error),
                    lower_CI95 = expm1(lower_CI95),
                    upperCI_CI95 = expm1(upperCI_CI95))
    }

  ## Spectral outlier screening


  ## Final results
  return(out)

}
