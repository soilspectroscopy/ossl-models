
## Packages
library("tidyverse")
library("mlr3")
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
                         models.dir = "~/mnt-ossl/ossl_models/"){

  ## Formatting and preprocessing
  mir.spectral.range <- paste0("scan_mir.", seq(600, 4000, by = 2), "_abs")
  visnir.spectral.range <- paste0("scan_visnir.", seq(400, 2500, by = 2), "_ref")
  nir.neospectra.spectral.range <- paste0("scan_nir.", seq(1350, 2550, by = 2), "_ref")
  ncomps = 120

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
  ## Values from "out/trustworthiness_q_critical_values.csv"
  ## Check for updates
  q.critical <- tibble(model_name = c("mir_mlr3..eml_kssl_v1.2", "mir_mlr3..eml_ossl_v1.2",
                                      "nir.neospectra_mlr3..eml_ossl_v1.2",
                                      "visnir_mlr3..eml_kssl_v1.2", "visnir_mlr3..eml_ossl_v1.2"),
                       q_critical = c(0.0245, 0.0452, 0.00000274, 0.00129, 0.000812))

  # Full and employed PC models (120 comps)
  pca.model.full <- pca.model
  class(pca.model.full) <- "prcomp"

  pca.model.employed <- pca.model
  class(pca.model.employed) <- "prcomp"
  pca.model.employed$sdev <- pca.model.employed$sdev[1:ncomps]
  pca.model.employed$rotation <- pca.model.employed$rotation[,1:ncomps]

  # Scores
  pca.scores.full <- dplyr::bind_cols(data.prep[,1], predict(pca.model.full, data.prep))
  pca.scores.employed <- dplyr::bind_cols(data.prep[,1], predict(pca.model.employed, data.prep))

  # Back-transformed spectra
  spectra.bt.full <- as.matrix(pca.scores.full[,-1]) %*% t(pca.model.full$rotation)
  spectra.bt.full <- t((t(spectra.bt.full) * pca.model.full$scale) + pca.model.full$center)
  spectra.bt.full <- dplyr::bind_cols(pca.scores.full[,1], spectra.bt.full)

  spectra.bt.employed <- as.matrix(pca.scores.employed[,-1]) %*% t(pca.model.employed$rotation)
  spectra.bt.employed <- t((t(spectra.bt.employed) * pca.model.employed$scale) + pca.model.employed$center)
  spectra.bt.employed <- dplyr::bind_cols(pca.scores.employed[,1], spectra.bt.employed)

  # Q stats - sum of squared differences between full and employed PCs back-transformed spectra
  q.stats <- apply((spectra.bt.employed[,-1]-spectra.bt.full[,-1])^2, MARGIN = 1, sum)

  q.critical <- q.critical %>%
    dplyr::filter(grepl(spectra.type, model_name)) %>%
    dplyr::filter(grepl(subset.type, model_name)) %>%
    dplyr::pull(q_critical)

  # Flagging spectral outlier
  q.stats <- tibble::tibble(q_stats = q.stats,
                            q_critical = q.critical) %>%
    dplyr::mutate(spectral_outlier = ifelse(q_stats >= q_critical, TRUE, FALSE)) %>%
    dplyr::select(spectral_outlier)

  ## Final results
  out <- bind_cols(out, q.stats)

  return(out)

}


## Testing

# spectra.file = "sample-data/235157XS01.0",
# spectra.file = "sample-data/icr056141.0",
# spectra.file = "sample-data/101453MD01.asd",
# spectra.file = "sample-data/235157MD01.asd",
# spectra.file = "sample-data/sample_mir_data.csv",
# spectra.file = "sample-data/sample_visnir_data.csv",
# spectra.type = "mir"
# spectra.type = "visnir",

## Neospectra
test1 <- predict.ossl(target = "clay.tot_usda.a334_w.pct",
                      spectra.file = "sample-data/sample_neospectra_data.csv",
                      spectra.type = "nir.neospectra",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test2 <- predict.ossl(target = "log..oc_usda.c729_w.pct",
                      spectra.file = "sample-data/sample_neospectra_data.csv",
                      spectra.type = "nir.neospectra",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test3 <- predict.ossl(target = "log..p.ext_usda.a1070_mg.kg",
                      spectra.file = "sample-data/sample_neospectra_data.csv",
                      spectra.type = "nir.neospectra",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")


## VisNIR
test4 <- predict.ossl(target = "clay.tot_usda.a334_w.pct",
                      spectra.file = "sample-data/235157MD01.asd",
                      spectra.type = "visnir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test5 <- predict.ossl(target = "clay.tot_usda.a334_w.pct",
                      spectra.file = "sample-data/sample_visnir_data.csv",
                      spectra.type = "visnir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test6 <- predict.ossl(target = "log..oc_usda.c729_w.pct",
                      spectra.file = "sample-data/235157MD01.asd",
                      spectra.type = "visnir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test7 <- predict.ossl(target = "log..oc_usda.c729_w.pct",
                      spectra.file = "sample-data/sample_visnir_data.csv",
                      spectra.type = "visnir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

## MIR
test8 <- predict.ossl(target = "clay.tot_usda.a334_w.pct",
                      spectra.file = "sample-data/icr056141.0",
                      spectra.type = "mir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test9 <- predict.ossl(target = "clay.tot_usda.a334_w.pct",
                      spectra.file = "sample-data/sample_mir_data.csv",
                      spectra.type = "mir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test10 <- predict.ossl(target = "log..oc_usda.c729_w.pct",
                      spectra.file = "sample-data/icr056141.0",
                      spectra.type = "mir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")

test11 <- predict.ossl(target = "log..oc_usda.c729_w.pct",
                      spectra.file = "sample-data/sample_mir_data.csv",
                      spectra.type = "mir",
                      subset.type = "ossl",
                      geo.type = "na",
                      models.dir = "~/mnt-ossl/ossl_models/")
