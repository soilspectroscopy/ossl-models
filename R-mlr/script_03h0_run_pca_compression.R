
# Packages
packages <- c("tidyverse", "qs", "parallel", "lubridate",
              "mlr3verse", "mlr3extralearners", "mlr3pipelines",
              "yardstick")

invisible(lapply(packages, library, character.only = TRUE))

# Dirs and reference ranges
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

mir.spectral.range <- paste0("scan_mir.", seq(600, 4000, by = 2), "_abs")
visnir.spectral.range <- paste0("scan_visnir.", seq(400, 2500, by = 2), "_ref")
nir.neospectra.spectral.range <- paste0("scan_nir.", seq(1350, 2550, by = 2), "_ref")

# MIR, VisNIR and NIR from Neospectra
spectra.type <- c("mir", "visnir", "nir.neospectra")

# Using KSSL only or the whole OSSL
subset <- c("kssl", "ossl")

# Adding (ll) or not adding (na) geocovariates to models
geocovariates <- c("na")

# Basic structure
modeling.combinations <- tibble(spectra_type = spectra.type) %>%
  crossing(subset = subset) %>%
  crossing(geo = geocovariates)

# Target models
modeling.combinations <- modeling.combinations %>%
  dplyr::filter(!(grepl("neospectra", spectra_type) & subset == "kssl")) %>%
  dplyr::filter(!(geo == "ll"))

# Model names
modeling.combinations <- modeling.combinations %>%
  mutate(model_name = paste0(spectra_type, "_cubist_",
                             subset, "_", geo, "_v1.2"), .before = 1)

modeling.combinations

# Soil properties
soil.properties <- read_csv("./out/ossl_models_soil_properties.csv")

soil.properties <- soil.properties %>%
  filter(include == TRUE) %>%
  mutate(export_name = ifelse(log == TRUE, paste0("log..", soil_property), soil_property)) %>%
  select(-include, -log)

soil.properties.names <- soil.properties %>%
  pull(soil_property)

## Reading OSSL level 1
rm.ossl <- qread(paste0(db.dir, "ossl_all_L1_v1.2.qs"))

# Preparing the bind of soil data level 1 for Neospectra
neospectra.soillab <- rm.ossl %>%
  dplyr::select(id.layer_uuid_txt, id.scan_local_c,
                all_of(soil.properties.names)) %>%
  filter(grepl("XS|XN", id.scan_local_c)) %>%
  mutate(id.scan_local_c = gsub("XS|XN", "", id.scan_local_c))

# Keeping only the important columns
rm.ossl <- rm.ossl %>%
  dplyr::select(id.layer_uuid_txt, dataset.code_ascii_txt,
                any_of(soil.properties.names),
                all_of(visnir.spectral.range),
                all_of(mir.spectral.range))

## Reading Neospectra datasets
neospectra.nir <- qread(paste0(db.dir, "neospectra_nir_v1.2.qs"))
neospectra.nir

# Spectra collected by multiple instruments
neospectra.nir.avg <- neospectra.nir %>%
  select(id.sample_local_c, all_of(nir.neospectra.spectral.range)) %>%
  group_by(id.sample_local_c) %>%
  summarise_all(mean)

neospectra.nir.reps <- neospectra.nir %>%
  select(id.sample_local_c, all_of(nir.neospectra.spectral.range))

rm.neospectra.avg <- inner_join(neospectra.soillab, neospectra.nir.avg,
                            by = c("id.scan_local_c" = "id.sample_local_c"))

rm.neospectra.reps <- inner_join(neospectra.soillab, neospectra.nir.reps,
                                by = c("id.scan_local_c" = "id.sample_local_c"))

# Selecting only important columns
rm.neospectra.avg <- rm.neospectra.avg %>%
  select(id.layer_uuid_txt,
         any_of(soil.properties.names), all_of(nir.neospectra.spectral.range))

rm.neospectra.reps <- rm.neospectra.reps %>%
  select(id.layer_uuid_txt,
         any_of(soil.properties.names), all_of(nir.neospectra.spectral.range))

# Preparing named list of datasets
# Selecting only important columns, spectra range and rows with available spectra
data.list <- list(
  "mir_cubist_kssl_v1.2" = {rm.ossl %>%
      dplyr::filter(dataset.code_ascii_txt == "KSSL.SSL") %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names),
                    all_of(mir.spectral.range)) %>%
      dplyr::filter(!is.na(scan_mir.1500_abs))},
  "mir_cubist_ossl_v1.2" = {rm.ossl %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names),
                    all_of(mir.spectral.range)) %>%
      dplyr::filter(!is.na(scan_mir.1500_abs))},
  "visnir_cubist_kssl_v1.2" = {rm.ossl %>%
      dplyr::filter(dataset.code_ascii_txt == "KSSL.SSL") %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names),
                    all_of(visnir.spectral.range)) %>%
      dplyr::filter(!is.na(scan_visnir.1500_ref))},
  "visnir_cubist_ossl_v1.2" = {rm.ossl %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names),
                    all_of(visnir.spectral.range)) %>%
      dplyr::filter(!is.na(scan_visnir.1500_ref))},
  "nir.neospectra_cubist_ossl_v1.2" = {rm.neospectra.avg %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names),
                    all_of(nir.neospectra.spectral.range)) %>%
      dplyr::filter(!is.na(scan_nir.1500_ref))},
  "nir.neospectra_cubist_ossl_v1.2_reps" = {rm.neospectra.reps %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names),
                    all_of(nir.neospectra.spectral.range)) %>%
      dplyr::filter(!is.na(scan_nir.1500_ref))}
)

# # Different sizes
# lapply(data.list, function(x) dim(x))

prep.list <- lapply(data.list, function(x){
  x %>%
    dplyr::select(-id.layer_uuid_txt, -any_of(soil.properties.names)) %>%
    as.matrix() %>%
    prospectr::standardNormalVariate(X = .) %>%
    as_tibble() %>%
    bind_cols({x %>%
        dplyr::select(id.layer_uuid_txt, any_of(soil.properties.names))}, .)
})

# # Checking preprocessing. Names are kept consistently across list objects and table columns
# lapply(data.list, function(x) dim(x))

pca.list <- mclapply(1:length(prep.list), function(i) {

  x <- prep.list[[i]] %>%
    dplyr::select(starts_with("scan_")) %>%
    as.data.frame()

  prcomp(x, center = T, scale = T)

}, mc.cores = length(prep.list))

names(pca.list) <- names(prep.list)

# # Checking the number of components of each spectra type
# lapply(pca.list, function(x) ncol(x$rotation))
#
# # Checking how many components explain 95% of the original variance
# lapply(pca.list, function(x) {which(cumsum(x$sdev/sum(x$sdev)) > 0.95)[1]})
#
# # Checking how many components explain 99% of the original variance
# lapply(pca.list, function(x) {which(cumsum(x$sdev/sum(x$sdev)) > 0.99)[1]})
#
# # Checking how many components explain 99.9%of the original variance
# lapply(pca.list, function(x) {which(cumsum(x$sdev/sum(x$sdev)) > 0.999)[1]})
#
# # Checking how many components explain 99.99% of the original variance
# lapply(pca.list, function(x) {which(cumsum(x$sdev/sum(x$sdev)) > 0.9999)[1]})

## Saving simplified pca models only for prediction purposes (m = modified):
# Omitting 5th object of prcomp, i.e x table (scores)
for(i in 1:length(pca.list)){
  export.name <- paste0(dir, "pca.ossl/mpca_", names(pca.list)[i], ".qs")
  if(file.exists(export.name)){file.remove(export.name)}
  qsave(pca.list[[i]][1:4], export.name)
}

## Saving the scores to be used as predictors
for(i in 1:length(pca.list)){

  idata <- prep.list[[i]] %>%
    dplyr::select(id.layer_uuid_txt, all_of(soil.properties.names))

  iscores <- pca.list[[i]]$x %>%
    as_tibble()

  idata.export <- bind_cols(idata, iscores)

  export.name <- paste0(dir, "pca.ossl/pca_scores_", names(pca.list)[i], ".qs")
  if(file.exists(export.name)){file.remove(export.name)}

  qsave(idata.export, export.name)

}
