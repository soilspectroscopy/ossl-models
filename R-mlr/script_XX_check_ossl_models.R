
library("tidyverse")
library("stringr")

dir.ossl <- "/mnt/soilspec4gg/ossl/ossl_models"

## PCA models
pca.models <- paste0(dir.ossl, "/pca.ossl")
pca.models <- list.files(pca.models)
pca.models <- grep("*v1.2.rds$", pca.models, value = T)
pca.models <- grep("^pca", pca.models, value = T)
pca.models

## Soil properties
soil.properties <- c("silt.tot_usda.c62_w.pct",
                     "clay.tot_usda.a334_w.pct",
                     "sand.tot_usda.c60_w.pct",
                     "log..oc_usda.c729_w.pct",
                     "log..n.tot_usda.a623_w.pct",
                     "log..caco3_usda.a54_w.pct",
                     "log..al.ox_usda.a59_w.pct",
                     "ph.h2o_usda.a268_index",
                     "ph.cacl2_usda.a481_index",
                     "log..k.ext_usda.a725_cmolc.kg",
                     "log..mg.ext_usda.a724_cmolc.kg",
                     "log..ca.ext_usda.a722_cmolc.kg",
                     "acidity_usda.a795_cmolc.kg",
                     "bd_usda.a4_g.cm3")

soillab <- read_csv("out/ossl_level1_names_soillab.csv")

soillab <- soillab %>%
  mutate(variable_name = paste0(analyte, " (", ossl_unit_level0, ")")) %>%
  rename(variable_join = ossl_name_level0) %>%
  select(variable_join, variable_name)

## Fitted models
models.dir <- paste(dir.ossl, soil.properties, sep = "/")

base.url <- "http://s3.us-east-1.wasabisys.com/soilspectroscopy/ossl_models/"

listed.models <- sapply(models.dir, function(x) {list.files(x, pattern = "*v1.2.rds$")}) %>%
  enframe() %>%
  unnest(value) %>%
  rename(variable = name, model_name = value) %>%
  mutate(model_path = variable, .before = 1) %>%
  mutate(variable = str_replace(variable, paste0(dir.ossl, "/"), "")) %>%
  mutate(variable_join = gsub("log..", "", variable)) %>%
  left_join(soillab, by = "variable_join") %>%
  select(-variable_join) %>%
  relocate(variable_name, .after = variable) %>%
  separate(model_name, into = c("spectra_type", "other"), sep = "_mlr..eml_", remove = F) %>%
  mutate(other = str_replace(other, ".rds", "")) %>%
  separate(other, into = c("subset", "geo", "version"), sep = "_") %>%
  mutate(variable_description = paste0("https://soilspectroscopy.github.io/ossl-manual/ossl-database-description.html#", str_replace(variable, "log..", "")),
         .before = model_name) %>%
  mutate(model_summary = paste0(base.url, variable, "/", model_name, "_resultsFit.txt"),
         model_scatterplot = paste0(base.url, variable, "/ap.", model_name, ".png"),
         model_url = paste0(base.url, variable, "/", model_name),
         pca_model_mir = case_when(spectra_type == "mir" & subset == "kssl" ~ paste0(base.url, "pca.ossl/pca_mir_kssl_v1.2.rds"),
                                   spectra_type == "mir" & subset == "ossl" ~ paste0(base.url, "pca.ossl/pca_mir_ossl_v1.2.rds"),
                                   spectra_type == "visnir.mir" & subset == "ossl" ~ paste0(base.url, "pca.ossl/pca_mir_ossl_v1.2.rds"),
                                   TRUE ~ ""),
         pca_model_visnir = case_when(spectra_type == "nir.neospectra" & subset == "ossl" ~ paste0(base.url, "pca.ossl/pca_nir_neospectra_v1.2.rds"),
                                      spectra_type == "visnir.mir" & subset == "kssl" ~ paste0(base.url, "pca.ossl/pca_visnir_kssl_v1.2.rds"),
                                      spectra_type == "visnir.mir" & subset == "ossl" ~ paste0(base.url, "pca.ossl/pca_visnir_ossl_v1.2.rds"),
                                      TRUE ~ "")) %>%
  mutate(model_description = paste0("Ensemble ML model to predict ", variable_name),
         .after = model_name) %>%
  relocate(model_path, .after = version)

listed.models %>%
  glimpse()

write_csv(listed.models, "out/list_ossl_models_v1.2.csv")

## Binding with soillab info

# listed.models <- read_csv("out/list_ossl_models_v1.2.csv")
# clipr::write_clip(listed.models)
