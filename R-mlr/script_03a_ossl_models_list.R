
# Packages
library("tidyverse")
library("qs")

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

# Model name
# modeling.combinations <- modeling.combinations %>%
#   mutate(model_name = paste0(spectra_type, "_mlr3..eml_", subset, "_", geo, "_v1.2"), .before = 1)

modeling.combinations <- modeling.combinations %>%
  mutate(model_name = paste0(spectra_type, "_cubist_", subset, "_", geo, "_v1.2"), .before = 1)

soil.properties <- read_csv("./out/ossl_models_soil_properties.csv")

soil.properties <- soil.properties %>%
  filter(include == TRUE) %>%
  mutate(export_name = ifelse(log == TRUE, paste0("log..", soil_property), soil_property)) %>%
  select(-include, -log)

soil.properties.names <- soil.properties %>%
  pull(soil_property)

# Final modeling combinations
modeling.combinations <- soil.properties %>%
  crossing(modeling.combinations)

write_csv(modeling.combinations, "./out/modeling_combinations_v1.2.csv")

# Modeling combinations
modeling.combinations <- read_csv("./out/modeling_combinations_v1.2.csv")
count.table <- read_csv("./out/tab_dataset_count.csv")

# Defining available data from ossl-import
count.table <- count.table %>%
  filter(dataset %in% c("KSSL.SSL", "OSSL")) %>%
  rename(subset = dataset) %>%
  mutate(subset = recode(subset,
                         "KSSL.SSL" = "kssl",
                         "OSSL" = "ossl")) %>%
  pivot_longer(-all_of(c("soil_property", "subset")),
               names_to = "spectra_type", values_to = "count") %>%
  mutate(spectra_type = recode(spectra_type,
                               "n_mir" = "mir",
                               "n_visnir" = "visnir",
                               "n_neospectra" = "nir.neospectra"))

# Defining models with at least 500 observations
modeling.combinations <- left_join(modeling.combinations,
                                   count.table,
                                   by = c("soil_property", "spectra_type", "subset"))

modeling.combinations <- modeling.combinations %>%
  filter(count > 500) %>%
  filter(!(soil_property == "efferv_usda.a479_class"))

write_csv(modeling.combinations, "out/fitted_modeling_combinations_v1.2.csv")

# Available soil properties
modeling.combinations %>%
  distinct(soil_property) %>%
  count()

# Final modeling combinations
modeling.combinations %>%
  count(spectra_type, subset)

modeling.combinations.overview <- modeling.combinations %>%
  select(export_name, description, spectra_type, subset, count) %>%
  pivot_wider(names_from = "spectra_type", values_from = "count")

# clipr::write_clip(modeling.combinations.overview)
