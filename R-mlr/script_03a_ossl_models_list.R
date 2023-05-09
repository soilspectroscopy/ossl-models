
# Packages
library("tidyverse")
library("qs")

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
