
## Packages
library("tidyverse")
library("qs")

## Folders
mnt.dir <- "~/projects/temp/ossl_models/comparison/"

## Importing data
important.columns <- c("id.layer_uuid_txt", "dataset.code_ascii_txt", "id.tile")
soil.properties <- c("c.tot_usda.a622_w.pct", "s.tot_usda.a624_w.pct")

mir.spectral.range <- paste0("scan_mir.", seq(600, 4000, by = 2), "_abs")
visnir.spectral.range <- paste0("scan_visnir.", seq(400, 2500, by = 2), "_ref")

ossl <- qread_url("https://storage.googleapis.com/soilspec4gg-public/ossl_all_L1_v1.2.qs")

oss.ovl <- qread_url("https://storage.googleapis.com/soilspec4gg-public/ossl_overlay_v1.2.qs")

# starts_with("clm") = geocovariates
# id.tile = block CV

oss.ovl <- oss.ovl %>%
  select(any_of(important.columns), starts_with("clm"))

## This will find the samples with both VisNIR and MIR data,
## With spatial locations and associated climate information,
## and complete cases for both total C and total S

ossl.subset <- ossl %>%
  select(any_of(important.columns), any_of(soil.properties),
         any_of(visnir.spectral.range), any_of(mir.spectral.range)) %>%
  filter(!(is.na(scan_mir.1500_abs))) %>%
  filter(!(is.na(scan_visnir.1500_ref))) %>%
  inner_join(oss.ovl, by = c("dataset.code_ascii_txt", "id.layer_uuid_txt")) %>%
  filter(!is.na(c.tot_usda.a622_w.pct)) %>%
  filter(!is.na(s.tot_usda.a624_w.pct)) %>%
  relocate(id.tile, .after = dataset.code_ascii_txt)

## Basically only samples from KSSL
ossl.subset %>%
  count(dataset.code_ascii_txt)

selected.ids <- ossl.subset %>%
  pull(id.layer_uuid_txt)

qsave(ossl.subset, paste0(mnt.dir, "ossl_comp_subset_v1.2.qs"))

## Reading PCA scores from KSSL
pca.kssl.mir <- qread(file.path(dirname(mnt.dir), "pca_scores_mir_mlr..eml_kssl_v1.2.qs"))

pca.kssl.mir.subset <- pca.kssl.mir %>%
  filter(id.layer_uuid_txt %in% selected.ids) %>%
  select(id.layer_uuid_txt, starts_with("PC"))

pca.kssl.mir.subset.export <- ossl.subset %>%
  select(all_of(important.columns), all_of(soil.properties), starts_with("clm")) %>%
  left_join(pca.kssl.mir.subset, by = "id.layer_uuid_txt")

qsave(pca.kssl.mir.subset.export, paste0(mnt.dir, "pca_scores_mir_kssl_subset_v1.2.qs"))

pca.kssl.visnir <- qread(file.path(dirname(mnt.dir), "pca_scores_visnir_mlr..eml_kssl_v1.2.qs"))

pca.kssl.visnir.subset <- pca.kssl.visnir %>%
  filter(id.layer_uuid_txt %in% selected.ids) %>%
  select(id.layer_uuid_txt, starts_with("PC"))

pca.kssl.visnir.subset.export <- ossl.subset %>%
  select(all_of(important.columns), all_of(soil.properties), starts_with("clm")) %>%
  left_join(pca.kssl.visnir.subset, by = "id.layer_uuid_txt")

qsave(pca.kssl.visnir.subset.export, paste0(mnt.dir, "pca_scores_visnir_kssl_subset_v1.2.qs"))

## Modeling combinations
modeling.combinations <- tibble(model_type = "ensemble") %>%
  crossing(base_predictions = c("cv5", "insample")) %>%
  bind_rows(tibble(model_type = "cubist", base_predictions = "insample")) %>%
  crossing(soil_property = soil.properties) %>%
  crossing(spectra = c("mir", "visnir")) %>%
  crossing(geocovariates = c("na", "ll")) %>%
  crossing(pca_compression = c("n120", "cumvar99perc"))

modeling.combinations

write_csv(modeling.combinations, paste0(mnt.dir, "modeling_combinations.csv"))
