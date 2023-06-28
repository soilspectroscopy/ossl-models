
## Packages
library("tidyverse")
library("qs")
library("fs")

## Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

## Fitted models
fitted.modeling.combinations <- read_csv("./out/fitted_modeling_combinations_v1.2.csv",
                                         show_col_types = FALSE)

fitted.modeling.combinations <- fitted.modeling.combinations %>%
  mutate(unit_transform = export_name) |>
  mutate(unit_transform = ifelse(grepl("log..", unit_transform), "log1p", "original")) %>%
  select(-count)

fitted.modeling.combinations

frontend.table <- tibble(file_description = c("train_data", "model", "performance_summary",
                                              "10cv_predictions", "validation_plot",
                                              "error_model", "error_predictions"),
                         file_code = c("task_", "model_", "perfmetrics_",
                                       "cvpred_", "valplot_",
                                       "error_model_", "error_pred_"),
                         file_extension = c(".qs", ".qs", ".csv",
                                            ".qs", ".png",
                                            ".qs", ".qs"))

frontend.table <- fitted.modeling.combinations %>%
  crossing(frontend.table)

base.url <- "https://storage.googleapis.com/soilspec4gg-public/models/"
manual.url <- "https://soilspectroscopy.github.io/ossl-manual/ossl-database-description.html#"

frontend.table <- frontend.table %>%
  mutate(file_url = paste0(base.url,
                           export_name,
                           "/",
                           file_code,
                           model_name,
                           file_extension)) %>%
  select(-file_code, -file_extension, -export_name)

frontend.table <- frontend.table %>%
  pivot_wider(names_from = "file_description", values_from = "file_url")

frontend.table <- frontend.table %>%
  mutate(description_url = paste0(manual.url, soil_property),
         .after = description) %>%
  mutate(pca_model_url = paste0(base.url, "pca.ossl/mpca_", gsub("_ll|_na", "", model_name), ".qs"),
         pca_scores_url = paste0(base.url, "pca.ossl/pca_scores_", gsub("_ll|_na", "", model_name), ".qs"),)

write_csv(frontend.table, "out/front_end_table.csv")

## Test

pca.scores <- qread_url("https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/pca_scores_nir.neospectra_cubist_ossl_v1.2.qs")
