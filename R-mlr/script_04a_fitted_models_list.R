
## Packages
library("tidyverse")
library("yardstick")
library("qs")
library("fs")

## Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

## Fitted models
fitted.modeling.combinations <- read_csv("./out/fitted_modeling_combinations_v1.2.csv",
                                         show_col_types = FALSE)

fitted.modeling.combinations <- fitted.modeling.combinations %>%
  select(soil_property, description, model_name, export_name) %>%
  mutate(unit_transform = export_name) |>
  mutate(unit_transform = ifelse(grepl("log..", unit_transform), "log1p", "original"))

fitted.modeling.combinations

files.table <- tibble(file_description = c("train data", "model", "performance",
                                           "10cv predictions", "validation plot",
                                           "error model", "error predictions"),
                      file_code = c("task_", "model_", "perfmetrics_",
                                    "cvpred_", "valplot_",
                                    "error_model_", "error_pred_"),
                      file_extension = c(".qs", ".qs", ".csv",
                                         ".qs", ".png",
                                         ".qs", ".qs"))

files.table <- fitted.modeling.combinations %>%
  crossing(files.table)

files.table <- files.table %>%
  mutate(file_url = paste0("https://storage.googleapis.com/soilspec4gg-public/models/",
                             export_name,
                             "/",
                             file_code,
                             model_name,
                             file_extension)) %>%
  select(-file_code, -file_extension, -export_name)

write_csv(files.table, "out/fitted_models_access.csv")

## Test

test.url <- files.table %>%
  filter(grepl("clay", soil_property)) %>%
  filter(grepl("nir.neospectra", model_name)) %>%
  filter(grepl("train data", file_description)) %>%
  pull(file_url)

train.clay.nir.neospectra <- qread_url(test.url)
