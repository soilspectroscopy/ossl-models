
## Packages
library("tidyverse")
library("qs")

## Dirs
dir <- "~/mnt-ossl/ossl_models/"

## Reference files
fitted.modeling.combinations <- read_csv("./out/fitted_modeling_combinations_v1.2.csv",
                                         show_col_types = FALSE)

reference.table <- tibble(file_description = c("train_data", "model", "performance_summary",
                                              "10cv_predictions", "validation_plot",
                                              "error_model", "error_predictions"),
                         file_code = c("task_", "model_", "perfmetrics_",
                                       "cvpred_", "valplot_",
                                       "error_model_", "error_pred_"),
                         file_extension = c(".qs", ".qs", ".csv",
                                            ".qs", ".png",
                                            ".qs", ".qs"))

reference.table <- fitted.modeling.combinations %>%
  crossing(reference.table)

reference.table <- reference.table %>%
  mutate(path = paste0(dir, export_name, "/", file_code, model_name, file_extension))

reference.table <- reference.table %>%
  mutate(exists = file.exists(path))

reference.table <- reference.table %>%
  mutate(size = file.size(path))

## Current files
ossl.models.dir <- list.files(dir, recursive = T, full.names = F)

ossl.models.dir <- tibble(file_path = paste0(dir, ossl.models.dir))

ossl.models.dir <- ossl.models.dir %>%
  filter(grepl("cubist", file_path))

reference.table.check <- reference.table %>%
  select(path)

check.table <- anti_join(ossl.models.dir, reference.table.check, by = c("file_path" = "path"))
check.table

# remove.table <- check.table %>%
#   filter(!grepl("pca", file_path))
#
# apply(remove.table, MARGIN=1, FUN = function(x) file.remove(x))

## Missing files

missing.table <- reference.table %>%
  filter(!exists)

missing.table.train <- missing.table %>%
  filter(grepl("model_|task_", file_code))

write_csv(missing.table.train, "R-mlr/sandbox/corrupted/redo_train.csv")

missing.table.evaluation <- missing.table %>%
  filter(grepl("cvpred_", file_code))

write_csv(missing.table.evaluation, "R-mlr/sandbox/corrupted/redo_evaluation.csv")

missing.table.error <- missing.table %>%
  filter(grepl("error_", file_code))

write_csv(missing.table.error, "R-mlr/sandbox/corrupted/redo_error.csv")
