
## Dirs
dir <- "~/mnt-ossl/ossl_models/"

## Copying to repository
files <- list.files(dir, pattern = "pca", recursive = T, full.names = T)
files

files <- tibble(path = files) %>%
  mutate(export_path = str_replace(path, "mlr3..eml", "cubist"))

file.copy(files$path, files$export_path, overwrite = T)
