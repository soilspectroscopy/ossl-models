
## Packages
library("tidyverse")
library("qs")

## Dirs
dir <- "~/mnt-ossl/ossl_models/"

## OSSL models dir
ossl.models.dir <- list.files(dir, recursive = T)

ossl.models.dir <- tibble(file_path = ossl.models.dir)

ossl.models.dir <- ossl.models.dir %>%
  mutate(file_path = str_replace(file_path, "/home/jsafanelli/mnt-ossl/ossl_models/", ""))

ossl.models.dir <- ossl.models.dir %>%
  mutate(local_path = "/mnt/soilspec4gg/ossl/ossl_models/mlr3",
         public_path = "https://storage.googleapis.com/soilspec4gg-public/models",
         .before = 1)

write_csv(ossl.models.dir, "out/ossl_models_directory_tree.csv")
