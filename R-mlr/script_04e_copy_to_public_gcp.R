
library("fs")

##  List of files of interest
ossl.models.files <- read_csv("out/ossl_models_directory_tree.csv")
ossl.models.files

## Temporary folder
rsync.dir <- "~/mnt-ossl/rsync/"
dir.exists(rsync.dir)
if(!dir.exists(rsync.dir)){dir.create(rsync.dir)}

## Copying to temporary folder to run rsync
ossl.models.files <- ossl.models.files %>%
  mutate(current = paste0(local_path, file_path)) %>%
  mutate(rsync = paste0(rsync.dir, file_path))

## Creating subfolders
folders <- unique(dirname(ossl.models.files$rsync))
sapply(folders, function(x) if(!dir.exists(x)){dir.create(x)})

## Copying files
file.copy(ossl.models.files$current, ossl.models.files$rsync, overwrite = T)

## Listing current files
system("gsutil ls gs://soilspec4gg-public/models")

## Deleting current files
# system("gsutil rm -r gs://soilspec4gg-public/models")

## Create a new folder using the GUI

## Check rsync without performing it
system("gsutil -m rsync -r -n ~/mnt-ossl/rsync gs://soilspec4gg-public/models")

## Execute rsync
system("gsutil -m rsync -r ~/mnt-ossl/rsync gs://soilspec4gg-public/models")

## Delete temporary folder
dir_delete(rsync.dir)
