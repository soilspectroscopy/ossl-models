
## Libraries
library("tidyverse")
library("qs")

## Folder organization
target.folder <- "~/mnt-ossl/ossl_models/"

## Files
ossl.files <- read_csv("out/ossl_models_directory_tree.csv")

ossl.files <- ossl.files %>%
  filter(grepl(".rds", file_path))

## Deleting old qs files
ossl.files <- ossl.files %>%
  mutate(old_file = file_path,
         old_file = gsub(".rds", ".qs", old_file))

error.log <- c()

i=1
for(i in 1:nrow(ossl.files)){

  ibase <- ossl.files[[i,"local_path"]]
  ifile <- ossl.files[[i,"old_file"]]

  full.path <- paste(ibase, ifile, sep = "/")

  if(file.exists(full.path)){
    file.remove(full.path)
  } else {
    error.log[i] <- i
  }

}

## Checking object size

i=1
for(i in 1:nrow(ossl.files)) {

  ibase <- ossl.files[[i,"local_path"]]
  iread <- ossl.files[[i,"file_path"]]

  read.path <- paste(ibase, iread, sep = "/")

  ossl.files[[i,"size"]] <- file.size(read.path)

}

ossl.files <- ossl.files %>%
  filter(size > 0)

## Reexporting

i=1
for(i in 1:nrow(ossl.files)) {

  ibase <- ossl.files[[i,"local_path"]]
  iread <- ossl.files[[i,"file_path"]]
  iexport <- ossl.files[[i,"old_file"]]

  read.path <- paste(ibase, iread, sep = "/")

  file <- readRDS(read.path)

  export.path <- paste(ibase, iexport, sep = "/")

  qsave(file, export.path)

}

## Deleting rds files
i=1
for(i in 1:nrow(ossl.files)) {

  ibase <- ossl.files[[i,"local_path"]]
  iread <- ossl.files[[i,"file_path"]]

  read.path <- paste(ibase, iread, sep = "/")

  file.remove(read.path)

}
