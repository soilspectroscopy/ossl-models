
# Libraries
library("tidyverse")

# Target folder for storing the local files
# The last folder must be named 'ossl_models'
# Parent folder can vary
target.folder <- "~/projects/temp/ossl_models"

# Repository tree and structure used for reference
ossl.files <- read_csv("out/ossl_models_directory_tree.csv")

# Specify your spectral range of interest
# Or omit this step if interested in all spectral ranges
visnir.range <- "visnir"
mir.range <- "mir"
nir.range <- "neospectra"

ossl.files <- ossl.files %>%
  filter(grepl(nir.range, file_path))

# Creating subfolders recursevely for storing the models and ancillary files
subfolders <- unique(dirname(ossl.files$file_path))

for(i in 1:length(subfolders)){
  isubfolder <- subfolders[i]
  subfolder.path <- paste(target.folder, isubfolder, sep = "/")
  if(!dir.exists(subfolder.path)){dir.create(subfolder.path)} else next
}

# Downloading the files
i=1
for(i in 1:nrow(ossl.files)) {

  ifile <- ossl.files[[i,"file_path"]]
  iremote.folder <- ossl.files[[i,"public_path"]]

  remote.file <- paste(iremote.folder, ifile, sep = "/")
  local.file <- paste(target.folder, ifile, sep = "/")

  if(file.exists(local.file)){file.remove(local.file)}

  tryCatch(
    expr = {download.file(remote.file, local.file, mode = "wb")},
    error = function(e){"Retry this file."})

}

# Test reading a downloaded file. For example
test <- qs::qread("~/projects/temp/ossl_models/clay.tot_usda.a334_w.pct/cvpred_nir.neospectra_cubist_ossl_na_v1.2.qs")
