
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5759693.svg)](https://doi.org/10.5281/zenodo.5759693)

# Fitting of global Soil Spectral models

[<img src="./img/soilspec4gg-logo_fc.png" alt="SoilSpec4GG logo" width="250"/>](https://soilspectroscopy.org/)

A repository for all model fitting development work for the [Soil Spectroscopy for Global
Good](https://soilspectroscopy.org) project.

Important files in this repository:

- The metadata table [OSSL_models_meta.csv](out/OSSL_models_meta.csv) contains total list of models currently available;  
- The table [global_layers1km.csv](out/global_layers1km.csv) contains URL of the global geographical covariates [WorldClim2.1](https://www.worldclim.org/data/worldclim21.html) layers and [MODIS LST](https://doi.org/10.5281/zenodo.1420114); 
- Folder [`R-mlr`](R-mlr/README.md) contains explanation of modeling steps used to produce default models used in the OSSL;
- Script [`R/model_accuracy_stats.R`](R-mlr/model_accuracy_stats.R) is used to generate accuracy plots;  

For more advanced uses of the soil spectral libraries **we advise to contact the original data producers** 
especially to get help with using, extending and improving the original SSL data.

To load the complete analysis-ready dataset (486MB) as a single table in R and run predictive modeling please use:

```
rep = "http://s3.us-east-1.wasabisys.com/soilspectroscopy/"
rm.ossl = readRDS(url(paste0(rep, "ossl_import/rm.ossl_v1.rds", "rb")))
dim(rm.ossl)
## 152,146 obs. of 2962 variables
```

If you fit your own models and/or if you are interested in contributing 
to this project please contact us and help us make better open soil data for global good!

Other tools and repositories of interest:

- OSSL documentation: <https://soilspectroscopy.github.io/ossl-manual/>;
- OSSL Explorer: <https://explorer.soilspectroscopy.org>;
- OSSL Engine: <https://engine.soilspectroscopy.org>;
- Data import repository: <https://github.com/soilspectroscopy/ossl-imports>;
