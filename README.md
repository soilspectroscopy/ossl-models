# Fitting of global Soil Spectral models

[<img src="./img/soilspec4gg-logo_fc.png" alt="SoilSpec4GG logo" width="250"/>](https://soilspectroscopy.org/)

A repository for all model fitting development work for the [Soil Spectroscopy for Global
Good](https://soilspectroscopy.org) project.

Type of models considered:

- The metadata table [OSSL_models_meta.csv](out/OSSL_models_meta.csv) contains total list of models currently available;  
- Folder `R` contains modeling steps used to produce default models used in the OSSL;
- Script `R/model_accuracy_stats.R` is used to generate accuracy plots;  

For more advanced uses of the soil spectral libraries **we advise to contact the original data producers** 
especially to get help with using, extending and improving the original SSL data.

Other tools and repositories of interest:

- OSSL documentation: <https://soilspectroscopy.github.io/ossl-manual/>;
- OSSL Explorer: <https://explorer.soilspectroscopy.org>;
- OSSL Engine: <https://engine.soilspectroscopy.org>;
- Data import repository: <https://github.com/soilspectroscopy/ossl-imports>;
