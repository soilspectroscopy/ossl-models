
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5759693.svg)](https://doi.org/10.5281/zenodo.5759693)

# Fitting OSSL models

[<img src="./img/soilspec4gg-logo_fc.png" alt="SoilSpec4GG logo" width="250"/>](https://soilspectroscopy.org/)

This is repository for all model fitting development for the [Soil
Spectroscopy for Global Good](https://soilspectroscopy.org) project and
based on the Open Soil Spectral Library
[(OSSL)](https://soilspectroscopy.github.io/ossl-manual/).

The `README` of folder [`R-mlr`](R-mlr/README.md) explains the steps
used to calibrate the models used in the OSSL.

In summary, we have provided 5 different model types depending on the
availability of samples in the OSSL database, without the use of
geocovariates (`na` code).

The model types are composed of two different subsets, i.e. using the
KSSL soil spectral library alone (`kssl` code) or the full OSSL database
(`ossl` code), in combinations with three spectral ranges: VisNIR
(`visnir` code), NIR from the Neospectra instrument (`nir.neospectra`
code), and MIR (`mir` code).

We have used the [MLR3 framework](https://mlr3book.mlr-org.com/) for
fitting machine learning ensemble models.

| spectra_type   | subset | geo | model_name                            |
|:---------------|:-------|:----|:--------------------------------------|
| mir            | kssl   | na  | mir_mlr3..eml_kssl_na_v1.2            |
| mir            | ossl   | na  | mir_mlr3..eml_ossl_na_v1.2            |
| nir.neospectra | ossl   | na  | nir.neospectra_mlr3..eml_ossl_na_v1.2 |
| visnir         | kssl   | na  | visnir_mlr3..eml_kssl_na_v1.2         |
| visnir         | ossl   | na  | visnir_mlr3..eml_ossl_na_v1.2         |

The ensemble model (coding name `mlr3..eml`) is composed of a linear
regression (meta-learner, `regr.lm`) of base learners.

Base learners employed were: Elastic net (`glmnet`), Random Forest
(`ranger`), XGBoost trees (`xgboost`), and Cubist (`cubist`).

MLR3 allows the base learners to feed the meta-learner with
cross-validated (`cv`) or plain calibration predictions (`insample`),
but we used `insample` to speed up processing.

Hyperparameter optimization was done with internal resampling (`inner`)
using 5-fold cross-validation and a smaller subset for speeding up this
operation. This task was performed with random search of the
hyperparameter space testing up to 20 configurations to find the lowest
RMSE. The final model with the best optimal hyperparameters was fitted
at the end with the full train data.

The search space was composed of crossed combinations, i.e. each model
is not individually tuned:

For predictors, we have use the first 120 PCs of the compressed spectra,
a threshold that considers the tradeoff between spectral representation
and compression magnitude.

The following soil properties have models fitted:

| soil_property                  | description                                         | unit_transform |
|:-------------------------------|:----------------------------------------------------|:---------------|
| acidity_usda.a795_cmolc.kg     | Acidity, BaCl2-TEA Extractable, pH 8.2              | log1p          |
| aggstb_usda.a1_w.pct           | Aggregate Stability                                 | log1p          |
| al.dith_usda.a65_w.pct         | Aluminum, Crystalline, Total Pedogenic Iron         | log1p          |
| al.ext_usda.a1056_mg.kg        | Aluminum, Extractable, Mehlich3                     | log1p          |
| al.ext_usda.a69_cmolc.kg       | Aluminum, Extractable, KCl                          | log1p          |
| al.ox_usda.a59_w.pct           | Aluminum, Amorphous, Total Non-Crystalline Iron     | log1p          |
| awc.33.1500kPa_usda.c80_w.frac | Available Water Content, Difference 33-1500 kPa     | log1p          |
| b.ext_mel3_mg.kg               | Boron, Extractable, Mehlich3                        | log1p          |
| bd_usda.a4_g.cm3               | Bulk Density, \<2mm fraction, Clod                  | log1p          |
| c.tot_usda.a622_w.pct          | Carbon, Total NCS                                   | log1p          |
| ca.ext_usda.a1059_mg.kg        | Calcium, Extractable, Mehlich3                      | log1p          |
| ca.ext_usda.a722_cmolc.kg      | Calcium, Extractable, NH4OAc                        | log1p          |
| caco3_usda.a54_w.pct           | Carbonate, \<2mm Fraction                           | log1p          |
| cec_usda.a723_cmolc.kg         | CEC, pH 7.0, NH4OAc, 2M KCl displacement            | log1p          |
| cf_usda.c236_w.pct             | Coarse Fragments, Greater 2mm                       | log1p          |
| clay.tot_usda.a334_w.pct       | Clay                                                | original       |
| cu.ext_usda.a1063_mg.kg        | Copper, Extractable, Mehlich3                       | log1p          |
| ec_usda.a364_ds.m              | Electrical Conductivity, (w/w)                      | log1p          |
| fe.dith_usda.a66_w.pct         | Iron, Crystalline, Total Pedogenic Iron             | log1p          |
| fe.ext_usda.a1064_mg.kg        | Iron, Extractable, Mehlich3                         | log1p          |
| fe.ox_usda.a60_w.pct           | Iron, Amorphous, Total Non-Crystalline Iron         | log1p          |
| k.ext_usda.a1065_mg.kg         | Potassium, Extractable, Mehlich3                    | log1p          |
| k.ext_usda.a725_cmolc.kg       | Potassium, Extractable, NH4OAc, 2M KCl displacement | log1p          |
| mg.ext_usda.a1066_mg.kg        | Magnesium, Extractable, Mehlich3                    | log1p          |
| mg.ext_usda.a724_cmolc.kg      | Magnesium, Extractable, NH4OAc, 2M KCl displacement | log1p          |
| mn.ext_usda.a1067_mg.kg        | Manganese, Extractable, Mehlich3                    | log1p          |
| mn.ext_usda.a70_mg.kg          | Manganese, Extractable, KCl                         | log1p          |
| n.tot_usda.a623_w.pct          | Nitrogen, Total NCS                                 | log1p          |
| na.ext_usda.a1068_mg.kg        | Sodium, Extractable, Mehlich3                       | log1p          |
| na.ext_usda.a726_cmolc.kg      | Sodium, Extractable, NH4OAc, 2M KCl displacement    | log1p          |
| oc_usda.c729_w.pct             | Organic Carbon, Total C without CaCO3, S prep       | log1p          |
| p.ext_usda.a1070_mg.kg         | Phosphorus, Extractable, Mehlich3                   | log1p          |
| p.ext_usda.a270_mg.kg          | Phosphorus, Extractable, Bray1                      | log1p          |
| p.ext_usda.a274_mg.kg          | Phosphorus, Extractable, Olsen                      | log1p          |
| ph.cacl2_usda.a481_index       | pH, 1:2 Soil-CaCl2 Suspension                       | original       |
| ph.h2o_usda.a268_index         | pH, 1:1 Soil-Water Suspension                       | original       |
| s.tot_usda.a624_w.pct          | Sulfur, Total NCS                                   | log1p          |
| sand.tot_usda.c60_w.pct        | Sand, Total, S prep                                 | original       |
| silt.tot_usda.c62_w.pct        | Silt, Total, S prep                                 | original       |
| wr.1500kPa_usda.a417_w.pct     | Water Retention, 15 Bar (1500 kPa)                  | log1p          |
| wr.33kPa_usda.a415_w.pct       | Water Retention, 1/3 Bar (33 kPa)                   | log1p          |
| zn.ext_usda.a1073_mg.kg        | Zinc, Extractable, Mehlich3                         | log1p          |

Final evaluation was performed with external (`outer`) 10-fold cross
validation of the final models using root mean square error (`rmse`),
mean error (`bias`),

Validation plots are available in the [`out/plots`](out/plots) folder.

Final fitted models along with their performance can be found in
[out/fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)

# Using OSSL models in your computer

To load the complete analysis-ready models, train data, cross-validated
predictions, validation performance metrics and validation plot in R,
please use the public urls described in
[out/fitted_models_access.csv](out/fitted_models_access.csv)

`qs` is a serialized and compressed file format that is faster than
native R `rds`. You need to have [qs
package](https://github.com/traversc/qs) version \>= 0.25.5 to load
files direct from the URLs.

Example:

``` r
library("tidyverse")
library("mlr3verse")
library("qs")

files.table <- readr::read_csv("out/fitted_models_access.csv", show_col_types = F)

files.table %>%
  dplyr::filter(grepl("clay", soil_property)) %>%
  dplyr::filter(grepl("nir.neospectra", model_name)) %>%
  dplyr::select(file_url)

trained.model <- qs::qread_url("https://storage.googleapis.com/soilspec4gg-public/models/clay.tot_usda.a334_w.pct/model_nir.neospectra_mlr3..eml_ossl_na_v1.2.qs")

train.data <- qs::qread_url("https://storage.googleapis.com/soilspec4gg-public/models/clay.tot_usda.a334_w.pct/task_nir.neospectra_mlr3..eml_ossl_na_v1.2.qs")

task <- mlr3::as_task_regr(train.data, id = "train", target = "clay.tot_usda.a334_w.pct", type = "regression")
task$set_col_roles("id.layer_uuid_txt", roles = "name")

predictions <- trained.model$predict(task) %>% data.table::as.data.table()

predictions <- dplyr::left_join(task$row_names, predictions, by = c("row_id" = "row_ids"))
```

> NOTE: For using the trained model with new predictions, the PC scores
> should be predicted from the [Standard Normal Variate (SNV)
> preprocessed](https://cran.r-project.org/web/packages/prospectr/vignettes/prospectr.html#scatter-and-baseline-corrections)
> spectra following the [same specifications and
> nomenclature](https://soilspectroscopy.github.io/ossl-manual/neospectra-database.html#database-description)
> of the OSSL database (in the previous example we used the
> `nir.neospectra` model).

PCA models are available in:  
-<https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_mir_mlr3..eml_kssl_v1.2.qs>.  
-<https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_mir_mlr3..eml_ossl_v1.2.qs>.  
-<https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_nir.neospectra_mlr3..eml_ossl_v1.2.qs>.  
-<https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_visnir_mlr3..eml_kssl_v1.2.qs>.  
-<https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_visnir_mlr3..eml_ossl_v1.2.qs>.

Example:

``` r
pca.model <- qread_url("https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_nir.neospectra_mlr3..eml_ossl_v1.2.qs")

class(pca.model) <- "prcomp"

new.scores <- predict(pca.model, new.spectra)
```

If you fit your own models and/or if you are interested in contributing
to this project, please contact us and help us make better open soil
data for global good!

# Other tools and repositories of interest

For more advanced uses of the soil spectral libraries **we advise to
contact the original data producers** especially to get help with using,
extending, and improving the original SSL data.

- OSSL Documentation: <https://soilspectroscopy.github.io/ossl-manual/>;
- OSSL Explorer: <https://explorer.soilspectroscopy.org>;
- OSSL Engine: <https://engine.soilspectroscopy.org>;
- OSSL datasets import repository:
  <https://github.com/soilspectroscopy/ossl-imports>;
