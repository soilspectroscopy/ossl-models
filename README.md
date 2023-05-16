
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5759693.svg)](https://doi.org/10.5281/zenodo.5759693)

# Fitting OSSL models

[<img src="./img/soilspec4gg-logo_fc.png" alt="SoilSpec4GG logo" width="250"/>](https://soilspectroscopy.org/)

This is the repository for all model calibration development for the
[Soil Spectroscopy for Global Good](https://soilspectroscopy.org)
project and based on the [Open Soil Spectral Library
(OSSL)](https://soilspectroscopy.github.io/ossl-manual/) database.

We have used the [MLR3 framework](https://mlr3book.mlr-org.com/) for
fitting machine learning ensemble models.

The `README` of folder [`R-mlr`](R-mlr/README.md) explains the steps
used to calibrate the models used in the OSSL.

In summary, we have provided 5 different model types depending on the
availability of samples in the database, without the use of
geocovariates (`na` code) as predictors.

The model types are composed of two different subsets, i.e. using the
KSSL soil spectral library alone (`kssl` code) or the full OSSL database
(`ossl` code), in combination with three spectral types: VisNIR
(`visnir` code), NIR from the Neospectra instrument (`nir.neospectra`
code), and MIR (`mir` code).

| spectra_type   | subset | geo | model_name                            |
|:---------------|:-------|:----|:--------------------------------------|
| mir            | kssl   | na  | mir_mlr3..eml_kssl_na_v1.2            |
| mir            | ossl   | na  | mir_mlr3..eml_ossl_na_v1.2            |
| nir.neospectra | ossl   | na  | nir.neospectra_mlr3..eml_ossl_na_v1.2 |
| visnir         | kssl   | na  | visnir_mlr3..eml_kssl_na_v1.2         |
| visnir         | ossl   | na  | visnir_mlr3..eml_ossl_na_v1.2         |

The ensemble model (coding name `mlr3..eml`) is composed of a linear
regression (meta-learner, `regr.lm`) of base learners.

Base learners were fitted using Elastic net (`glmnet`), Random Forest
(`ranger`), XGBoost trees (`xgboost`), and Cubist (`cubist`).

MLR3 allows the base learners to feed the meta-learner with
cross-validated (`cv`) or plain calibration predictions (`insample`),
but we used `insample` to speed up processing.

Hyperparameter optimization was done with internal resampling (`inner`)
using 5-fold cross-validation and a smaller subset for speeding up this
operation. This task was performed with a random search of the
hyperparameter space testing up to 20 configurations to find the lowest
RMSE. The final model with the best optimal hyperparameters was fitted
at the end with the full train data.

The search space was composed of crossed combinations, i.e. each model
is not individually tuned. Check [`R-mlr/README.md`](R-mlr/README.md)
for more details.

As predictors, we have use the first 120 PCs of the compressed spectra,
a threshold that considers the tradeoff between spectral representation
and compression magnitude.

The following soil properties have models fitted:

| soil_property                  | description                                         | unit_transform | export_name                         |
|:-------------------------------|:----------------------------------------------------|:---------------|:------------------------------------|
| acidity_usda.a795_cmolc.kg     | Acidity, BaCl2-TEA Extractable, pH 8.2              | log1p          | log..acidity_usda.a795_cmolc.kg     |
| aggstb_usda.a1_w.pct           | Aggregate Stability                                 | log1p          | log..aggstb_usda.a1_w.pct           |
| al.dith_usda.a65_w.pct         | Aluminum, Crystalline, Total Pedogenic Iron         | log1p          | log..al.dith_usda.a65_w.pct         |
| al.ext_usda.a1056_mg.kg        | Aluminum, Extractable, Mehlich3                     | log1p          | log..al.ext_usda.a1056_mg.kg        |
| al.ext_usda.a69_cmolc.kg       | Aluminum, Extractable, KCl                          | log1p          | log..al.ext_usda.a69_cmolc.kg       |
| al.ox_usda.a59_w.pct           | Aluminum, Amorphous, Total Non-Crystalline Iron     | log1p          | log..al.ox_usda.a59_w.pct           |
| awc.33.1500kPa_usda.c80_w.frac | Available Water Content, Difference 33-1500 kPa     | log1p          | log..awc.33.1500kPa_usda.c80_w.frac |
| b.ext_mel3_mg.kg               | Boron, Extractable, Mehlich3                        | log1p          | log..b.ext_mel3_mg.kg               |
| bd_usda.a4_g.cm3               | Bulk Density, \<2mm fraction, Clod                  | log1p          | log..bd_usda.a4_g.cm3               |
| c.tot_usda.a622_w.pct          | Carbon, Total NCS                                   | log1p          | log..c.tot_usda.a622_w.pct          |
| ca.ext_usda.a1059_mg.kg        | Calcium, Extractable, Mehlich3                      | log1p          | log..ca.ext_usda.a1059_mg.kg        |
| ca.ext_usda.a722_cmolc.kg      | Calcium, Extractable, NH4OAc                        | log1p          | log..ca.ext_usda.a722_cmolc.kg      |
| caco3_usda.a54_w.pct           | Carbonate, \<2mm Fraction                           | log1p          | log..caco3_usda.a54_w.pct           |
| cec_usda.a723_cmolc.kg         | CEC, pH 7.0, NH4OAc, 2M KCl displacement            | log1p          | log..cec_usda.a723_cmolc.kg         |
| cf_usda.c236_w.pct             | Coarse Fragments, Greater 2mm                       | log1p          | log..cf_usda.c236_w.pct             |
| clay.tot_usda.a334_w.pct       | Clay                                                | original       | clay.tot_usda.a334_w.pct            |
| cu.ext_usda.a1063_mg.kg        | Copper, Extractable, Mehlich3                       | log1p          | log..cu.ext_usda.a1063_mg.kg        |
| ec_usda.a364_ds.m              | Electrical Conductivity, (w/w)                      | log1p          | log..ec_usda.a364_ds.m              |
| fe.dith_usda.a66_w.pct         | Iron, Crystalline, Total Pedogenic Iron             | log1p          | log..fe.dith_usda.a66_w.pct         |
| fe.ext_usda.a1064_mg.kg        | Iron, Extractable, Mehlich3                         | log1p          | log..fe.ext_usda.a1064_mg.kg        |
| fe.ox_usda.a60_w.pct           | Iron, Amorphous, Total Non-Crystalline Iron         | log1p          | log..fe.ox_usda.a60_w.pct           |
| k.ext_usda.a1065_mg.kg         | Potassium, Extractable, Mehlich3                    | log1p          | log..k.ext_usda.a1065_mg.kg         |
| k.ext_usda.a725_cmolc.kg       | Potassium, Extractable, NH4OAc, 2M KCl displacement | log1p          | log..k.ext_usda.a725_cmolc.kg       |
| mg.ext_usda.a1066_mg.kg        | Magnesium, Extractable, Mehlich3                    | log1p          | log..mg.ext_usda.a1066_mg.kg        |
| mg.ext_usda.a724_cmolc.kg      | Magnesium, Extractable, NH4OAc, 2M KCl displacement | log1p          | log..mg.ext_usda.a724_cmolc.kg      |
| mn.ext_usda.a1067_mg.kg        | Manganese, Extractable, Mehlich3                    | log1p          | log..mn.ext_usda.a1067_mg.kg        |
| mn.ext_usda.a70_mg.kg          | Manganese, Extractable, KCl                         | log1p          | log..mn.ext_usda.a70_mg.kg          |
| n.tot_usda.a623_w.pct          | Nitrogen, Total NCS                                 | log1p          | log..n.tot_usda.a623_w.pct          |
| na.ext_usda.a1068_mg.kg        | Sodium, Extractable, Mehlich3                       | log1p          | log..na.ext_usda.a1068_mg.kg        |
| na.ext_usda.a726_cmolc.kg      | Sodium, Extractable, NH4OAc, 2M KCl displacement    | log1p          | log..na.ext_usda.a726_cmolc.kg      |
| oc_usda.c729_w.pct             | Organic Carbon, Total C without CaCO3, S prep       | log1p          | log..oc_usda.c729_w.pct             |
| p.ext_usda.a1070_mg.kg         | Phosphorus, Extractable, Mehlich3                   | log1p          | log..p.ext_usda.a1070_mg.kg         |
| p.ext_usda.a270_mg.kg          | Phosphorus, Extractable, Bray1                      | log1p          | log..p.ext_usda.a270_mg.kg          |
| p.ext_usda.a274_mg.kg          | Phosphorus, Extractable, Olsen                      | log1p          | log..p.ext_usda.a274_mg.kg          |
| ph.cacl2_usda.a481_index       | pH, 1:2 Soil-CaCl2 Suspension                       | original       | ph.cacl2_usda.a481_index            |
| ph.h2o_usda.a268_index         | pH, 1:1 Soil-Water Suspension                       | original       | ph.h2o_usda.a268_index              |
| s.tot_usda.a624_w.pct          | Sulfur, Total NCS                                   | log1p          | log..s.tot_usda.a624_w.pct          |
| sand.tot_usda.c60_w.pct        | Sand, Total, S prep                                 | original       | sand.tot_usda.c60_w.pct             |
| silt.tot_usda.c62_w.pct        | Silt, Total, S prep                                 | original       | silt.tot_usda.c62_w.pct             |
| wr.1500kPa_usda.a417_w.pct     | Water Retention, 15 Bar (1500 kPa)                  | log1p          | log..wr.1500kPa_usda.a417_w.pct     |
| wr.33kPa_usda.a415_w.pct       | Water Retention, 1/3 Bar (33 kPa)                   | log1p          | log..wr.33kPa_usda.a415_w.pct       |
| zn.ext_usda.a1073_mg.kg        | Zinc, Extractable, Mehlich3                         | log1p          | log..zn.ext_usda.a1073_mg.kg        |

Final evaluation was performed with external (`outer`) 10-fold cross
validation of the tuned models using root mean square error (`rmse`),
mean error (`bias`), R squared (`rsq`), Lin’s concordance correlation
coefficient (`ccc`), and ratio of performance to the interquartile range
(`rpiq`).

Validation plots are available in the [`out/plots`](out/plots) folder.

The final fitted models description along with their performance can be
found in
[out/fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)

# Using OSSL models on your computer

To load the complete analysis-ready models, train data, cross-validated
predictions, validation performance metrics, and validation plot in R,
please use the public URLs described in
[out/fitted_models_access.csv](out/fitted_models_access.csv)

`qs` is a serialized and compressed file format that is faster than
native R `rds`. You need to have [qs
package](https://github.com/traversc/qs) version \>= 0.25.5 to load
files direct from the URLs.

> NOTE: For using the trained models with new spectra, the PC scores
> should be predicted from the [Standard Normal Variate (SNV)
> preprocessed](https://cran.r-project.org/web/packages/prospectr/vignettes/prospectr.html#scatter-and-baseline-corrections)
> spectra following the [same specifications and
> nomenclature](https://soilspectroscopy.github.io/ossl-manual/neospectra-database.html#database-description)
> of the OSSL database (in the previous example we used the
> `nir.neospectra` model).

We provided a prediction in this repository both [examples of
datasets](sample-data) and a [prediction
function](R-mlr/OSSL_functions.R) that incorporates all of these
operations.

Please, check the example datasets for formatting your dataset to the
minimum level required for the prediction function. You can provide
either `csv` files or directly `asd` or opus (`.0`) for VisNIR and MIR
scans, respectively.

With it, the user can get table results for the soil property of
interest (with confidence interval) and a flag (spectral_outlier) which
screens a **potential spectral misrepresentation**.

The prediction functions requires the
[tidyverse](https://www.tidyverse.org/),
[mlr3](https://mlr3book.mlr-org.com/intro.html),
[qs](https://github.com/traversc/qs),
[asdreader](https://github.com/pierreroudier/asdreader),
[opusreader2](https://github.com/spectral-cockpit/opusreader2), and
[matrixStats](https://cran.rstudio.com/web/packages/matrixStats/index.html)
packages.

Parameters:

- `target`: the soil property of interest. Log-transformed properties
  must have `log..` appended at the beginning of the name as indicated
  in the `export_name` column.  
- `spectra.file`: the path for the spectral measurements, either a `csv`
  table (following the sample-data specifications), `asd`, or opus
  (`.0`) file.  
- `spectra.type`: the spectral range of interest. Either `visnir`,
  `nir.neospectra`, or `mir`.  
- `subset.type`: the subset of interest, either the whole `ossl` or the
  `kssl` alone.  
- `geo.type`: only available for `na`.  
- `models.dir`: the path for the `ossl_models` folder. Should follow the
  same structure and naming code as represented in
  [out/fitted_models_access.csv](out/fitted_models_access.csv).

> Please note that the soil properties indication follows the export
> name. In addition, check
> [out/fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)
> for the complete list of models of each soil property as some spectra
> types are not completely available.

``` r
source("R-mlr/OSSL_functions.R")

list.files("sample-data")
```

    ## [1] "101453MD01.asd"             "235157MD01.asd"            
    ## [3] "235157XS01.0"               "icr056141.0"               
    ## [5] "sample_mir_data.csv"        "sample_neospectra_data.csv"
    ## [7] "sample_visnir_data.csv"

``` r
clay.visnir.ossl <- predict.ossl(target = "clay.tot_usda.a334_w.pct",
                                 spectra.file = "sample-data/101453MD01.asd",
                                 spectra.type = "visnir",
                                 subset.type = "ossl",
                                 geo.type = "na",
                                 models.dir = "~/mnt-ossl/ossl_models/")

clay.visnir.ossl |>
  knitr::kable()
```

| sample_id | clay.tot_usda.a334_w.pct | std_error | lower_CI95 | upper_CI95 | spectral_outlier |
|----------:|-------------------------:|----------:|-----------:|-----------:|:-----------------|
|         1 |                 29.18769 | 0.0017156 |   29.18293 |   29.19245 | FALSE            |

``` r
oc.mir.kssl <- predict.ossl(target = "log..oc_usda.c729_w.pct",
                            spectra.file = "sample-data/sample_mir_data.csv",
                            spectra.type = "mir",
                            subset.type = "kssl",
                            geo.type = "na",
                            models.dir = "~/mnt-ossl/ossl_models/")

oc.mir.kssl |>
  knitr::kable()
```

| sample_id | log..oc_usda.c729_w.pct | std_error | lower_CI95 | upper_CI95 | spectral_outlier |
|----------:|------------------------:|----------:|-----------:|-----------:|:-----------------|
|         1 |               0.1005917 | 0.0007044 |  0.0984421 |  0.1027455 | FALSE            |
|         2 |               0.1653977 | 0.0007029 |  0.1631263 |  0.1676737 | FALSE            |
|         3 |               0.5200828 | 0.0005102 |  0.5179318 |  0.5222369 | FALSE            |
|         4 |               1.3324552 | 0.0005617 |  1.3288214 |  1.3360947 | TRUE             |
|         5 |               0.2194784 | 0.0008252 |  0.2166886 |  0.2222746 | FALSE            |
|         6 |               0.8639680 | 0.0005760 |  0.8609902 |  0.8669505 | TRUE             |
|         7 |               1.3676510 | 0.0005056 |  1.3643307 |  1.3709760 | FALSE            |
|         8 |               0.5309905 | 0.0008374 |  0.5274366 |  0.5345526 | FALSE            |
|         9 |               2.0153683 | 0.0003929 |  2.0120813 |  2.0186589 | FALSE            |
|        10 |               1.3325443 | 0.0006622 |  1.3282613 |  1.3368351 | FALSE            |
|        11 |               0.4262675 | 0.0007427 |  0.4233308 |  0.4292104 | FALSE            |
|        12 |              14.5348483 | 0.0007157 | 14.5040220 | 14.5657360 | TRUE             |
|        13 |               1.3003476 | 0.0005039 |  1.2971324 |  1.3035674 | FALSE            |
|        14 |              23.6976542 | 0.0017670 | 23.5768878 | 23.8190140 | TRUE             |
|        15 |               4.7397521 | 0.0017102 |  4.7125853 |  4.7670482 | TRUE             |
|        16 |               0.6108798 | 0.0009207 |  0.6067693 |  0.6150009 | FALSE            |
|        17 |               0.4825308 | 0.0004385 |  0.4807274 |  0.4843363 | FALSE            |
|        18 |              34.0265527 | 0.0008844 | 33.9406932 | 34.1126233 | FALSE            |
|        19 |               0.0926400 | 0.0008598 |  0.0900358 |  0.0952504 | FALSE            |
|        20 |               0.1802776 | 0.0007664 |  0.1777696 |  0.1827909 | FALSE            |

PCA models are available in `ossl_models/pca.ossl` folder:

- <https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_mir_mlr3..eml_kssl_v1.2.qs>.  
- <https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_mir_mlr3..eml_ossl_v1.2.qs>.  
- <https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_nir.neospectra_mlr3..eml_ossl_v1.2.qs>.  
- <https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_visnir_mlr3..eml_kssl_v1.2.qs>.  
- <https://storage.googleapis.com/soilspec4gg-public/models/pca.ossl/mpca_visnir_mlr3..eml_ossl_v1.2.qs>.

If you fit your own models and/or if you are interested in contributing
to this project, please contact us and help us make better open soil
data for global good!

# Other tools and repositories of interest

For more advanced uses of the soil spectral libraries, **we advise to
contacting the original data producers** especially to get help with
using, extending, and improving the original SSL data.

- OSSL Documentation: <https://soilspectroscopy.github.io/ossl-manual/>;
- OSSL Explorer: <https://explorer.soilspectroscopy.org>;
- OSSL Engine: <https://engine.soilspectroscopy.org>;
- OSSL datasets import repository:
  <https://github.com/soilspectroscopy/ossl-imports>;
