
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5759693.svg)](https://doi.org/10.5281/zenodo.5759693)

# Fitting OSSL models

[<img src="./img/soilspec4gg-logo_fc.png" alt="SoilSpec4GG logo" width="250"/>](https://soilspectroscopy.org/)

This is the repository for all model calibration development for the
[Soil Spectroscopy for Global Good](https://soilspectroscopy.org)
project and based on the [Open Soil Spectral Library
(OSSL)](https://soilspectroscopy.github.io/ossl-manual/) database.

We have used the [MLR3 framework](https://mlr3book.mlr-org.com/) for
fitting machine learning models, specifically with the [Cubist
algorithm](https://cran.r-project.org/web/packages/Cubist/vignettes/cubist.html).

The `README` of folder [`R-mlr`](R-mlr/README.md) explains the steps
used to calibrate the models, while this page provides an overview of
the modeling and results.

In summary, we have provided 5 different model types depending on the
availability of samples in the database, without the use of
geocovariates (`na` code), i.e., site information or environmental
layers are not used, only the spectra data.

The model types are composed of two different subsets, i.e. using the
KSSL soil spectral library alone (`kssl` code) or the full OSSL database
(`ossl` code), in combination with three spectral types: VisNIR
(`visnir` code), NIR from the Neospectra instrument (`nir.neospectra`
code), and MIR (`mir` code).

After running some internal evaluations, we recommend using the KSSL
models only when i) there is not any model available using the full OSSL
database for a given combination of soil property and spectral region of
interest; ii) the spectra to be predicted has the same instrument
manufacturer/model as the spectrometers used to build the KSSL VisNIR
and MIR libraries; iii) the KSSL library is representative for the new
spectra, based both on spectral similarity and range of soil properties
of interest of the target samples. **Otherwise, use the models with
`ossl` acronym**.

| spectra_type   | subset | geo | model_name                         |
|:---------------|:-------|:----|:-----------------------------------|
| mir            | kssl   | na  | mir_cubist_kssl_na_v1.2            |
| mir            | ossl   | na  | mir_cubist_ossl_na_v1.2            |
| nir.neospectra | ossl   | na  | nir.neospectra_cubist_ossl_na_v1.2 |
| visnir         | kssl   | na  | visnir_cubist_kssl_na_v1.2         |
| visnir         | ossl   | na  | visnir_cubist_ossl_na_v1.2         |

The machine learning algorithm Cubist (coding name `cubist`) takes
advantage of a decision-tree splitting method but fits linear regression
models at each terminal leaf. It also uses a boosting mechanism
(sequential trees adjusted by weights) that allows the growth of a
forest by tuning the number of committees. We haven’t used the
correction of final predictions by the nearest neighbors’ influence due
to the lack of this feature in the MLR3 framework.

Hyperparameter optimization was done with internal resampling (`inner`)
using 5-fold cross-validation and a smaller subset for speeding up this
operation. This task was performed with a grid search of the
hyperparameter space testing up to 5 configurations to find the lowest
RMSE. The final model with the best optimal hyperparameters was fitted
at the end with the full train data.

As predictors, we have used the first 120 PCs of the compressed spectra,
a threshold that considers the trade-off between spectral representation
and compression magnitude. The remaining farther components were used
for the trustworthiness test, i.e., if a sample is underrepresented in
respect of the feature space from the calibration set because of their
unique features presented in those farther components, then a flag is
raised and indicated in the results.

Soil properties with available models can be found in
[fitted_modeling_combinations_v1.2.csv](out/fitted_modeling_combinations_v1.2.csv).

Final evaluation was performed with external (`outer`) 10-fold cross
validation of the tuned models using root mean square error (`rmse`),
mean error (`bias`), R squared (`rsq`), Lin’s concordance correlation
coefficient (`ccc`), and the ratio of performance to the interquartile
range (`rpiq`).

The cross-validated predictions were used to estimate the unbiased error
(absolute residual), which were further employed to calibrate
uncertainty models via [conformal prediction
intervals](https://en.wikipedia.org/wiki/Conformal_prediction). The
error model is calibrated using the same fine-tuned structure of the
respective response model. Conformity scores are estimated for a defined
confidence level, which was defined to 68% to approximate one standard
deviation.

Some soil properties were natural-log transformed (with offset = 1, that
is why we use `log1p` function) to improve the prediction performance of
highly skewed distributed soil properties. They were back-transformed
only at the end after running all modeling steps, including performance
estimation and the definition of the uncertainty intervals.

The final fitted models along with their performance metrics can be
found in
[fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)

Validation plots are available in the [`out/plots`](out/plots) folder.

# OSSL Engine

The [OSSL Engine](https://engine.soilspectroscopy.org/) is a web
platform where the user can upload spectra from the VisNIR (400-2500
nm), NIR (1350-2550 nm), or MIR (600-4000 cm<sup>-1</sup>) ranges and
get predictions back with uncertainty estimation and the
representativeness flag. Please, check it out!

Please, check the example datasets for formatting your spectra to the
minimum required level. You can provide either csv files or directly asd
or opus (.0) for VisNIR and MIR scans, respectively.

We recommend using the OSSL model type when getting predictions. KSSL is
recommended when the spectra to be predicted have the same instrument
manufacturer/model as the KSSL library and are represented by the range
of soil properties of interest. Otherwise, use the OSSL subset.

# Using OSSL models on your computer

To load the complete analysis-ready models, training data,
cross-validated predictions, validation performance metrics, and
validation plots in R, please download them to the same directory
structure as specified in
[fitted_models_access.csv](out/fitted_models_access.csv).

Script
**[script_04f_download_ossl_models.R](R-mlr/script_04f_download_ossl_models.R)**
describes the automated steps to grab all the required files to your
computer.

`qs` is a serialized and compressed file format that is faster than
native R `rds`. You need to have [qs
package](https://github.com/traversc/qs) version \>= 0.25.5 to load
files direct from the URLs.

> NOTE: For using the trained models on new spectra, the spectra must
> have the [same
> range](https://soilspectroscopy.github.io/ossl-manual/neospectra-database.html#database-description)
> of the OSSL models, i.e. 400-2500 nm for VisNIR, 600-4000
> cm<sup>-1</sup> for MIR, and 1350-2550 for NIR (Neospectra)..

We provided in this repository both [examples of datasets](sample-data)
and a [prediction function](R-mlr/OSSL_functions.R) that preprocess and
provide all outputs.

Please, check the example datasets for formatting your spectra to the
minimum required level of the prediction function. You can provide
either `csv` files or directly `asd` or opus (`.0`) for VisNIR and MIR
scans, respectively.

The results table has the prediction value (already back-transformed if
log transformation was used) for the soil property of interest, standard
deviation and uncertainty band, and a flag column for **potential
underrepresented** samples given the OSSL calibration data, which is
calculated based on principal components and Q statistics.

The prediction function requires the
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
  [fitted_models_access.csv](out/fitted_models_access.csv).

All files that represent the **ossl_models** directory tree for local
run or online access are described in
[ossl_models_directory_tree.csv](out/ossl_models_directory_tree.csv).

> Please note that the soil properties indication follows the export
> name. Check
> [fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)
> for the complete list of models of each soil property as some spectral
> types may not be available. More importantly, for natural-log soil
> properties, the upper and lower bands were estimated before the
> back-transformation.

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

| sample_id | clay.tot_usda.a334_w.pct |  std_dev |    lower |    upper | underrepresented |
|----------:|-------------------------:|---------:|---------:|---------:|:-----------------|
|         1 |                 33.35406 | 9.900726 | 23.45334 | 43.25479 | FALSE            |

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

| sample_id | oc_usda.c729_w.pct |   std_dev |      lower |      upper | underrepresented |
|----------:|-------------------:|----------:|-----------:|-----------:|:-----------------|
|         1 |          0.0232091 | 0.0387254 | -0.0149378 |  0.0628333 | FALSE            |
|         2 |          0.0438117 | 0.0774065 | -0.0311812 |  0.1246095 | FALSE            |
|         3 |          0.4533481 | 0.0906541 |  0.3325472 |  0.5851001 | FALSE            |
|         4 |          1.0724492 | 0.0625031 |  0.9505346 |  1.2019838 | TRUE             |
|         5 |          0.0490533 | 0.0355205 |  0.0130686 |  0.0863161 | FALSE            |
|         6 |          0.8021375 | 0.0527328 |  0.7118660 |  0.8971693 | TRUE             |
|         7 |          1.2250589 | 0.0807371 |  1.0588345 |  1.4047039 | FALSE            |
|         8 |          0.3362929 | 0.1103006 |  0.2035415 |  0.4836867 | FALSE            |
|         9 |          1.9155340 | 0.0538360 |  1.7665917 |  2.0724947 | FALSE            |
|        10 |          1.3158559 | 0.0914151 |  1.1218836 |  1.5275601 | FALSE            |
|        11 |          0.4234255 | 0.0786920 |  0.3195847 |  0.5354376 | FALSE            |
|        12 |         13.3578608 | 0.2007397 | 10.9575133 | 16.2400532 | TRUE             |
|        13 |          1.4059354 | 0.1456846 |  1.0999981 |  1.7564430 | FALSE            |
|        14 |         33.1607521 | 0.0746672 | 30.7872847 | 35.7114397 | TRUE             |
|        15 |          6.1175378 | 0.1536174 |  5.1697557 |  7.2109157 | TRUE             |
|        16 |          0.5058015 | 0.0990479 |  0.3700963 |  0.6549481 | FALSE            |
|        17 |          0.3885227 | 0.0553760 |  0.3156664 |  0.4654135 | FALSE            |
|        18 |         34.2093266 | 0.0610237 | 32.1842980 | 36.3579298 | FALSE            |
|        19 |          0.2090189 | 0.0927047 |  0.1064461 |  0.3211006 | FALSE            |
|        20 |          0.2358387 | 0.0734969 |  0.1512271 |  0.3266690 | FALSE            |

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
