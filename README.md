
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

Soil properties with available models can be found in
[fitted_modeling_combinations_v1.2.csv](out/fitted_modeling_combinations_v1.2.csv).

Final evaluation was performed with external (`outer`) 10-fold cross
validation of the tuned models using root mean square error (`rmse`),
mean error (`bias`), R squared (`rsq`), Lin’s concordance correlation
coefficient (`ccc`), and ratio of performance to the interquartile range
(`rpiq`).

The final fitted models description along with their performance can be
found in
[fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)

Validation plots are available in the [`out/plots`](out/plots) folder.

# Using OSSL models on your computer

To load the complete analysis-ready models, train data, cross-validated
predictions, validation performance metrics, and validation plot in R,
please use the public URLs described in
[fitted_models_access.csv](out/fitted_models_access.csv)

`qs` is a serialized and compressed file format that is faster than
native R `rds`. You need to have [qs
package](https://github.com/traversc/qs) version \>= 0.25.5 to load
files direct from the URLs.

> NOTE: For using the trained models with new spectra, the PC scores
> should be predicted from the [Standard Normal Variate (SNV)
> preprocessed](https://cran.r-project.org/web/packages/prospectr/vignettes/prospectr.html#scatter-and-baseline-corrections)
> spectra following the [same specifications and
> nomenclature](https://soilspectroscopy.github.io/ossl-manual/neospectra-database.html#database-description)
> of the OSSL database.

We provided in this repository both [examples of datasets](sample-data)
and a [prediction function](R-mlr/OSSL_functions.R) that incorporates
all of above operations.

Please, check the example datasets for formatting your spectra to the
minimum level required for the prediction function. You can provide
either `csv` files or directly `asd` or opus (`.0`) for VisNIR and MIR
scans, respectively.

With it, the user can get table results for the soil property of
interest (with uncertainty) and a flag column for a **potential spectral
misrepresentation** of the OSSL models to your measurements. The
uncertainty was derived from the standard error of the predictions
obtained from the linear model (meta-learner) of the ensemble models.

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
> name. In addition, check
> [fitted_models_performance_v1.2.csv](out/fitted_models_performance_v1.2.csv)
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
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(target)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(target))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
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

| sample_id | oc_usda.c729_w.pct | std_error | lower_CI95 | upper_CI95 | spectral_outlier |
|----------:|-------------------:|----------:|-----------:|-----------:|:-----------------|
|         1 |          0.1005917 | 0.0007044 |  0.0984421 |  0.1027455 | FALSE            |
|         2 |          0.1653977 | 0.0007029 |  0.1631263 |  0.1676737 | FALSE            |
|         3 |          0.5200828 | 0.0005102 |  0.5179318 |  0.5222369 | FALSE            |
|         4 |          1.3324552 | 0.0005617 |  1.3288214 |  1.3360947 | TRUE             |
|         5 |          0.2194784 | 0.0008252 |  0.2166886 |  0.2222746 | FALSE            |
|         6 |          0.8639680 | 0.0005760 |  0.8609902 |  0.8669505 | TRUE             |
|         7 |          1.3676510 | 0.0005056 |  1.3643307 |  1.3709760 | FALSE            |
|         8 |          0.5309905 | 0.0008374 |  0.5274366 |  0.5345526 | FALSE            |
|         9 |          2.0153683 | 0.0003929 |  2.0120813 |  2.0186589 | FALSE            |
|        10 |          1.3325443 | 0.0006622 |  1.3282613 |  1.3368351 | FALSE            |
|        11 |          0.4262675 | 0.0007427 |  0.4233308 |  0.4292104 | FALSE            |
|        12 |         14.5348483 | 0.0007157 | 14.5040220 | 14.5657360 | TRUE             |
|        13 |          1.3003476 | 0.0005039 |  1.2971324 |  1.3035674 | FALSE            |
|        14 |         23.6976542 | 0.0017670 | 23.5768878 | 23.8190140 | TRUE             |
|        15 |          4.7397521 | 0.0017102 |  4.7125853 |  4.7670482 | TRUE             |
|        16 |          0.6108798 | 0.0009207 |  0.6067693 |  0.6150009 | FALSE            |
|        17 |          0.4825308 | 0.0004385 |  0.4807274 |  0.4843363 | FALSE            |
|        18 |         34.0265527 | 0.0008844 | 33.9406932 | 34.1126233 | FALSE            |
|        19 |          0.0926400 | 0.0008598 |  0.0900358 |  0.0952504 | FALSE            |
|        20 |          0.1802776 | 0.0007664 |  0.1777696 |  0.1827909 | FALSE            |

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
