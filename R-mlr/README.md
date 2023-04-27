OSSL: Global Soil Spectral Calibration Models
================
Jose L. Safanelli (<jsafanelli@woodwellclimate.org>), Tomislav Hengl
(<tom.hengl@opengeohub.org>), Leandro Parente
(<leandro.parente@opengeohub.org>), and Jonathan Sanderman
(<jsanderman@woodwellclimate.org>)
27 April, 2023



[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5759693.svg)](https://doi.org/10.5281/zenodo.5759693)

[<img src="../img/soilspec4gg-logo_fc.png" alt="SoilSpec4GG logo" width="250"/>](https://soilspectroscopy.org/)

[<img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" />](http://creativecommons.org/licenses/by-sa/4.0/)

This work is licensed under a [Creative Commons Attribution-ShareAlike
4.0 International
License](http://creativecommons.org/licenses/by-sa/4.0/).

## Open Soil Spectral Library

Part of: <https://github.com/soilspectroscopy>  
Project: [Soil Spectroscopy for Global
Good](https://soilspectroscopy.org)  
Last update: 2023-04-27  
Dataset:
[OSSL](https://soilspectroscopy.github.io/ossl-manual/ossl-database-description.html)

    ## Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
    ## had status 1

The directory/folder path:

``` r
dir <- "/mnt/soilspec4gg/ossl/ossl_models/"
db.dir <- "/mnt/soilspec4gg/ossl/ossl_import/"
```

### Overview

This tutorial explains the steps required to fit Global Soil Spectral
Calibration Models for the purpose of the [Soil Spectroscopy for Global
Good project](https://soilspectroscopy.org).

We use Ensemble Machine Learning with meta-learner/stacking where the
meta-learner is fitted using 5-fold cross-validation with model
optimization ([Hengl & MacMillan, 2019](#ref-hengl2019predictive)).

Modeling framework is based on the [mlr
package](https://mlr.mlr-org.com/) ([Bischl et al.,
2016](#ref-bischl2016mlr)), which is currently not developed further.
Migration to the newer modeling framework
[mlr3](https://mlr3book.mlr-org.com/) is pending.

The following code is implemented in a high-performance computing
environment with many steps fully parallelized and enough RAM (&gt;
450GB). Running this code on a standard computer without subsetting data
is not recommended.

For each soil property of interest, the possible model combinations are:

``` r
# MIR, VisNIR and NIR from Neospectra
spectra.type <- c("mir", "visnir", "nir.neospectra")

# Using KSSL only or the whole OSSL
subset <- c("kssl", "ossl")

# Adding (ll) or not (na) geocovariates to models
geocovariates <- c("ll", "na")

# Basic structure
modeling.combinations <- tibble(spectra_type = spectra.type) %>%
  crossing(subset = subset) %>%
  crossing(geo = geocovariates)
```

At this stage, Neospectra models are fitted using only the OSSL subset,
as they have minimal differences (OSSL \~ KSSL). In addition, models
with geocovariates are omitted as additional tests are required.

``` r
# Target models
modeling.combinations <- modeling.combinations %>%
  dplyr::filter(!(grepl("neospectra", spectra_type) & subset == "kssl")) %>%
  dplyr::filter(!(geo == "ll"))

# Model name
modeling.combinations <- modeling.combinations %>%
  mutate(model_name = paste0(spectra_type, "_mlr..eml_", subset, "_", geo, "_v1.2"), .before = 1)

modeling.combinations %>%
  knitr::kable()
```

| model\_name                              | spectra\_type  | subset | geo |
|:-----------------------------------------|:---------------|:-------|:----|
| mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |

After database update v1.2, the following properties of interest have
models fitted:

``` r
soil.properties <- c("silt.tot_usda.c62_w.pct",
                     "clay.tot_usda.a334_w.pct",
                     "sand.tot_usda.c60_w.pct",
                     "log..oc_usda.c729_w.pct",
                     "log..n.tot_usda.a623_w.pct",
                     "log..caco3_usda.a54_w.pct",
                     "log..al.ox_usda.a59_w.pct",
                     "ph.h2o_usda.a268_index",
                     "ph.cacl2_usda.a481_index",
                     "log..k.ext_usda.a725_cmolc.kg",
                     "log..mg.ext_usda.a724_cmolc.kg",
                     "log..ca.ext_usda.a722_cmolc.kg",
                     "acidity_usda.a795_cmolc.kg",
                     "bd_usda.a4_g.cm3")

soil.properties.original <- gsub("log..", "", soil.properties)
```

**Note: some soil properties were natural-log transformed to help
models.**

Final modeling combinations:

``` r
modeling.combinations <- tibble(soil_property = soil.properties) %>%
  crossing(modeling.combinations)

modeling.combinations %>%
  knitr::kable()
```

| soil\_property                   | model\_name                              | spectra\_type  | subset | geo |
|:---------------------------------|:-----------------------------------------|:---------------|:-------|:----|
| acidity\_usda.a795\_cmolc.kg     | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| acidity\_usda.a795\_cmolc.kg     | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| acidity\_usda.a795\_cmolc.kg     | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| acidity\_usda.a795\_cmolc.kg     | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| acidity\_usda.a795\_cmolc.kg     | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| bd\_usda.a4\_g.cm3               | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| bd\_usda.a4\_g.cm3               | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| bd\_usda.a4\_g.cm3               | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| bd\_usda.a4\_g.cm3               | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| bd\_usda.a4\_g.cm3               | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| clay.tot\_usda.a334\_w.pct       | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| clay.tot\_usda.a334\_w.pct       | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| clay.tot\_usda.a334\_w.pct       | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| clay.tot\_usda.a334\_w.pct       | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| clay.tot\_usda.a334\_w.pct       | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..al.ox\_usda.a59\_w.pct      | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..al.ox\_usda.a59\_w.pct      | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..al.ox\_usda.a59\_w.pct      | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..al.ox\_usda.a59\_w.pct      | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..al.ox\_usda.a59\_w.pct      | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..ca.ext\_usda.a722\_cmolc.kg | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..ca.ext\_usda.a722\_cmolc.kg | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..ca.ext\_usda.a722\_cmolc.kg | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..ca.ext\_usda.a722\_cmolc.kg | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..ca.ext\_usda.a722\_cmolc.kg | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..caco3\_usda.a54\_w.pct      | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..caco3\_usda.a54\_w.pct      | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..caco3\_usda.a54\_w.pct      | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..caco3\_usda.a54\_w.pct      | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..caco3\_usda.a54\_w.pct      | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..k.ext\_usda.a725\_cmolc.kg  | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..k.ext\_usda.a725\_cmolc.kg  | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..k.ext\_usda.a725\_cmolc.kg  | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..k.ext\_usda.a725\_cmolc.kg  | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..k.ext\_usda.a725\_cmolc.kg  | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..mg.ext\_usda.a724\_cmolc.kg | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..mg.ext\_usda.a724\_cmolc.kg | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..mg.ext\_usda.a724\_cmolc.kg | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..mg.ext\_usda.a724\_cmolc.kg | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..mg.ext\_usda.a724\_cmolc.kg | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..n.tot\_usda.a623\_w.pct     | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..n.tot\_usda.a623\_w.pct     | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..n.tot\_usda.a623\_w.pct     | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..n.tot\_usda.a623\_w.pct     | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..n.tot\_usda.a623\_w.pct     | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| log..oc\_usda.c729\_w.pct        | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| log..oc\_usda.c729\_w.pct        | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| log..oc\_usda.c729\_w.pct        | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| log..oc\_usda.c729\_w.pct        | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| log..oc\_usda.c729\_w.pct        | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| ph.cacl2\_usda.a481\_index       | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| ph.cacl2\_usda.a481\_index       | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| ph.cacl2\_usda.a481\_index       | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| ph.cacl2\_usda.a481\_index       | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| ph.cacl2\_usda.a481\_index       | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| ph.h2o\_usda.a268\_index         | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| ph.h2o\_usda.a268\_index         | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| ph.h2o\_usda.a268\_index         | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| ph.h2o\_usda.a268\_index         | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| ph.h2o\_usda.a268\_index         | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| sand.tot\_usda.c60\_w.pct        | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| sand.tot\_usda.c60\_w.pct        | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| sand.tot\_usda.c60\_w.pct        | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| sand.tot\_usda.c60\_w.pct        | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| sand.tot\_usda.c60\_w.pct        | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |
| silt.tot\_usda.c62\_w.pct        | mir\_mlr..eml\_kssl\_na\_v1.2            | mir            | kssl   | na  |
| silt.tot\_usda.c62\_w.pct        | mir\_mlr..eml\_ossl\_na\_v1.2            | mir            | ossl   | na  |
| silt.tot\_usda.c62\_w.pct        | nir.neospectra\_mlr..eml\_ossl\_na\_v1.2 | nir.neospectra | ossl   | na  |
| silt.tot\_usda.c62\_w.pct        | visnir\_mlr..eml\_kssl\_na\_v1.2         | visnir         | kssl   | na  |
| silt.tot\_usda.c62\_w.pct        | visnir\_mlr..eml\_ossl\_na\_v1.2         | visnir         | ossl   | na  |

``` r
write_csv(modeling.combinations, "../out/modeling_combinations_v1.2.csv")
```

## Loading regression matrix

The regression-matrices were produced in the [OSSL-imports
repository](https://github.com/soilspectroscopy/ossl-imports).

``` r
## Reading OSSL level 1
rm.ossl <- qread(paste0(db.dir, "ossl_all_L1_v1.2.qs"))

rm.ossl <- rm.ossl %>%
  dplyr::select(id.layer_uuid_txt, dataset.code_ascii_txt,
                any_of(soil.properties.original), starts_with("scan_visnir"), starts_with("scan_mir"))

# head(names(rm.ossl), 20)
# tail(names(rm.ossl), 20)

## Reading Neospectra datasets
neospectra.soilsite <- qread(paste0(db.dir, "neospectra_soilsite_v1.2.qs"))
neospectra.soillab <- qread(paste0(db.dir, "neospectra_soillab_v1.2.qs"))
neospectra.nir <- qread(paste0(db.dir, "neospectra_nir_v1.2.qs"))

# Averaging spectra collected by multiple instruments
neospectra.nir <- neospectra.nir %>%
  dplyr::select(id.sample_local_c, starts_with("scan_nir")) %>%
  group_by(id.sample_local_c) %>%
  summarise_all(mean)

rm.neospectra <- left_join(neospectra.soilsite, neospectra.soillab, by = "id.sample_local_c") %>%
  left_join(neospectra.nir, by = "id.sample_local_c")

# Selecting only important columns
rm.neospectra <- rm.neospectra %>%
  dplyr::select(id.sample_local_c, dataset.code_ascii_txt,
                any_of(soil.properties.original), starts_with("scan_nir"))

# Preparing named list of datasets
# Selecting only important columns, spectra range and rows with available spectra
data.list <- list(
  "mir_mlr..eml_kssl" = {rm.ossl %>%
      dplyr::filter(dataset.code_ascii_txt == "KSSL.SSL") %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.original), starts_with("scan_mir")) %>%
      dplyr::filter(!is.na(scan_mir.1500_abs))},
  "mir_mlr..eml_ossl" = {rm.ossl %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.original), starts_with("scan_mir")) %>%
      dplyr::filter(!is.na(scan_mir.1500_abs))},
  "visnir_mlr..eml_kssl" = {rm.ossl %>%
      dplyr::filter(dataset.code_ascii_txt == "KSSL.SSL") %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.original), starts_with("scan_visnir")) %>%
      dplyr::filter(!is.na(scan_visnir.1500_ref))},
  "visnir_mlr..eml_ossl" = {rm.ossl %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.original), starts_with("scan_visnir")) %>%
      dplyr::filter(!is.na(scan_visnir.1500_ref))},
  "nir.neospectra_mlr..eml_ossl" = {rm.neospectra %>%
      mutate(id.layer_uuid_txt = id.sample_local_c) %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.original), starts_with("scan_nir")) %>%
      dplyr::filter(!is.na(scan_nir.1500_ref))})

# Different sizes
lapply(data.list, function(x) dim(x))
```

### Preprocessing

To remove baseline offset, additive/multiplicative scatter, and
multicollinearity in spectra from multiple sources,
[SNV](https://cran.r-project.org/web/packages/prospectr/vignettes/prospectr.html#scatter-and-baseline-corrections)
preprocessing is used before PCA compression and model calibration:

``` r
prep.list <- lapply(data.list, function(x){
  x %>%
    dplyr::select(-id.layer_uuid_txt, -any_of(soil.properties.original)) %>%
    as.matrix() %>%
    prospectr::standardNormalVariate(X = .) %>%
    as_tibble() %>%
    bind_cols({x %>%
      dplyr::select(id.layer_uuid_txt, any_of(soil.properties.original))}, .)
  })

# lapply(data.list, function(x) dim(x))
```

### PCA compression

Compress and save
[PCA](http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/)
models for the combinations of spectra. This is used to reduce the
number of dimensions and also control multicollinearity ([Chang, Laird,
Mausbach, & Hurburgh Jr, 2001](#ref-chang2001near)).

Fit PCA models in parallel:

``` r
registerDoMC(length(prep.list))

pca.list <- foreach(i=1:length(prep.list)) %dopar% {
  x <- prep.list[[i]] %>%
    dplyr::select(starts_with("scan")) %>%
    as.data.frame()
  prcomp(x)
}

try(stopCluster())

gc()

# names(pca.list) = names(prep.list)
# lapply(pca.list, function(x) dim(x))
# lapply(pca.list, function(x) dim(x))
# str(pca.lst[[1]]$x[,1:20])
```

<!-- Save as RDS files: -->
<!-- ```{r, eval=FALSE} -->
<!-- ## save complete objects: -->
<!-- for(i in 1:length(pca.lst)){  saveRDS.gz(pca.lst[[i]], paste0("pca.ossl/", pca.names[i]))  } -->
<!-- pcam.names = paste0("m", c("pca_mir_kssl_v1.rds", "pca_visnir_kssl_v1.rds", "pca_mir_ossl_v1.rds", "pca_visnir_ossl_v1.rds")) -->
<!-- ## save models only / for prediction purposes: -->
<!-- for(i in 1:length(pca.lst)){  saveRDS.gz(pca.lst[[i]][1:4], paste0("pca.ossl/", pcam.names[i]))  } -->
<!-- ``` -->
<!-- ## Model fitting -->
<!-- For soil-calibration model fitting we use first 120 PCs. This method was first time introduced by @chang2001near,  -->
<!-- although at this stage the authors used only top 10 PCA components, here we use 120. -->
<!-- We focus on target soil variables for which there is enough training data. -->
<!-- ```{r} -->
<!-- n.spc = 120 -->
<!-- t.vars = c("log..oc_usda.calc_wpct", "log..n.tot_usda.4h2_wpct", "silt.tot_usda.3a1_wpct", -->
<!--            "clay.tot_usda.3a1_wpct", "sand.tot_usda.3a1_wpct", "log..ecec_usda.4b4_cmolkg", -->
<!--            "ph.h2o_usda.4c1_index", "ph.cacl2_usda.4c1_index", "log..al.kcl_usda.4b3_cmolkg", -->
<!--            "log..k.ext_usda.4b1_cmolkg", "log..caco3_usda.4e1_wpct", "log..mg.ext_usda.4b1_cmolkg", -->
<!--            "log..ca.ext_usda.4b1_cmolkg", "log..gyp_usda.4e2_wpct", -->
<!--            "log..cec.ext_usda.4b1_cmolkg", "bsat_usda.4b4_wpct", "bd.od_usda.3b2_gcm3") -->
<!-- which(!gsub("log..", "", t.vars) %in% names(rm.ossl)) -->
<!-- x <- lapply(as.vector(unlist(paste0("./models/", t.vars, "/"))), dir.create, recursive=TRUE, showWarnings=FALSE) -->
<!-- ``` -->
<!-- In addition to PC's we will also use a number of additional covariates: -->
<!-- ```{r} -->
<!-- geo.sel = names(rm.ossl)[grep("clm_", names(rm.ossl))] -->
<!-- str(geo.sel) -->
<!-- harm.sel = names(rm.ossl)[grep("dataset.code_ascii_c_", names(rm.ossl))] -->
<!-- ``` -->
<!-- For each soil property, we basically fit the following six type of models: -->
<!-- 1. MIR only KSSL subset `mir_mlr..eml_kssl_na_v1.rds`; -->
<!-- 2. MIR only KSSL subset with geographic covariates `mir_mlr..eml_kssl_ll_v1.rds`; -->
<!-- 3. VisNIR only ALL data `visnir_mlr..eml_ossl_na_v1.rds`; -->
<!-- 4. VisNIR only ALL data with geographic covariates `visnir_mlr..eml_ossl_ll_v1.rds`; -->
<!-- 5. VisNIR and MIR fusion `visnir.mir_mlr..eml_ossl_na_v1.rds`; -->
<!-- 6. VisNIR and MIR fusion with geographic covariates `visnir.mir_mlr..eml_ossl_ll_v1.rds`; -->
<!-- ```{r} -->
<!-- mn.lst = c("mir_mlr..eml_kssl_na_v1", "mir_mlr..eml_kssl_ll_v1", "visnir_mlr..eml_ossl_na_v1", -->
<!--            "mir_mlr..eml_ossl_ll_v1", "visnir.mir_mlr..eml_ossl_na_v1", "visnir.mir_mlr..eml_ossl_ll_v1") -->
<!-- ``` -->
<!-- Using geographic covariates implies that lon-lat coordinates are available at relatively high -->
<!-- location accuracy and not higher location error than 1000 m. For each of the models -->
<!-- we specify a list of covariates: -->
<!-- ```{r} -->
<!-- ## Target 6 models per variable: -->
<!-- pr.lst0 = list(c(paste0("mir.PC", 1:n.spc)), -->
<!--                c(paste0("mir.PC", 1:n.spc), "hzn_depth", geo.sel), -->
<!--                c(paste0("visnir.PC", 1:n.spc)), -->
<!--                c(paste0("mir.PC", 1:n.spc), "hzn_depth", geo.sel, harm.sel), -->
<!--                c(paste0("mir.PC", 1:n.spc), paste0("visnir.PC", 1:n.spc), harm.sel), -->
<!--                c(paste0("mir.PC", 1:n.spc), paste0("visnir.PC", 1:n.spc), "hzn_depth", geo.sel, harm.sel)) -->
<!-- ``` -->
<!-- Training takes at least 24hrs on 80t server with 450GB RAM: -->
<!-- ```{r, eval=FALSE} -->
<!-- ## Run train ---- -->
<!-- max.loc.accuracy=1000 -->
<!-- ## run in parallel for 6 combinations -->
<!-- source('functions_ossl.models.R') -->
<!-- for(tv in t.vars){ -->
<!--   for(k in 1:length(pr.lst0)){ -->
<!--     if(k==1|k==2){ -->
<!--       X.pc = as.data.frame(pca.lst[[1]]$x[,1:n.spc]) -->
<!--       colnames(X.pc) = paste0("mir.PC", 1:n.spc) -->
<!--     } -->
<!--     if(k==3){ -->
<!--       X.pc = as.data.frame(pca.lst[[2]]$x[,1:n.spc]) -->
<!--       colnames(X.pc) = paste0("visnir.PC", 1:n.spc) -->
<!--     } -->
<!--     if(k==4){ -->
<!--       X.pc = as.data.frame(pca.lst[[3]]$x[,1:n.spc]) -->
<!--       colnames(X.pc) = paste0("mir.PC", 1:n.spc) -->
<!--     } -->
<!--     if(k==5|k==6){ -->
<!--       X1.pc = as.data.frame(pca.lst[[3]]$x[,1:n.spc]) -->
<!--       colnames(X1.pc) = paste0("mir.PC", 1:n.spc) -->
<!--       X2.pc = as.data.frame(pca.lst[[4]]$x[,1:n.spc]) -->
<!--       colnames(X2.pc) = paste0("visnir.PC", 1:n.spc) -->
<!--       ov.r = intersect(row.names(X1.pc), row.names(X2.pc)) -->
<!--       X.pc = cbind(X1.pc[which(row.names(X1.pc) %in% ov.r),], X2.pc[which(row.names(X2.pc) %in% ov.r),]) -->
<!--     } -->
<!--     if(nrow(X.pc)>0){ -->
<!--       if(k==1|k==3){ -->
<!--         X = cbind(rm.ossl[as.integer(row.names(X.pc)), c(gsub("log..", "", tv), "ID")], X.pc) -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = FALSE) ) -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = FALSE, rf.feature=FALSE) ) -->
<!--       } -->
<!--       if(k==2){ -->
<!--         X = cbind(rm.ossl[as.integer(row.names(X.pc)), c(gsub("log..", "", tv), "location.error_any_m", "ID", "hzn_depth", geo.sel)], X.pc) -->
<!--         X = X[!is.na(X$location.error_any_m) & X$location.error_any_m < max.loc.accuracy,] -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = FALSE) ) -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = FALSE, rf.feature=FALSE) ) -->
<!--       } -->
<!--       if(k==4|k==6){ -->
<!--         X = cbind(rm.ossl[as.integer(row.names(X.pc)), c(gsub("log..", "", tv), "location.error_any_m", "ID", "hzn_depth", geo.sel, harm.sel)], X.pc) -->
<!--         X = X[!is.na(X$location.error_any_m) & X$location.error_any_m < max.loc.accuracy,] -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = TRUE) ) -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = TRUE, rf.feature=FALSE) ) -->
<!--       } -->
<!--       if(k==5){ -->
<!--         X = cbind(rm.ossl[as.integer(row.names(X.pc)), c(gsub("log..", "", tv), "ID", harm.sel)], X.pc) -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = FALSE) ) -->
<!--         try( i <- train.ossl(tv, pr.var=pr.lst0[[k]], X, model.name=mn.lst[k], hzn_depth = FALSE, rf.feature=FALSE) ) -->
<!--       } -->
<!--     } -->
<!--     try( cat_eml(tv, model.name=mn.lst[k]) ) -->
<!--   } -->
<!-- } -->
<!-- ``` -->
<!-- ## Model evaluation -->
<!-- We can export the summary accuracy statistics for all soil variables into a single table: -->
<!-- ```{r, eval=FALSE} -->
<!-- p.lst = c("RMSE", "R.square", "N.tot", "N.outliers") -->
<!-- acc.mat = data.frame(matrix(nrow=length(t.vars), ncol=2+length(mn.lst)*4)) -->
<!-- colnames(acc.mat) = c("variable", "std", sapply(mn.lst, function(i){paste0(i, "_", p.lst)})) -->
<!-- acc.mat$variable = t.vars -->
<!-- acc.mat$std = sapply(t.vars, function(i){ if(length(grep("log..", i))>0) { sd(log1p(rm.ossl[,gsub("log..","",i)]), na.rm = TRUE) } else { sd(rm.ossl[,i], na.rm = TRUE) } }) -->
<!-- for(i in 1:nrow(acc.mat)){ -->
<!--   for(j in 1:length(mn.lst)){ -->
<!--     in.rds = paste0(dir, t.vars[i], "/", mn.lst[j], ".rds") -->
<!--     if(file.exists(in.rds)){ -->
<!--       t.m = readRDS.gz(in.rds) -->
<!--       x.s = summary(t.m$learner.model$super.model$learner.model) -->
<!--       RMSE = signif(sqrt(sum(t.m$learner.model$super.model$learner.model$residuals^2) / t.m$learner.model$super.model$learner.model$df.residual), 3) -->
<!--       acc.mat[i,paste0(mn.lst[j], "_RMSE")] = RMSE -->
<!--       acc.mat[i,paste0(mn.lst[j], "_R.square")] = round(x.s$adj.r.squared, 3) -->
<!--       acc.mat[i,paste0(mn.lst[j], "_N.tot")] = x.s$df[2] -->
<!--       acc.mat[i,paste0(mn.lst[j], "_N.outliers")] = sum(abs(x.s$residuals) > 3*RMSE) -->
<!--       gc() -->
<!--     } -->
<!--   } -->
<!-- } -->
<!-- write.csv(acc.mat, "../out/accuracy_matrix_ossl_models.csv") -->
<!-- ``` -->
<!-- We also produce standard accuracy plots based on the 5-fold cross-validation by using: -->
<!-- ```{r, eval=FALSE} -->
<!-- ## Accuracy plots ---- -->
<!-- library(hexbin) -->
<!-- library(plotKML) -->
<!-- library(lattice) -->
<!-- par(family = "sans") -->
<!-- for(i in 1:nrow(acc.mat)){ -->
<!--   for(j in 1:length(mn.lst)){ -->
<!--     in.rds = paste0(dir, t.vars[i], "/", mn.lst[j], ".rds") -->
<!--     out.file = paste0(dir, t.vars[i], "/ap.", mn.lst[j], ".rds.png") -->
<!--     if(file.exists(in.rds) & !file.exists(out.file)){ -->
<!--       t.m = readRDS.gz(in.rds) -->
<!--       yh = t.m$learner.model$super.model$learner.model$fitted.values -->
<!--       meas = t.m$learner.model$super.model$learner.model$model[,t.vars[i]] -->
<!--       t.var.breaks = quantile(meas, c(0.001, 0.01, 0.999), na.rm=TRUE) -->
<!--       plot_hexbin(varn=t.vars[i], breaks=c(t.var.breaks[1], seq(t.var.breaks[2], t.var.breaks[3], length=25)), meas=ifelse(meas<0, 0, meas), pred=ifelse(yh<0, 0, yh), main=t.vars[i], out.file=out.file, log.plot=FALSE, colorcut=c(0,0.01,0.02,0.03,0.06,0.12,0.20,0.35,1.0)) -->
<!--       #gc() -->
<!--     } -->
<!--   } -->
<!-- } -->
<!-- ``` -->
<!-- Example of an accuracy plot: -->
<!-- ```{r ac-soc1, echo=FALSE, fig.cap="Accuracy plot for `log..oc_usda.calc_wpct/mir_mlr..eml_kssl_na_v1.rds`.", out.width="60%"} -->
<!-- knitr::include_graphics("./models/log..oc_usda.calc_wpct/ap.mir_mlr..eml_kssl_na_v1.rds.png") -->
<!-- ``` -->
<!-- ```{r, eval=FALSE} -->
<!-- #save.image.pigz(file=paste0(dir, "ossl.models.RData"), n.cores=80) -->
<!-- #rmarkdown::render("R-mlr/README.Rmd") -->
<!-- ``` -->

## References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-bischl2016mlr" class="csl-entry">

Bischl, B., Lang, M., Kotthoff, L., Schiffner, J., Richter, J.,
Studerus, E., … Jones, Z. M. (2016). <span class="nocase">mlr: Machine
Learning in R</span>. *The Journal of Machine Learning Research*,
*17*(1), 5938–5942. Retrieved from
<https://www.jmlr.org/papers/volume17/15-066/15-066.pdf>

</div>

<div id="ref-chang2001near" class="csl-entry">

Chang, C.-W., Laird, D., Mausbach, M. J., & Hurburgh Jr, C. R. (2001).
Near-infrared reflectance spectroscopy–principal components regression
analyses of soil properties. *Soil Science Society of America Journal*,
*65*(2), 480.
doi:[10.2136/sssaj2001.652480x](https://doi.org/10.2136/sssaj2001.652480x)

</div>

<div id="ref-hengl2019predictive" class="csl-entry">

Hengl, T., & MacMillan, R. A. (2019). *<span class="nocase">Predictive
Soil Mapping with R</span>* (p. 370). Wageningen, the Netherlands: Lulu.
Retrieved from <https://soilmapper.org>

</div>

</div>
