
## Packages
library("tidyverse")
library("mlr3")
library("qs")

options(scipen = 999)

## Neospectra model files
task <- qread("~/mnt-ossl/ossl_models/clay.tot_usda.a334_w.pct/task_nir.neospectra_mlr3..eml_ossl_na_v1.2.qs")
model <- qread("~/mnt-ossl/ossl_models/clay.tot_usda.a334_w.pct/model_nir.neospectra_mlr3..eml_ossl_na_v1.2.qs")
pca.model <- qread("~/mnt-ossl/ossl_models/pca.ossl/mpca_nir.neospectra_mlr3..eml_ossl_v1.2.qs")
pca.scores <- qread("~/mnt-ossl/ossl_models/pca.ossl/pca_scores_nir.neospectra_mlr3..eml_ossl_v1.2.qs")

## Important columns
nir.neospectra.spectral.range <- paste0("scan_nir.", seq(1350, 2550, by = 2), "_ref")
soil.properties.names <- c("clay.tot_usda.a334_w.pct")
ncomps <- 120

# Q statistics from training set
names(pca.model)
lapply(pca.model, class)
lapply(pca.model, length)

pca.model.full <- pca.model

pca.model.employed <- pca.model
pca.model.employed$sdev <- pca.model.employed$sdev[1:ncomps]
pca.model.employed$rotation <- pca.model.employed$rotation[,1:ncomps]

class(pca.model.full) <- "prcomp"
class(pca.model.employed) <- "prcomp"

pca.scores.full <- pca.scores %>%
  select(starts_with("id"), starts_with("PC"))

pca.bt.full <- as.matrix(pca.scores.full[,-1]) %*% t(pca.model.full$rotation)
pca.bt.full <- t((t(pca.bt.full) * pca.model.full$scale) + pca.model.full$center)
pca.bt.full <- bind_cols(pca.scores.full[,1], pca.bt.full)
pca.bt.full

pca.scores.employed <- pca.scores %>%
  select(starts_with("id"), all_of(paste0("PC", seq(1, ncomps, by = 1))))

pca.bt.employed <- as.matrix(pca.scores.employed[,-1]) %*% t(pca.model.employed$rotation)
pca.bt.employed <- t((t(pca.bt.employed) * pca.model.employed$scale) + pca.model.employed$center)
pca.bt.employed <- bind_cols(pca.scores.employed[,1], pca.bt.employed)
pca.bt.employed

q.stats.original <- apply((pca.bt.employed[,-1]-pca.bt.full[,-1])^2, MARGIN = 1, sum)

q.stats.original <- tibble(id.layer_uuid_txt = pca.bt.employed[,1],
                           q_stats = q.stats.original)

ggplot(q.stats.original) +
  geom_histogram(aes(x = q_stats)) +
  theme_light()

## Q critical
# https://doi.org/10.1016/j.geoderma.2023.116491
# https://doi.org/10.1016/j.chemolab.2011.04.002
# https://sci-hub.se/https://doi.org/10.1080/00401706.1979.10489779

E <- cov(pca.bt.employed[,-1]-pca.bt.full[,-1])
teta1 <- sum(diag(E))^1
teta2 <- sum(diag(E))^2
teta3 <- sum(diag(E))^3
h0 <- 1-((2*teta1*teta3)/(3*teta2^2))
Ca <- 2.57
Qa <- teta1*(1-(teta2*h0*((1-h0)/teta1^2))+((sqrt(Ca*(2*teta2*h0^2)))/teta1))^(1/h0)
Qa

## Test set

# Selected ids
neospectra.site <- qread("~/mnt-ossl/ossl_import/neospectra_soilsite_v1.2.qs")

selected.ids <- neospectra.site %>%
  filter(!(location.country_iso.3166_txt == "USA")) %>%
  pull(id.sample_local_c)

neospectra.site <- neospectra.site %>%
  filter(id.sample_local_c %in% selected.ids)

# NIR
neospectra.nir <- qread("~/mnt-ossl/ossl_import/neospectra_nir_v1.2.qs")

neospectra.nir <- neospectra.nir %>%
  dplyr::select(id.sample_local_c, all_of(nir.neospectra.spectral.range)) %>%
  group_by(id.sample_local_c) %>%
  summarise_all(mean)

neospectra.nir <- neospectra.nir %>%
  filter(id.sample_local_c %in% selected.ids)

set.seed(1993)
neospectra.nir.outliers <- neospectra.nir %>%
  sample_n(10)

new.ids <- paste("OUT", seq(1, 10, 1))
set.seed(1993)
error <- rnorm(601, mean = 0, sd = 0.05)
error.table <- tibble(error = error, name = nir.neospectra.spectral.range)

neospectra.nir.outliers <- neospectra.nir.outliers %>%
  mutate(id.sample_local_c = new.ids) %>%
  pivot_longer(-id.sample_local_c) %>%
  left_join(error.table) %>%
  mutate(value = value+error) %>%
  select(-error) %>%
  pivot_wider()

neospectra.nir <- bind_rows(neospectra.nir.outliers, neospectra.nir)

## Preprocessing

neospectra.nir.prep <- neospectra.nir %>%
  dplyr::select(-id.sample_local_c, -any_of(soil.properties.names)) %>%
  as.matrix() %>%
  prospectr::standardNormalVariate(X = .) %>%
  as_tibble() %>%
  bind_cols({neospectra.nir %>%
      dplyr::select(id.sample_local_c, any_of(soil.properties.names))}, .)

## PCA compression
neospectra.nir.pca.full <- predict(pca.model.full, neospectra.nir.prep[,-1]) %>%
  bind_cols(neospectra.nir.prep[,1], .)

neospectra.nir.pca.employed <- predict(pca.model.employed, neospectra.nir.prep[,-1]) %>%
  bind_cols(neospectra.nir.prep[1], .)

# Visualization
# Only neospectra.nir.pca.full as both employed and full have the same response
# for the first 120 PCs

neospectra.nir.pca <- neospectra.nir.pca.full %>%
  mutate(flag = ifelse(grepl("OUT", id.sample_local_c), "outlier", "normal"))

ggplot(data = pca.scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_point(data = neospectra.nir.pca, aes(x = PC1, y = PC2, color = flag), size = 1.5, alpha = 1) +
  scale_color_manual(values = c("gold", "red")) +
  labs(title = "Neospectra calibration (black), test from Africa (gold), synthetic outliers (red)") +
  theme_light() + theme(legend.position = "bottom")

# Reconstruction
neospectra.nir.bt.full <- as.matrix(neospectra.nir.pca.full[,-1]) %*% t(pca.model.full$rotation)
neospectra.nir.bt.full <- t((t(neospectra.nir.bt.full) * pca.model.full$scale) + pca.model.full$center)
neospectra.nir.bt.full <- bind_cols(neospectra.nir.pca.full[,1], neospectra.nir.bt.full)
neospectra.nir.bt.full

neospectra.nir %>%
  pivot_longer(-id.sample_local_c, names_to = "wavelength", values_to = "reflectance") %>%
  filter(!grepl("OUT", id.sample_local_c)) %>%
  mutate(wavelength = as.numeric(gsub("scan_nir.|_ref", "", wavelength))) %>%
  ggplot() +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            linewidth = 0.5, alpha = 0.5) +
  labs(title = "Original test spectra") +
  theme_light()

neospectra.nir.bt.full %>%
  pivot_longer(-id.sample_local_c, names_to = "wavelength", values_to = "reflectance") %>%
  filter(!grepl("OUT", id.sample_local_c)) %>%
  mutate(wavelength = as.numeric(gsub("scan_nir.|_ref", "", wavelength))) %>%
  ggplot() +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            linewidth = 0.5, alpha = 0.5) +
  labs(title = "SNV test spectra") +
  theme_light()

neospectra.nir.bt.full %>%
  pivot_longer(-id.sample_local_c, names_to = "wavelength", values_to = "reflectance") %>%
  filter(grepl("OUT", id.sample_local_c)) %>%
  mutate(wavelength = as.numeric(gsub("scan_nir.|_ref", "", wavelength))) %>%
  ggplot() +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            linewidth = 0.5, alpha = 0.5) +
  labs(title = "SNV outlier spectra") +
  theme_light()

neospectra.nir.bt.employed <- as.matrix(neospectra.nir.pca.employed[,-1]) %*% t(pca.model.employed$rotation)
neospectra.nir.bt.employed <- t((t(neospectra.nir.bt.employed) * pca.model.employed$scale) + pca.model.employed$center)
neospectra.nir.bt.employed <- bind_cols(neospectra.nir.pca.employed[,1], neospectra.nir.bt.employed)
neospectra.nir.bt.employed

neospectra.nir.bt.employed %>%
  pivot_longer(-id.sample_local_c, names_to = "wavelength", values_to = "reflectance") %>%
  filter(!grepl("OUT", id.sample_local_c)) %>%
  mutate(wavelength = as.numeric(gsub("scan_nir.|_ref", "", wavelength))) %>%
  ggplot() +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            linewidth = 0.5, alpha = 0.5) +
  labs(title = "Back-transformed SNV test spectra") +
  theme_light()

neospectra.nir.bt.employed %>%
  pivot_longer(-id.sample_local_c, names_to = "wavelength", values_to = "reflectance") %>%
  filter(grepl("OUT", id.sample_local_c)) %>%
  mutate(wavelength = as.numeric(gsub("scan_nir.|_ref", "", wavelength))) %>%
  ggplot() +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            linewidth = 0.5, alpha = 0.5) +
  labs(title = "Back-transformed SNV outlier spectra") +
  theme_light()

## Q statistics of selected neospectra
## Sum of squared differences between reconstructed full and reconstructed employed in the models
## Highlights if important features are being missed from the models
## The predicted Q stats will be compared to the Q stats distribution from the training set

q.stats <- apply((neospectra.nir.bt.employed[,-1]-neospectra.nir.bt.full[,-1])^2, MARGIN = 1, sum)

q.stats <- tibble(neospectra.nir.bt.employed[,1],
                  q_stats = q.stats,
                  q_critical_99perc = Qa) %>%
  mutate(flag = q_stats >= q_critical_99perc)

q.stats %>%
  filter(!grepl("OUT", id.sample_local_c)) %>%
  summarise(q_stats_mean = mean(q_stats),
            q_stats_min = min(q_stats),
            q_stats_max = max(q_stats))

q.stats %>%
  filter(grepl("OUT", id.sample_local_c)) %>%
  summarise(q_stats_mean = mean(q_stats),
            q_stats_min = min(q_stats),
            q_stats_max = max(q_stats))

# Visualization
ggplot(q.stats) +
  geom_histogram(aes(x = log10(q_stats))) +
  geom_vline(aes(xintercept = Qa)) +
  labs(title = "Q-stats distribution for test set") +
  theme_light()

