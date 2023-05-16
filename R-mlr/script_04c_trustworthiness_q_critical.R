
## Packages
library("tidyverse")
library("mlr3")
library("qs")

options(scipen = 999)

mnt.dir <- "~/mnt-ossl/ossl_models/"

pca.list <- list.files(paste0(mnt.dir, "pca.ossl/"), pattern = "mpca_")
pca.list <- gsub("mpca_|\\.qs", "", pca.list)
pca.list

q.list <- list()

i=1
for(i in 1:length(pca.list)) {

  iname <- pca.list[i]

  ## PCA files
  pca.model <- qread(paste0(mnt.dir, "pca.ossl/mpca_", iname, ".qs"))
  pca.scores <- qread(paste0(mnt.dir, "pca.ossl/pca_scores_", iname, ".qs"))

  ## Spectral ranges
  if(grepl("mir", iname)) {

    spectral.range <- paste0("scan_mir.", seq(600, 4000, by = 2), "_abs")

  } else if(grepl("visnir", iname)) {

    spectral.range <- paste0("scan_visnir.", seq(400, 2500, by = 2), "_ref")

  } else if(grepl("neospectra", iname)) {

    spectral.range <- paste0("scan_nir.", seq(1350, 2550, by = 2), "_ref")

  }

  ncomps <- 120

  ## Q statistics from training set

  # Full PCs
  pca.model.full <- pca.model
  class(pca.model.full) <- "prcomp"

  # Scores
  pca.scores.full <- pca.scores %>%
    select(starts_with("id"), starts_with("PC"))

  # Back-transformed spectra
  pca.bt.full <- as.matrix(pca.scores.full[,-1]) %*% t(pca.model.full$rotation)
  pca.bt.full <- t((t(pca.bt.full) * pca.model.full$scale) + pca.model.full$center)
  pca.bt.full <- bind_cols(pca.scores.full[,1], pca.bt.full)

  # Employed PC (120 PCs)
  pca.model.employed <- pca.model
  class(pca.model.employed) <- "prcomp"
  pca.model.employed$sdev <- pca.model.employed$sdev[1:ncomps]
  pca.model.employed$rotation <- pca.model.employed$rotation[,1:ncomps]

  # Scores
  pca.scores.employed <- pca.scores %>%
    select(starts_with("id"), all_of(paste0("PC", seq(1, ncomps, by = 1))))

  # Back-transformed spectra
  pca.bt.employed <- as.matrix(pca.scores.employed[,-1]) %*% t(pca.model.employed$rotation)
  pca.bt.employed <- t((t(pca.bt.employed) * pca.model.employed$scale) + pca.model.employed$center)
  pca.bt.employed <- bind_cols(pca.scores.employed[,1], pca.bt.employed)

  # Q stats - sum of squared differences between full and employed PCs back-transformed spectra
  q.stats <- apply((pca.bt.employed[,-1]-pca.bt.full[,-1])^2, MARGIN = 1, sum)

  q.stats <- tibble(id.layer_uuid_txt = pca.bt.employed[,1],
                    q_stats = q.stats)

  # Visualization
  # ggplot(q.stats) +
  #   geom_histogram(aes(x = q_stats)) +
  #   theme_light()

  ## Q critical for new prediction
  # 1% significance level

  # https://doi.org/10.1016/j.geoderma.2023.116491
  # https://doi.org/10.1016/j.chemolab.2011.04.002
  # https://sci-hub.se/https://doi.org/10.1080/00401706.1979.10489779

  E <- cov(pca.bt.employed[,-1]-pca.bt.full[,-1])
  teta1 <- sum(diag(E))^1
  teta2 <- sum(diag(E))^2
  teta3 <- sum(diag(E))^3
  h0 <- 1-((2*teta1*teta3)/(3*teta2^2))
  Ca <- 2.57 # 1% significance level
  Qa <- teta1*(1-(teta2*h0*((1-h0)/teta1^2))+((sqrt(Ca*(2*teta2*h0^2)))/teta1))^(1/h0)
  Qa

  # Saving result
  out <- tibble(model_name = iname, q_critical = Qa)
  q.list[[i]] <- out

  cat(paste0("Run ", i, "/", length(pca.list), "\n"))

}

q.critical.table <- Reduce(bind_rows, q.list)
q.critical.table

write_csv(q.critical.table, "out/trustworthiness_q_critical_values.csv")
