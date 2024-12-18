
## Packages
library("tidyverse")
library("yardstick")
library("ggpattern")
library("qs")
library("fs")

## Dirs
dir <- "~/mnt-ossl/ossl_models/"
db.dir <- "~/mnt-ossl/ossl_import/"

## Fitted models
fitted.modeling.combinations <- read_csv("./out/fitted_modeling_combinations_v1.2.csv",
                                         show_col_types = FALSE)

## Custom functions

kge_r <- function(truth,
                  estimate,
                  ...) {
  cor(truth, estimate, method = "pearson")
}

kge_alpha <- function(truth,
                      estimate,
                      ...) {

  var(estimate)/var(truth)

}

kge_beta <-function(truth,
                    estimate,
                    ...) {

  mean(estimate)/mean(truth)

}

kge_vec <- function(truth,
                    estimate,
                    ...) {
  calc_r <- kge_r(truth, estimate)
  calc_alpha <- kge_alpha(truth, estimate)
  calc_beta <- kge_beta(truth, estimate)

  ed = sqrt((calc_r-1)^2 + (calc_alpha-1)^2 + (calc_beta-1)^2)
  1-ed

}

## Full table with final performance

perf.results <- list()

i=1
for(i in 1:nrow(fitted.modeling.combinations)) {

  # Parameters
  isoil_property = fitted.modeling.combinations[[i,"soil_property"]]
  imodel_name = fitted.modeling.combinations[[i,"model_name"]]
  iexport_name = fitted.modeling.combinations[[i,"export_name"]]
  ispectra_type = fitted.modeling.combinations[[i,"spectra_type"]]
  isubset = fitted.modeling.combinations[[i,"subset"]]
  igeo = fitted.modeling.combinations[[i,"geo"]]

  cat(paste0("Run ", i, "/", nrow(fitted.modeling.combinations), " at ", lubridate::now(), "\n"))

  # Predictions
  cv.results <- qread(paste0(dir,
                             iexport_name,
                             "/cvpred_",
                             imodel_name,
                             ".qs"))

  performance.metrics <- cv.results %>%
    summarise(n = n(),
              rmse = rmse_vec(truth = truth, estimate = response),
              bias = msd_vec(truth = truth, estimate = response),
              rsq = rsq_vec(truth = truth, estimate = response),
              nse = rsq_trad_vec(truth = truth, estimate = response),
              kge = kge_vec(truth = truth, estimate = response),
              kge_r = kge_r(truth = truth, estimate = response),
              kge_alpha = kge_alpha(truth = truth, estimate = response),
              kge_beta = kge_beta(truth = truth, estimate = response),
              ccc = ccc_vec(truth = truth, estimate = response, bias = T),
              rpiq = rpiq_vec(truth = truth, estimate = response))

  performance.metrics <- performance.metrics %>%
    mutate(soil_property = isoil_property,
           model_name = imodel_name,
           unit_transform = ifelse(grepl("log..", iexport_name), "log1p", "original"),
           .before = 1)

  perf.results[[i]] <- performance.metrics

}

performance.extra <- Reduce(bind_rows, perf.results) %>%
  mutate_if(is.numeric, round, 3)

write_csv(performance.extra, "out/fitted_models_performance_v1.2_extra.csv")

## Visualizations
library("corrplot")

M <- performance.extra %>%
  as_tibble() %>%
  select(-soil_property, -model_name, -unit_transform) %>%
  cor()

corrplot(M, method = 'number', type = 'lower', diag = FALSE) # colorful number

## Bar plot

performance <- performance.extra %>%
  select(-n, -unit_transform)

performance <- performance %>%
  mutate(info = model_name, .after = model_name) %>%
  mutate(info = gsub("cubist_|_na_v1.2", "", info)) %>%
  separate(info, into = c("spectra_type", "library"), sep = "_")

performance <- performance %>%
  filter(library == "ossl") %>%
  mutate(trust = ifelse(rpiq >= 2, "yes", "no"))

rank <- performance %>%
  filter(spectra_type == "mir") %>%
  arrange(ccc) %>%
  pull(soil_property)

labels <- performance %>%
  filter(spectra_type == "mir") %>%
  arrange(ccc) %>%
  pull(soil_property)

p.rank <- performance %>%
  mutate(soil_property = factor(soil_property, levels = rank)) %>%
  mutate(spectra_type = factor(spectra_type,
                               levels = rev(c("mir", "visnir", "nir.neospectra")))) %>%
  ggplot(aes(x = soil_property, y = nse,
             fill = spectra_type, pattern = trust)) +
  geom_col_pattern(position = position_dodge(preserve = "single"),
                   pattern_color = 'gray70',
                   size = 0.1,
                   pattern_angle = 45,
                   pattern_density = 0.075,
                   pattern_spacing = 0.015,
                   pattern_key_scale_factor = 1,
                   show.legend = T) +
  scale_fill_manual(values = rev(c("darkblue", "salmon", "gold"))) +
  scale_pattern_manual(values = c(yes = "none", no = "stripe")) +
  scale_x_discrete(labels = labels) +
  labs(x = NULL, y = "NSE", fill = NULL) +
  coord_flip() +
  theme_light() +
  theme(legend.position = "bottom") +
  guides(pattern = "none")

ggsave("~/plot_rank_10CV_nse.png", p.rank,
       width = 5, height = 5, units = "in", dpi = 300, scale = 1.5)
