
## Packages
library("tidyverse")
library("lubridate")
library("mlr3verse")
library("yardstick")
library("qs")

## Folders
mnt.dir <- "~/projects/temp/ossl_models/comparison/"

## Importing data
modeling.combinations <- read_csv(file.path(mnt.dir, "modeling_combinations.csv"))

## Evaluation pipeline
lgr::get_logger("mlr3")$set_threshold("warn")
future::plan("multisession")

performance.list <- list()

i=1
for(i in 1:nrow(modeling.combinations)) {

  # Parameters
  imodel_type <- modeling.combinations[[i,"model_type"]]
  ibase_predictions <- modeling.combinations[[i,"base_predictions"]]
  isoil_property <- modeling.combinations[[i,"soil_property"]]
  ispectra <- modeling.combinations[[i,"spectra"]]
  igeocovariates <- modeling.combinations[[i,"geocovariates"]]
  ipca_compression <- modeling.combinations[[i,"pca_compression"]]

  # Task
  sel.data <- qread(paste0(mnt.dir,
                           "autotune/",
                           paste("task",
                                 imodel_type,
                                 ibase_predictions,
                                 isoil_property,
                                 ispectra,
                                 igeocovariates,
                                 ipca_compression,
                                 sep = "_"),
                           ".qs"))

  # Create regression task
  task <- as_task_regr(sel.data, id = "train", target = isoil_property, type = "regression")

  # Defining id column
  task$set_col_roles("id.layer_uuid_txt", roles = "name")

  # For block cv. If 'id.tile' not in the data.frame, default to random CV
  if("id.tile" %in% colnames(sel.data)) {
    task$set_col_roles("id.tile", roles = "group")
  }

  # Autotuned model
  at <- qread(paste0(mnt.dir,
                     "autotune/",
                     paste("model",
                           imodel_type,
                           ibase_predictions,
                           isoil_property,
                           ispectra,
                           igeocovariates,
                           ipca_compression,
                           sep = "_"),
                     ".qs"))

  # Outer 10-CV evaluation
  tuned_learner <- at$learner

  set.seed(1993)
  rr = mlr3::resample(task = task,
                      learner = tuned_learner,
                      resampling = rsmp("cv", folds = 10))

  cv.results <- lapply(rr$predictions("test"), function(x) as.data.table(x))
  cv.results <- Reduce(rbind, cv.results)

  # Metrics
  performance.metrics <- cv.results %>%
    summarise(n = n(),
              rmse = rmse_vec(truth = truth, estimate = response),
              bias = msd_vec(truth = truth, estimate = response),
              rsq = rsq_vec(truth = truth, estimate = response),
              ccc = ccc_vec(truth = truth, estimate = response, bias = T),
              rpiq = rpiq_vec(truth = truth, estimate = response))

  performance.list[[i]] <- performance.metrics

  perfomance.annotation <- paste0("Lin's CCC = ", round(performance.metrics[[1,"ccc"]], 3),
                                  "\nRMSE = ", round(performance.metrics[[1,"rmse"]], 3))


  # Plot
  p.hex <- ggplot(cv.results, aes(x = truth, y = response)) +
    geom_hex(bins = 30, alpha = 0.75) +
    geom_abline(intercept = 0, slope = 1) +
    geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2),
              label = perfomance.annotation) +
    scale_fill_viridis_c(trans = "log10") +
    labs(x = "log(observed)", y = "log(predicted)", fill = bquote(log[10](count)),
         title = isoil_property) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key.size = unit(0.35, "cm"))

  r.max <- max(layer_scales(p.hex)$x$range$range)
  r.min <- min(layer_scales(p.hex)$x$range$range)

  s.max <-max(layer_scales(p.hex)$y$range$range)
  s.min <-min(layer_scales(p.hex)$y$range$range)

  t.max <-round(max(r.max,s.max),1)
  t.min <-round(min(r.min,s.min),1)

  p.hex <- p.hex + coord_equal(xlim=c(t.min,t.max),ylim=c(t.min,t.max))

  ggsave(paste0(mnt.dir,
                "evaluation/",
                paste("plot",
                      imodel_type,
                      ibase_predictions,
                      isoil_property,
                      ispectra,
                      igeocovariates,
                      ipca_compression,
                      sep = "_"),
                ".png"),
         p.hex, width = 6, height = 5, units = "in", scale = 1)

  cat(paste0("Run ", i, "/", nrow(modeling.combinations), "\n"))

}

final.performance <- Reduce(bind_rows, performance.list)
final.performance

final.performance <- bind_cols(modeling.combinations, final.performance)
final.performance

# Arrange by lowest RMSE
final.performance <- final.performance %>%
  select(soil_property, spectra, everything()) %>%
  arrange(soil_property, spectra, rmse)

write_csv(final.performance, "out/comparison_modeling_configurations.csv")

# All combinations plot
final.performance <- read_csv("out/comparison_modeling_configurations.csv")

p.side <- ggplot(final.performance, aes(x = interaction(base_predictions, pca_compression, geocovariates),
                              y = ccc,
                              color = model_type)) +
  geom_boxplot() +
  facet_grid(soil_property~spectra, scales = "free_y") +
  labs(x = NULL, y = "Lin's CCC", color = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45)); p.side

ggsave(paste0(mnt.dir, "plot_all_configurations.png"), p.side, width = 6, height = 5, units = "in", scale = 1)

# Best combinations

# Model type
final.performance %>%
  group_by(model_type) %>%
  summarise(ccc = median(ccc))

ggplot(final.performance, aes(x = model_type, y = ccc, color = soil_property)) +
  geom_boxplot() + labs(color = NULL) +
  theme_light() + theme(legend.position = "bottom")

# Base predictions
final.performance %>%
  group_by(base_predictions) %>%
  summarise(ccc = median(ccc))

ggplot(final.performance, aes(x = base_predictions, y = ccc, color = soil_property)) +
  geom_boxplot() + labs(color = NULL) +
  theme_light() + theme(legend.position = "bottom")

# Spectra
final.performance %>%
  group_by(spectra) %>%
  summarise(ccc = median(ccc))

ggplot(final.performance, aes(x = spectra, y = ccc, color = soil_property)) +
  geom_boxplot() + labs(color = NULL) +
  theme_light() + theme(legend.position = "bottom")

# Geocovariates
final.performance %>%
  group_by(geocovariates) %>%
  summarise(ccc = median(ccc))

ggplot(final.performance, aes(x = geocovariates, y = ccc, color = soil_property)) +
  geom_boxplot() + labs(color = NULL) +
  theme_light() + theme(legend.position = "bottom")

# PCA compression
final.performance %>%
  group_by(pca_compression) %>%
  summarise(ccc = median(ccc))

ggplot(final.performance, aes(x = pca_compression, y = ccc, color = soil_property)) +
  geom_boxplot() + labs(color = NULL) +
  theme_light() + theme(legend.position = "bottom")

# Best combination
final.performance %>%
  filter(geocovariates == "na") %>%
  filter(pca_compression == "n120") %>%
  filter(base_predictions == "cv5") %>%
  filter(model_type == "ensemble")

p.best <- final.performance %>%
  filter(geocovariates == "na") %>%
  filter(pca_compression == "n120") %>%
  filter(!(model_type == "ensemble" & base_predictions == "insample")) %>%
  ggplot(aes(x = model_type,
             y = ccc,
             color = model_type)) +
  geom_boxplot() +
  facet_grid(soil_property~spectra, scales = "free_y") +
  labs(x = NULL, y = "Lin's CCC", color = NULL,
       title = paste0("Best configuration: model ensemble, without geocovariates,\n",
                      "CV predictions as input, and PCA compression with 120 comps")) +
  theme_light() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
        plot.title = element_text(size = 10)); p.best

ggsave(paste0(mnt.dir, "plot_best_configuration.png"), p.best, width = 6, height = 5, units = "in", scale = 1)
