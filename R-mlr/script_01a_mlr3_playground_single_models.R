
## Packages
library("tidyverse")
library("yardstick")
library("mlr3verse")
library("qs")

## Folders
dir <- "/mnt/soilspec4gg/ossl/ossl_models/"
db.dir <- "/mnt/soilspec4gg/ossl/ossl_import/"

## Modeling combinations
modeling.combinations <- read_csv("./out/modeling_combinations_v1.2.csv")
modeling.combinations

## Parameters and regression matrix
i=1
isoil_property = modeling.combinations[[i,"soil_property"]]
imodel_name = modeling.combinations[[i,"model_name"]]
igeo = modeling.combinations[[i,"geo"]]
imodel_name.pca <- str_replace(imodel_name, paste0("_", igeo), "")
ipca <- qread(paste0(dir, "pca.ossl/pca_scores_", imodel_name.pca, ".qs"))

## Modeling pipeline
ncomps = 120
block = FALSE
block_var = NULL

# OTher basic parameters
predictors <- paste0("PC", seq(1, ncomps))
regression.matrix <- ipca

# Regression matrix
regression.matrix <- regression.matrix %>%
  select(id.layer_uuid_txt, all_of(block_var), all_of(isoil_property), all_of(predictors)) %>%
  filter(!is.na(!!as.name(isoil_property))) %>%
  as.data.frame()

# Modeling pipeline
set.seed(1993)

# Create regression task
task <- as_task_regr(regression.matrix, id = "train", target = isoil_property, type = "regression")

# Defining id column
task$set_col_roles("id.layer_uuid_txt", roles = "name")

# For block cv. If not in the data.frame, default to random CV
if(block == TRUE) {
  task$set_col_roles(block_var, roles = "group")
}

# Define learners

# Using some defaults
# Tune: s (lambda), and alpha (0 ridge, 1 lasso)
lrn_glmnet <- lrn("regr.glmnet",
                  s = to_tune(0.0001, 1),
                  alpha = to_tune(0, 1))

# # Using some defaults, num.trees=100, replace=TRUE
# # Tune mtry, min.bucket, min.node.size
# lrn_ranger <- lrn("regr.ranger", num.trees = 100, replace=TRUE, num.threads = 1L,
#                   seed = 1993,
#                   mtry = to_tune(ceiling(sqrt(ncomps)), ceiling(ncomps/2)),
#                   max.depth = to_tune(5, 50),
#                   min.node.size = to_tune(5, 50))

# # Using some defaults, booster = "gbtree", subsample = 0.67, colsample_bytree=0.67
# # Tune eta (learning rate), max_depth, min_child_weight
# lrn_xgboost <- lrn("regr.xgboost", booster = "gbtree",
#                    subsample = 0.67,
#                    colsample_bytree = 0.67,
#                    nthread = 1L,
#                    eta = to_tune(0.3, 0.5),
#                    max_depth = to_tune(3, 10),
#                    min_child_weight = to_tune(1, 10))

# # Using some defaults, neighbors = 0
# # Tune committees
# lrn_cubist <- lrn("regr.cubist",
#                   committees = to_tune(5, 20),
#                   neighbors = 0) # Tune not work for this HP

# Two basic HPO options

# Grid search
# Resolution 5 means that within the range, 5 equidistant values will be tested
# If we have 3 HP to be tuned, a grid of 5^3 will be searched: 125 combinations
at_glmnet = auto_tuner(tuner = tnr("grid_search", resolution = 5, batch_size = 32), # batch_size ~ ncores
                       learner = lrn_glmnet,
                       resampling = rsmp("cv", folds = 5),
                       measure = msr("regr.rmse"),
                       terminator = trm("none"),
                       store_models = TRUE)

# Random search
# Random values are picked within the range of HP and combined together
# n_evals means that 20 different combinations will be tested
at_glmnet = auto_tuner(tuner = tnr("random_search"),
                       learner = lrn_glmnet,
                       resampling = rsmp("cv", folds = 5),
                       measure = msr("regr.rmse"),
                       terminator = trm("evals", n_evals = 20),
                       store_models = TRUE)

# Fit
at_glmnet$train(task)

# Final model
at_glmnet$tuning_result

# Setting otimal HP to the learner
at_glmnet$tuning_instance

# Calibration results after HPO, not CV predictions
results <- as.data.table(at_glmnet$predict(task))
ggplot(results, aes(x = truth, y = response)) + geom_point(size = 0.1)

# If we are interested in getting unbiased validation from CV, then we need to
# use nested resampling. However, as we focus on HPO and will just report the
# metrics from it, the true validation will be performed with external
# datasets

# Outer subsetting for evaluation
# 20 combinations x 5 outer folds
outer_resampling = rsmp("cv", folds = 5)
rr = resample(task, at_glmnet, outer_resampling, store_models = TRUE)

# tabular results
as.data.table(rr)

# Each outer fold has its owns hyperparameters with internal performance
extract_inner_tuning_results(rr)

# These are the outer score after HPO
rr$score()

# This is the final score, averaged
rr$aggregate()

# Each fold has its own models and HP
rr$learners[[1]]$tuning_result
rr$learners[[2]]$tuning_result
rr$learners[[3]]$tuning_result
rr$learners[[4]]$tuning_result
rr$learners[[5]]$tuning_result

# Getting outer predictions
cv.result <- as.data.table(rr, reassemble_learners = TRUE, convert_predictions = TRUE, predict_sets = "test")
cv.result$prediction

full.cv.results <- lapply(cv.result$prediction, function(x) as.data.table(x))
full.cv.results <- Reduce(rbind, full.cv.results)

ggplot(full.cv.results, aes(x = truth, y = response)) + geom_point(size = 0.1)

# Comparison
results %>%
  summarise(rmse = rmse_vec(truth = truth, estimate = response),
            bias = msd_vec(truth = truth, estimate = response),
            ccc = ccc_vec(truth = truth, estimate = response, bias = T),
            rsq = rsq_vec(truth = truth, estimate = response))

full.cv.results %>%
  summarise(rmse = rmse_vec(truth = truth, estimate = response),
            bias = msd_vec(truth = truth, estimate = response),
            ccc = ccc_vec(truth = truth, estimate = response, bias = T),
            rsq = rsq_vec(truth = truth, estimate = response))

