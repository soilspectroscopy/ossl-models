
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

# Other basic parameters
predictors <- paste0("PC", seq(1, ncomps))
regression.matrix <- ipca

# Regression matrix
set.seed(1993)
regression.matrix <- regression.matrix %>%
  select(id.layer_uuid_txt, all_of(block_var), all_of(isoil_property), all_of(predictors)) %>%
  filter(!is.na(!!as.name(isoil_property))) %>%
  sample_n(1000) %>%
  as.data.frame()

# Modeling pipeline
set.seed(1993)

# Task
train.task <- as_task_regr(regression.matrix, target = isoil_property, id = "train")

train.task$set_col_roles("id.layer_uuid_txt", roles = "name")

# Learners
learner_glmnet = lrn("regr.glmnet", predict_type = "response")

learner_ranger = lrn("regr.ranger", predict_type = "response",
                     replace = TRUE, num.threads = 1)

learner_xgboost = lrn("regr.xgboost", predict_type = "response",
                      booster = "gbtree", nthread = 1,
                      subsample = 0.67)

learner_cubist = lrn("regr.cubist", predict_type = "response",
                     neighbors = 0, committees = 1, unbiased = FALSE, seed = 1993)

# Internal CV
base_glmnet = po("learner_cv", learner_glmnet, id = "base_glmnet",
                 param_vals = list(resampling.method = "cv",
                                   resampling.folds = 5))

base_ranger = po("learner_cv", learner_ranger, id = "base_ranger",
                 param_vals = list(resampling.method = "cv",
                                   resampling.folds = 5))

base_xgboost = po("learner_cv", learner_xgboost, id = "base_xgboost",
                  param_vals = list(resampling.method = "cv",
                                    resampling.folds = 5))

base_cubist = po("learner_cv", learner_cubist, id = "base_cubist",
                 param_vals = list(resampling.method = "cv",
                                   resampling.folds = 5))

# Stacking in the pipeline
base_learners = gunion(list(
  base_glmnet,
  base_ranger,
  base_xgboost,
  base_cubist)) %>>%
  po("featureunion", id = "predictions")

# Meta-learner: linear model of base learners' cross-predictions
learner_lm = lrn("regr.lm", predict_type = "response")

meta_learner = po("learner", learner_lm, id = "meta_lm")

stack = gunion(list(base_learners %>>% meta_learner))

# Visualizing model graph
stack$plot(html = FALSE)

# Setting as a learner
learner_ensemble = as_learner(stack)
learner_ensemble$id = "ensemble"
learner_ensemble$predict_type = "response"

# Fitting without HPO, default HPs from mlr3
# learner_ensemble$train(train.task)
# results <- learner_ensemble$predict(train.task)
# ggplot(results, aes(x = truth, y = response)) + geom_point()

# Hyperparameters space, all crossed, i.e. not tuned separately
search_space_ensemble = ps(
  base_glmnet.alpha = p_dbl(0, 1),
  base_glmnet.lambda = p_dbl(0.0001, 1.0),
  base_ranger.num.trees = p_int(20, 100),
  base_ranger.min.node.size = p_int(5, 20),
  base_xgboost.nrounds = p_int(20, 100),
  base_xgboost.eta = p_dbl(0.3, 0.5),
  base_xgboost.max_depth = p_int(5, 20)
)

# Inner resampling for HPO
inner_resampling = rsmp("cv", folds = 5)

# Parallelization
future::plan("multisession")

# Auto tuner
at = auto_tuner(tuner = tnr("random_search"),
                learner = learner_ensemble,
                resampling = inner_resampling,
                measure = msr("regr.rmse"),
                search_space = search_space_ensemble,
                terminator = trm("evals", n_evals = 20),
                store_models = TRUE)

at$train(train.task)

# Overview of HPO
at$tuning_instance

# Overview of best HPs
at$tuning_result

# CV results
cv.results <- tibble(response = at$model$learner$model$meta_lm$model$fitted.values,
                     residual = at$model$learner$model$meta_lm$model$residuals) %>%
  mutate(truth = response+residual)

ggplot(cv.results, aes(x = truth, y = response)) + geom_point()

# Metrics
cv.results %>%
  summarise(rmse = rmse_vec(truth = truth, estimate = response),
            bias = msd_vec(truth = truth, estimate = response),
            ccc = ccc_vec(truth = truth, estimate = response, bias = T),
            rsq = rsq_vec(truth = truth, estimate = response))
