saveRDS.gz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

## remove extreme residuals and retrain models
retrain.ossl <- function(t.var, model.name, out.m.rds, t=3, SL.library = c("regr.ranger", "regr.xgboost", "regr.cvglmnet", "regr.cubist", "regr.plsr")){
  in.dir=paste0("./models/", t.var, "/")
  out.dir=paste0("./modelsF/", t.var, "/")
  t.m = readRDS.gz(paste0(in.dir, model.name, ".rds"))
  res = t.m$learner.model$super.model$learner.model$residuals
  q <- quantile(res, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(res)
  sel.rm <- !(res > (q[1] - t*iqr) & res < (q[2] + t*iqr))
  #summary(sel.rm)
  id.rm = as.integer(attr(res, "names")[which(sel.rm)])
  #t.var = all.vars(t.m$learner.model$super.model$learner.model$terms)[1]
  if(missing(out.m.rds)){ out.m.rds <- paste0(out.dir, model.name, ".rds") }
  out.t.mrf = paste0(in.dir, model.name, "_t.mrf.rds")
  out.t.xgb = paste0(in.dir, model.name, "_t.xgb.rds")
  if(!file.exists(out.m.rds) & file.exists(out.t.mrf) & file.exists(out.t.xgb)){
    rm.x = readRDS.gz(paste0(in.dir, model.name, "_rm.rds"))
    rm.x = rm.x[-id.rm,]
    parallelMap::parallelStartSocket(parallel::detectCores())
    tskF <- mlr::makeRegrTask(data = rm.x, target = t.var, blocking = attr(rm.x, "ID")[-id.rm])
    var.mod1 <- readRDS.gz(out.t.mrf)
    var.mod2 <- readRDS.gz(out.t.xgb)
    lrn.rf = mlr::setHyperPars(mlr::makeLearner(SL.library[1]), par.vals = getHyperPars(var.mod1$learner))
    lrn.xg = mlr::setHyperPars(mlr::makeLearner(SL.library[2]), par.vals = getHyperPars(var.mod2$learner))
    lrnsE <- list(lrn.rf, lrn.xg, mlr::makeLearner(SL.library[3]), mlr::makeLearner(SL.library[4]), mlr::makeLearner(SL.library[5]))
    init.m <- mlr::makeStackedLearner(base.learners = lrnsE, predict.type = "response", method = "stack.cv", super.learner = "regr.lm", resampling=makeResampleDesc(method = "CV", blocking.cv=TRUE))
    t.m <- mlr::train(init.m, tskF)
    #summary(t.m$learner.model$super.model$learner.model)
    saveRDS.gz(t.m, out.m.rds)
    parallelMap::parallelStop()
  }
}

## Model fine-tuning wrapper
## by default run spatial CV
train.ossl <- function(t.var, pr.var, X, model.name, out.dir=paste0("./models/", t.var, "/"), SL.library = c("regr.ranger", "regr.xgboost", "regr.cvglmnet", "regr.cubist"), discrete_ps = makeParamSet(makeDiscreteParam("mtry", values = seq(5, round(length(pr.var)*.7), by=5))), ctrl = mlr::makeTuneControlGrid(), rdesc = mlr::makeResampleDesc("CV", iters = 2L), outer = mlr::makeResampleDesc("CV", iters = 2L), inner = mlr::makeResampleDesc("Holdout"), ctrlF = mlr::makeFeatSelControlRandom(maxit = 20), xg.model_Params, hzn_depth=TRUE, out.m.rds, save.rm=TRUE, rf.feature=TRUE, xg.size=2e3){
  require(mlr)
  ## https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
  if(missing(xg.model_Params)){
    xg.model_Params <- makeParamSet(
      #nrounds=50; max_depth=4; eta=0.2; subsample=1; min_child_weight=4; colsample_bytree=0.6
      makeDiscreteParam("nrounds", value=c(20,50)),
      makeDiscreteParam("max_depth", value=c(3,4,5,6)),
      makeDiscreteParam("eta", value=c(0.3,0.4)),
      makeDiscreteParam("subsample", value=c(1)),
      makeDiscreteParam("min_child_weight", value=c(1,4)),
      makeDiscreteParam("colsample_bytree", value=c(0.6))
    )
  }
  if(substr(t.var, 1, 5)=="log.."){
    X[,t.var] = log1p(X[,gsub("log..", "", t.var)])
  }
  ## regression matrix
  rm.x = X[,c(t.var, pr.var)]
  r.sel = complete.cases(rm.x)
  rm.x = rm.x[r.sel,]
  attr(rm.x, "ID") = as.factor(X$ID[r.sel])
  if(save.rm==TRUE){ saveRDS.gz(rm.x, paste0(out.dir, model.name, "_rm.rds")) }
  parallelMap::parallelStartSocket(parallel::detectCores())
  ## subset to speed up fine-tuning
  sel.rf = sample.int(nrow(rm.x), size=xg.size)
  tsk0 <- mlr::makeRegrTask(data = rm.x[sel.rf,], target = t.var, blocking = as.factor(X$ID[r.sel][sel.rf]))
  out.t.mrf = paste0(out.dir, model.name, "_t.mrf.rds")
  if(!file.exists(out.t.mrf)){
    ## fine-tune mtry
    resR.lst = tuneParams(mlr::makeLearner("regr.ranger", num.threads = round(parallel::detectCores()/length(discrete_ps$pars$mtry$values)), num.trees=85), task = tsk0, resampling = rdesc, par.set = discrete_ps, control = ctrl)
    ## feature selection
    lrn.rf = mlr::makeLearner("regr.ranger", num.threads = parallel::detectCores(), mtry=resR.lst$x$mtry, num.trees=85, importance="impurity")
    if(rf.feature==TRUE){
      lrn1 = mlr::makeFeatSelWrapper(lrn.rf, resampling = inner, control = makeFeatSelControlSequential(method="sfbs", maxit = 20), show.info=TRUE) ##
    } else {
      lrn1 = lrn.rf
    }
    var.mod1 = mlr::train(lrn1, task = tsk0)
    saveRDS.gz(var.mod1, out.t.mrf)
  } else {
    var.mod1 = readRDS.gz(out.t.mrf)
  }
  if(rf.feature==TRUE){
    var.sfeats1 = mlr::getFeatSelResult(var.mod1)
  } else {
    var.sfeats1 = data.frame(x=pr.var)
  }
  out.t.xgb = paste0(out.dir, model.name, "_t.xgb.rds")
  if(!file.exists(out.t.xgb)){
    ## fine-tune xgboost
    sel.xg = sample.int(nrow(rm.x), size=xg.size)
    tsk0s <- mlr::makeRegrTask(data = rm.x[sel.xg,], target = t.var, blocking = as.factor(X$ID[r.sel][sel.xg]))
    resX.lst = mlr::tuneParams(mlr::makeLearner("regr.xgboost"), task = tsk0s, resampling = rdesc, par.set = xg.model_Params, control = ctrl)
    lrn.xg = mlr::makeLearner("regr.xgboost", par.vals = list(objective ='reg:squarederror'))
    lrn.xg = mlr::setHyperPars(lrn.xg, par.vals = resX.lst$x)
    lrn2 = mlr::makeFeatSelWrapper(lrn.xg, resampling = inner, control = ctrlF, show.info=TRUE)
    var.mod2 = mlr::train(lrn2, task = tsk0s)
    saveRDS.gz(var.mod2, out.t.xgb)
  } else {
    var.mod2 = readRDS.gz(out.t.xgb)
  }
  var.sfeats2 = mlr::getFeatSelResult(var.mod2)
  ## new shorter formula
  ## we add depth otherwise not a 3D model
  if(hzn_depth==TRUE){
    formulaString.y = as.formula(paste(t.var, ' ~ ', paste(c("hzn_depth", union(var.sfeats1$x, var.sfeats2$x)), collapse="+")))
  } else {
    formulaString.y = as.formula(paste(t.var, ' ~ ', paste(union(var.sfeats1$x, var.sfeats2$x), collapse="+")))
  }
  ## save formula:
  ## remove all covariates without enough variation
  #s.lst = sapply(all.vars(formulaString.y), function(i){sd(rm.x[,i], na.rm=TRUE)})
  #all.vars(formulaString.y)[which(s.lst<5)]
  ## final EML model
  if(missing(out.m.rds)){ out.m.rds <- paste0(out.dir, model.name, ".rds") }
  if(!file.exists(out.m.rds)){
    tskF <- mlr::makeRegrTask(data = rm.x[,all.vars(formulaString.y)], target = t.var, blocking = as.factor(X$ID[r.sel]))
    var.mod1 <- readRDS.gz(out.t.mrf)
    var.mod2 <- readRDS.gz(out.t.xgb)
    lrn.rf = mlr::setHyperPars(mlr::makeLearner(SL.library[1]), par.vals = getHyperPars(var.mod1$learner))
    lrn.xg = mlr::setHyperPars(mlr::makeLearner(SL.library[2]), par.vals = getHyperPars(var.mod2$learner))
    lrnsE <- list(lrn.rf, lrn.xg, mlr::makeLearner(SL.library[3]), mlr::makeLearner(SL.library[4]))
    init.m <- mlr::makeStackedLearner(base.learners = lrnsE, predict.type = "response", method = "stack.cv", super.learner = "regr.lm", resampling=makeResampleDesc(method = "CV", blocking.cv=TRUE))
    t.m <- mlr::train(init.m, tskF)
    #t.m$learner.model$super.model$learner.model
    saveRDS.gz(t.m, out.m.rds)
  }
  parallelMap::parallelStop()
}

cat_eml = function(t.var, model.name, in.dir=paste0("./models/", t.var, "/"), out.dir=paste0("./models/", t.var, "/"), n.max=50){
  r.file = paste0(out.dir, model.name, "_resultsFit.txt")
  if(!file.exists(r.file)){
    out.m.rds = paste0(out.dir, model.name, ".rds")
    t.m = readRDS.gz(out.m.rds)
    x.s = summary(t.m$learner.model$super.model$learner.model)
    cat("Results of ensemble model fitting 'ranger', 'xgboost', 'Cubist', 'cvglmnet' ...:\n", file=r.file)
    cat("\n", file=r.file, append=TRUE)
    cat(paste("Variable:", t.var, "\n"), file=r.file, append=TRUE)
    cat(paste("R-square:", round(x.s$adj.r.squared, 3), "\n"), file=r.file, append=TRUE)
    cat(paste("Fitted values sd:", signif(sd(t.m$learner.model$super.model$learner.model$fitted.values), 3), "\n"), file=r.file, append=TRUE)
    cat(paste("RMSE:", signif(sqrt(sum(t.m$learner.model$super.model$learner.model$residuals^2) / t.m$learner.model$super.model$learner.model$df.residual), 3), "\n\n"), file=r.file, append=TRUE)
    sink(file=r.file, append=TRUE, type="output")
    cat("EML model summary:", file=r.file, append=TRUE)
    print(x.s)
    cat("\n", file=r.file, append=TRUE)
    imp.rds = paste0(in.dir, model.name, "_t.mrf.rds")
    if(file.exists(imp.rds)){
      var.mod1 <- readRDS.gz(imp.rds)
      cat("Variable importance:\n", file=r.file, append=TRUE)
      xl <- as.data.frame(mlr::getFeatureImportance(var.mod1)$res)
      write.csv(xl[order(xl$importance, decreasing=TRUE),], paste0(out.dir, model.name, "_rf_varImportance.csv"))
      print(xl[order(xl$importance, decreasing=TRUE),][1:n.max,])
    }
    sink()
  }
}


## confidence limits based on RMSE:
pfunL <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)
  panel.abline(0,1,lty=1,lw=2,col="black")
  panel.abline(0+m$Summary$logRMSE,1,lty=3,lw=2,col="black")
  panel.abline(0-m$Summary$logRMSE,1,lty=3,lw=2,col="black")
}

pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)
  panel.abline(0,1,lty=1,lw=2,col="black")
}

plot_hexbin <- function(varn, breaks, main, meas, pred, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1), pal = openair::openColours("increment", 18)[-18], in.file, log.plot, out.file){ ## pal=R_pal[["bpy_colors"]][1:18]
  require("hexbin"); require("plotKML"); require("latticeExtra"); require("openair")
  if(missing(out.file)){ out.file = paste0("./outputs/plot_CV_", varn, ".png") }
  if(!file.exists(out.file)){
    if(missing(pred)){
      m <- readRDS.gz(in.file)
      #pred <- t.m$learner.model$super.model$learner.model$fitted.values
      pred <- m$predictions
    }
    if(missing(meas)){
      meas <- t.m$learner.model$super.model$learner.model$model[,1]
    }
    if(log.plot==TRUE){
      R.squared = yardstick::ccc(data.frame(pred, meas), truth="meas", estimate="pred")
      #pred <- 10^pred-1
      pred = expm1(pred)
      #meas <- 10^meas-1
      meas = expm1(meas)
      d.meas <- min(meas, na.rm=TRUE)
    } else {
      d.meas <- min(meas, na.rm=TRUE)
      R.squared = yardstick::ccc(data.frame(pred, meas), truth="meas", estimate="pred")
    }
    main.txt = paste0(main, "  (CCC: ", signif(R.squared$.estimate, 3), ")")
    png(file = out.file, res = 150, width=850, height=850, type="cairo")
    if(log.plot==TRUE){
      pred <- pred+ifelse(d.meas==0, 1, d.meas)
      meas <- meas+ifelse(d.meas==0, 1, d.meas)
      lim <- range(breaks)+ifelse(d.meas==0, 1, d.meas)
      meas <- ifelse(meas<lim[1], lim[1], ifelse(meas>lim[2], lim[2], meas))
      plt <- hexbinplot(meas~pred, colramp=colorRampPalette(pal), main=main.txt, ylab="measured", xlab="predicted", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, scales=list(x = list(log = 2, equispaced.log = FALSE), y = list(log = 2, equispaced.log = FALSE)), asp=1, xbins=50, ybins=50, xlim=lim, ylim=lim, panel=pfun, colorcut=colorcut)
    } else {
      lim <- range(breaks)
      meas <- ifelse(meas<lim[1], lim[1], ifelse(meas>lim[2], lim[2], meas))
      plt <- hexbinplot(meas~pred, colramp=colorRampPalette(pal), main=main.txt, ylab="measured", xlab="predicted", type="g", lwd=1, lcex=8, inner=.2, cex.labels=.8, xlim=lim, ylim=lim, asp=1, xbins=50, ybins=50, panel=pfun, colorcut=colorcut)
    }
    print(plt)
    dev.off()
  }
}
