## OSSL functions

subset.ossl <- function(ossl.x, tvar.y="ph.h2o_usda.4c1_index", dataset.y, harm.x="dataset.code_ascii_c",
                        visnir.x=paste0("scan_visnir.", seq(350, 2500, by=2), "_pcnt"),
                        mir.x=paste0("scan_mir.", seq(600, 4000, by=2), "_abs"),
                        geo.x="clm_", ID="ID"){
  if(missing(ossl.x)){
    ossl.x = readRDS(url("http://s3.us-east-1.wasabisys.com/soilspectroscopy/ossl_import/rm.ossl_v1", "rb"))
  }
  if(missing(dataset.y)){
    dataset.y = levels(as.factor(ossl.x$dataset.code_ascii_c))
  }
  if(!is.null(geo.x)){
    geo.sel = names(ossl.x)[grep(geo.x, names(ossl.x))]
  } else {
    geo.sel = NULL
  }
  x.lst = c(tvar.y, harm.x, visnir.x, mir.x, geo.sel, ID)
  sel.r = complete.cases(ossl.x[,x.lst]) & ossl.x$dataset.code_ascii_c %in% dataset.y
  out = ossl.x[sel.r,x.lst]
  if(!is.null(harm.x)){
    out <- fastDummies::dummy_cols(out, select_columns = harm.x)[,-which(names(out)==harm.x)]
  }
  ## remove layers that have no variation
  c.sd = sapply(out[,-which(names(out)==ID)], function(i){var(i, na.rm=TRUE)})
  c.r = which(c.sd == 0)
  if(length(c.r)>0){
    out = out[,-c.r]
  }
  return(out)
}

## Target variables
site.name = c("id.location_olc_c", "id.layer_uuid_c", "observation.ogc.schema.title_ogc_txt",
              "observation.ogc.schema_idn_url", "observation.date.begin_iso.8601_yyyy.mm.dd",
              "observation.date.end_iso.8601_yyyy.mm.dd", "location.address_utf8_txt", "location.country_iso.3166_c",
              "location.method_any_c", "surveyor.title_utf8_txt", "surveyor.contact_ietf_email",
              "surveyor.address_utf8_txt", "longitude_wgs84_dd", "latitude_wgs84_dd",
              "location.error_any_m", "dataset.title_utf8_txt", "dataset.owner_utf8_txt",
              "dataset.code_ascii_c", "dataset.address_idn_url",
              "dataset.license.title_ascii_txt", "dataset.license.address_idn_url", "dataset.doi_idf_c",
              "dataset.contact.name_utf8_txt", "dataset.contact_ietf_email", "id.project_ascii_c",
              "id.user.site_ascii_c", "pedon.taxa_usda_c", "pedon.completeness_usda_uint8",
              "layer.sequence_usda_uint16", "layer.type_usda_c", "layer.field.label_any_c",
              "layer.upper.depth_usda_cm", "layer.lower.depth_usda_cm", "horizon.designation_usda_c",
              "horizon.designation.discontinuity_usda_c", "layer.structure.type_usda_c",
              "layer.structure.grade_usda_c", "layer.texture_usda_c")
soilab.name = c("id.layer_uuid_c", "id.layer_local_c", "sample.doi_idf_c", "sample.contact.name_utf8_txt",
                "sample.contact.email_ietf_email", "acid.tea_usda4b2_cmolkg", "al.dith_usda.4g1_wpct",
                "al.kcl_usda.4b3_cmolkg", "al.ox_usda.4g2_wpct", "bsat_usda.4b4_wpct",
                "bd.clod_usda.3b1_gcm3", "bd.od_usda.3b2_gcm3", "ca.ext_usda.4b1_cmolkg", "c.tot_usda.4h2_wpct",
                "caco3_usda.4e1_wpct", "cec.ext_usda.4b1_cmolkg", "gyp_usda.4e2_wpct",
                "ecec_usda.4b4_cmolkg", "ec.w_usda.4f1_dsm", "oc_usda.calc_wpct", "fe.dith_usda.4g1_wpct",
                "fe.kcl_usda.4b3_mgkg", "fe.ox_usda.4g2_wpct", "mg.ext_usda.4b1_cmolkg", "n.tot_usda.4h2_wpct",
                "ph.kcl_usda.4c1_index", "ph.h2o_usda.4c1_index", "ph.cacl2_usda.4c1_index",
                "ph.naf_usda.4c1_index", "p.ext_usda.4d6_mgkg", "p.olsn_usda.4d5_mgkg",
                "k.ext_usda.4b1_cmolkg", "sand.tot_usda.3a1_wpct", "wpg2_usda.3a2_wpct",
                "silt.tot_usda.3a1_wpct", "clay.tot_usda.3a1_wpct", "na.ext_usda.4b1_cmolkg",
                "s.tot_usda.4h2_wpct", "sum.bases_4b4b2a_cmolkg", "wr.33kbar_usda.3c1_wpct",
                "wr.1500kbar_usda.3c2_wpct", "al.meh3_usda.4d6_wpct", "as.meh3_usda.4d6_mgkg",
                "ba.meh3_usda.4d6_mgkg", "ca.meh3_usda.4d6_mgkg", "cd.meh3_usda.4d6_wpct",
                "co.meh3_usda.4d6_mgkg", "cr.meh3_usda.4d6_mgkg", "cu.meh3_usda.4d6_mgkg",
                "p.meh3_usda.4d6_mgkg", "k.meh3_usda.4d6_mgkg", "na.meh3_usda.4d6_mgkg",
                "mg.meh3_usda.4d6_mgkg", "fe.meh3_usda.4d6_mgkg", "pb.meh3_usda.4d6_mgkg",
                "zn.meh3_usda.4d6_mgkg", "mo.meh3_usda.4d6_mgkg", "si.meh3_usda.4d6_mgkg",
                "sr.meh3_usda.4d6_mgkg")
mir.name = c("id.scan_uuid_c", "id.scan_local_c", "id.layer_uuid_c", "id.layer_local_c", "model.name_utf8_txt",
             "model.code_any_c", "method.light.source_any_c",
             "method.preparation_any_c", "scan.file_any_c", "scan.date.begin_iso.8601_yyyy.mm.dd",
             "scan.date.end_iso.8601_yyyy.mm.dd", "scan.license.title_ascii_txt", "scan.license.address_idn_url",
             "scan.doi_idf_c", "scan.contact.name_utf8_txt", "scan.contact.email_ietf_email",
             paste0("scan_mir.", seq(600, 4000, by=2), "_abs"))
visnir.name = c("id.scan_uuid_c", "id.scan_local_c", "id.layer_uuid_c", "id.layer_local_c", "model.name_utf8_txt",
                "model.code_any_c", "method.light.source_any_c", "method.preparation_any_c",
                "scan.file_any_c", "scan.date.begin_iso.8601_yyyy.mm.dd", "scan.date.end_iso.8601_yyyy.mm.dd",
                "scan.license.title_ascii_txt", "scan.license.address_idn_url", "scan.doi_idf_c",
                "scan.contact.name_utf8_txt", "scan.contact.email_ietf_email",
                paste0("scan_visnir.", seq(350, 2500, by=2), "_pcnt"))

## Function to find a mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## translate names
transvalues = function(df, out.name, in.name, fun.lst){
  if(!length(out.name)==length(in.name)){
    stop("Arguments 'out.name' and 'in.name' not equal length")
  }
  if(missing(fun.lst)){
    fun.lst = as.list(rep("x*1", length(out.name)))
  }
  ## https://stackoverflow.com/questions/61094854/storing-functions-in-an-r-list
  utility.fns = lapply(1:length(fun.lst), function(i){function(x){eval(parse(text = fun.lst[[i]]) )}})
  out <- as.data.frame(lapply(1:length(in.name), function(i){utility.fns[[i]](df[,in.name[i]])}))
  names(out) = out.name
  return(out)
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

hor2xyd <- function(x, U="UHDICM", L="LHDICM", treshold.T=15){
  x$DEPTH <- x[,U] + (x[,L] - x[,U])/2
  x$THICK <- x[,L] - x[,U]
  sel <- x$THICK < treshold.T
  ## begin and end of the horizon:
  x1 <- x[!sel,]; x1$DEPTH = x1[,L]
  x2 <- x[!sel,]; x2$DEPTH = x1[,U]
  y <- do.call(rbind, list(x, x1, x2))
  return(y)
}

plot_gh <- function(pnts, out.pdf, world, lats, longs, crs_goode = "+proj=igh", fill.col="yellow"){
  # https://wilkelab.org/practicalgg/articles/goode.html
  require(cowplot)
  require(sf)
  require(rworldmap)
  require(ggplot2)
  if(missing(world)){ world <- sf::st_as_sf(rworldmap::getMap(resolution = "low")) }
  if(missing(lats)){
    lats <- c(
      90:-90, # right side down
      -90:0, 0:-90, # third cut bottom
      -90:0, 0:-90, # second cut bottom
      -90:0, 0:-90, # first cut bottom
      -90:90, # left side up
      90:0, 0:90, # cut top
      90 # close
    )
  }
  if(missing(longs)){
    longs <- c(
      rep(180, 181), # right side down
      rep(c(80.01, 79.99), each = 91), # third cut bottom
      rep(c(-19.99, -20.01), each = 91), # second cut bottom
      rep(c(-99.99, -100.01), each = 91), # first cut bottom
      rep(-180, 181), # left side up
      rep(c(-40.01, -39.99), each = 91), # cut top
      180 # close
    )
  }
  goode_outline <-
    list(cbind(longs, lats)) %>%
    st_polygon() %>%
    st_sfc(
      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    )
  # now we need to work in transformed coordinates, not in long-lat coordinates
  goode_outline <- st_transform(goode_outline, crs = crs_goode)
  # get the bounding box in transformed coordinates and expand by 10%
  xlim <- st_bbox(goode_outline)[c("xmin", "xmax")]*1.1
  ylim <- st_bbox(goode_outline)[c("ymin", "ymax")]*1.1
  # turn into enclosing rectangle
  goode_encl_rect <-
    list(
      cbind(
        c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
      )
    ) %>%
    st_polygon() %>%
    st_sfc(crs = crs_goode)
  # calculate the area outside the earth outline as the difference
  # between the enclosing rectangle and the earth outline
  goode_without <- st_difference(goode_encl_rect, goode_outline)
  m <- ggplot(world) + geom_sf(fill = "gray80", color = "black", size = 0.5/.pt) +
    geom_sf(data = goode_without, fill = "white", color = "NA") +
    geom_sf(data = goode_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
    cowplot::theme_minimal_grid() + theme(panel.background = element_rect(fill = "#56B4E950", color = "white", size = 1),  panel.grid.major = element_line(color = "gray30", size = 0.25)) +
    geom_sf(data = pnts, size = 0.8, shape = 21, fill = fill.col, color="black") +
    #geom_sf(data = pnts, size = 1, pch="+", color="black") +
    coord_sf(crs = crs_goode, xlim = 0.95*xlim, ylim = 0.95*ylim, expand = FALSE)
  ggsave(out.pdf, m, dpi=150, height = 5.35, width = 9)
}

extract.cog <- function(pnts, cog.lst, url="http://s3.us-east-1.wasabisys.com/soilspectroscopy/layers1km/", mc.cores = 16, local=TRUE){
  if(local==TRUE){
    ov.tmp = parallel::mclapply(1:length(cog.lst), function(j){ terra::extract(terra::rast(cog.lst[j]), terra::vect(pnts)) }, mc.cores = mc.cores)
  } else {
    if(mc.cores>1){
      ov.tmp = parallel::mclapply(1:length(cog.lst), function(j){ terra::extract(terra::rast(paste0("/vsicurl/", url, cog.lst[j])), terra::vect(pnts)) }, mc.cores = mc.cores)
    } else {
      ov.tmp = lapply(1:length(cog.lst), function(j){ terra::extract(terra::rast(paste0("/vsicurl/", url, cog.lst[j])), terra::vect(pnts)) })
    }
  }
  suppressMessages( ov.tmp <- dplyr::bind_cols(lapply(ov.tmp, function(i){i[,2]})) )
  names(ov.tmp) = tools::file_path_sans_ext(basename(cog.lst))
  return(ov.tmp)
}

predict.ossl <- function(t.var, mir.raw, visnir.raw, lon, lat, hzn_depth=10, ossl.model, ossl.pca.mir, ossl.pca.visnir, spc.type="mir", subset.type="ossl", geo.type="na", n.spc=60, sd=TRUE, cog.dir="/data/WORLDCLIM/", ylim=NULL, dataset.code_ascii_c="KSSL.SSL"){ ## =c(0,100)
  ## check that input scans pass some minimum checks
  if(spc.type == "mir" | spc.type == "visnir.mir"){
    if(!class(mir.raw)=="data.frame"){
      stop("Input dataset '*.raw' not a correctly formated scan file. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
    if(nrow(mir.raw)>1000 | ncol(mir.raw)<1700){
      stop("Input dataset '*.raw' dimensions invalid. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
  }
  if(spc.type == "visnir" | spc.type == "visnir.mir"){
    if(!class(visnir.raw)=="data.frame"){
      stop("Input dataset '*.raw' not a correctly formated scan file. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
    if(nrow(visnir.raw)>1000 | ncol(visnir.raw)<1200){
      stop("Input dataset '*.raw' dimensions invalid. See https://soilspectroscopy.github.io/ossl-manual/ for examples.")
    }
  }
  if(missing(ossl.model)){
    model.rds = paste0("http://s3.us-east-1.wasabisys.com/soilspectroscopy/ossl_models/", t.var, "/", spc.type, "_mlr..eml_", subset.type, "_", geo.type, "_v1.rds")
    ossl.model = readRDS(url(model.rds, "rb"))
  }
  ## convert to PCs
  if(spc.type == "mir" | spc.type == "visnir.mir"){
    wn = as.numeric(gsub("X", "", names(mir.raw)))
    spc = as.matrix(mir.raw)
    #colnames(spc) = paste(wn)
    spc = as.data.frame(prospectr::resample(spc, wn, seq(600, 4000, by=2), interpol = "spline"))
    spc = lapply(spc, function(j){ round(ifelse(j<0, NA, ifelse(j>3, NA, j))*1000) })
    spc = as.data.frame(do.call(cbind, spc))
    names(spc) = paste0("scan_mir.", seq(600, 4000, by=2), "_abs")
    class(ossl.pca.mir) = "prcomp"
    X1.pc = as.data.frame(predict(ossl.pca.mir, newdata=spc))[,1:n.spc]
    colnames(X1.pc) = paste0("mir.PC", 1:n.spc)
  } else {
    X1.pc = NA
  }
  if(spc.type == "visnir" | spc.type == "visnir.mir"){
    wn = as.numeric(gsub("X", "", names(visnir.raw)))
    spc = as.matrix(visnir.raw)
    #colnames(spc) = paste(wn)
    spc = as.data.frame(prospectr::resample(spc, wn, seq(350, 2500, by=2), interpol = "spline"))
    spc = lapply(spc, function(j){ round(ifelse(j<0, NA, ifelse(j>1, NA, j))*100, 1) })
    spc = as.data.frame(do.call(cbind, spc))
    names(spc) = paste0("scan_visnir.", seq(350, 2500, by=2), "_pcnt")
    class(ossl.pca.visnir) = "prcomp"
    X2.pc = as.data.frame(predict(ossl.pca.visnir, newdata=spc))[,1:n.spc]
    colnames(X2.pc) = paste0("visnir.PC", 1:n.spc)
  } else {
    X2.pc = NA
  }
  ## obtain GeoTIFF values
  if(geo.type=="ll"){
    pnts = SpatialPoints(data.frame(lon, lat), proj4string = CRS("EPSG:4326"))
    cog.lst = paste0(cog.dir, ossl.model$features[grep("clm_", ossl.model$features)], ".tif")
    ov.tmp = extract.cog(pnts, cog.lst)
  } else {
    ov.tmp = NA
  }
  ## Bind all covariates together
  X = do.call(cbind, list(X1.pc, X2.pc, ov.tmp, data.frame(hzn_depth=hzn_depth)))
  X = X[,which(unlist(lapply(X, function(x) !all(is.na(x)))))]
  X$dataset.code_ascii_c = factor(rep(dataset.code_ascii_c, nrow(X)), levels = c("NEON.SSL", "KSSL.SSL", "CAF.SSL", "AFSIS1.SSL", "ICRAF.ISRIC", "LUCAS.SSL"))
  X <- fastDummies::dummy_cols(X, select_columns = "dataset.code_ascii_c")
  ## predict
  pred = predict(ossl.model, newdata=X[,ossl.model$features])
  ## uncertainty
  if(sd==TRUE){
    out.c <- as.matrix(as.data.frame(mlr::getStackedBaseLearnerPredictions(ossl.model, newdata=X[,ossl.model$features])))
    cf = eml.cf(ossl.model)
    model.error <- sqrt(matrixStats::rowSds(out.c, na.rm=TRUE)^2 * cf)
  }
  ## Return result as a data.frame:
  out = data.frame(pred.mean=pred$data$response, pred.error=model.error)
  ## back-transform
  if(length(grep("log..", t.var))>0){
    out$tpred.mean = expm1(out$pred.mean)
    if(!is.null(ylim)) {
      out$tpred.mean = ifelse(out$tpred.mean< ylim[1], ylim[1], ifelse(out$tpred.mean > ylim[2], ylim[2], out$tpred.mean))
    }
    out$lower.1std = expm1(out$pred.mean - out$pred.error)
    out$lower.1std = ifelse(out$lower.1std<0, NA, out$lower.1std)
    out$upper.1std = expm1(out$pred.mean + out$pred.error)
    out$upper.1std = ifelse(out$upper.1std<0, NA, out$upper.1std)
  } else {
    out$lower.1std = out$pred.mean - out$pred.error
    out$upper.1std = out$pred.mean + out$pred.error
  }
  if(!is.null(ylim)) {
    out$lower.1std = ifelse(out$lower.1std < ylim[1], ylim[1], ifelse(out$lower.1std > ylim[2], ylim[2], out$lower.1std))
    out$upper.1std = ifelse(out$upper.1std < ylim[1], ylim[1], ifelse(out$upper.1std > ylim[2], ylim[2], out$upper.1std))
  }
  return(list(pred=out, x=X, model=ossl.model$learner.model$super.model$learner.model$model, cf=cf))
}

eml.cf = function(t.m){
  m.train = t.m$learner.model$super.model$learner.model$model
  m.terms = all.vars(t.m$learner.model$super.model$learner.model$terms)
  eml.MSE0 = matrixStats::rowSds(as.matrix(m.train[,m.terms[-1]]), na.rm=TRUE)^2
  eml.MSE = deviance(t.m$learner.model$super.model$learner.model)/df.residual(t.m$learner.model$super.model$learner.model)
  ## correction factor:
  eml.cf = eml.MSE/mean(eml.MSE0, na.rm = TRUE)
  return(eml.cf)
}

## Model fine-tuning wrapper
## by default run spatial CV
train.ossl <- function(t.var, pr.var, X, model.name, out.dir=paste0("./models/", t.var, "/"), SL.library = c("regr.ranger", "regr.xgboost", "regr.cvglmnet", "regr.cubist"), discrete_ps = makeParamSet(makeDiscreteParam("mtry", values = seq(5, round(length(pr.var)*.6), by=5))), ctrl = mlr::makeTuneControlGrid(), rdesc = mlr::makeResampleDesc("CV", iters = 2L), outer = mlr::makeResampleDesc("CV", iters = 2L), inner = mlr::makeResampleDesc("Holdout"), ctrlF = mlr::makeFeatSelControlRandom(maxit = 20), xg.model_Params, hzn_depth=TRUE, out.m.rds, save.rm=FALSE, rf.feature=TRUE, xg.size=2e3){
  require(mlr)
  ## "regr.glmboost",
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

cat_eml = function(t.var, model.name, out.dir=paste0("./models/", t.var, "/"), n.max=50){
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
    imp.rds = paste0(out.dir, model.name, "_t.mrf.rds")
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

