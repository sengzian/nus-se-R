
SL.xgb_oob = function(Y, X, newX, oob_X = NULL,  family, obsWeights, id, ntrees = 1000,
                      max_depth = 4, shrinkage = 0.1, minobspernode = 10,
                      params = list(),
                      nthread = 1,
                      verbose = 0,
                      save_period = NULL,
                      ...) {
  .SL.require("xgboost")
  if(packageVersion("xgboost") < 0.6) stop("SL.xgboost requires xgboost version >= 0.6, try help(\'SL.xgboost\') for details")

  if (!is.matrix(X)) {
    X = model.matrix(~ . - 1, X)
  }

  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)

  if (family$family == "gaussian") {
    model = xgboost::xgboost(data = xgmat, objective="reg:linear", nrounds = ntrees,
                max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                verbose = verbose, nthread = nthread, params = params,
                save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                verbose = verbose, nthread = nthread, params = params,
                save_period = save_period)
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective="multi:softmax", nrounds = ntrees,
                max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                verbose = verbose, num_class = length(unique(Y)), nthread = nthread,
                params = params,
                save_period = save_period)
  }

  if (!is.matrix(newX)) {
    newX = model.matrix(~ . - 1, newX)
  }

  pred = predict(model, newdata = newX)
  
  fit = list(object = model)
  class(fit) = c("SL.xgb_oob")

  # Start - OOB handling
  pred_oob = NULL
  not_null_oob = !is.null(oob_X)
  if(not_null_oob) 
  {
    if (!is.matrix(oob_X)) {oob_X = model.matrix(~ . - 1, oob_X)}
    pred_oob = predict(model, newdata = oob_X)
  }
  # End - OOB handling

  out = list(pred = pred, fit = fit, pred_oob = pred_oob)
  return(out)
}


predict.SL.xgb_oob <- function(object, newdata, family, ...) {
  .SL.require("xgboost")
  if(packageVersion("xgboost") < 0.6) stop("SL.xgboost requires xgboost version >= 0.6, try help(\'SL.xgboost\') for details")

  if (!is.matrix(newdata)) {
    newdata = model.matrix(~ . - 1, newdata)
  }
  pred = predict(object$object, newdata = newdata)
  return(pred)
}


create.SL.xgboost = function(tune = list(ntrees = c(1000), max_depth = c(4), shrinkage = c(0.1),
                             minobspernode = c(10)), detailed_names = F, env = .GlobalEnv,
                             name_prefix = "SL.xgb") {

  tuneGrid = expand.grid(tune, stringsAsFactors=F)

  names = rep("", nrow(tuneGrid))

  for (i in seq(nrow(tuneGrid))) {
    g = tuneGrid[i,]
    if (detailed_names) {
      name = paste(name_prefix, g$ntrees, g$max_depth, g$shrinkage, g$minobspernode, sep=".")
    } else {
      name = paste(name_prefix, i, sep=".")
    }
    names[i] = name
    eval(parse(text = paste0(name, "= function(..., ntrees = ", g$ntrees, ", max_depth = ", g$max_depth, ", shrinkage=", g$shrinkage, ", minobspernode=", g$minobspernode, ") SL.xgboost(..., ntrees = ntrees, max_depth = max_depth, shrinkage=shrinkage, minobspernode=minobspernode)")), envir = env)
  }
  results = list(grid = tuneGrid, names = names)
  invisible(results)
}
