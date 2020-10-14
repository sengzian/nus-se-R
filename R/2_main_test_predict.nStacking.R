
predict.nStacking <- function(object, newdata, onlySL = FALSE, ...) {
  if (missing(newdata)) {
    stop("Please provide testing dataset for NUS-SE")
  }
  if (!object$control$saveFitLibrary) {
    stop("This NUS-SE fit was created using control$saveFitLibrary = FALSE, so new predictions cannot be made.")
  }

  k <- length(object$libraryNames)
  predY <- matrix(NA, nrow = nrow(newdata), ncol = k)
  colnames(predY) <- object$libraryNames
  if (onlySL) {
    
    whichLibrary <- which(object$coef > 0)
    predY <- matrix(0, nrow = nrow(newdata), ncol = k)
    for (mm in whichLibrary) {
      family <- object$family
      predY[, mm] <- do.call('predict', list(object = object$fitLibrary[[mm]],
                                             newdata = newdata,
                                             family = family))
    }
    getPred <- object$method$computePred(predY = predY, coef = object$coef, control = object$control)
    out <- list(pred = getPred, library.predict = predY)
  } else {
    for (mm in seq(k)) {
      family <- object$family
      predY[, mm] <- do.call('predict', list(object = object$fitLibrary[[mm]],
                                             newdata = newdata,
                                             family = family))
    }
    getPred <- object$method$computePred(predY = predY, coef = object$coef, control = object$control)
    out <- list(pred = getPred, library.predict = predY)
  }
  return(out)
}
