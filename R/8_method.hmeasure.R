method.hmeasure <- function(nlopt_method = "NLOPT_LN_SBPLX", bounds = c(0, Inf), normalize = TRUE, nloptr_iter=100) {

  # Contributed by marcus, sz
  nlopt_global <- c("NLOPT_GN_DIRECT",
                    "NLOPT_GN_DIRECT_L",
                    "NLOPT_GN_DIRECT_L_RAND",
                    "NLOPT_GN_DIRECT_NOSCAL",
                    "NLOPT_GN_DIRECT_L_NOSCAL",
                    "NLOPT_GN_DIRECT_L_RAND_NOSCAL",
                    "NLOPT_GN_ORIG_DIRECT",
                    "NLOPT_GN_ORIG_DIRECT_L",
                    "NLOPT_GN_CRS2_LM",
                    "NLOPT_GN_ISRES")
  nlopt_local <- c("NLOPT_LN_PRAXIS",
                   "NLOPT_LN_COBYLA",
                   "NLOPT_LN_NEWUOA_BOUND",
                   "NLOPT_LN_NELDERMEAD",
                   "NLOPT_LN_SBPLX",
                   "NLOPT_LN_BOBYQA")

  if (!(nlopt_method %in% c(nlopt_global, nlopt_local))) {
    stop("supplied 'nlopt_method' value not supported")
  }


  require = c('hmeasure', 'nloptr')
  
  computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
    
    .cvRisk_hmeasure <- function(par, Z, Y) {
      predictions <- crossprod(t(Z), par)
      H_object    <- hmeasure::HMeasure(true.class = Y, scores = predictions)
      h_measure   <- H_object$metrics$H
      cvRisk <- 1 - h_measure
      return(cvRisk)
    }

    # Step 1 - Coef Initialization
    coef_init <- rep(1/ncol(Z), ncol(Z))
    names(coef_init) <- libraryNames

    res <- nloptr::nloptr(x0 = coef_init,
                  eval_f = .cvRisk_hmeasure,
                  lb = rep(bounds[1], ncol(Z)),
                  ub = rep(bounds[2], ncol(Z)),
                  opts = list(algorithm = nlopt_method, xtol_rel = 1e-08, maxeval=nloptr_iter),
                  Z = Z,
                  Y = Y)


    if (res$status < 1 || res$status > 4) {
      warning(res$message)
    }
    coef <- res$solution
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)
    if (normalize) coef <- coef/sum(coef)

    hm <- apply(Z, 2, function(x) hmeasure::HMeasure(true.class = Y, scores = x)$metrics$H)
    cvRisk    <- 1 - hm
    names(coef) <- libraryNames

    out <- list(cvRisk = cvRisk, coef = coef, optimizer = res)
    return(out)
  }



  computePred = function(predY, coef, control, ...) {
    out <- crossprod(t(predY), coef)
    return(out)
  }


  out =list(require=require, computeCoef=computeCoef, computePred=computePred)
  invisible(out)
}




