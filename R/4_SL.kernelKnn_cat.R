SL.kernelKnn_cat <- function (Y, X, newX, oob_X=NULL, family, k = 5, method = "euclidean", weights_function = NULL, extrema = F, h = 1, ...) 
{
    # if the weights_function is NULL then a simple k-nearest-neighbor is performed
    # https://cran.r-project.org/web/packages/KernelKnn/KernelKnn.pdf

    #Start: to handle categorical variable for X and newX
    c_fac_X  <- get_df_factor_col(X)
    c_fac_newX  <- get_df_factor_col(newX)

    if(length(c_fac_X) > 0) {
      X[,c_fac_X] <- sapply(X[,c_fac_X], as.numeric)
    }
    if(length(c_fac_newX) > 0) {
      newX[,c_fac_newX] <- sapply(newX[,c_fac_newX], as.numeric)
    }
    #End: to handle categorical variable for X and newX

    nusser:::.SL.require("KernelKnn")
    if (family$family != "gaussian" && min(Y) == 0) {
        Y = Y + 1
    }
    if (family$family == "gaussian") {
        levels = NULL
    }
    else {
        levels = unique(Y)
    }
    pred = KernelKnn::KernelKnn(data = X, TEST_data = newX, y = Y, 
        k = k, method = method, h = h, weights_function = weights_function, 
        regression = family$family == "gaussian", Levels = levels)

    if (family$family == "binomial") {
        pred = pred[, 2]
    }

    fit = list(k = k, method = method, weights_function = weights_function, 
        extrema = extrema, h = h, X = X, Y = Y, family = family)

    # Start - OOB handling
    pred_oob = NULL
    not_null_oob = !is.null(oob_X)
    if(not_null_oob) 
    {
        c_fac_oobX  <- get_df_factor_col(oob_X)
        if(length(c_fac_oobX) > 0) {
          oob_X[,c_fac_oobX] <- sapply(oob_X[,c_fac_oobX], as.numeric)
        }

        pred_oob = KernelKnn::KernelKnn(data = X, TEST_data = oob_X, y = Y, 
        k = k, method = method, h = h, weights_function = weights_function, 
        regression = family$family == "gaussian", Levels = levels)

        if (family$family == "binomial") {pred_oob = pred_oob[, 2]}
    }
    # End - OOB handling

    out = list(pred = pred, fit = fit, pred_oob = pred_oob)

    class(out$fit) = "SL.kernelKnn_cat"
    return(out)
}

predict.SL.kernelKnn_cat <- function (object, newdata, ...) 
{
    #Start: to handle categorical variable for newdata
    c_fac_newdata  <- get_df_factor_col(newdata)
    if(length(c_fac_newdata) > 0) {
      newdata[,c_fac_newdata] <- sapply(newdata[,c_fac_newdata], as.numeric)
    }
    #End: to handle categorical variable for newdata
    
    nusser:::.SL.require("KernelKnn")
    if (object$family$family == "gaussian") {
        levels = NULL
    }
    else {
        levels = unique(object$Y)
    }
    pred = KernelKnn::KernelKnn(data = object$X, TEST_data = newdata, 
        y = object$Y, k = object$k, method = object$method, h = object$h, 
        weights_function = object$weights_function, regression = object$family$family == 
            "gaussian", Levels = levels)
    if (object$family$family == "binomial") {
        pred = pred[, 2]
    }
    return(pred)
}