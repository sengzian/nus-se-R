SL.nnet_cat <- function (Y, X, newX, oob_X = NULL, family, obsWeights, size = 2, ...) 
{

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


    nusser:::.SL.require("nnet")
    if (family$family == "gaussian") {
        fit.nnet <- nnet::nnet(x = X, y = Y, size = size, linout = TRUE, trace = FALSE, maxit = 500, weights = obsWeights)
    }
    if (family$family == "binomial") {
        fit.nnet <- nnet::nnet(x = X, y = Y, size = size, trace = FALSE, maxit = 500, linout = FALSE, weights = obsWeights)
    }
    pred <- predict(fit.nnet, newdata = newX, type = "raw")

    fit <- list(object = fit.nnet)
    
    # Start - OOB handling
    pred_oob = NULL
    not_null_oob = !is.null(oob_X)
    if(not_null_oob) 
    {
        c_fac_oobX  <- get_df_factor_col(oob_X)
        if(length(c_fac_oobX) > 0) {
          oob_X[,c_fac_oobX] <- sapply(oob_X[,c_fac_oobX], as.numeric)
        }

        pred_oob <- predict(fit.nnet, newdata = oob_X, type = "raw")
    }
    # End - OOB handling
    
    out = list(pred = pred, fit = fit, pred_oob = pred_oob)
    class(out$fit) <- c("SL.nnet_cat")
    return(out)
}
predict.SL.nnet_cat <- function (object, newdata, ...) 
{
    #Start: to handle categorical variable for newdata
    c_fac_newdata  <- get_df_factor_col(newdata)
    if(length(c_fac_newdata) > 0) {
      newdata[,c_fac_newdata] <- sapply(newdata[,c_fac_newdata], as.numeric)
    }
    #End: to handle categorical variable for newdata

    nusser:::.SL.require("nnet")
    pred <- predict(object$object, newdata = newdata, type = "raw")
    return(pred)
}

