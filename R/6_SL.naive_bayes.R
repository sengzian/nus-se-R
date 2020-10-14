

# Rweka

SL.nb <- function (Y, X, newX, oob_X = NULL, family, obsWeights, ...) 
{

	require("e1071")
	fit_nb <- naiveBayes(X, Y, laplace = 0)
	pred <- predict(fit_nb, newdata = newX, type = "raw")[,2]
	fit <- list(object = fit_nb)

    # Start - OOB handling
    pred_oob = NULL
    not_null_oob = !is.null(oob_X)
    if(not_null_oob) 
    {
        pred_oob <- predict(fit_nb, newdata = oob_X, type = "raw")[,2]
    }
    # End - OOB handling
    out = list(pred = pred, fit = fit, pred_oob = pred_oob)
    class(out$fit) <- c("SL.nb")
    return(out)
}


predict.SL.nb <- function (object, newdata, ...) 
{
    require("e1071")
    pred <- predict(object$object, newdata = newdata, type = "raw")[,2]
    return(pred)
}