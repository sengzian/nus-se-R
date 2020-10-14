# Rweka

SL.J48 <- function (Y, X, newX, oob_X = NULL, family, obsWeights, ...) 
{

	require("RWeka")
	fit_J48 <- RWeka:::J48(formula=as.formula("model_target ~ ."), data = cbind(X, model_target=as.factor(Y)), control = Weka_control(), options = NULL)
	pred <- predict(fit_J48, newdata = newX, type = "probability")[,2]
	fit <- list(object = fit_J48)
	rJava:::.jcache(fit$object$classifier)
    # Start - OOB handling
    pred_oob = NULL
    not_null_oob = !is.null(oob_X)
    if(not_null_oob) 
    {
        pred_oob <- predict(fit_J48, newdata = oob_X, type = "probability")[,2]
    }
    # End - OOB handling
    out = list(pred = pred, fit = fit, pred_oob = pred_oob)
    class(out$fit) <- c("SL.J48")
    return(out)
}


predict.SL.J48 <- function (object, newdata, ...) 
{
    require("RWeka")
    pred <- predict(object$object, newdata = newdata, type = "probability")[,2]
    return(pred)
}