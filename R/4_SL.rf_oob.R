# randomForest{randomForest}

SL.rf_oob <- function(Y, X, newX, oob_X = NULL, family, mtry = ifelse(family$family == "gaussian",
                            max(floor(ncol(X)/3), 1), floor(sqrt(ncol(X)))), ntree = 1000,
                            nodesize = ifelse(family$family == "gaussian", 5, 1),
                            maxnodes = NULL,importance = FALSE, ...) {
	.SL.require('randomForest')
	if (family$family == "gaussian") {
		fit.rf <- randomForest::randomForest(Y ~ ., data = X, ntree = ntree, xtest = newX, keep.forest = TRUE, mtry = mtry, nodesize = nodesize, maxnodes = maxnodes, importance = importance)
		pred <- fit.rf$test$predicted
		fit <- list(object = fit.rf)
	}
	if (family$family == "binomial") {
		fit.rf <- randomForest::randomForest(y = as.factor(Y), x = X, ntree = ntree, xtest = newX, keep.forest = TRUE, mtry = mtry, nodesize = nodesize, maxnodes = maxnodes, importance = importance)
		pred <- fit.rf$test$votes[, 2]
		fit <- list(object = fit.rf)
	}

	# Start - OOB handling
	pred_oob = NULL
	not_null_oob = !is.null(oob_X)
	if(not_null_oob) 
	{
		if (family$family == "binomial") {pred_oob <- predict(fit.rf, newdata = oob_X, type = 'vote')[,2]}
	}
	# End - OOB handling

	out = list(pred = pred, fit = fit, pred_oob = pred_oob)
	class(out$fit) <- c("SL.rf_oob")
	return(out)
}

predict.SL.rf_oob <- function(object, newdata, family, ...) {
	.SL.require('randomForest')
	if (family$family == "gaussian") {
		pred <- predict(object$object, newdata = newdata, type = 'response')
	}
	if (family$family == "binomial") {
		pred <- predict(object$object, newdata = newdata, type = 'vote')[,2]
	}
	pred
}
