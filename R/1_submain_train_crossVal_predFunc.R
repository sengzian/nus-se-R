
.crossValFunc <- function(valid, Y, dataX, id, obsWeights, library, 
                         kScreen=NULL, k, p, libraryNames, saveCVFitLibrary, OOB_X, verbose, family) {
    tempLearn <- dataX[-valid, , drop = FALSE]
    tempOutcome <- Y[-valid]
    tempValid <- dataX[valid, , drop = FALSE]
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]

    out <- matrix(NA, nrow = nrow(tempValid), ncol = k)
    if(!is.null(OOB_X)) {oob_out <- matrix(NA, nrow = nrow(OOB_X), ncol = k)
        } else {oob_out <- NULL}
        
    if(saveCVFitLibrary){
        model_out <- vector(mode = "list", length = k)
    }else{
        model_out <- NULL
    }

    for(s in seq(k)) {
        pred_fn = get(library$library$predAlgorithm[s])
        testAlg <- try(do.call(pred_fn, list(Y = tempOutcome, X = tempLearn, newX = tempValid, oob_X = OOB_X, family = family, id = tempId, obsWeights = tempObsWeights)))
        if(inherits(testAlg, "try-error")) {
            warning(paste("Error in algorithm", library$library$predAlgorithm[s], "\n  The Algorithm will be removed (i.e. given weight 0) \n" ))
        } else {
            out[, s] <- testAlg$pred
            if(!is.null(OOB_X)) {oob_out[, s] <- testAlg$pred_oob}
            if(saveCVFitLibrary){
                model_out[[s]] <- testAlg$fit
            }
        }
        if (verbose) message(paste("CV", libraryNames[s]))
    } #end library
    if(saveCVFitLibrary){
        names(model_out) <- libraryNames
    }
    invisible(list(out = out, model_out = model_out, oob_out = oob_out))
} 


# for base learner prediction
.predFunc <- function(subset_iter, index_k, lib, Y, X, newX, OOB_X, family, id, obsWeights, verbose, control, libraryNames, fitLibEnv, algos) {
    pred_fn = get(lib$predAlgorithm[index_k])
    testAlg <- try(do.call(pred_fn, list(Y = Y, X = X, newX = newX, oob_X = OOB_X, family = family, id = id, obsWeights = obsWeights)))
    if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", lib$predAlgorithm[index_k], " on full data", "\n  The Algorithm will be removed (i.e. given weight 0) \n" ))
        out <- rep.int(NA, times = nrow(newX))
    } else {
        out     <- testAlg$pred
        out_oob <- testAlg$pred_oob
        if (control$saveFitLibrary) {
            previous_classifiers <- (subset_iter-1)*algos
            eval(bquote(fitLibrary[[.(index_k+previous_classifiers)]] <- .(testAlg$fit)), envir = fitLibEnv)
        }
    }

    invisible(list(out=out, out_oob=out_oob))
}
