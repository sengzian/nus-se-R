nStacking_nus_doob <- function(Y, X, newX = NULL, family = binomial(), SL.library,
                         method = 'method.NNLS', id = NULL, verbose = FALSE, control = list(),
                         cvControl = list(), obsWeights = NULL,
                         num_resampling = 2, data_this_iter, LN_n_weight, nus) {

    time_start = proc.time()
    env_nStacking <- environment()
    # Prep - Method Object, Control
    if (is.character(method)) {
        if (exists(method, mode = 'list')) {
            method <- get(method, mode = 'list')
        } else if (exists(method, mode = 'function')) {
            method <- get(method, mode = 'function')()
        }
    } else if (is.function(method)) {
        method <- method()
    }
    if(!is.list(method)) {
        stop("method is not in the appropriate format. Check out help('method.template')")
    }
    if(!is.null(method$require)) {
        sapply(method$require, function(x) require(force(x), character.only = TRUE))
    }
    # get defaults for controls and make sure in correct format
    control <- do.call('nStacking.control', control)
    cvControl <- do.call('nStacking.CV.control', cvControl)

# Prep - Lib, Z-meta, Dimensions, N and Z careful
    library <- nusser:::.createLibrary(SL.library)
    call <- match.call(expand.dots = TRUE)
    if(!inherits(X, 'data.frame')) message('X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs')
    varNames <- colnames(X)
    p <- dim(X)[2L]
    k <- nrow(library$library)


    if(p < 2 & !identical(library$screenAlgorithm, "All")) {
        warning('Screening algorithms specified in combination with single-column X.')
    }

# Prep - Lib in Environment - put fitLibrary in it's own environment to locate later
    fitLibEnv <- new.env()
    libraryNames <- expand.grid(library$library$predAlgorithm, 1:num_resampling)
    libraryNames <- paste0(libraryNames[,1], "_", libraryNames[,2])
    assign('fitLibrary', vector('list', length = num_resampling*k), envir = fitLibEnv)
    assign('libraryNames', libraryNames, envir = fitLibEnv)
    evalq(names(fitLibrary) <- libraryNames, envir = fitLibEnv)


# Prep
    if (!is.numeric(Y)) {
        stop("the outcome Y must be a numeric vector")
    }
    if(is.character(family))
        family <- get(family, mode="function", envir=parent.frame())
    if(is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (family$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
        stop("'method.AUC' is designed for the 'binomial' family only")
    }

    resampling_subsets_index <- nus_list(num_resampling, data_this_iter, 
                                        ir = nus$ir, index = nus$index, 
                                        replacement = nus$replacement, 
                                        LN_n_weight_vec=LN_n_weight)

    list_z              <- vector(mode = "list", length = length(resampling_subsets_index))
    list_predY          <- vector(mode = "list", length = length(resampling_subsets_index))
    names(list_z)       <- paste0("Classifiers_by_subset",1:num_resampling)
    names(list_predY)   <- paste0("Classifiers_by_subset",1:num_resampling)
    list_z_cv           <- vector(mode = "list", length = length(resampling_subsets_index))
    list_z_oob          <- vector(mode = "list", length = length(resampling_subsets_index))
    list_validRows      <- vector(mode = "list", length = length(resampling_subsets_index))
    complete_L1         <- matrix(NA, nrow(data_this_iter$train), k*num_resampling)
    complete_L0_predY   <- matrix(NA, nrow(data_this_iter$train), k*num_resampling)
    errorsInCVLibrary   <- vector(mode = "list", length = length(resampling_subsets_index))
    errorsInLibrary     <- vector(mode = "list", length = length(resampling_subsets_index))
    cvFitLibrary        <- vector(mode = "list", length = length(resampling_subsets_index))

time_train_start = proc.time()
for (subset_iter in 1:length(resampling_subsets_index)) 
{ # Start - forloop-subset
    # - L0 subsets Classifiers - fitLibEnv$fitLibrary
    # - L1 meta data           - list_z   
    # - L0 prediction          - list_predY # for training performance

    L1_predY_subset_iter <- get_L1_predY_doob(subset_iter, resampling_subsets_index, 
                                             X, Y, data_this_iter,
                                             library, libraryNames, 
                                             k, p, control, cvControl, family,
                                             obsWeights = NULL, id = NULL, newX = NULL,
                                             env_nStacking, fitLibEnv, verbose) 
    list_z[[subset_iter]]              <- L1_predY_subset_iter$Z_subset_iter
    list_predY[[subset_iter]]          <- L1_predY_subset_iter$predY_subset_iter
    errorsInCVLibrary[[subset_iter]]   <- L1_predY_subset_iter$errorsInCVLibrary
    errorsInLibrary[[subset_iter]]     <- L1_predY_subset_iter$errorsInLibrary
    cvFitLibrary[[subset_iter]]        <- L1_predY_subset_iter$cvFitLibrary_subset
    list_z_cv[[subset_iter]]           <- L1_predY_subset_iter$Z_cv
    list_z_oob[[subset_iter]]          <- L1_predY_subset_iter$Z_oob
    list_validRows[[subset_iter]]      <- L1_predY_subset_iter$validRows
} # End - forloop-subset

    errorsInCVLibrary       <- unlist(errorsInCVLibrary) # in .crossValFunc
    errorsInLibrary         <- unlist(errorsInLibrary) # in .predFunc
    complete_L1             <- do.call(cbind, list_z)
    complete_L0_predY       <- do.call(cbind, list_predY)
    colnames(complete_L1)       <- fitLibEnv$libraryNames
    colnames(complete_L0_predY) <- fitLibEnv$libraryNames
    names(errorsInCVLibrary)    <- fitLibEnv$libraryNames
    names(errorsInLibrary)      <- fitLibEnv$libraryNames
    
# Train - get Coef using Z - 
    # Compute weights for each algorithm in library.
    obsWeights <- rep(1, nrow(complete_L1))
    # meta learner
    getCoef <- method$computeCoef(Z = complete_L1, Y = Y, libraryNames = libraryNames,
                                  obsWeights = obsWeights, control = control,
                                  verbose = verbose,
                                  errorsInLibrary = errorsInCVLibrary)
    coef <- getCoef$coef
    names(coef) <- libraryNames

    time_train = proc.time() - time_train_start

    # Set a default in case the method does not return the optimizer result.
    if (!("optimizer" %in% names(getCoef))) {
        getCoef["optimizer"] <- NA
    }

# Train - whichScreen
    time_predict_start = proc.time()

# Train - SL.predict
    # NUS-SE prediction on training dataset for training performance
    getPred <- method$computePred(predY = complete_L0_predY, coef = coef, control = control)
    time_predict = proc.time() - time_predict_start

    # Clean up when errors in library.
    if(sum(errorsInCVLibrary) > 0) {
        getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
    }

    time_end = proc.time()

    # Compile execution times.
    times = list(everything = time_end - time_start,
                 train = time_train,
                 predict = time_predict)

    # Put everything together in a list.
    out <- list(
        call = call,
        libraryNames = libraryNames,
        SL.library = library,
        SL.predict = getPred,
        coef = coef,
        library.predict = complete_L0_predY,
        validRows = list_validRows,
        z_oob = list_z_oob,
        z_cv = list_z_cv,
        Z = complete_L1,
        cvRisk = getCoef$cvRisk,
        family = family,
        fitLibrary = get('fitLibrary', envir = fitLibEnv),
        cvFitLibrary = cvFitLibrary,
        varNames = varNames,
        method = method,
        control = control,
        cvControl = cvControl,
        errorsInCVLibrary = errorsInCVLibrary,
        errorsInLibrary = errorsInLibrary,
        metaOptimizer = getCoef$optimizer,
        LN_n_weight = LN_n_weight,
        resampling_subsets_index = resampling_subsets_index,
        times = times
    )
    class(out) <- c("nStacking")
    return(out)
}



get_L1_predY_doob <- function(subset_iter, resampling_subsets_index, 
                        X, Y, data_this_iter,
                        library, libraryNames, 
                        k, p, control, cvControl, family,
                        obsWeights = NULL, id = NULL, newX = NULL,
                        env_nStacking, fitLibEnv, verbose
                        ) 
{
    temp_predY_predOOB  <- vector(mode = "list", length = k)
    temp_predY_cbind    <- matrix(NA, nrow(data_this_iter$train), k)

    subset_index    <- resampling_subsets_index[[subset_iter]]
    subset_X        <- X[subset_index,]
    subset_Y        <- Y[subset_index]
    OOB_X           <- X[-subset_index,]

    N       <- dim(subset_X)[1L]
    Z_cv    <- matrix(NA, N, k)
    rownames(Z_cv) <- as.numeric(rownames(subset_X))

    validRows <- CVFolds(N = N, id = id, Y = subset_Y, cvControl = cvControl, envir_nStacking = env_nStacking)

    if(is.null(newX)) {newX <- subset_X}
    if(is.null(id)) {id <- seq(N)}
    if(is.null(obsWeights)) {obsWeights <- rep(1, N)} # test observation weights
    if(!identical(length(obsWeights), N)) {stop("obsWeights vector must have the same dimension as subset_X")}
    if(!identical(length(id), N)) {stop("id vector must have the same dimension as subset_X")} 


    # -----------------------------  Start - GET L1     ----------------------------- 
    # Train - .crossValFUN
    # OOB_X = NULL, no prediction on OOB is perform.
    crossValFUN_out <- lapply(validRows, FUN = nusser:::.crossValFunc, 
                              Y = subset_Y, dataX = subset_X, 
                              id = id, obsWeights = obsWeights, 
                              library = library, 
                              k = k, p = p, libraryNames = libraryNames,
                              saveCVFitLibrary = control$saveCVFitLibrary,
                              OOB_X = NULL, verbose= verbose, family= family
                              )
    # Train - L1 Data - Z
    Z_cv[unlist(validRows, use.names = FALSE), ] <- do.call('rbind', lapply(crossValFUN_out, "[[", "out"))
    
    # Train - save CV classifiers
    if(control$saveCVFitLibrary){
        cvFitLibrary_subset <- lapply(crossValFUN_out, "[[", "model_out")
    }else{
        cvFitLibrary_subset <- NULL
    }

    # -----------------------------  Start - GET predY     ----------------------------- 
    # Train - library.predict | .predFun
    for(index_k in 1:k) {
    temp_predY_predOOB[[index_k]] <- nusser:::.predFunc(subset_iter=subset_iter, index_k=index_k, lib = library$library, 
                            Y = subset_Y, X = subset_X, newX = subset_X, OOB_X = OOB_X,
                            family = family, id = id, obsWeights = obsWeights, verbose = verbose,
                            control = control, libraryNames = libraryNames, fitLibEnv=fitLibEnv, algos =k)
    }
    temp_predY      <- lapply(temp_predY_predOOB, "[[", "out")
    temp_predOOB    <- lapply(temp_predY_predOOB, "[[", "out_oob")
    temp_predY_cbind    <- do.call('cbind', temp_predY)
    temp_predOOB_cbind  <- do.call('cbind', temp_predOOB)

    if(sum(rownames(OOB_X)!=rownames(temp_predOOB_cbind)) >0) {
      stop("Check rownames for temp_predOOB_cbind")
    } else {
        Z <- rbind(temp_predOOB_cbind, Z_cv)
        Z <- Z[order(as.numeric(rownames(Z))),]
    }

    # -----------------------------  End - GET L1     -----------------------------

    # Check for errors. If any algorithms had errors, replace entire column with
    # 0 even if error is only in one fold.
    errorsInCVLibrary <- apply(Z, 2, function(x) anyNA(x))
    if (sum(errorsInCVLibrary) > 0) {
        Z[, as.logical(errorsInCVLibrary)] <- 0
    }
    if (all(Z == 0)) {
        stop("All algorithms dropped from library (CV step for meta data)")
    }

    # check for errors
    errorsInLibrary <- apply(temp_predY_cbind, 2, function(algorithm) anyNA(algorithm))
    if (sum(errorsInLibrary) > 0) {
        if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
            stop(paste0("Failed (prediction) algorithm(s)\n"))
        }
    }
    # -----------------------------  End - GET predY     ----------------------------- 
    list(Z_subset_iter=Z, predY_subset_iter=temp_predY_cbind, 
        Z_cv = Z_cv, Z_oob = temp_predOOB_cbind, validRows=validRows,
        errorsInCVLibrary= errorsInCVLibrary, errorsInLibrary=errorsInLibrary, cvFitLibrary_subset=cvFitLibrary_subset)
}



    

