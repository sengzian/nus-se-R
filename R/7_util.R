.onAttach <- function(...) {
  packageStartupMessage('nusser loaded')
  packageStartupMessage('Version: ', utils::packageDescription('nusser')$Version)
  packageStartupMessage('Package created on ', utils::packageDescription('nusser')$Date, '\n')
}

get_df_singular_col <- function(df.x) {which(sapply(df.x, function(x) length(unique(x)))==1)}
get_df_factor_col <- function(df.x) {which(sapply(df.x, class)=="factor")}


nStacking.control <- function(saveFitLibrary = TRUE, 
                                 saveCVFitLibrary = FALSE,
                                 trimLogit = 0.001) {
  if(trimLogit > 0.5) {
    warning('trimLogit must be less than 0.5, will replace with trimLogit = 0.001')
    trimLogit <- 0.001
  }
  list(saveFitLibrary = saveFitLibrary, trimLogit = trimLogit, saveCVFitLibrary = saveCVFitLibrary)
}

nStacking.CV.control <- function(V = 10L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL){

  V <- as.integer(V)
  if(!is.null(validRows)) {
    if(!is.list(validRows)) {
      stop('validRows must be a list of length V containing the row numbers for the corresponding validation set')
    }
    if(!identical(V, length(validRows))) {
      stop('V and length(validRows) must be identical')
    }
  }
  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}

trimLogit <- function(x, trim=0.00001) {
	x[x < trim] <- trim
	x[x > (1-trim)] <- (1-trim)
	foo <- log(x/(1-x))
	return(foo)
}

.SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}


.createLibrary <- function(SL.library) {
	if (is.character(SL.library)) { 
		k <- length(SL.library)
		whichScreen <- matrix(1, nrow = 1, ncol = k)
		screenAlgorithm <- "All"
		library <- data.frame(predAlgorithm = SL.library, rowScreen = 1, stringsAsFactors=FALSE)
	} else if (is.list(SL.library)) {
		predNames <- sapply(SL.library, FUN = "[", 1)
		NumberScreen <- (sapply(SL.library, FUN = length) - 1)
		if (sum(NumberScreen == 0) > 0) {
			for(ii in which(NumberScreen == 0)) {
				SL.library[[ii]] <- c(SL.library[[ii]], "All")
				NumberScreen[ii] <- 1
			}
		}
		screenAlgorithmFull <- unlist(lapply(SL.library, FUN="[", -1))
		screenAlgorithm <- unique(screenAlgorithmFull)
		
		library <- data.frame(predAlgorithm = rep(predNames, times=NumberScreen), rowScreen = match(screenAlgorithmFull, screenAlgorithm), stringsAsFactors = FALSE)
	} else {
	  stop('format for SL.library is not recognized')
	}
	
	out <- list(library = library, screenAlgorithm = screenAlgorithm)
	return(out)
}
