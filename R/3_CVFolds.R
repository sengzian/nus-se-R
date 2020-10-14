

CVFolds <- function(N, id, Y, cvControl, envir_nStacking=NULL){

	if(!is.null(cvControl$validRows)) {
	  return(cvControl$validRows)
	}
	stratifyCV <- cvControl$stratifyCV
	shuffle <- cvControl$shuffle
	V <- cvControl$V
	
	if(stratifyCV=="CV") {
		if(shuffle) {
			if(is.null(id)) {
				validRows <- split(sample(1:N), rep(1:V, length=N))
			} else {
				n.id <- length(unique(id))
				id.split <- split(sample(1:n.id), rep(1:V, length=n.id))
				validRows <- vector("list", V)
				for(v in seq(V)) {
					validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
				}
			}
		} else {
			if(is.null(id)) {
				validRows <- split(1:N, rep(1:V, length=N))
			} else {
				n.id <- length(unique(id))
				id.split <- split(1:n.id, rep(1:V, length=n.id))
				validRows <- vector("list", V)
				for(v in seq(V)) {
					validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
				}
			}
		}
	} else if (stratifyCV=="SCV") {
		if(length(unique(Y)) != 2) {
			stop("stratifyCV only implemented for binary Y")
		}
		if(sum(Y) < V) {
			cat(paste0("kfold original= ", V,"\n"))
			V <- sum(Y)
			cvControl$V <- sum(Y)
			assign('cvControl', cvControl, envir = envir_nStacking)
			print(paste0("kfold readjust to ", V))
		}
		if(shuffle) {
			if(is.null(id)) {
			  wiY0 <- which(Y == 0)
			  wiY1 <- which(Y == 1)
			  rowsY0 <- split(sample(wiY0), rep(1:V, length=length(wiY0)))
			  rowsY1 <- split(sample(wiY1), rep(1:V, length=length(wiY1)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for(vv in seq(V)) {
          validRows[[vv]] <- c(rowsY0[[vv]], rowsY1[[vv]])
        }
			} else {
				stop("stratified sampling with id not currently implemented")
			}
		} else {
			if(is.null(id)) {
				within.split <- suppressWarnings(tapply(1:N, INDEX = Y, FUN = split, rep(1:V)))
				validRows <- vector("list", length = V)
				names(validRows) <- paste(seq(V))
				for(vv in seq(V)) {
					validRows[[vv]] <- c(within.split[[1]][[vv]], within.split[[2]][[vv]])
				}
			} else {
				stop("stratified sampling with id not currently implemented")
			}
		}
	}

	return(validRows)
}
