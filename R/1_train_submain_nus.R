
nus <- function(data_this_iter, ir = 1, index = T, replacement = F, LN_n_weight_vec) {
  # ir = Imbalance Ratio. (how many times majority instances are over minority instances)

  target <- data_this_iter$position_class
  data <- data_this_iter$train

  p <- data[which(data[ ,target] == "pos"), ]
  n <- data[which(data[ ,target] == "neg"), ]

  # Step 3 - injection into resampling weight
  n <- n[sample(nrow(n), nrow(p) * ir, replace = replacement, prob= LN_n_weight_vec), ]
  result <- rbind(p, n)

  if(index==T) {result <- round(as.numeric(rownames(result)))}

  return(result)
}


nus_list <- function(num_resampling, data_this_iter, ir = 1, index = T, replacement = F, LN_n_weight_vec) {
  res <- vector(mode = "list", length = num_resampling)
  for (i in 1:num_resampling) {res[[i]] <- nus(data_this_iter, ir, index, replacement, LN_n_weight_vec)}
  
  res
}


gen_rds_LN <- function(data_this_iter, k_value, dist="HVDM") {
  target_pos        <- data_this_iter$position_class
  trn_data          <- data_this_iter$train
  trn_target        <- as.character(trn_data[[target_pos]])
  LN_all_hard       <- vector(length = nrow(trn_data))
  trn_distances     <- nusser:::distances(tgt = target_pos, dat = trn_data, dist = "HVDM")
  trn_distances     <- trn_distances$distm
  row_n             <- which(trn_data[ ,target_pos] == "neg")

  for (i in 1:nrow(trn_data)) {
    instance_class_distance <- data.frame(class=trn_data[[target_pos]], distance=trn_distances[,i])
    instance_class_distance <- instance_class_distance[-i,]
    instance_class_distance <- dplyr::arrange(instance_class_distance, distance, desc(class))
    instance_class_distance <- instance_class_distance[1:k_value,]
    LN_all_hard[i] <- sum(instance_class_distance$class=="pos")/k_value  # hard majority, has higher weight
  }

  LN_df <- data.frame(trn_target, num_pos = LN_all_hard)
  LN_n_hard = LN_all_hard[row_n]

  list(distances=trn_distances, LN_n_hard=LN_n_hard, k_value=k_value, dist= dist)
}



# get LN_n and raw distances data
gen_rds_LN_hard_higher <- function(LN_n_hard) {

  LN_n = LN_n_hard
  unique_LN_n_Ori <- sort(unique(LN_n))

  # Step 1.1 - outlier_reset_weight
  LN_n[LN_n==1] <- 0
  unique_LN_n_reset <- sort(unique(LN_n))

  list(LN_n=LN_n, unique_LN_n_Ori= unique_LN_n_Ori, unique_LN_n_reset=unique_LN_n_reset)
}


# get LN_n_weight
assign_LN_n_weight <- function (LN_object, exp_coef=1) {

  LN_n <- LN_object$LN_n

  # Step 2 - exp_coef - degree of increase weight for LN exponentially        # Function for resampling weight
  LN_n_weight <- sapply(LN_n, function(x, coef) exp(coef*x), coef=exp_coef)   # exp format
  LN_n_weight <- LN_n_weight/sum(LN_n_weight)                                 # normalized (prob)
                              
  LN_n_weight
}


distances <- function(tgt, dat, dist, p=2) {
  k <- 1
  if(p<1) stop("The parameter p must be >=1!")
  p<- switch(dist,
             "Chebyshev"=0,
             "Manhattan"=1,
             "Euclidean"=2,
             "Canberra"=-1,
             "Overlap"=-2,
             "HEOM"=-3,
             "HVDM"=-4,
             "p-norm"=p,
             stop("Distance measure not available!"))
  
  
  if (class(dat[,tgt]) == "numeric" & p <= -4) stop("distance measure selected only available for classification tasks")
  
  nomatr <- c() 
  for (col in seq.int(dim(dat)[2])) {
    if (class(dat[,col]) %in% c('factor','character')) {
      nomatr <- c(nomatr, col)
    }
  }
  
  nomatr <- setdiff(nomatr, tgt)
  numatr <- setdiff(seq.int(dim(dat)[2]), c(nomatr,tgt))
  
  nomData <- t(sapply(subset(dat, select = nomatr), as.integer))
  numData <- t(subset(dat, select = numatr))
  
  if (length(numatr) & p == -2) {
    stop("Can not compute Overlap metric with numeric attributes!")
  }
  if (length(nomatr) & p >= -1) {
    stop("Can not compute ", dist ," distance with nominal attributes!")
  }
  
  tgtData <- dat[, tgt]
  n <- length(tgtData)
  res <- matrix(0.0, nrow = k, ncol = n)
  if (class(tgtData) != "numeric") {tgtData <- as.integer(tgtData)}
  
  Cl <- length(unique(tgtData))
  nnom <- length(nomatr)
  nnum <- length(numatr)
  
  distm <- matrix(0.0, nrow = n, ncol = n)
  numD <- matrix(0.0, nrow = nnum, ncol = n)
  
  
  storage.mode(numData) <- "double"
  storage.mode(nomData) <- "integer"
  storage.mode(res) <- "integer"
  storage.mode(tgtData) <- "double"
  storage.mode(distm) <- "double"
  storage.mode(numD) <- "double"
  
  neig <- .Fortran("F_neighbours", 
                   tgtData = tgtData,  # tgt data
                   numData = numData, #numeric data
                   nomData = nomData, #nominal data
                   p = as.integer(p), # code for distance metric
                   k = as.integer(k), # nr of neighbours
                   n = as.integer(n), # nr of examples in the data
                   nnum = as.integer(nnum), # nr of numeric attributes
                   nnom = as.integer(nnom), # nr of nominal attributes
                   Cl = as.integer(Cl), # number of different classes in the target variable
                   distm = distm,
                   numD = numD,
                   res = res) # output  
  neig
}