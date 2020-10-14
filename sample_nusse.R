
library(nusser)
setwd(paste0(.libPaths()[1],"/nusser/data"))

set.seed(4)
data_this_iter = list()
data_this_iter[["train"]]       <- read.csv(file = "yeast6-5-5tra.dat",header = T,sep = ',')
data_this_iter[["test"]]        <- read.csv(file = "yeast6-5-5tst.dat",header = T,sep = ',')
data_this_iter[["position_class"]]  <- ncol(data_this_iter$train)
data_this_iter[["model_target"]]    <- colnames(data_this_iter$train)[data_this_iter$position_class]
count.tab <- table(data_this_iter$train[,data_this_iter$position_class])
data_this_iter[["original_majority_label"]] <- names(count.tab)[count.tab==max(count.tab)]
data_this_iter[["original_minority_label"]] <- names(count.tab)[count.tab==min(count.tab)]
data_this_iter[["majority_count"]] <- count.tab[count.tab==max(count.tab)]
data_this_iter[["minority_count"]] <- count.tab[count.tab==min(count.tab)]
data_this_iter$train[,data_this_iter$position_class] <- factor(ifelse(data_this_iter$train[,data_this_iter$position_class]==data_this_iter$original_minority_label,"pos","neg"), levels = c("neg","pos"))
data_this_iter$test[,data_this_iter$position_class] <- factor(ifelse(data_this_iter$test[,data_this_iter$position_class]==data_this_iter$original_minority_label,"pos","neg"), levels = c("neg","pos"))
factor_pos <- which(sapply(data_this_iter$train, class)=="factor")
train.x.y <- data_this_iter$train
train.x <- data_this_iter$train[,-data_this_iter$position_class]
train.y <- data_this_iter$train[,data_this_iter$position_class]
test.x  <- data_this_iter$test[,-data_this_iter$position_class]
test.y  <- data_this_iter$test[,data_this_iter$position_class]
train.y.num <- as.numeric(train.y)-1
eval_y_true <- as.numeric(test.y) -1 


sl_choice = c("SL.kernelKnn_cat"
              , "SL.nnet_cat"
              , "SL.rf_oob"
              , "SL.xgb_oob" 
              , "SL.J48" 
              , "SL.nb" )
method.AUC_SBPLX =  nusser::method.AUC(nlopt_method = "NLOPT_LN_SBPLX", optim_method = NULL)
sl_method_vec   = c("method.AUC_SBPLX", "method.hmeasure", "method.CC_nloglik")             
sl_method       = sl_method_vec[2]
clt             = list() 
cvclt           = list(stratifyCV="SCV", V=10)
num_resampling  = 10
nus = list(mode             = "nus"
          ,ir                   = 1
          ,index                = TRUE
          ,replacement          = FALSE
          ,k_value              = 5
          ,exp_coef             = 1
          ,outlier_reset_weight = TRUE
          ,weight               = "HARD_T")


# Get local neighbourhood information
LN_information <- gen_rds_LN(data_this_iter, k_value=nus$k_value, dist="HVDM")
LN_information <- gen_rds_LN_hard_higher(LN_information$LN_n_hard)
LN_n_weight    <- assign_LN_n_weight(LN_information, exp_coef=nus$exp_coef) 

# Get NUS-SE model
nusse_model <-  nStacking_nus_doob(Y = train.y.num, X = train.x, newX = NULL, family = binomial(), SL.library=sl_choice,
                                  method = sl_method, id = NULL, verbose = FALSE, control = clt,
                                  cvControl = cvclt, obsWeights = NULL, num_resampling = num_resampling, 
                                  data_this_iter=data_this_iter, LN_n_weight=LN_n_weight, nus=nus)

# Prediction of NUS-SE model
nusse_pred <- predict.nStacking(object = nusse_model, newdata=test.x, onlySL = T)
pred_prob  <- as.numeric(nusse_pred$pred)

performance <- hmeasure::HMeasure(true.class = eval_y_true, scores = pred_prob)
performance$metrics$H
performance$metrics$AUC