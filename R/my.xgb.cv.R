#' Title Cross-validation for evaluating the performance of subgrouping using win-ratio based value function
#'
#' @param params the list of parameters for xgboost
#' @param data the data of predictors
#' @param label treatment label of subjects
#' @param comp.ind The vector of comparision results
#' @param nrounds the maximum number of boosting rounds.
#' @param nfold number of folds in cross validation
#' @param maximize a logical value, if \code{TRUE}, maximize the loss function, here we set it to \code{FALSE}
#' @param verbose a logical value, whether print out the messages of fitting the model
#' @param early_stopping_rounds the number of rounds before stopping if no improvement, defualt is 5
#'
#' @return a vector of errors
#' @import xgboost
#' @import caret
#' @export
#'
#' @examples NULL
my.xgb.cv<-function(    params,
                        data,
                        label,
                        comp.ind,
                        nrounds,
                        nfold,
                        maximize,
                        verbose,               # silent,
                        early_stopping_rounds # stop if no improvement for 10 consecutive trees
){
  n=dim(data)[1]
  set.seed(123)
  folds<-caret::createFolds(1:n, k=nfold, list=T, returnTrain=TRUE)
  cv.error<-list()
  for(f in 1:nfold){
    #cat(paste0("Cross-validation fold ",f, "\n"))
    data.train <- data[folds[[f]],]
    data.test <- data[-folds[[f]],]
    label.train <- label[folds[[f]]]
    label.test <- label[-folds[[f]]]
    id.train <<- expand.grid(i=1:length(intersect(which(label==0),folds[[f]])), j=1:length(intersect(which(label==1),folds[[f]])))
    id.test <<- expand.grid(i=1:length(setdiff(which(label==0),folds[[f]])), j=1:length(setdiff(which(label==1),folds[[f]])))
    comp.ind.train <<- as.numeric(matrix(comp.ind, nrow = sum(label==0), ncol = sum(label==1))[which(label==0) %in% folds[[f]],which(label==1) %in% folds[[f]]])
    comp.ind.test <<- as.numeric(matrix(comp.ind, nrow = sum(label==0), ncol = sum(label==1))[!(which(label==0) %in% folds[[f]]),!(which(label==1) %in% folds[[f]])])


    dtrain <- xgb.DMatrix(as.matrix(data.train), label = label.train)
    dtest <- xgb.DMatrix(as.matrix(data.test), label = label.test)
    #param <- list(max_depth = 2, eta = 0.001, silent = 1,
    #               objective = Myloss, eval_metric = evalerror,verbose = 1,lambda=1,base_score=0,colsample_bytree=1,min_child_weight=0)
    watchlist <- list(test = dtest)

    model <- xgb.train(params, dtrain, nrounds = nrounds, watchlist,
                       early_stopping_rounds = early_stopping_rounds,
                       verbose = verbose,
                       maximize = maximize)

    # record cv error
    cv.error[[f]]  <- model$evaluation_log$test_OTR_error

  }

  #cv.error.matrix<-do.call("rbind",cv.error)
  #test_OTR_error_mean<-c()
  cv.rounds<-min(sapply(cv.error,function(x) length(x)))
  cv.error.matrix<-matrix(nrow=nfold,ncol=cv.rounds)
  for(f in 1:nfold){
    for(m in 1:cv.rounds){
      cv.error.matrix[f,m]<-cv.error[[f]][m]
    }
  }
  test_OTR_error_mean<-colMeans(cv.error.matrix, na.rm = T)
  return(test_OTR_error_mean)
}

Myloss.train <- function(preds, dtrain) {
  trt01p <- getinfo(dtrain, "label")
  id=id.train
  comp.ind=comp.ind.train


  arm.val <- c(1,0)
  ## (1) Get Time to event Data Ready ##
  dat <- data.frame(trt01p)
  n<-dim(dat)[1]
  dat$id<-c(1:n)
  dat$f <- preds
  dat$pred <- 1/(1+exp(-preds))
  dat$predg <- exp(preds)/(1+exp(preds))^2

  gH1.r1.dat <- sum((dat$trt01p==arm.val[1])*dat$pred) #sum of prob for treatment arm in subgroup 1
  gH0.r1.dat <- sum((dat$trt01p==arm.val[2])*dat$pred) # sum of prob for control arm in subgroup 1

  gH1.r2.dat <- sum((dat$trt01p==arm.val[1])*(1-dat$pred)) # sum of prob for treatment arm in subgroup 2
  gH0.r2.dat <- sum((dat$trt01p==arm.val[2])*(1-dat$pred)) # sum of prob for control arm in subgroup 2

  ## subgroup 1 --- r1 ##
  if(gH1.r1.dat > 0 & gH0.r1.dat > 0){ # both arms have subjects in subgroup 1
    # Num win in subgroup 1
    Nw.r1 <- sum(dat$pred[id[,"i"]] * dat$pred[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==1))
    # Num lose in subgroup 1
    Nl.r1 <- sum(dat$pred[id[,"i"]] * dat$pred[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==-1))

    if(Nw.r1==0) Nw.r1 <- 0.5 # give it a small value if treatment arm always loss
    if(Nl.r1==0) Nl.r1 <- 0.5 # give it a small value if treatment arm always win

    # win ratio of subgroup 1
    Rw.r1 <-  Nw.r1/ Nl.r1
    #Rw.r1.g <-(gNw.r1*Nl.r1-gNl.r1*Nw.r1)/Nl.r1^2

    # gradient of log num win in subgroup 1
    log.Nw.r1.g <- sapply(1:n, function(i) sum(dat$pred[which(dat$trt01p==1)] * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==1)) +
                            sum(dat$pred[which(dat$trt01p==0)] * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==1))  ) /Nw.r1
    log.Nl.r1.g <- sapply(1:n, function(i) sum(dat$pred[which(dat$trt01p==1)] * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==-1)) +
                            sum(dat$pred[which(dat$trt01p==0)] * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==-1))  ) /Nl.r1
    log.Rw.r1.g <- log.Nw.r1.g - log.Nl.r1.g

  } else{
    Rw.r1 = 1 # do not contribute to the value function
    log.Rw.r1.g = 0 # do not contribute to the gradient
  }

  ## subgroup 2 --- r2 ##
  if(gH1.r2.dat > 0 & gH0.r2.dat > 0){ # both arms have subjects in subgroup 2
    # Num  win in subgroup 2
    Nw.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind== 1))
    # Num loss in subgroup 2
    Nl.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind== -1))

    if(Nw.r2==0) Nw.r2 <- 0.5 # give it a small value if equals 0
    if(Nl.r2==0) Nl.r2 <- 0.5 # give it a small value if equals 0

    # win ratio of subgroup 2
    Rw.r2 <-  Nw.r2/ Nl.r2

    # gradient of log num win in subgroup2
    #Rw.r2.g <-(gNw.r2*Nl.r2-gNl.r2*Nw.r2)/Nl.r2^2
    log.Nw.r2.g <- - sapply(1:n, function(i) sum((1-dat$pred[which(dat$trt01p==1)]) * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==1)) +
                              sum((1-dat$pred[which(dat$trt01p==0)]) * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==1))  ) /Nw.r2
    log.Nl.r2.g <- - sapply(1:n, function(i) sum((1-dat$pred[which(dat$trt01p==1)]) * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==-1)) +
                              sum((1-dat$pred[which(dat$trt01p==0)]) * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==-1))  ) /Nl.r2
    log.Rw.r2.g <- log.Nw.r2.g - log.Nl.r2.g

  } else {
    Rw.r2 = 1 # do not contribute to the value function
    log.Rw.r2.g = 0 # do not contribute to the gradient
  }

  g.p <- (sum(dat$pred)*log.Rw.r1.g + log(Rw.r1) - sum(1-dat$pred)*log.Rw.r2.g + log(Rw.r2))
  #h.p <- (2*rmst.diff.r1.g + sum(dat$pred)*rmst.diff.r1.h + 2*rmst.diff.r2.g - sum(1-dat$pred)*rmst.diff.r2.h)
  g <-  dat$predg*(-1)*g.p
  #h <- (-1)*( (dat$predg)^2 * h.p  + g.p*dat$predh)
  g <- g[order(dat$id)]
  #h <- h[order(dat$id)]
  h<-rep(0.00001,n)


  return(list(grad = g, hess = h))

}



evalerror.test <- function(preds, dtrain) {
  trt01p <- getinfo(dtrain, "label")

  id=id.test
  comp.ind=comp.ind.test

  arm.val <- c(1,0)

  ## (1) Get Time to event Data Ready ##
  dat <- data.frame(trt01p)
  dat$pred <- 1/(1+exp(-preds))
  #dat<-dat[order(dat$aval),]
  n<-dim(dat)[1]

  ## (2) value function ##

  gH1.r1.dat <- sum((dat$trt01p==arm.val[1])*dat$pred) #sum of prob for treatment arm in subgroup 1
  gH0.r1.dat <- sum((dat$trt01p==arm.val[2])*dat$pred) # sum of prob for control arm in subgroup 1

  gH1.r2.dat <- sum((dat$trt01p==arm.val[1])*(1-dat$pred)) # sum of prob for treatment arm in subgroup 2
  gH0.r2.dat <- sum((dat$trt01p==arm.val[2])*(1-dat$pred)) # sum of prob for control arm in subgroup 2

  ## subgroup 1 --- r1 ##
  if(gH1.r1.dat > 0 & gH0.r1.dat > 0){  # both arms have subjects in subgroup 1
    # Num win in subgroup 1
    Nw.r1 <- sum(dat$pred[id[,"i"]] * dat$pred[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==1))
    # Num lose in subgroup 1
    Nl.r1 <- sum(dat$pred[id[,"i"]] * dat$pred[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==-1))

    if(Nw.r1==0) Nw.r1 <- 0.5 # give it a small value if equals 0
    if(Nl.r1==0) Nl.r1 <- 0.5 # give it a small value if equals 0

    # win ratio of subgroup 1
    Rw.r1 <-  Nw.r1/ Nl.r1
  } else{
    Rw.r1 = 1 # do not contribute to the value function
  }

  ## subgroup 2 --- r2 ##
  if(gH1.r2.dat > 0 & gH0.r2.dat > 0){  # both arms have subjects in subgroup 2
    # Num win in subgroup 2
    Nw.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==1))
    # Num loss in subgroup 2
    Nl.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==-1))

    if(Nw.r2==0) Nw.r2 <- 0.5 # give it a small value if equals 0
    if(Nl.r2==0) Nl.r2 <- 0.5 # give it a small value if equals 0

    # win ratio of subgroup 2
    Rw.r2 <-  Nw.r2/ Nl.r2

  } else {
    Rw.r2 = 1 # do not contribute to the value function
  }

  # err is negative of value function
  err <- (-1)*( sum(dat$pred)*log(Rw.r1) - sum(1-dat$pred)*log(Rw.r2) )

  return(list(metric = "OTR_error", value = err))
}
