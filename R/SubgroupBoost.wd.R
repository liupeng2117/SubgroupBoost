#' Title Win-difference based value function guided subgrouping
#'
#' @param dat a dataframe, rows are subject and columns are predictors,
#' (censored) or 1 (event).
#' @param info a dataframe of treatment indicator variable \code{trt01p}, event indicator variable \code{evnt},
#' and time variable \code{aval}. \code{trt01p} is coded as -1 (control) or 1 (treatment),and \code{evnt} is coded as 0
#' (censored) or 1(event). If there are two survival outcomes, use \code{evnt1},\code{evnt2} and \code{aval1},\code{aval2}
#' to distinguish the two outcomes.
#' @param comparison The types of outcome to compare in win-difference based value function. The possible choices are
#' \itemize{
#' \item{\code{survival survival}}{Two time-to-event outcomes}
#' \item{\code{survival}}{one time-to-event outcome}
#' \item{\code{continous}}{one continous outcome, larger value means more benefit from treatment}
#' \item{\code{binary}}{one binary outcome, larger value means more benefit from treatment}
#' \item{\code{ordinal}}{one ordinal discrete outcome, larger value means more benefit from treatment}
#' }
#'
#'
#' @details add details later
#' @return a fitted gradient boosting model
#' @export SubgroupBoost.wd
#' @import xgboost
#'
#' @examples NULL
SubgroupBoost.wd <- function(dat, info, comparison){
  #dat contains only predictive variables, info contains treatment and outcomes

  ## ----- Define comparison rule ------ ##
  #compare two survival observations
  comp.surv<-function(Y1a, Y1b){
    # if do not belong to any of the cases below(i.e missing), set it to unkown
    out=0

    if(anyNA(Y1a)==F & anyNA(Y1b)==F){
      if(Y1a[2] != Y1b[2]){
        ind<-which.min(c(Y1a[2], Y1b[2]))
        if(c(Y1a[1],Y1b[1])[ind] == 1 & Y1a[2] > Y1b[2]){
          out = -1
        } else if(c(Y1a[1],Y1b[1])[ind] == 1 & Y1a[2] < Y1b[2]){
          out = 1
        } else if(c(Y1a[1],Y1b[1])[ind] == 0){
          out=0
        }
      }

      if(Y1a[2] == Y1b[2]){
        if(Y1a[1] == Y1b[1]) {
          out=0
        } else if(Y1a[1] == 0 & Y1b[1] == 1){
          out=-1
        } else if(Y1a[1] == 1 & Y1b[1] == 0){
          out=1
        }
      }
    }

    return(out)
  }

  # get final comparison output from the two comparisons
  get.final.out<-function(out1, out2){
    if(out1 == 1) {
      final.out=1
    } else if(out1 == -1){
      final.out=-1
    } else if(out1 == 0 & out2 == 1){
      final.out=1
    } else if(out1 == 0 & out2 == -1){
      final.out=-1
    } else if(out1 == 0 & out2 == 0){
      final.out=0
    }
    return(final.out)
  }

  switch(comparison,
         "survival survival"={
           compare <- function(a, b){
             #outcome Y1
             Y1a <- as.numeric(a[c("evnt1","aval1")])
             Y1b <- as.numeric(b[c("evnt1","aval1")])

             #outcome Y2
             Y2a <- as.numeric(a[c("evnt2","aval2")])
             Y2b <- as.numeric(b[c("evnt2","aval2")])

             #compare Y1 and Y2 seperately
             out1<- comp.surv(Y1a, Y1b)
             out2<- comp.surv(Y2a, Y2b)

             #get final output
             final.out<-get.final.out(out1, out2)

             return(final.out)
           }
         },
         "survival"={
           compare <- function(a, b){
             #outcome Y1
             Y1a <- as.numeric(a[c("evnt","aval")])
             Y1b <- as.numeric(b[c("evnt","aval")])
             #compare Y1
             out1<- comp.surv(Y1a, Y1b)

             return(out1)
           }
         },
         "continous"={
           compare <- function(a, b){
             #outcome Y1
             Y1a <- as.numeric(a[c("aval")])
             Y1b <- as.numeric(b[c("aval")])
             #compare Y1
             out1<- comp.conts(Y1a, Y1b)

             return(out1)
           }
         },
         "binary"={
           compare <- function(a, b){
             #outcome Y1
             Y1a <- as.numeric(a[c("aval")])
             Y1b <- as.numeric(b[c("aval")])
             #compare Y1
             out1<- comp.conts(Y1a, Y1b)

             return(out1)
           }
         },
         "ordinal"={
           compare <- function(a, b){
             #outcome Y1
             Y1a <- as.numeric(a[c("aval")])
             Y1b <- as.numeric(b[c("aval")])
             #compare Y1
             out1<- comp.conts(Y1a, Y1b)

             return(out1)
           }
         }
  )


  if(is.function(comparison)){
    compare <- comparison
  }

  n<-nrow(dat)
  #index matrix for comparison
  id<-expand.grid(i=which(info$trt01p==0),j=which(info$trt01p==1))
  n1=length(which(info$trt01p==0))
  n2=length(which(info$trt01p==1))
  # indicator values of pairwise comparison <== this step is computing intensive **
  comp.ind<-sapply(1:(n1*n2), function(x) compare(info[id[x,"i"],], info[id[x,"j"],]))
  #comp.temp<-function(i,j) compare(info[i,], info[j,])
  #comp.ind<-as.numeric(outer(1:n,1:n, "comp.temp"))

  ##----  Customized Loss and Error Function ----##

  Myloss <- function(preds, dtrain) {
    trt01p <- getinfo(dtrain, "label")

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
      # Num win in subgroup 1, a
      Nw.r1 <- sum(dat$pred[id[,"i"]] * dat$pred[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==1))
      # Num lose in subgroup 1, b
      Nl.r1 <- sum(dat$pred[id[,"i"]] * dat$pred[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==-1))

      if(Nw.r1==0) Nw.r1 <- 1 # give it a small value if treatment arm always loss
      if(Nl.r1==0) Nl.r1 <- 1 # give it a small value if treatment arm always win

      #e
      NN.r1 <- gH1.r1.dat*gH0.r1.dat
      # win diff of subgroup 1
      Rw.r1 <-  Nw.r1/NN.r1 - Nl.r1/NN.r1

      #Rw.r1.g <-(gNw.r1*Nl.r1-gNl.r1*Nw.r1)/Nl.r1^2

      # gradient of log num win in subgroup 1
      Nw.r1.g <- sapply(1:n, function(i) sum(dat$pred[which(dat$trt01p==1)] * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==1)) +
                              sum(dat$pred[which(dat$trt01p==0)] * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==1))  )
      Nl.r1.g <- sapply(1:n, function(i) sum(dat$pred[which(dat$trt01p==1)] * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==-1)) +
                              sum(dat$pred[which(dat$trt01p==0)] * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==-1))  )

      #g
      delta.e <- (dat$trt01p==arm.val[1])*sum((dat$trt01p==arm.val[2])*dat$pred)+sum((dat$trt01p==arm.val[1])*dat$pred)*(dat$trt01p==arm.val[2])


      Rw.r1.g <-  (Nw.r1.g*NN.r1 - Nw.r1* delta.e)/NN.r1^2 - (Nl.r1.g*NN.r1 -Nl.r1*delta.e)/NN.r1^2

    } else{
      Rw.r1 = 0 # do not contribute to the value function
      Rw.r1.g = 0 # do not contribute to the gradient
    }

    ## subgroup 2 --- r2 ##
    if(gH1.r2.dat > 0 & gH0.r2.dat > 0){ # both arms have subjects in subgroup 2
      # Num  win in subgroup 2
      Nw.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind== 1))
      # Num loss in subgroup 2
      Nl.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind== -1))

      if(Nw.r2==0) Nw.r2 <- 1 # give it a small value if equals 0
      if(Nl.r2==0) Nl.r2 <- 1 # give it a small value if equals 0

      #f
      NN.r2 <- gH1.r2.dat*gH0.r2.dat
      # win diff of subgroup 2
      Rw.r2 <-  Nw.r2/NN.r2 - Nl.r2/NN.r2

      # gradient of log num win in subgroup2
      #Rw.r2.g <-(gNw.r2*Nl.r2-gNl.r2*Nw.r2)/Nl.r2^2
      Nw.r2.g <- - sapply(1:n, function(i) sum((1-dat$pred[which(dat$trt01p==1)]) * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==1)) +
                                sum((1-dat$pred[which(dat$trt01p==0)]) * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==1))  )
      Nl.r2.g <- - sapply(1:n, function(i) sum((1-dat$pred[which(dat$trt01p==1)]) * (dat$trt01p[i]==arm.val[2]) * (dat$trt01p[which(dat$trt01p==1)]==arm.val[1]) * (comp.ind[id[,"i"]==i]==-1)) +
                                sum((1-dat$pred[which(dat$trt01p==0)]) * (dat$trt01p[which(dat$trt01p==0)]==arm.val[2]) * (dat$trt01p[i]==arm.val[1]) * (comp.ind[id[,"j"]==i]==-1))  )


      #h
      delta.f <- -(dat$trt01p==arm.val[1])*sum((dat$trt01p==arm.val[2])*(1-dat$pred)) - (dat$trt01p==arm.val[2])*sum((dat$trt01p==arm.val[1])*(1-dat$pred))


      Rw.r2.g <- (Nw.r2.g*NN.r2 - Nw.r2* delta.f)/NN.r2^2 - (Nl.r2.g*NN.r2 -Nl.r2*delta.f)/NN.r2^2

    } else {
      Rw.r2 = 0 # do not contribute to the value function
      Rw.r2.g = 0 # do not contribute to the gradient
    }

    g.p <- (sum(dat$pred)*Rw.r1.g + Rw.r1 - sum(1-dat$pred)*Rw.r2.g + Rw.r2)
    #h.p <- (2*rmst.diff.r1.g + sum(dat$pred)*rmst.diff.r1.h + 2*rmst.diff.r2.g - sum(1-dat$pred)*rmst.diff.r2.h)
    g <-  dat$predg*(-1)*g.p
    #h <- (-1)*( (dat$predg)^2 * h.p  + g.p*dat$predh)
    g <- g[order(dat$id)]
    #h <- h[order(dat$id)]
    h<-rep(0.00001,n)


    return(list(grad = g, hess = h))

  }



  evalerror <- function(preds, dtrain) {
    trt01p <- getinfo(dtrain, "label")

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

      #e
      NN.r1 <- gH1.r1.dat*gH0.r1.dat
      # win diff of subgroup 1
      Rw.r1 <-  Nw.r1/NN.r1 - Nl.r1/NN.r1
    } else{
      Rw.r1 = 0 # do not contribute to the value function
    }

    ## subgroup 2 --- r2 ##
    if(gH1.r2.dat > 0 & gH0.r2.dat > 0){  # both arms have subjects in subgroup 2
      # Num win in subgroup 2
      Nw.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==1))
      # Num loss in subgroup 2
      Nl.r2 <- sum((1-dat$pred)[id[,"i"]] * (1-dat$pred)[id[,"j"]] * (dat$trt01p[id[,"i"]]==arm.val[2]) * (dat$trt01p[id[,"j"]]==arm.val[1]) * (comp.ind==-1))

      if(Nw.r2==0) Nw.r2 <- 0.5 # give it a small value if equals 0
      if(Nl.r2==0) Nl.r2 <- 0.5 # give it a small value if equals 0

      #f
      NN.r2 <- gH1.r2.dat*gH0.r2.dat
      # win diff of subgroup 2
      Rw.r2 <-  Nw.r2/NN.r2 - Nl.r2/NN.r2

    } else {
      Rw.r2 = 0 # do not contribute to the value function
    }

    # err is negative of value function
    err <- (-1)*( sum(dat$pred)*Rw.r1 - sum(1-dat$pred)*Rw.r2 )

    return(list(metric = "OTR_error", value = err))
  }


  ##---- Let's boost ----##

  # Grid Search #

  hyper_grid <- expand.grid(
    #eta = c(0.001,.005, .01, .05, .1),
    eta = c(.005, .01),
    max_depth = c(2,4),
    #subsample = c(.65, .8, 1),
    #colsample_bytree = c(.8, 1),
    #lambda =c(1,3,5),
    optimal_trees = 0,
    min_error = 0
  )

  if(1){
    cat("CV to find the optimal parameter setting \n")
    for(i in 1:nrow(hyper_grid)) {
      #cat(paste0(i," eta=",hyper_grid$eta[i],", max_depth=",hyper_grid$max_depth[i],"\n"))
      print(i)
      #source("C:/Users/liupen10/OneDrive - Merck Sharp & Dohme, Corp/desktop/project/my.xgb.cv.R")
      # create parameter list #
      params <- list(
        objective = Myloss.train.diff,
        eval_metric = evalerror.test.diff,
        eta = hyper_grid$eta[i],
        max_depth = hyper_grid$max_depth[i],
        lambda = 1,
        min_child_weight = 0,
        base_score=0, #equivalent to inital prob=0.5
        colsample_bytree = 1,
        subsample=1
      )


      xgb.tune.error<-my.xgb.cv.diff(
        params=params,
        data = as.matrix(dat),
        label = info$trt01p,
        comp.ind = comp.ind,
        nrounds = 500,
        nfold = 5,
        #objective,
        #eval_metric,
        maximize = F,
        verbose = 0,               # silent,
        early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
      )
      #cat("test_error:\n")
      #print(xgb.tune.error)
      # add min training error and trees to grid
      hyper_grid$optimal_trees[i] <- which.min(xgb.tune.error)
      hyper_grid$min_error[i] <- min(xgb.tune.error)
      #print(hyper_grid[i,])

      # add min training error and trees to grid
      #hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_OTR_error_mean)
      #hyper_grid$min_error[i] <- min(xgb.tune$evaluation_log$test_OTR_error_mean)
    }

    hyper_grid <- hyper_grid[order(hyper_grid$min_error),]
    cat(" Top Five Fitting per CV \n")
    print (hyper_grid[1:5,])
  }


  ##----- Train a model based on the best CV parameter set -----##
  cat(" Train a model based on the best CV parameter set \n")
  dtrain <- xgb.DMatrix(as.matrix(dat),label = info$trt01p)
  param <- list(max_depth = hyper_grid$max_depth[1], eta = hyper_grid$eta[1], silent = 1,
                objective = Myloss, eval_metric = evalerror,verbose = 1,lambda=1,base_score=0,colsample_bytree=1,min_child_weight=0)
  watchlist <- list(train = dtrain)

  cat("Train Model based on Optimal Parameter Setting from CV \n")
  model <- xgb.train(param, dtrain, nrounds = hyper_grid$optimal_trees[1],watchlist)

  cat('Model Fitting Finished \n')

  return(model)

}
