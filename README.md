# SubgroupBoost

This tutorial guide the reader through the analysis of treatment performing and non-performing subgroup identification using our package SubgroupBoost. 

## Installation

You can install the released version of SubgroupBoost from GitHub with:

``` r
devtools::install_github("liupeng2117/SubgroupBoost")
```

## Data input and preprocessing

```r
library(SubgroupBoost)
```

* Load data. In this tutorial, we will use a simulated data that comes with SubgroupBoost package to demonstrate. 

```r
data("simdata")
names(simdata)
str(simdata$train)
str(simdata$test)
```

* The data is a list containing 4 elements: training data, testing data, true subgroup label for training data, and true subgroup label for testing data. The training data contains 600 samples, and testing data contains 5000 samples. The data columns has 26 variables: `trt01p` is treatment indicator. `evnt1` and `aval1` are for the primary outcome overall survival (OS). `evnt1` is the indicator of event (death) for OS, and `aval1` is the time value. `evnt2` and `aval2` are for the secondary outcome time-to-progression (TTP), where `evnt2` is the indicator of event (disease progression), and `aval2` is the time value. Variable `s1` is the predictive variable that interact with treatment, while variables z1 to z20 are not related to treatment effect.

* Preprocess the training and testing data
```r
#----- data manipulation-----#
#training data
dat.train=simdata[[1]][,-c(1:5)]
info.train=simdata[[1]][,1:5]
dtrain <- xgb.DMatrix(as.matrix(dat.train),label = info.train$trt01p)

#testing data
dat.test=simdata[[2]][,-c(1:5)]
info.test=simdata[[2]][,1:5]
dtest <- xgb.DMatrix(as.matrix(dat.test),label = info.test$trt01p)
```

* Before model fitting, we define a function to predict the subgroup assignment. In this function, we first calculate soft assignment probabilities by a sigmoid transformation of predicted value, then seperate the samples into two groups by 0.5 cutoff.  
```r
SubgroupBoost.predict<-function(model, data){
  # obtain subgroup assignment for dtrain
  pred<-predict(model, data)
  prob<-1/(1+exp(-pred))
  assign<-ifelse(prob>0.5,1,0)
  return(group_assign=assign)
}

```

## RMST
* Use RMST based gradient tree boosting method to identify subgroups of patients. Here we choose OS as input, as RMST based method can only take one survival outcome.
```r
#----- RMST death -----#
dat1=data.frame(dat.train,trt01p=info.train$trt01p, evnt=info.train$evnt1, aval=info.train$aval1)
model_RMST<-SubgroupBoost.RMST(dat1)

#check variables selected
importance=xgb.importance(model_RMST$feature_names, model_RMST)$Feature
importance

# [1] "s1"  "z13" "z20" "z18" "z6"  "z15" "z11" "z7"  "z9"  "z4"  "z3"  "z2"  "z10" "z8"  "z19" "z1" 

#calculate prediction accuracy for training data and testing data
RMST.predict.train<-SubgroupBoost.predict(model_RMST, dtrain)
table(RMST.predict.train, simdata[[3]])

#RMST.predict.train   0   1
#                 0 260  79
#                 1  30 231
                 
RMST.predict.test<-SubgroupBoost.predict(model_RMST, dtest)
table(RMST.predict.test, simdata[[4]])

#RMST.predict.test    0    1
#                0 2212  504
#                1  269 2015
```

## Win-difference
* Use win-difference based gradient tree boosting method to identify subgroups of patients.
```r
#----- win-ratio ----#
comparison="survival survival"
model_wd<-SubgroupBoost.wd(dat, info, comparison)

#check variables selected
importance=xgb.importance(model_wd$feature_names, model_wd)$Feature
importance

#calculate prediction accuracy for training data and testing data
wd.predict.train<-SubgroupBoost.predict(model_wd, dtrain)
table(wd.predict.train, simdata[[3]])

wd.predict.test<-SubgroupBoost.predict(model_wd, dtest)
table(wd.predict.test, simdata[[4]])

```

