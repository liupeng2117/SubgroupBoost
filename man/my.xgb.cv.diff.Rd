% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my.xgb.cv.diff.R
\name{my.xgb.cv.diff}
\alias{my.xgb.cv.diff}
\title{Title Cross-validation for evaluating the performance of subgrouping using win-difference based value function}
\usage{
my.xgb.cv.diff(
  params,
  data,
  label,
  comp.ind,
  nrounds,
  nfold,
  maximize,
  verbose,
  early_stopping_rounds
)
}
\arguments{
\item{params}{the list of parameters for xgboost}

\item{data}{the data of predictors}

\item{label}{treatment label of subjects}

\item{comp.ind}{The vector of comparision results}

\item{nrounds}{the maximum number of boosting rounds.}

\item{nfold}{number of folds in cross validation}

\item{maximize}{a logical value, if \code{TRUE}, maximize the loss function, here we set it to \code{FALSE}}

\item{verbose}{a logical value, whether print out the messages of fitting the model}

\item{early_stopping_rounds}{the number of rounds before stopping if no improvement, defualt is 5}
}
\value{
a vector of errors
}
\description{
Title Cross-validation for evaluating the performance of subgrouping using win-difference based value function
}
\examples{
NULL
}
