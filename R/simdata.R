#' Simulated training and testing data of two survival outcomes
#'
#' A list of simulated data, which include 600 patients for training in \code{train}, 5000 patients for testing in {test},
#' and their true subgroup membership labels in \code{label.train} and \code{label.test}. In training and testing set, predictive variable
#' 's1' controls the treatment benefit, patients with positive 's1' benefit from treatment, and patients with negative 's1' do not benefit
#' from treatment.Variables 'z1' to 'z20' are confounding baseline variables, 'trt01p' denote treatment assignment, 'evnt1' and 'aval1'
#' are the event indicator(0 censor, 1 event) and survival outcome for the primary outcome overal survival (OS). 'evnt2' and 'aval2' are
#' the event indicator(0 censor, 1 event) and survival outcome for the secondary outcome time to prgression (TTP). \code{label.train} and
#' \code{label.test} denotes the true subgroup membership (0 tretment nonperforming subgroup, 1 treatment performing subgoup).
#'
#' @docType data
#'
#' @usage data(simdata)
#'
#' @format a list containing training data \code{train}, testing data \code{test}, subgroup label \code{label.train}, and subgroup label \code{label.test}
#'
#' @keywords datasets
#'
#' @examples NULL
"simdata"
