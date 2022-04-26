// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppSMC.h which pulls Rcpp.h in for us
#include "RcppSMC.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppSMC so that the build process will know what to do
//
// [[Rcpp::depends(RcppSMC)]]

// Simple example of computing stable calculations of the logarithm of the sum
// of weights which is of interest for computing e.g. normalising constants. The
// function returns the log-sum as well as the normalized logarithmic weights.
//
// via the exports attribute we tell Rcpp to make this function available from R
//
// [[Rcpp::export]]
Rcpp::List rcppsmc_logNormWeightsCpp(arma::vec logWeights){
  double logSum = smc::stableLogSumWeights(logWeights);
  logWeights -= logSum;
  return Rcpp::List::create(Rcpp::Named("logSum") = logSum,
                            Rcpp::Named("logNormWeights") = logWeights);
}
