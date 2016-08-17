#ifndef _spSurv_MPT_tools_super_h
#define _spSurv_MPT_tools_super_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

// from conditional Ys to cumulative probs
void Ys_to_probs(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& probs, int maxL);
void Ys_to_probs2(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& probs, int maxL);

// loglogistic baseline survival functions
double S0MPT(double t, double th1, double th2, Rcpp::NumericVector probs, int maxL, bool MPT, int dist);
double logf0MPT(double t, double th1, double th2, Rcpp::NumericVector probs, int maxL, bool MPT, int dist);
double logh0MPT(double t, double th1, double th2, Rcpp::NumericVector probs, int maxL, bool MPT, int dist);

#endif
