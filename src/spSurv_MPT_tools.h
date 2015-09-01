#ifndef _spSurv_MPT_tools_h
#define _spSurv_MPT_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

// from conditional Ys to cumulative probs
void Ys_to_probs(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& probs, int maxL);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AFTlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log Survival function of t given xi
double AFTSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log likelihood given data, frailties and parameters 
void AFTloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll);

// log likelihood given frailties, parameters and data of block i
void AFTloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi);

// Calculate 1.0/likelihood for CPO
arma::vec AFTinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2);


/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double POlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log Survival function of t given xi
double POSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log likelihood given data, frailties and parameters 
void POloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll);

// log likelihood given frailties, parameters and data of block i
void POloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi);

// Calculate 1.0/likelihood for CPO
arma::vec POinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PHlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log Survival function of t given xi
double PHSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log likelihood given data, frailties and parameters 
void PHloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll);

// log likelihood given frailties, parameters and data of block i
void PHloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi);

// Calculate 1.0/likelihood for CPO
arma::vec PHinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AHlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log Survival function of t given xi
double AHSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2);

// log likelihood given data, frailties and parameters 
void AHloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll);

// log likelihood given frailties, parameters and data of block i
void AHloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi);

// Calculate 1.0/likelihood for CPO
arma::vec AHinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2);


#endif
