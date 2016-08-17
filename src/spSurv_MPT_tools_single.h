#ifndef _spSurv_MPT_tools_single_h
#define _spSurv_MPT_tools_single_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"
#include "spSurv_MPT_tools_super.h"

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PHlogpdf(double t, double th1, double th2, Rcpp::NumericVector probs,
                int maxL, bool MPT, int dist, double xibeta);
// log survival function of t given xi
double PHlogsurv(double t, double th1, double th2, Rcpp::NumericVector probs, 
                 int maxL, bool MPT, int dist, double xibeta);
// log cdf of t given xi
double PHlogcdf(double t, double th1, double th2, Rcpp::NumericVector probs, 
                int maxL, bool MPT, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double PHlogsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector probs, 
                     int maxL, bool MPT, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void PHloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
              const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
              int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void PHloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                    const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
                    int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll, 
                    int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec PHinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs, 
                   int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double POlogpdf(double t, double th1, double th2, Rcpp::NumericVector probs,
                int maxL, bool MPT, int dist, double xibeta);
// log survival function of t given xi
double POlogsurv(double t, double th1, double th2, Rcpp::NumericVector probs, 
                 int maxL, bool MPT, int dist, double xibeta);
// log cdf of t given xi
double POlogcdf(double t, double th1, double th2, Rcpp::NumericVector probs, 
                int maxL, bool MPT, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double POlogsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector probs, 
                     int maxL, bool MPT, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void POloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
              const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
              int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void POloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                    const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
                    int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll, 
                    int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec POinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs, 
                   int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AFTlogpdf(double t, double th1, double th2, Rcpp::NumericVector probs,
                 int maxL, bool MPT, int dist, double xibeta);
// log survival function of t given xi
double AFTlogsurv(double t, double th1, double th2, Rcpp::NumericVector probs, 
                  int maxL, bool MPT, int dist, double xibeta);
// log cdf of t given xi
double AFTlogcdf(double t, double th1, double th2, Rcpp::NumericVector probs, 
                 int maxL, bool MPT, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double AFTlogsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector probs, 
                      int maxL, bool MPT, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void AFTloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
               const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
               int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void AFTloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                     const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
                     int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll, 
                     int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec AFTinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                    const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs, 
                    int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta);

#endif
