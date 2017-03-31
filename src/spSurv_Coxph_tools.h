#ifndef _spSurv_Coxph_tools_h
#define _spSurv_Coxph_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"
#include "spSurv_spatialtools.h"

//////////////////////////////////////////////////////////////////////
// Cox PH common
/////////////////////////////////////////////////////////////////////

// Get cutpoints from Exp(hcen)
Rcpp::NumericVector Cutpoints(double hcen, int M1);

// Calculate picewise constant baseline cumulative hazard function
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
double Lambda0t(double t, Rcpp::NumericVector h, Rcpp::NumericVector d);

// Calculate picewise constant baseline cumulative hazard function at t=(t_1, ..., t_n)
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
arma::vec Lambda0tvec(Rcpp::NumericVector t, Rcpp::NumericVector h, Rcpp::NumericVector d);

// Claculate cdf of survival time given xi
double Foft(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta);

// Claculate pdf of survival time given xi
double foft(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta);

// Claculate hazard function of survival time given xi
double lambdat(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta);

//calculate Fxi^{-1}(u)
double Finvofu(double u, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta, double lower, double upper);

// Calculate  M(ti), i=1, ..., n;
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
void GetMt(Rcpp::IntegerVector& Mt, const Rcpp::NumericVector& t, Rcpp::NumericVector d);

// Calculate mk = sum_i I(M(ti)=k), k=1, ..., M with m0=0;
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
void Getmk(Rcpp::IntegerVector& mk, const Rcpp::IntegerVector& Mt);

// Calculate lk, k=1, ..., M with m0=0;
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
void Getlk(Rcpp::NumericVector& lk, const Rcpp::IntegerVector& Mt, int M1, Rcpp::NumericVector d, 
           const Rcpp::NumericVector& t, const Rcpp::NumericVector& Xbeta);

//////////////////////////////////////////////////////////////////////
// Independent Cox PH
/////////////////////////////////////////////////////////////////////
// Calculate CPO for Independent Cox PH
arma::vec LinvIndeptCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, Rcpp::NumericVector Xbeta,
                        Rcpp::NumericVector h, Rcpp::NumericVector d);

//////////////////////////////////////////////////////////////////////
// spatial Copula Cox PH
/////////////////////////////////////////////////////////////////////
// Claculate transformed survival time vector z
void Getz(arma::vec& z, const Rcpp::NumericVector& t, Rcpp::NumericVector h, Rcpp::NumericVector d, 
          const Rcpp::NumericVector& Xbeta, int n);

// Calculate CPO for spatial Copula Cox PH
arma::vec LinvSpCopulaCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, Rcpp::NumericVector Xbeta, 
                          Rcpp::NumericVector h, Rcpp::NumericVector d, arma::mat Cinv, arma::vec z);

#endif
