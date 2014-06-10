#ifndef _spSurv_Coxph_tools_h
#define _spSurv_Coxph_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

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
           const Rcpp::NumericVector& t, const arma::vec& Xbeta);

//Sample hcen using adaptive M-H for Cox PH with random cutpoints;
void sample_hcen(double& hcen, int& rejhcen, double& hSnew, double& hbarnew, Rcpp::NumericVector h, double r0, double h0, double V0,
                 double hl0, double hs0, double hadapter, int iscan);

//////////////////////////////////////////////////////////////////////
// Independent Cox PH
/////////////////////////////////////////////////////////////////////
//Sample beta using adaptive M-H for independent Cox PH;
void indept_sample_beta(arma::vec& beta, int& rejbeta, arma::mat& Snew, arma::vec& betabarnew, const arma::mat& X, 
                 const arma::vec& Lamb0, arma::vec mu0, arma::mat Sig0, int p, double l0, arma::mat S0, 
                 double adapter, int iscan);

// Calculate CPO for Independent Cox PH
arma::vec LinvIndeptCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, arma::vec Xbeta, Rcpp::NumericVector h, Rcpp::NumericVector d);

//////////////////////////////////////////////////////////////////////
// spatial Copula Cox PH
/////////////////////////////////////////////////////////////////////
// Claculate transformed survival time vector z
void Getz(arma::vec& z, const Rcpp::NumericVector& t, Rcpp::NumericVector h, Rcpp::NumericVector d, const arma::vec& Xbeta, int n);

//Sample beta using adaptive M-H for spatial Copula Cox PH;
void spCopula_sample_beta(arma::vec& beta, int& rejbeta, arma::mat& Snew, arma::vec& betabarnew, const arma::mat& X, 
                 const arma::vec& Lamb0, arma::vec mu0, arma::mat Sig0, int p, double l0, arma::mat S0, 
                 double adapter, int iscan, arma::vec& z, const arma::mat& Cinv, 
                 const Rcpp::NumericVector& t, Rcpp::NumericVector h, Rcpp::NumericVector d, int n);

// Calculate CPO for spatial Copula Cox PH
arma::vec LinvSpCopulaCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, arma::vec Xbeta, Rcpp::NumericVector h, Rcpp::NumericVector d, 
                          arma::mat Cinv, arma::vec z);

// Get density or survival Plots for Cox PH
RcppExport SEXP CoxPHplots(SEXP xpred_, SEXP tgrid_, SEXP beta_, SEXP h_, SEXP d_, SEXP probs_);

#endif
