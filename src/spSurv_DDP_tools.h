#ifndef _spSurv_DDP_tools_h
#define _spSurv_DDP_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"
#include "spSurv_spatialtools.h"

//////////////////////////////////////////////////////////////////////
// DDP common
/////////////////////////////////////////////////////////////////////
// Caculate CDF at points ygrid
arma::vec Fmix(arma::vec ygrid, arma::vec mu, arma::vec sig, arma::vec w);

// Caculate pdf at points ygrid
arma::vec fmix(arma::vec ygrid, arma::vec mu, arma::vec sig, arma::vec w);

// calculate marginal CDF of y: a mixture of normals
double Fofy(double y, arma::vec w, arma::vec mu, arma::vec sig);

// calculate marginal pdf of y: a mixture of normals
double fofy(double y, arma::vec w, arma::vec mu, arma::vec sig);

// calculate F^{-1}(u)
double DDP_Finvofu(double u, arma::vec wma, arma::vec mu, arma::vec sig, double lower, double upper);

// From V to w
void DDP_Vtow(arma::vec& w, Rcpp::NumericVector V, int N);

// sample(1:N, prob=w), where w.size()=N
int DDP_sample(arma::vec w);

//Sample K;
void DDP_sample_K(Rcpp::IntegerVector& K, const Rcpp::NumericVector& y, const arma::mat& Xbeta, 
                  arma::vec w, Rcpp::NumericVector tau, int n, int N);

// calculate pnorm((y-mu)/sig) for vectors of mu and sig
arma::vec Phivec(double y, arma::vec mu, NumericVector sig);

//////////////////////////////////////////////////////////////////////
// ANOVA DDP
/////////////////////////////////////////////////////////////////////
// Sample beta;
void anovaDDP_sample_beta(arma::mat& beta, const Rcpp::NumericVector& y, const arma::mat& X, 
                     const Rcpp::NumericVector& tau2, const Rcpp::IntegerVector& nK, 
                     const Rcpp::IntegerMatrix& Kind, arma::vec mu, arma::mat Sig, 
                     arma::mat invSig, int N, int p);

//Sample simga2;
void anovaDDP_sample_sigma2(Rcpp::NumericVector& tau2, const Rcpp::NumericVector& y, const arma::mat& Xbeta, 
                     const Rcpp::IntegerVector& nK, const Rcpp::IntegerMatrix& Kind, double nua, double nub, int N);

// Calculate CPO for anovaDDP
arma::vec anovaDDP_Linv(Rcpp::NumericVector yobs, Rcpp::IntegerVector delta, arma::mat X, arma::mat beta, arma::vec sig, arma::vec wei );

//////////////////////////////////////////////////////////////////////
// spatial Copula DDP
/////////////////////////////////////////////////////////////////////
// Sample y_i when delta_i=0
void spCopula_sample_y(Rcpp::NumericVector& y, Rcpp::NumericVector& rejy, arma::mat& zPhi, arma::vec& z, arma::vec w, 
                       const Rcpp::NumericVector& yobs, const Rcpp::IntegerVector& delta, const arma::mat& Xbeta, Rcpp::NumericVector tau, 
                       Rcpp::IntegerVector K, const arma::mat& Cinv, int n, int N, int iscan, int nburn);

// Sample beta;
void spCopula_sample_beta(arma::mat& beta, Rcpp::NumericVector& rejbeta, arma::mat& zPhi, arma::vec& z, arma::vec w, 
                          const Rcpp::NumericVector& y, const arma::mat& X, Rcpp::NumericVector tau2, 
                          const Rcpp::IntegerVector& nK, const Rcpp::IntegerMatrix& Kind, arma::vec mu, 
                          arma::mat Sig, arma::mat invSig, const arma::mat& Cinv, int n, int N, int p, int iscan, int nburn);

// Sample beta blockwise;
void spCopula_sample_beta_block(arma::mat& beta, Rcpp::NumericVector& rejbeta, arma::mat& zPhi, arma::vec& z, arma::vec w, 
      const Rcpp::NumericVector& y, const arma::mat& X, Rcpp::NumericVector tau2, const Rcpp::IntegerVector& nK, 
      const Rcpp::IntegerMatrix& Kind, arma::vec mu, arma::mat Sig, arma::mat invSig, const arma::mat& Cinv, int n, int N, int p);

//Sample simga2;
void spCopula_sample_sigma2(Rcpp::NumericVector& tau2, Rcpp::NumericVector& rejsigma, arma::mat& zPhi, arma::vec& z, 
                            arma::vec w, const Rcpp::NumericVector& y, const arma::mat& Xbeta, const Rcpp::IntegerVector& nK, 
                            const Rcpp::IntegerMatrix& Kind, double nua, double nub, const arma::mat& Cinv, int n, int N, 
                            int iscan, int nburn);

//Sample simga2 blockwise;
void spCopula_sample_sigma2_block(Rcpp::NumericVector& tau2, Rcpp::NumericVector& rejsigma, arma::mat& zPhi, arma::vec& z, 
      arma::vec w, const Rcpp::NumericVector& y, const arma::mat& Xbeta, const Rcpp::IntegerVector& nK, 
      const Rcpp::IntegerMatrix& Kind, double nua, double nub, const arma::mat& Cinv, int n, int N);

//Sample V;
void spCopula_sample_V(Rcpp::NumericVector& V, Rcpp::NumericVector& rejV, arma::mat& zPhi, arma::vec& z, arma::vec& w, 
                       const Rcpp::IntegerVector& nK, double alpha, const arma::mat& Cinv, int N, int iscan, int nburn);

//Sample V blockwise;
void spCopula_sample_V_block(Rcpp::NumericVector& V, Rcpp::NumericVector& rejV, arma::mat& zPhi, arma::vec& z, arma::vec& w, 
      const Rcpp::IntegerVector& nK, double alpha, const arma::mat& Cinv, int N);

// Calculate CPO for spatial Copula DDP
arma::vec spCopula_Linv(Rcpp::NumericVector yobs, Rcpp::IntegerVector delta, arma::mat X, arma::mat beta, arma::vec sig, arma::vec wei, 
                    arma::mat Cinv, arma::vec zj);

// using anovaDDP to get initial chain
void spCopulaInitial(Rcpp::IntegerVector& K, Rcpp::NumericVector& y, arma::mat& beta, Rcpp::NumericVector& tau2, Rcpp::NumericVector& V,
                  arma::vec& w, double& alpha, arma::vec& mu, arma::mat& Sig, arma::mat& invSig, const Rcpp::NumericVector& yobs,
                  const Rcpp::IntegerVector& delta, const arma::mat& X, arma::vec m0, arma::mat S0, arma::mat Sig0, int k0, 
                  double a0, double b0, double nua, double nub, arma::mat invS0, arma::vec invS0m0 );

// Get density or survival Plots for DDP
RcppExport SEXP DDPplots(SEXP xpred_, SEXP tgrid_, SEXP beta_, SEXP sigma2_, SEXP w_, SEXP CI_);

// Spatial Maps for transformed spatial process using LDDPM-spatial model
RcppExport SEXP PredMapsZ(SEXP ds0n_, SEXP dnn_, SEXP beta_, SEXP sigma2_, SEXP w_, SEXP theta1_, SEXP theta2_, SEXP z_);

#endif
