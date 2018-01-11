#ifndef _spSurv_spatialtools_h
#define _spSurv_spatialtools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

// inverse of (ZK^{-1}Z' + In/lambda2) and its determinant using low-rank kriging
void Inverse_kriging(const arma::mat& Z, const arma::mat& K, double lambda2, const arma::mat& In,
                     arma::mat& Cinv, double& logdetC);

// Get distance matrix for s1 s2 ... sn and s1 s2 ... sm, where s1 is column vector with coordinates
Rcpp::NumericMatrix Dist(Rcpp::NumericMatrix si, Rcpp::NumericMatrix sj);

// Get distance matrix for s1 s2 ... sn and s1 s2 ... sm, where s1 is column vector with coordinates
RcppExport SEXP DistMat(SEXP si_, SEXP sj_);

//////////////////////////////////////////////////////////////////////
// spatial density things
/////////////////////////////////////////////////////////////////////
// M-distance function
double Mdist(arma::vec x1, arma::vec x2, const arma::mat& Sinv, double phi);

// log density for new x0 given the data (X,y);
void logf_spatdens(double y0, arma::vec x0, const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                   double th1, double exp_th2, double phi, const arma::mat& Sinv, 
                   Rcpp::IntegerMatrix& kyindex, double& logf);

// log likelihood given the data (X,y);
void loglik_spatdens(const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                     double phi, const arma::mat& Sinv, Rcpp::NumericVector& lognormy, 
                     Rcpp::IntegerMatrix& kyindex, double& logf);

// log likelihood given the data (X,y) with the parametric parts removed;
void loglik_spatdens_q(const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                       double phi, const arma::mat& Sinv, Rcpp::IntegerMatrix& kyindex, double& logf);

// log posterior of y_i with the parametric part removed; used for updating censored outcomes
void logq_yi_spatdens(double y0, arma::vec x0, int indx, const Rcpp::NumericVector& y, const arma::mat& X, 
                      int J, double cpar, double th1, double exp_th2, double phi, const arma::mat& Sinv, 
                      Rcpp::IntegerMatrix& kyindex, double& logf);

// posterior prob of phi=0; used for updating phi
void prob_phi0_spatdens(const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                        const Rcpp::NumericVector& phiseq, const arma::mat& Sinv,
                        Rcpp::IntegerMatrix& kyindex, double q0phi, double a0phi, double b0phi,
                        double& ratio);

//////////////////////////////////////////////////////////////////////
// spatial Copula things
/////////////////////////////////////////////////////////////////////
// make transformation on theta: log(theta1/(1.0-theta1)); log(theta2);
arma::vec trans_theta(arma::vec theta);

// transform back to get theta;
arma::vec trans_theta_inv(arma::vec trans);

//Sample theta using adaptive M-H for spatial Copula Model using FSA;
void spCopula_sample_theta_FSA(arma::vec& theta, double& rejtheta, arma::mat& spSnew, arma::vec& thetabarnew, 
                               arma::mat& Cinv, double& logdetC, double theta1a, double theta1b, 
                               double theta2a, double theta2b, double spl0, arma::mat spShat, const arma::mat& dnn, 
                               double spadapter, int iscan, int nburn, const arma::vec& z, int n, const arma::mat& dnm, 
                               const arma::mat& dmm, const arma::mat clustindx);

#endif
