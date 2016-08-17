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
// spatial Copula things
/////////////////////////////////////////////////////////////////////
// Preprocess C^{-1} to get Cinv using FSA
void GetCinv_FSA(int n, double theta1, double theta2, const arma::mat& dnn, const arma::mat& dnm, const arma::mat& dmm, 
                const Rcpp::IntegerVector& blocki, arma::mat& Cinv, double& logdetC);

// Preprocess C^{-1} to get Cinv directly
void GetCinv(int n, double theta1, double theta2, const arma::mat& dnn, arma::mat& Cinv, double& logdetC);

// make transformation on theta: log(theta1/(1.0-theta1)); log(theta2);
arma::vec trans_theta(arma::vec theta);

// transform back to get theta;
arma::vec trans_theta_inv(arma::vec trans);

//Sample theta using adaptive M-H for spatial Copula Model;
void spCopula_sample_theta(arma::vec& theta, int& rejtheta, arma::mat& spSnew, arma::vec& thetabarnew, arma::mat& Cinv, double& logdetC, 
                 double theta1a, double theta1b, double theta2a, double theta2b, double spl0, arma::mat spS0, const arma::mat& dnn, 
                 double spadapter, int iscan, const arma::vec& z, int n);

//Sample theta using adaptive M-H for spatial Copula Model using FSA;
void spCopula_sample_theta_FSA(arma::vec& theta, int& rejtheta, arma::mat& spSnew, arma::vec& thetabarnew, arma::mat& Cinv, double& logdetC, 
                 double theta1a, double theta1b, double theta2a, double theta2b, double spl0, arma::mat spS0, const arma::mat& dnn, 
                 double spadapter, int iscan, const arma::vec& z, int n, const arma::mat& dnm, const arma::mat& dmm, 
                const Rcpp::IntegerVector& blocki);

#endif
