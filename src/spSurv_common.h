#ifndef _spSurv_common_h
#define _spSurv_common_h

#include <RcppArmadillo.h>

#define ESMALL 1e-10  /* small number */
#define ELARGE 1e+10 /* large number */
typedef Rcpp::NumericMatrix::iterator mat_iterator;
using namespace Rcpp;

//Truncated normal N(y;mu,s)I(a<y<b)
double rtexp(double a, double b); // Rejection algorithm with a truncated expoential proposal for N(0,1)I(a<x<b) when a is very large: |a| < |b|
double trun_rnorm(const double mu, const double s, double a, double b);

// generate multivariate normal (mu, sigma)
arma::vec mvrnorm(arma::vec mu, arma::mat sigma);

// generate Wishart random matrices
arma::mat rwish(arma::mat Sig, int n);

// sample(Nseq, prob=w), where Nseq is n-dim vector, and w.siz()=n
int sample(Rcpp::IntegerVector Nseq, Rcpp::NumericVector w);

// Get distance matrix for s1 s2 ... sn and s1 s2 ... sm, where s1 is column vector with coordinates
Rcpp::NumericMatrix Dist(Rcpp::NumericMatrix si, Rcpp::NumericMatrix sj);

// calculate qnorm(x) for a vector of x
arma::vec qnormvec(arma::vec x);

// calculate pnorm(x) for a vector of x
arma::vec pnormvec(arma::vec x);

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

// Get distance matrix for s1 s2 ... sn and s1 s2 ... sm, where s1 is column vector with coordinates
RcppExport SEXP DistMat(SEXP si_, SEXP sj_);

#endif
