#ifndef _spSurv_common_h
#define _spSurv_common_h

#include <RcppArmadillo.h>

#define ESMALL 1e-10  /* small number */
#define ELARGE 1e+10 /* large number */
#define SYSMIN 1e-305  /* small number */
#define SYSMAX 1e+305 /* large number */
#define LOGSYSMAX 702.28845336318397585 /* large number */
#define LOGSYSMIN -702.28845336318397585 /* large number */
typedef Rcpp::NumericMatrix::iterator mat_iterator;
using namespace Rcpp;

//Truncated normal N(y;mu,s)I(a<y<b)
double rtexp(double a, double b); // Rejection algorithm with a truncated expoential proposal for N(0,1)I(a<x<b) when a is very large: |a| < |b|
double trun_rnorm(const double mu, const double s, double a, double b);

// generate multivariate normal (mu, sigma)
arma::vec mvrnorm(arma::vec mu, arma::mat sigma);

// density of multivariate normal (mu, sigma)
double mvdnorm(arma::vec x, arma::vec mu, arma::mat sigma, bool logd);

// generate Wishart random matrices
arma::mat rwish(arma::mat Sig, int n);

// sample(Nseq, prob=w), where Nseq is n-dim vector, and w.siz()=n
int sample(Rcpp::IntegerVector Nseq, Rcpp::NumericVector w);

// calculate qnorm(x) for a vector of x
arma::vec qnormvec(arma::vec x);

// calculate pnorm(x) for a vector of x
arma::vec pnormvec(arma::vec x);

// Matern correlation function
double rho_Matern(double dis, double nu, double phi);

// process convolution bivariate Gaussian kernel
double kernel_G(double dis, double phi);

// Powered Exponential
double pow_exp(double dis, double phi, double nu);

#endif
