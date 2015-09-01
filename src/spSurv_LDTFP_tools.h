#ifndef _spSurv_LDTFP_tools_h
#define _spSurv_LDTFP_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

// logconditional density of y given xij, vi
void ldensldtfp(double y, double xbetavi, const arma::vec& xbetatfi, double sigma2, double& loglik, int maxL);

// conditional cdf of y given xij, vi
void cdfldtfp(double y, double xbetavi, const arma::vec& xbetatfi, double sigma2, double& loglik, int maxL);

// log likelihood of v_i
void loglikldtfpvi(double vi, double meanvi, double varvi, int ind1, int ind2, const Rcpp::NumericVector& y, 
      const arma::vec& Xbeta, const arma::mat& xbetatf, double sigma2, double& loglik, int maxL, 
      double vm, int indm1, int indm2);

// log likelihood given data, frailties and parameters 
void loglikldtfp(const Rcpp::NumericVector& y, const arma::vec& Xbetav, const arma::mat& Xbetatf, double sigma2, 
      Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, double& loglik, int maxL);

// initial update 
void startlrcoefldtfp(int nri, int kk, int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      arma::mat& betatf, const arma::mat& xtf, const arma::mat c0);

// function of beta
void logposldtfp(const arma::vec& betace, const arma::mat& betatf, const Rcpp::NumericVector& y, const arma::mat& xce, 
        const arma::vec& vn, const arma::mat& xtf, double sigma2, const arma::vec& betacepm, const arma::mat& betaceprec,
        Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, double& logpost, int maxL);

// update baseline reg coefficients using adaptive M-H
void update_regcoeff_adaptiveMH(arma::vec& betace, const arma::mat& betatf, const Rcpp::NumericVector& y, const arma::mat& xce, 
        const arma::vec& vn, const arma::mat& xtf, double sigma2, const arma::vec& betacepm, const arma::mat& betaceprec,
        Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, int maxL, double& rej, arma::mat& Snew, arma::vec& betabarnew, 
        int pce, int adpl0, arma::mat adpS0, double adapter, int iscan);

// update baseline reg coefficients using IWSL based M-H sampling
void update_regcoeff_iwls(arma::vec& betace, const arma::mat& betatf, const Rcpp::NumericVector& y, const arma::mat& xce, 
        const arma::vec& vn, const arma::mat& xtf, double sigma2, const arma::vec& betacepm, const arma::mat& betaceprec,
        Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, int maxL, double& rej);

// function of sigma2 
void loglikldtfpsig(const Rcpp::NumericVector& y, const arma::vec& Xbetav, const arma::mat& Xbetatf, double sigma2, 
      Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, double& loglik, int maxL, 
      double a0sig, double b0sig);

// log likelihood for updating tf logistic regressions
void compullldtfp(int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      const arma::vec& gamma, const arma::vec& xtx, const arma::mat& xtf, double& logp);

// find the mean vector m and covariance matrix C for the normal proposal used in IWLS
void tfcoeffproposal(int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      const arma::mat& xtf, const arma::mat invc0, const arma::vec& betatfkk, arma::vec& mbetatf, arma::mat& Cbetatf );

// update tf logistic regressions using IWLS based M-H
void update_tfcoeff_iwls(int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      const arma::mat& xtf, const arma::mat invc0, arma::vec& betatfkk, double& rej);

// update tf logistic regressions using slice sampling
void updatelrcoefldtfpss(int kk, int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      arma::mat& betatf, const arma::mat& xtf, const arma::mat c0);

// Calculate CPO for spatial LDTFP
arma::vec spldtfp_Linv(const Rcpp::NumericMatrix& tobs, const Rcpp::IntegerVector& type, const arma::vec& xbetav, 
      const arma::mat& xbetatf, double sigma2, int maxL);

#endif
