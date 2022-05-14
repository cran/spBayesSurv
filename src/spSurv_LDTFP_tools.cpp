#include "spSurv_LDTFP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// log conditional density of y given xij, vi
void ldensldtfp(double y, double xbetavi, const arma::vec& xbetatfi, double sigma2, double& loglik, int maxL){
  //temp variables
  int maxL1 = maxL+1;
  Rcpp::IntegerVector K(maxL1);
  // find where the data point is located at each level
  double sigma = std::sqrt(sigma2);
  double tmp2 = xbetavi;
  loglik = Rf_dnorm4(y, tmp2, sigma, true);
  double tmp1 = (y-tmp2)/sigma;
  if (tmp1>4.0) {
    tmp2 = 0.999968;
  }else if(tmp1<-4.0){
    tmp2 = 0.000032;
  }else {
    tmp2 = Rf_pnorm5(y, tmp2, sigma, true, false);
  }
  for(int j=0; j<maxL1; ++j){
    K[j] = (int)(std::pow(2,j)*tmp2)+1;
  }
  
  // ldtfp part
  int m=0; 
  for(int j1=1; j1<maxL1; ++j1){ 
    int k2 = m + K[j1-1];
    m += std::pow(2,j1-1);
    int ll = (K[j1-1]-1)*2 +1;
    tmp2 = xbetatfi[k2-1];
    tmp1 = std::exp(tmp2)/(1.0+std::exp(tmp2));
    if(K[j1]==ll) {
      loglik += std::log(tmp1);
    } else {
      tmp1 = 1.0 - tmp1; loglik += std::log(tmp1);
    }
  }
  loglik += (double)(maxL+0.0)*std::log(2.0);
}

// conditional cdf of y given xij, vi
void cdfldtfp(double y, double xbetavi, const arma::vec& xbetatfi, double sigma2, double& cdfval, int maxL){
  //R_CheckUserInterrupt();
  //temp variables
  int maxL1 = maxL+1;
  int nsets1 = std::pow(2, maxL)+1;
  Rcpp::IntegerVector K(maxL1);
  Rcpp::NumericVector prob(nsets1);
  Rcpp::NumericVector probc(nsets1);
  for(int i=0; i<maxL1; ++i) K[i]=1;
  for(int i=0; i<nsets1; ++i) {
    prob[i] = 0;
    probc[i] = 0;
  }
  double loglik=0;

  // find where the data point is located at each level
  double sigma = std::sqrt(sigma2);
  double tmp2 = xbetavi;
  double tmp1 = (y-tmp2)/sigma;
  if (tmp1>4.0) {
    tmp2 = 0.999968;
  }else if(tmp1<-4.0){
    tmp2 = 0.000032;
  }else {
    tmp2 = Rf_pnorm5(y, tmp2, sigma, true, false);
  }
  double cdfw = tmp2;
  double possi = (int)(std::pow(2,maxL)*tmp2)+1;
  
  // ldtfp part
  double accsum = 0;
  for(int i=1; i<nsets1; ++i){
    tmp2 = (i-1.0)/(nsets1-1.0);
    for(int j1=1; j1<=maxL; ++j1) K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    loglik=0;
    int m=0; 
    for(int j1=1; j1<=maxL; ++j1){ 
      int k2 = m + K[j1-1];
      m += std::pow(2,j1-1);
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = xbetatfi[k2-1];
      tmp1 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglik += std::log(tmp1);
      } else {
        tmp1 = 1.0 - tmp1; loglik += std::log(tmp1);
      }
    }
    prob[i] = std::exp(loglik);
    accsum += prob[i];
  }
  for(int i=1; i<nsets1; ++i){
    prob[i] = prob[i]/accsum;
    probc[i] = prob[i] + probc[i-1];
  }
  cdfval = probc[possi-1] + prob[possi]*((nsets1-1.0)*cdfw-(possi-1.0));
  if(cdfval>1.0) cdfval=1.0;
  if(cdfval<0.0) cdfval=0.0;
}

// log likelihood of v_i
void loglikldtfpvi(double vi, double meanvi, double varvi, int ind1, int ind2, const Rcpp::NumericVector& y, 
      const arma::vec& Xbeta, const arma::mat& xbetatf, double sigma2, double& loglik, int maxL, 
      double vm, int indm1, int indm2){
  // temp variables
  int maxL1 = maxL+1;
  Rcpp::IntegerVector K(maxL1);
  double sigma = std::sqrt(sigma2);
  double tmp1, tmp2, tmp3;
  loglik=0;
  
  // find where the data point is located at each level
  for(int i=ind1; i<=ind2; ++i){
    //R_CheckUserInterrupt();
    tmp1 = Xbeta[i] + vi;
    tmp3 = (y[i]-tmp1)/sigma;
    loglik += Rf_dnorm4(y[i], tmp1, sigma, true);
    if (tmp3>4.0) {
      tmp2 = 0.999968;
    } else if(tmp3<-4.0) {
      tmp2 = 0.000032;
    } else {
      tmp2 = Rf_pnorm5(y[i], tmp1, sigma, true, false);
    }
    for(int j1=0; j1<maxL1; ++j1){
      K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    }
  
    // ldtfp part
    int m=0; 
    for(int j1=1; j1<maxL1; ++j1){
      int k2 = m + K[j1-1];
      m = m + std::pow(2,j1-1);
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = xbetatf(k2-1,i);
      tmp3 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglik += std::log(tmp3);
      } else {
        tmp3 = 1.0 - tmp3; 
        loglik += std::log(tmp3);
      }
    }
    loglik += (maxL1-1.0)*std::log(2.0);
  }
  for(int i=indm1; i<=indm2; ++i){
    tmp1 = Xbeta[i] + vm;
    tmp3 = (y[i]-tmp1)/sigma;
    loglik += Rf_dnorm4(y[i], tmp1, sigma, true);
    if (tmp3>4.0) {
      tmp2 = 0.999968;
    } else if(tmp3<-4.0) {
      tmp2 = 0.000032;
    } else {
      tmp2 = Rf_pnorm5(y[i], tmp1, sigma, true, false);
    }
    for(int j1=0; j1<maxL1; ++j1){
      K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    }
  
    // ldtfp part
    int m=0; 
    for(int j1=1; j1<maxL1; ++j1){
      int k2 = m + K[j1-1];
      m = m + std::pow(2,j1-1);
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = xbetatf(k2-1,i);
      tmp3 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglik += std::log(tmp3);
      } else {
        tmp3 = 1.0 - tmp3; 
        loglik += std::log(tmp3);
      }
    }
    loglik += (maxL1-1.0)*std::log(2.0);
  }
  loglik += -0.5*std::pow((vi-meanvi), 2)/varvi;
}

// log likelihood of for v_i
void loglikldtfpvi2(double vi, int ind1, int ind2, const Rcpp::NumericVector& y, 
                   const arma::vec& Xbeta, const arma::mat& xbetatf, double sigma2, 
                   double& loglik, int maxL){
  // temp variables
  int maxL1 = maxL+1;
  Rcpp::IntegerVector K(maxL1);
  double sigma = std::sqrt(sigma2);
  double tmp1, tmp2, tmp3;
  loglik=0;
  
  // find where the data point is located at each level
  for(int i=ind1; i<=ind2; ++i){
    tmp1 = Xbeta[i] + vi;
    tmp3 = (y[i]-tmp1)/sigma;
    loglik += Rf_dnorm4(y[i], tmp1, sigma, true);
    if (tmp3>4.0) {
      tmp2 = 0.999968;
    } else if(tmp3<-4.0) {
      tmp2 = 0.000032;
    } else {
      tmp2 = Rf_pnorm5(y[i], tmp1, sigma, true, false);
    }
    for(int j1=0; j1<maxL1; ++j1){
      K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    }
    
    // ldtfp part
    int m=0; 
    for(int j1=1; j1<maxL1; ++j1){
      int k2 = m + K[j1-1];
      m = m + std::pow(2,j1-1);
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = xbetatf(k2-1,i);
      tmp3 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglik += std::log(tmp3);
      } else {
        tmp3 = 1.0 - tmp3; 
        loglik += std::log(tmp3);
      }
    }
    loglik += (maxL1-1.0)*std::log(2.0);
  }
}

// log likelihood given data, frailties and parameters 
void loglikldtfp(const Rcpp::NumericVector& y, const arma::vec& Xbetav, const arma::mat& xbetatf, double sigma2, 
      Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, double& loglik, int maxL){
  
  // temp variables
  int maxL1 = maxL+1;
  Rcpp::IntegerVector K(maxL1);
  int nrec = y.size();
  double sigma = std::sqrt(sigma2);
  double tmp1, tmp2, tmp3;
  for(int i=0; i<nobsbc.size(); ++i) nobsbc(i) = 0;
  loglik=0;
  
  // find where the data point is located at each level
  for(int i=0; i<nrec; ++i){
    //R_CheckUserInterrupt();
    tmp1 = Xbetav[i];
    tmp3 = (y[i]-tmp1)/sigma;
    loglik += Rf_dnorm4(y[i], tmp1, sigma, true);
    if (tmp3>4.0) {
      tmp2 = 0.999968;
    } else if(tmp3<-4.0) {
      tmp2 = 0.000032;
    } else {
      tmp2 = Rf_pnorm5(y[i], tmp1, sigma, true, false);
    }
    for(int j1=0; j1<maxL1; ++j1){
      K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    }
  
    // ldtfp part
    int j=0; int m=0; 
    for(int j1=1; j1<maxL1; ++j1){
      int j2 = j + K[j1];
      int k2 = m + K[j1-1];
      j = j + std::pow(2,j1);
      m = m + std::pow(2,j1-1);
      obsbc(j2-1, nobsbc[j2-1])=i;
      nobsbc[j2-1] += 1;
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = xbetatf(k2-1,i);
      tmp3 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglik += std::log(tmp3);
      } else {
        tmp3 = 1.0 - tmp3; 
        loglik += std::log(tmp3);
      }
    }
    loglik += (maxL1-1.0)*std::log(2.0);
  }
}

// initial update 
void startlrcoefldtfp(int nri, int kk, int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      arma::mat& betatf, const arma::mat& xtf, const arma::mat c0){
  //temp variables
  int ptf = c0.n_cols;
  arma::mat xtx(ptf,ptf);
  arma::vec xty(ptf); 
  arma::mat c0inv = arma::inv_sympd(c0);
  arma::vec beta = betatf.col(kk);
  
  // MLE with nri step of N-R
  for(int nhr=0; nhr<nri; ++nhr){
    xtx = c0inv;
    xty.fill(0.0);
    
    if(n1>0){
      for(int i=0; i<n1; ++i){
        int ll = obsbc(ii,i);
        //R_CheckUserInterrupt();
        arma::vec xtfll = xtf.col(ll);
        double eta = arma::dot(xtfll, beta);
        double mu = std::exp(eta)/(1+std::exp(eta));
        double tmp1 = std::exp(eta)/std::pow((1+std::exp(eta)),2);
        double ytilde = eta+(1-mu)/tmp1;
        xtx += xtfll*(arma::trans(xtfll))*tmp1;
        xty += xtfll*ytilde*tmp1;
      }
    }
    if(n2>0){
      for(int i=0; i<n2; ++i){
        int ll = obsbc(jj,i);
        //R_CheckUserInterrupt();
        arma::vec xtfll = xtf.col(ll);
        double eta = arma::dot(xtfll, beta);
        double mu = std::exp(eta)/(1+std::exp(eta));
        double tmp1 = std::exp(eta)/std::pow((1+std::exp(eta)),2);
        double ytilde = eta+(0-mu)/tmp1;
        xtx += xtfll*(arma::trans(xtfll))*tmp1;
        xty += xtfll*ytilde*tmp1;
      }
    }
    xtx = arma::inv_sympd(xtx);
    beta = xtx*xty;
  }
  betatf.col(kk) = beta;
}

// function of beta
void logposldtfp(const arma::vec& betace, const arma::mat& betatf, const Rcpp::NumericVector& y, const arma::mat& xce, 
        const arma::vec& vn, const arma::mat& xtf, double sigma2, const arma::vec& betacepm, const arma::mat& betaceprec,
        Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, double& logpost, int maxL){
  // temp variables
  int maxL1 = maxL+1;
  Rcpp::IntegerVector K(maxL1);
  int nrec = y.size();
  double sigma = std::sqrt(sigma2);
  double tmp1, tmp2, tmp3;
  for(int i=0; i<nobsbc.size(); ++i) nobsbc(i) = 0;
  double loglikn = 0;
  for(int i=0; i<nrec; ++i){
    //R_CheckUserInterrupt();
    tmp1=arma::dot(xce.col(i),betace) + vn[i];
    tmp3=(y[i]-tmp1)/sigma;
    loglikn += Rf_dnorm4(y[i], tmp1, sigma, true);
    if(tmp3>4.0){
      tmp2 = 0.999968;
    } else if(tmp3<-4.0) {
      tmp2 = 0.000032;
    } else {
      tmp2 = Rf_pnorm5(y[i], tmp1, sigma, true, false);
    }
    for(int j1=0; j1<maxL1; ++j1) K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    
    int j=0; int m=0;
    for(int j1=1; j1<maxL1; ++j1){
      int j2 = j + K[j1];
      int k2 = m + K[j1-1];
      j = j + std::pow(2,j1);
      m = m + std::pow(2,j1-1);
      obsbc(j2-1, nobsbc[j2-1])=i;
      nobsbc[j2-1] += 1;
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = arma::dot(xtf.col(i),betatf.col(k2-1));
      tmp3 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglikn += std::log(tmp3);
      } else {
        tmp3 = 1.0 - tmp3; loglikn += std::log(tmp3);
      }
    }
    loglikn += (maxL1-1.0)*std::log(2.0);
  }
  logpost = loglikn -0.5*arma::dot((betace-betacepm), betaceprec*(betace-betacepm));;
}

// update baseline reg coefficients using adaptive M-H
void update_regcoeff_adaptiveMH(arma::vec& betace, const arma::mat& betatf, const Rcpp::NumericVector& y, const arma::mat& xce, 
        const arma::vec& vn, const arma::mat& xtf, double sigma2, const arma::vec& betacepm, const arma::mat& betaceprec,
        Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, int maxL, double& rej, arma::mat& Snew, arma::vec& betabarnew, 
        int pce, int adpl0, arma::mat adpS0, double adapter, int iscan){
  // start
  const arma::mat Ip = arma::eye(pce, pce);
  arma::vec betaceold = betace;
  // denominator
  double tmpold = 0;
  logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, betacepm, betaceprec, nobsbc, obsbc, tmpold, maxL);
  if(iscan>adpl0){
    betace = mvrnorm(betaceold, Snew);
  }else{
    betace = mvrnorm(betaceold, adpS0);
  }
  // numerator
  double tmpnew = 0;
  logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, betacepm, betaceprec, nobsbc, obsbc, tmpnew, maxL);
  // decision
  double ratio = std::exp(tmpnew-tmpold);
  double uu = unif_rand();
  if(uu>ratio) {betace = betaceold; rej=1;}
  double nn = iscan+1;
  arma::vec betabarold = betabarnew;
  betabarnew = (nn)/(nn+1.0)*betabarold + betace/(nn+1.0);
  arma::mat Sold = Snew;
  Snew = (nn-1.0)/nn*Sold + adapter/nn*(nn*betabarold*betabarold.t() - (nn+1.0)*betabarnew*betabarnew.t() + betace*betace.t() + 0.01*Ip );
  //Rprintf( "tmpnew %f\n", ratio );
  //Rprintf( "tmpold %f\n", tmpold );
  //Rprintf( "S2%f\n", Snew(1,1) );
  //Rprintf( "S3%f\n", Snew(2,2) );
}

// update baseline reg coefficients using IWSL based M-H sampling
void update_regcoeff_iwls(arma::vec& betace, const arma::mat& betatf, const Rcpp::NumericVector& y, const arma::mat& xce, 
        const arma::vec& vn, const arma::mat& xtf, double sigma2, const arma::vec& betacepm, const arma::mat& betaceprec,
        Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, int maxL, double& rej){
  // start
  arma::vec yvec = as<vec>(y);
  arma::mat C = arma::inv_sympd(betaceprec + xce*xce.t()/sigma2);
  arma::vec m = C*(betaceprec*betacepm + xce*(yvec-vn)/sigma2);
  arma::vec betaceold = betace;
  // denominator
  double tmpold = 0;
  logposldtfp(betaceold, betatf, y, xce, vn, xtf, sigma2, betacepm, betaceprec, nobsbc, obsbc, tmpold, maxL);
  betace = mvrnorm(m, C);
  tmpold += -0.5*arma::dot( (betace-m), arma::solve(C, (betace-m)) );
  // numerator
  double tmpnew = 0;
  logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, betacepm, betaceprec, nobsbc, obsbc, tmpnew, maxL);
  tmpnew += -0.5*arma::dot( (betaceold-m), arma::solve(C, (betaceold-m)) );
  // decision
  double ratio = std::exp(tmpnew-tmpold);
  double uu = unif_rand();
  if(uu>ratio) {betace = betaceold; rej=1;}
  Rprintf( "ratio %f\n", ratio );
}

// function of sigma2 
void loglikldtfpsig(const Rcpp::NumericVector& y, const arma::vec& Xbetav, const arma::mat& xbetatf, double sigma2, 
      Rcpp::IntegerVector& nobsbc, Rcpp::IntegerMatrix& obsbc, double& loglik, int maxL, 
      double a0sig, double b0sig){
  
  // temp variables
  int maxL1 = maxL+1;
  Rcpp::IntegerVector K(maxL1);
  int nrec = y.size();
  double sigma = std::sqrt(sigma2);
  double tmp1, tmp2, tmp3;
  double loglikn = 0;
  for(int i=0; i<nobsbc.size(); ++i) nobsbc(i) = 0;
  
  // find where the data point is located at each level
  for(int i=0; i<nrec; ++i){
    //R_CheckUserInterrupt();
    tmp1 = Xbetav[i];
    tmp3 = (y[i]-tmp1)/sigma;
    loglikn += Rf_dnorm4(y[i], tmp1, sigma, true);
    if (tmp3>4.0) {
      tmp2 = 0.999968;
    } else if(tmp3<-4.0) {
      tmp2 = 0.000032;
    } else {
      tmp2 = Rf_pnorm5(y[i], tmp1, sigma, true, false);
    }
    for(int j1=0; j1<maxL1; ++j1){
      K[j1] = (int)(std::pow(2,j1)*tmp2)+1;
    }
  
    // ldtfp part
    int j=0; int m=0; 
    for(int j1=1; j1<maxL1; ++j1){
      int j2 = j + K[j1];
      int k2 = m + K[j1-1];
      j = j + std::pow(2,j1);
      m = m + std::pow(2,j1-1);
      obsbc(j2-1, nobsbc[j2-1])=i;
      nobsbc[j2-1] += 1;
      int ll = (K[j1-1]-1)*2 +1;
      tmp2 = xbetatf(k2-1,i);
      tmp3 = std::exp(tmp2)/(1.0+std::exp(tmp2));
      if(K[j1]==ll) {
        loglikn += std::log(tmp3);
      } else {
        tmp3 = 1.0 - tmp3; 
        loglikn += std::log(tmp3);
      }
    }
    loglikn += (maxL1-1.0)*std::log(2.0);
  }
  if(a0sig>0){
    loglik = loglikn - (a0sig+1.0)*std::log(sigma2) - b0sig/sigma2;
  } else {
    loglik = loglikn -std::log(sigma2);
  }
}

// log likelihood for updating tf logistic regressions
void compullldtfp(int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      const arma::vec& gamma, const arma::mat& xtx, const arma::mat& xtf, double& logp){
  //temp variables
  double loglikn = 0;
  double logpriorn = -0.5*arma::dot(gamma, xtx*gamma);
  if(n1>0){
    for(int i=0; i<n1; ++i){
      int ll = obsbc(ii,i);
      //R_CheckUserInterrupt();
      double eta = arma::dot(xtf.col(ll), gamma);
      loglikn += eta-std::log(1.0+std::exp(eta));
    }
  }
  if(n2>0){
    for(int i=0; i<n2; ++i){
      int ll = obsbc(jj,i);
      //R_CheckUserInterrupt();
      double eta = arma::dot(xtf.col(ll), gamma);
      loglikn += -std::log(1.0+std::exp(eta));
    }
  }
  logp = loglikn + logpriorn;
}

// find the mean vector m and covariance matrix C for the normal proposal used in IWLS
void tfcoeffproposal(int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      const arma::mat& xtf, const arma::mat invc0, const arma::vec& betatfkk, arma::vec& mbetatf, arma::mat& Cbetatf ){
  //temp variables
  int ptf = invc0.n_cols;
  arma::vec xty(ptf); xty.fill(0.0);
  arma::mat xtx = invc0;
  
  // IWLS step
	if(n1>0){
	  for(int i=0; i<n1; ++i){
  		int ll = obsbc(ii,i);
  		//R_CheckUserInterrupt();
  		arma::vec xtfll = xtf.col(ll);
  		double eta = arma::dot(xtfll, betatfkk);
  		double mu = std::exp(eta)/(1.0+std::exp(eta));
  		double tmp1 = std::exp(eta)/std::pow((1.0+std::exp(eta)),2);
  		double ytilde = eta+(1.0-mu)/tmp1;
  		xtx += xtfll*(arma::trans(xtfll))*tmp1;
  		xty += xtfll*ytilde*tmp1;
	  }
	}
	if(n2>0){
	  for(int i=0; i<n2; ++i){
		int ll = obsbc(jj,i);
		//R_CheckUserInterrupt();
		arma::vec xtfll = xtf.col(ll);
		double eta = arma::dot(xtfll, betatfkk);
		double mu = std::exp(eta)/(1.0+std::exp(eta));
		double tmp1 = std::exp(eta)/std::pow((1.0+std::exp(eta)),2);
		double ytilde = eta+(0.0-mu)/tmp1;
		xtx += xtfll*(arma::trans(xtfll))*tmp1;
		xty += xtfll*ytilde*tmp1;
	  }
	}
	Cbetatf = arma::inv_sympd(xtx);
	mbetatf = Cbetatf*xty;
}

// update tf logistic regressions using IWLS based M-H
void update_tfcoeff_iwls(int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      const arma::mat& xtf, const arma::mat invc0, arma::vec& betatfkk, double& rej){
  // start
  int ptf = invc0.n_cols;
  arma::vec m(ptf); m.fill(0.0);
  arma::mat C = invc0;
  arma::vec betatfkkold = betatfkk;
  // denominator
  tfcoeffproposal(ii, jj, n1, n2, obsbc, xtf, invc0, betatfkkold, m, C);
  double tmpold = 0;
  compullldtfp(ii, jj, n1, n2, obsbc, betatfkkold, invc0, xtf, tmpold);
  betatfkk = mvrnorm(m, C);
  tmpold += mvdnorm(betatfkk, m, C, true);
  // numerator
  tfcoeffproposal(ii, jj, n1, n2, obsbc, xtf, invc0, betatfkk, m, C);
  double tmpnew = 0;
  compullldtfp(ii, jj, n1, n2, obsbc, betatfkk, invc0, xtf, tmpnew);
  tmpnew += mvdnorm(betatfkkold, m, C, true);
  // decision
  double ratio = std::exp(tmpnew-tmpold);
  double uu = unif_rand();
  if(uu>ratio) {betatfkk = betatfkkold; rej=1;}
}

// update tf logistic regressions using slice sampling
void updatelrcoefldtfpss(int kk, int ii, int jj, int n1, int n2, const Rcpp::IntegerMatrix& obsbc, 
      arma::mat& betatf, const arma::mat& xtf, const arma::mat c0){
  //temp variables
  double llim, rlim, gllim, grlim;
  double xx0, xx1, gxx0, gxx1;
  double logy, uwork;
  int JJ, KK, mm;
  double win=0.4; mm=10;
  int ptf = c0.n_cols;
  arma::mat xtx = arma::inv_sympd(c0);
  arma::vec beta = betatf.col(kk);
  
  for(int i=0; i<ptf; ++i){
    //R_CheckUserInterrupt();
    int evali=1;
    xx0 = beta(i);
    compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, gxx0);
    logy = gxx0-exp_rand();
    uwork=unif_rand();
    llim=xx0-win*uwork;
    rlim=llim+win;
    uwork=unif_rand();
    JJ = (int)(mm*uwork);
    KK = (mm-1)-JJ;
    // repeat 
    ++evali;
    beta(i)=llim;
    compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, gllim);
    ++evali;
    beta(i)=rlim;
    compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, grlim);
    while((JJ>0)&(gllim>logy)){
      llim -= win;
      JJ -= 1;
      ++evali;
      beta(i)=llim;
      compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, gllim);
    }
    while((KK>0)&(grlim>logy)){
      rlim += win;
      KK -= 1;
      ++evali;
      beta(i)=rlim;
      compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, grlim);
    }
    xx1 = llim +(rlim-llim)*unif_rand();
    ++evali;
    beta(i)=xx1;
    compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, gxx1);
    while(gxx1<logy){
      if(xx1>xx0) rlim=xx1;
      if(xx1<xx0) llim=xx1;
      xx1 = llim +(rlim-llim)*unif_rand();
      ++evali;
      beta(i)=xx1;
      compullldtfp(ii, jj, n1, n2, obsbc, beta, xtx, xtf, gxx1);
    }
    betatf(i,kk)=xx1;
  }
}

// Calculate CPO for spatial LDTFP
arma::vec spldtfp_Linv(const Rcpp::NumericMatrix& tobs, const Rcpp::IntegerVector& type, const arma::vec& xbetav, 
      const arma::mat& xbetatf, double sigma2, int maxL){
  int nrec = type.size();
  arma::vec res(nrec);
  double tmp1, tmp2;
  for(int i=0; i<nrec; ++i){
    double y1 = std::log(tobs(i,0));
    double y2 = std::log(tobs(i,1));
    if(type[i]==2){
      cdfldtfp(y2, xbetav[i], xbetatf.col(i), sigma2, tmp1, maxL);
      res(i) = 1.0/tmp1;
    } else if (type[i]==3){
      cdfldtfp(y1, xbetav[i], xbetatf.col(i), sigma2, tmp1, maxL);
      cdfldtfp(y2, xbetav[i], xbetatf.col(i), sigma2, tmp2, maxL);
      res(i) = 1.0/(tmp2-tmp1);
    } else if (type[i]==0){
      cdfldtfp(y1, xbetav[i], xbetatf.col(i), sigma2, tmp1, maxL);
      res(i) = 1.0/(1.0-tmp1);
    } else {
      ldensldtfp(y1, xbetav[i], xbetatf.col(i), sigma2, tmp1, maxL);
      res(i) = std::exp(y1-tmp1);
    }
  }
  return(res);
}

// Get density or survival Plots for frailty LDTFP AFT
RcppExport SEXP frailtyGAFTplots(SEXP tgrid_, SEXP xcepred_, SEXP xtfpred_, SEXP betace_, 
                                 SEXP betatf_, SEXP v_, SEXP sigma2_, SEXP maxL_, SEXP CI_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  arma::vec tgrid = as<vec>(tgrid_);
  arma::mat xcepred = as<mat>(xcepred_); // npred by pce;
  arma::mat xtfpred = as<mat>(xtfpred_); // npred by ptf;
  arma::mat betace = as<mat>(betace_); // pce by nsave;
  arma::mat v = as<arma::mat>(v_); // npred by nsave;
  Rcpp::NumericVector vecArray(betatf_);
  Rcpp::IntegerVector arrayDims = vecArray.attr("dim");
  arma::cube betatf(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  //arma::mat v = as<mat>(v_); // npred by nsave;
  Rcpp::NumericVector sigma2(sigma2_); // nsave by 1;
  
  double CI = as<double>(CI_);
  int maxL = as<int>(maxL_);
  int nsave = sigma2.size();
  int ngrid = tgrid.size();
  int npred = xcepred.n_rows;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  
  // Temp variables
  arma::vec ygrid = arma::log(tgrid);
  Rcpp::NumericVector estfArray(nsave*ngrid*npred);
  arma::cube estf(estfArray.begin(), ngrid, nsave, npred, false);
  Rcpp::NumericVector estSArray(nsave*ngrid*npred);
  arma::cube estS(estSArray.begin(), ngrid, nsave, npred, false);
  Rcpp::NumericVector esthArray(nsave*ngrid*npred);
  arma::cube esth(esthArray.begin(), ngrid, nsave, npred, false);
  
  // things to save;
  arma::mat fhat(ngrid, npred);
  arma::mat fhatup(ngrid, npred);
  arma::mat fhatlow(ngrid, npred);
  arma::mat Shat(ngrid, npred);
  arma::mat Shatup(ngrid, npred);
  arma::mat Shatlow(ngrid, npred);
  arma::mat hhat(ngrid, npred);
  arma::mat hhatup(ngrid, npred);
  arma::mat hhatlow(ngrid, npred);
  
  // temp variables;
  double tmp1=0; 
  double tmp2=0;
  
  for(int i=0; i<nsave; ++i){
    arma::vec xibetav = xcepred*betace.col(i)+v.col(i);
    arma::mat xibetatf = arma::trans( xtfpred*betatf.slice(i) );
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        ldensldtfp(ygrid[k], xibetav(j), xibetatf.col(j), sigma2[i], tmp1, maxL);
        estf(k, i, j) = std::exp(tmp1)/tgrid[k];
        cdfldtfp(ygrid[k], xibetav(j), xibetatf.col(j), sigma2[i], tmp2, maxL);
        estS(k, i, j) = 1.0 - tmp2;
        esth(k, i, j) = std::exp(tmp1)/tgrid[k]/(1.0 - tmp2);
      }
    }
  }
  for(int j=0; j<npred; ++j){
    fhat.col(j) = arma::mean(estf.slice(j), 1);
    Shat.col(j) = arma::mean(estS.slice(j), 1);
    hhat.col(j) = arma::mean(esth.slice(j), 1);
    arma::mat temp = arma::sort(estf.slice(j),"ascend", 1);
    fhatlow.col(j) = temp.col(low);
    fhatup.col(j) = temp.col(up);
    temp = arma::sort(estS.slice(j),"ascend", 1);
    Shatlow.col(j) = temp.col(low);
    Shatup.col(j) = temp.col(up);
    temp = arma::sort(esth.slice(j),"ascend", 1);
    hhatlow.col(j) = temp.col(low);
    hhatup.col(j) = temp.col(up);
  }
  return List::create(Named("fhat")=fhat,
                      Named("Shat")=Shat,
                      Named("hhat")=hhat,
                      Named("fhatlow")=fhatlow,
                      Named("fhatup")=fhatup,
                      Named("Shatlow")=Shatlow,
                      Named("Shatup")=Shatup,
                      Named("hhatlow")=hhatlow,
                      Named("hhatup")=hhatup);
  END_RCPP
}

// BayesFactors 
RcppExport SEXP BayesFactor(SEXP betatf_, SEXP maxL_, SEXP gprior_, SEXP alpha_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  const arma::mat betatf = as<mat>(betatf_); // (2**maxL-2)*ptf by nsave;
  const int maxL = as<int>(maxL_);
  const arma::mat gprior = as<mat>(gprior_);
  const double alpha = as<double>(alpha_);
  int ntf = std::pow(2, maxL)-2;
  int nbetatf = betatf.n_rows;
  int ptf = nbetatf/ntf;
  int nsave = betatf.n_cols;
  
  
  // Temp variables
  double tmp2=0;
  Rcpp::NumericVector BF1(ptf, 0.0);
  Rcpp::NumericVector BF0(ptf, 0.0);
  
  // Test for individual covariate effect
  for(int i=0; i<ptf; ++i){
    arma::mat gammai(ntf, nsave);
    arma::vec veczero(ntf); veczero.fill(0.0);
    for(int j=0; j<ntf; ++j){
      gammai.row(j) = betatf.row(i+j*ptf);
    }
    arma::vec mean1 = arma::mean(gammai, 1);
    arma::mat cov1 = arma::cov(gammai.t());
    BF1[i] = mvdnorm(veczero, mean1, cov1, false);
  }
  for(int i=0; i<ptf; ++i){
    tmp2=0;
    for(int j=0; j<ntf; ++j){
      tmp2 += Rf_dnorm4(0, 0, std::sqrt( gprior(i,i)/alpha )/((int)(log(j+2)/log(2))+1.0), true);
    }
    BF0[i] = std::exp(tmp2);
  }
  NumericVector BFindividual = BF0/BF1;
  
  // overall test for covariates dependency of LDTFP (i.e. all coeffs=0 execpt intercept)
  double BF1overall, BF0overall;
  if(ptf>1){
    arma::mat gammai(ntf*(ptf-1), nsave);
    arma::vec zero1(ntf*(ptf-1)); zero1.fill(0.0);
    for(int j=0; j<ntf; ++j){
      gammai.rows(j*(ptf-1), (j+1)*(ptf-1)-1) = betatf.rows(j*ptf+1, (j+1)*ptf-1);
    }
    arma::vec mean1 = arma::mean(gammai, 1);
    arma::mat cov1 = arma::cov(gammai.t());
    BF1overall = mvdnorm(zero1, mean1, cov1, false);
    arma::vec zero0(ptf-1); zero0.fill(0.0);
    tmp2=0;
    for(int j=0; j<ntf; ++j){
      tmp2 += mvdnorm(zero0, zero0, gprior.submat(1,1,ptf-1,ptf-1)/(alpha*std::pow(((int)(log(j+2)/log(2))+1.0),2)), true);
    }
    BF0overall = std::exp(tmp2);
  }else{
    BF1overall = -1;
    BF0overall = 100;
  }
  double BFoverallLDTFP = BF0overall/BF1overall;
  
  // overall test for Normality of LDTFP (i.e. all coeffs=0)
  if(ptf>0){
    arma::vec zero1(nbetatf); zero1.fill(0.0);
    arma::vec mean1 = arma::mean(betatf, 1);
    arma::mat cov1 = arma::cov(betatf.t());
    BF1overall = mvdnorm(zero1, mean1, cov1, false);
    arma::vec zero0(ptf); zero0.fill(0.0);
    tmp2=0;
    for(int j=0; j<ntf; ++j){
      tmp2 += mvdnorm(zero0, zero0, gprior/(alpha*std::pow(((int)(log(j+2)/log(2))+1.0),2)), true);
    }
    BF0overall = std::exp(tmp2);
  }else{
    BF1overall = -1;
    BF0overall = 100;
  }
  double BFoverallParam = BF0overall/BF1overall;
  
  return List::create(Named("BFindividual")=BFindividual, 
                      Named("BFoverallLDTFP")=BFoverallLDTFP,
                      Named("BFoverallParam")=BFoverallParam);
  END_RCPP
}
