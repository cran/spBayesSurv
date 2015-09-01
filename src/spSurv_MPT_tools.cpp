#include "spSurv_MPT_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// from conditional Ys to cumulative probs
void Ys_to_probs(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& probs, int maxL){
  int nprob = probs.size();
  Rcpp::NumericVector temp(nprob,1.0);
  for(int j=0; j<maxL; ++j){
    int jk = std::pow(2, j);
    int j0 = jk-1;
    for(int k=0; k<jk; ++k){
      probs[2*k] = temp[k]*Ys[j0+k];
      probs[2*k+1] = temp[k]*(1.0-Ys[j0+k]);
    }
    for( int i=0; i<(2*jk); ++i) temp[i] = probs[i];
  }
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AFTlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(R_NegInf);
  if(t>SYSMAX) return(R_NegInf);
  
  // log density of log-logistic(exp(x'beta)*t);
  double tmp1 = 1.0+std::pow(exp(th1+xbetai)*t,exp(th2));
  double ll = th1+th2+(exp(th2)-1.0)*(th1+xbetai+log(t)) - 2.0*log(tmp1);
  
  // find where t is located at level maxL
  double St = 1.0/std::min(tmp1,SYSMAX);
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  
  ll += xbetai + maxL*log(2.0) + log(probs[kt]);
  return(ll);
}

// log Survival function of t given xi
double AFTSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(1.0);
  if(t>SYSMAX) return(0.0);
  
  // find where t is located at level maxL
  double tmp1 = 1.0+std::pow(exp(th1+xbetai)*t,exp(th2));
  double St = 1.0/std::min(tmp1,SYSMAX);
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }
  return(surv);
}

// log likelihood given data, frailties and parameters 
void AFTloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( AFTSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += AFTlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - AFTSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      ll += log( AFTSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - AFTSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  }
}

// log likelihood given frailties, parameters and data of block i
void AFTloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( AFTSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += AFTlogpdf(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - AFTSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else{
      ll += log( AFTSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) - AFTSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }
  }
}

// Calculate 1.0/likelihood for CPO
arma::vec AFTinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      res[i] = 1.0/AFTSurv(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==1){
      res[i] = 1.0/exp(AFTlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2));
    }else if(type[i]==2){
      res[i] = 1.0/( 1.0 - AFTSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      res[i] = 1.0/( AFTSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - AFTSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  } 
  return(res);
}


/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double POlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(R_NegInf);
  if(t>SYSMAX) return(R_NegInf);
  
  // log density of log-logistic(t);
  double tmp1 = 1.0+std::pow(exp(th1)*t,exp(th2));
  double ll = th1+th2+(exp(th2)-1.0)*(th1+log(t)) - 2.0*log(tmp1);
  
  // find where t is located at level maxL
  double St = 1.0/std::min(tmp1,SYSMAX);
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  
  // log f0t given Ys, theta
  ll += maxL*log(2.0) + log(probs[kt]);
  // surv given Ys, theta
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }

  ll += -xbetai - 2.0*log(1.0+(exp(-xbetai)-1.0)*surv);
  return(ll);
}

// log Survival function of t given xi
double POSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(1.0);
  if(t>SYSMAX) return(0.0);
  
  // find where t is located at level maxL
  double tmp1 = 1.0+std::pow(exp(th1)*t,exp(th2));
  double St = 1.0/std::min(tmp1,SYSMAX);
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }
  
  surv = exp( -xbetai+log(surv)-log(1.0+(exp(-xbetai)-1.0)*surv) );
  return(surv);
}

// log likelihood given data, frailties and parameters 
void POloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( POSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += POlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - POSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      ll += log( POSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - POSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  }
}

// log likelihood given frailties, parameters and data of block i
void POloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( POSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += POlogpdf(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - POSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else{
      ll += log( POSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) - POSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }
  }
}

// Calculate 1.0/likelihood for CPO
arma::vec POinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      res[i] = 1.0/POSurv(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==1){
      res[i] = 1.0/exp(POlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2));
    }else if(type[i]==2){
      res[i] = 1.0/( 1.0 - POSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      res[i] = 1.0/( POSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - POSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  } 
  return(res);
}
      
/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PHlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(R_NegInf);
  if(t>SYSMAX) return(R_NegInf);
  
  // log density of Weibull(t);
  double tmp1 = std::pow(exp(th1)*t,exp(th2));
  double ll = th1+th2+(exp(th2)-1.0)*(th1+log(t)) - tmp1;
  
  // find where t is located at level maxL
  double St = exp(-std::min(tmp1,LOGSYSMAX));
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  
  // log f0t given Ys, theta
  ll += maxL*log(2.0) + log(probs[kt]);
  // surv given Ys, theta
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }

  ll += xbetai + (exp(xbetai)-1.0)*log(surv);
  return(ll);
}

// log Survival function of t given xi
double PHSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(1.0);
  if(t>SYSMAX) return(0.0);
  
  // find where t is located at level maxL
  double tmp1 = std::pow(exp(th1)*t,exp(th2));
  double St = exp(-std::min(tmp1,LOGSYSMAX));
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }
  
  surv = exp( exp(xbetai)*log(surv) );
  return(surv);
}

// log likelihood given data, frailties and parameters 
void PHloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( PHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += PHlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - PHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      ll += log( PHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - PHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  }
}

// log likelihood given frailties, parameters and data of block i
void PHloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( PHSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += PHlogpdf(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - PHSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else{
      ll += log( PHSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) - PHSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }
  }
}

// Calculate 1.0/likelihood for CPO
arma::vec PHinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      res[i] = 1.0/PHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==1){
      res[i] = 1.0/exp(PHlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2));
    }else if(type[i]==2){
      res[i] = 1.0/( 1.0 - PHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      res[i] = 1.0/( PHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - PHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  } 
  return(res);
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AHlogpdf(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(R_NegInf);
  if(t>SYSMAX) return(R_NegInf);
  
  // log density of Weibull(t);
  double tmp1 = std::pow(exp(th1+xbetai)*t,exp(th2));
  double ll = th1+th2+(exp(th2)-1.0)*(th1+xbetai+log(t)) - tmp1;
  
  // find where t is located at level maxL
  double St = exp(-std::min(tmp1,LOGSYSMAX));
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  
  // log f0t given Ys, theta
  ll += maxL*log(2.0) + log(probs[kt]);
  // surv given Ys, theta
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }

  ll += (exp(-xbetai)-1.0)*log(surv);
  return(ll);
}

// log Survival function of t given xi
double AHSurv(double t, double xbetai, Rcpp::NumericVector probs, int maxL, double th1, double th2){
  //temp variables
  int nprob = probs.size();
  R_CheckUserInterrupt();
  if(t<SYSMIN) return(1.0);
  if(t>SYSMAX) return(0.0);
  
  // find where t is located at level maxL
  double tmp1 = std::pow(exp(th1+xbetai)*t,exp(th2));
  double St = exp(-std::min(tmp1,LOGSYSMAX));
  int kt = (int)(nprob-(double)(nprob)*St);
  if(kt==nprob) --kt;
  double surv = probs[kt]*( (double)(nprob)*St -(nprob-kt-1) );
  if(kt<nprob-1){
    for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
  }
  
  surv = exp( exp(-xbetai)*log(surv) );
  return(surv);
}

// log likelihood given data, frailties and parameters 
void AHloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( AHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += AHlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - AHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      ll += log( AHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - AHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  }
}

// log likelihood given frailties, parameters and data of block i
void AHloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta, 
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2, double& ll, int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      ll += log( AHSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else if(type[i]==1){
      ll += AHlogpdf(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2);
    }else if(type[i]==2){
      ll += log( 1.0 - AHSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }else{
      ll += log( AHSurv(t1[i], Xbeta[i]+vi, probs, maxL, th1, th2) - AHSurv(t2[i], Xbeta[i]+vi, probs, maxL, th1, th2) );
    }
  }
}

// Calculate 1.0/likelihood for CPO
arma::vec AHinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::IntegerVector& type, const Rcpp::NumericVector& Xbeta,
      const Rcpp::NumericVector& probs, int maxL, double th1, double th2){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    R_CheckUserInterrupt();
    if(type[i]==0){
      res[i] = 1.0/AHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2);
    }else if(type[i]==1){
      res[i] = 1.0/exp(AHlogpdf(t1[i], Xbeta[i], probs, maxL, th1, th2));
    }else if(type[i]==2){
      res[i] = 1.0/( 1.0 - AHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }else{
      res[i] = 1.0/( AHSurv(t1[i], Xbeta[i], probs, maxL, th1, th2) - AHSurv(t2[i], Xbeta[i], probs, maxL, th1, th2) );
    }
  } 
  return(res);
}
