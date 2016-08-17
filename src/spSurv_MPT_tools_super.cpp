#include "spSurv_MPT_tools_super.h"

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
// from conditional Ys to cumulative probs, when Ys[0] is fixed at 0.5
void Ys_to_probs2(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& probs, int maxL){
  int nprob = probs.size();
  Rcpp::NumericVector temp(nprob,0.5);
  for(int j=1; j<maxL; ++j){
    int jk = std::pow(2, j);
    int j0 = jk-2;
    for(int k=0; k<jk; ++k){
      probs[2*k] = temp[k]*Ys[j0+k];
      probs[2*k+1] = temp[k]*(1.0-Ys[j0+k]);
    }
    for( int i=0; i<(2*jk); ++i) temp[i] = probs[i];
  }
}

/////////////////////////////////////////////////////////////////////////
/////////////////// baseline suvival functions///////////////////////////
/////////////////////////////////////////////////////////////////////////
double S0MPT(double t, double th1, double th2, Rcpp::NumericVector probs, int maxL, bool MPT, int dist){
  int nprob = probs.size();
  if(t<0) t=0;
  double tmp1 = (log(t)+th1)*exp(th2);
  if(tmp1<LOGSYSMIN) tmp1=LOGSYSMIN;
  if(tmp1>LOGSYSMAX) tmp1=LOGSYSMAX;
  double surv;
  if(MPT){
    if(dist==1){
      surv = 1.0/(1.0+exp(tmp1));
    }else if (dist==2){
      surv = Rf_plnorm(t, -th1, exp(-th2), false, false);
    }else{
      surv = exp(-exp(tmp1));
    }
    int kt = (int)(nprob-(double)(nprob)*surv);
    if(kt==nprob) --kt;
    surv = probs[kt]*( (double)(nprob)*surv -(nprob-kt-1) );
    if(kt<nprob-1){
      for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
    }
  }else{
    if(dist==1){
      surv = 1.0/(1.0+exp(tmp1));
    }else if (dist==2){
      surv = Rf_plnorm(t, -th1, exp(-th2), false, false);
    }else{
      surv = exp(-exp(tmp1));
    }
  }
  if(surv<SYSMIN) surv=SYSMIN;
  return(surv);
}
double logf0MPT(double t, double th1, double th2, Rcpp::NumericVector probs, int maxL, bool MPT, int dist){
  int nprob = probs.size();
  if(t<0) t=0;
  double tmp1 = (log(t)+th1)*exp(th2);
  if(tmp1<LOGSYSMIN) tmp1=LOGSYSMIN;
  if(tmp1>LOGSYSMAX) tmp1=LOGSYSMAX;
  double ll, surv;
  if(MPT){
    if(dist==1){
      surv = 1.0/(1.0+exp(tmp1));
      int kt = (int)(nprob-(double)(nprob)*surv);
      if(kt==nprob) --kt;
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - 2.0*log(1.0+exp(tmp1)) + maxL*log(2.0) + log(probs[kt]);
    }else if (dist==2){
      surv = Rf_plnorm(t, -th1, exp(-th2), false, false);
      int kt = (int)(nprob-(double)(nprob)*surv);
      if(kt==nprob) --kt;
      ll = Rf_dlnorm(t, -th1, exp(-th2), true) + maxL*log(2.0) + log(probs[kt]);
    }else{
      surv = exp(-exp(tmp1));
      int kt = (int)(nprob-(double)(nprob)*surv);
      if(kt==nprob) --kt;
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - exp(tmp1) + maxL*log(2.0) + log(probs[kt]);
    }
  }else{
    if(dist==1){
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - 2.0*log(1.0+exp(tmp1));
    }else if (dist==2){
      ll = Rf_dlnorm(t, -th1, exp(-th2), true);
    }else{
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - exp(tmp1);
    }
  }
  return(ll);
}
double logh0MPT(double t, double th1, double th2, Rcpp::NumericVector probs, int maxL, bool MPT, int dist){
  int nprob = probs.size();
  if(t<0) t=0;
  double tmp1 = (log(t)+th1)*exp(th2);
  if(tmp1<LOGSYSMIN) tmp1=LOGSYSMIN;
  if(tmp1>LOGSYSMAX) tmp1=LOGSYSMAX;
  double ll, surv;
  if(MPT){
    int kt=0;
    if(dist==1){
      surv = 1.0/(1.0+exp(tmp1));
      kt = (int)(nprob-(double)(nprob)*surv);
      if(kt==nprob) --kt;
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - 2.0*log(1.0+exp(tmp1)) + maxL*log(2.0) + log(probs[kt]);
    }else if (dist==2){
      surv = Rf_plnorm(t, -th1, exp(-th2), false, false);
      kt = (int)(nprob-(double)(nprob)*surv);
      if(kt==nprob) --kt;
      ll = Rf_dlnorm(t, -th1, exp(-th2), true) + maxL*log(2.0) + log(probs[kt]);
    }else{
      surv = exp(-exp(tmp1));
      kt = (int)(nprob-(double)(nprob)*surv);
      if(kt==nprob) --kt;
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - exp(tmp1) + maxL*log(2.0) + log(probs[kt]);
    }
    surv = probs[kt]*( (double)(nprob)*surv -(nprob-kt-1) );
    if(kt<nprob-1){
      for(int j=kt+1; j<nprob; ++j ) surv += probs[j];
    }
    if(surv<SYSMIN) surv=SYSMIN;
    return(ll-log(surv));
  }else{
    if(dist==1){
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1 - log(1.0+exp(tmp1));
    }else if (dist==2){
      ll = Rf_dlnorm(t, -th1, exp(-th2), true)-Rf_plnorm(t, -th1, exp(-th2), false, true);;
    }else{
      ll = th2 + tmp1*(1.0-exp(-th2)) + th1;
    }
    return(ll);
  }
}
