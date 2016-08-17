#include "spSurv_BP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// from conditional Ys to cumulative probs
void Ys_to_weight(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& weight){
  int nYs = Ys.size();
  Rcpp::NumericVector newYs(nYs+1, 1.0);
  for(int j=0; j<nYs; ++j) newYs[j] = exp(Ys[j]);
  double den = Rcpp::sum(newYs);
  for(int j=0; j<=nYs; ++j) weight[j] = newYs[j]/den; 
}

/////////////////////////////////////////////////////////////////////////
/////////////////// baseline suvival functions///////////////////////////
/////////////////////////////////////////////////////////////////////////
double S0BP(double t, double th1, double th2, Rcpp::NumericVector w, bool BP, int dist){
  if(t<SYSMIN) return(1.0);
  double tmp1 = (std::log(t)+th1)*std::exp(th2);
  int J = w.size();
  if(J==1) BP=false;
  double surv=0, logitSt=0, logtemp=0, Ixprev=0, Ft=0;
  if(BP){
    if(dist==1){
      Ft = std::exp(tmp1)/(1.0+std::exp(tmp1));
    }else if (dist==2){
      Ft = Rf_pnorm5(tmp1, 0, 1, true, false);
    }else{
      Ft = 1.0-std::exp(-std::exp(tmp1));
    }
    if(Ft<SYSMIN) Ft=SYSMIN;
    logitSt = std::log(1.0-Ft)-std::log(Ft); //log(St/(1-St));
    if(logitSt<LOGSYSMIN) return(SYSMIN);
    logtemp = J*std::log(Ft); // J*log(1-St);
    Ixprev = 1.0 - std::exp(logtemp); 
    surv = w[0]*Ixprev;
    for(int j=1; j<J; ++j){
      logtemp += logitSt + std::log((J-j+1.0)/(j+0.0));
      Ixprev -= std::exp(logtemp);
      surv += w[j]*Ixprev;
    }
  }else{
    if(dist==1){
      surv = 1.0/(1.0+std::exp(tmp1));
    }else if (dist==2){
      surv = Rf_pnorm5(tmp1, 0, 1, false, false);
    }else{
      surv = std::exp(-std::exp(tmp1));
    }
  }
  if(surv>SYSMIN){
    return(surv);
  }else{
    return(SYSMIN);
  }
}
double F0BP(double t, double th1, double th2, Rcpp::NumericVector w, bool BP, int dist){
  if(t<SYSMIN) return(SYSMIN);
  double tmp1 = (std::log(t)+th1)*std::exp(th2);
  int J = w.size();
  if(J==1) BP=false;
  double surv=0, logitSt=0, logtemp=0, Ixprev=0, Ft=0;
  if(BP){
    if(dist==1){
      Ft = std::exp(tmp1)/(1.0+std::exp(tmp1));
    }else if (dist==2){
      Ft = Rf_pnorm5(tmp1, 0, 1, true, false);
    }else{
      Ft = 1.0-std::exp(-std::exp(tmp1));
    }
    if(Ft<SYSMIN) Ft=SYSMIN;
    logitSt = std::log(1.0-Ft)-std::log(Ft); //log(St/(1-St));
    if(logitSt<LOGSYSMIN) return(SYSMIN);
    logtemp = J*std::log(Ft); // J*log(1-St);
    Ixprev = 1.0 - std::exp(logtemp); 
    surv = w[0]*Ixprev;
    for(int j=1; j<J; ++j){
      logtemp += logitSt + std::log((J-j+1.0)/(j+0.0));
      Ixprev -= std::exp(logtemp);
      surv += w[j]*Ixprev;
    }
  }else{
    if(dist==1){
      surv = std::exp(tmp1)/(1.0+std::exp(tmp1));
    }else if (dist==2){
      surv = Rf_pnorm5(tmp1, 0, 1, false, false);
    }else{
      surv = std::exp(-std::exp(tmp1));
    }
  }
  return(surv);
}
double logf0BP(double t, double th1, double th2, Rcpp::NumericVector w, bool BP, int dist){
  if(t<SYSMIN) return(LOGSYSMIN);
  double tmp1 = (std::log(t)+th1)*std::exp(th2);
  if(tmp1>LOGSYSMAX) return(LOGSYSMIN);
  int J = w.size();
  if(J==1) BP=false;
  double ll=0, logft=0, logitSt=0, logtemp=0, Ft=0;
  if(BP){
    if(dist==1){
      Ft = std::exp(tmp1)/(1.0+std::exp(tmp1));
      logft = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - 2.0*std::log(1.0+std::exp(tmp1));
    }else if (dist==2){
      Ft = Rf_pnorm5(tmp1, 0, 1, true, false);
      logft = Rf_dlnorm(t, -th1, std::exp(-th2), true);
    }else{
      Ft = 1.0-std::exp(-std::exp(tmp1));
      logft = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - std::exp(tmp1);
    }
    if(Ft<SYSMIN) return(LOGSYSMIN);
    logitSt = std::log(1.0-Ft)-std::log(Ft); //log(St/(1-St));
    if(logitSt<LOGSYSMIN) return(LOGSYSMIN);
    logtemp = std::log(J) + (J-1.0)*std::log(Ft); // log(J)+(J-1)*log(1-St);
    ll = w[0]*std::exp(logtemp+logft);
    for(int j=1; j<J; ++j){
      logtemp += logitSt + std::log((J-j+0.0)/(j+0.0));
      ll += w[j]*std::exp(logtemp+logft);
    }
    ll = std::log(ll);
  }else{
    if(dist==1){
      ll = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - 2.0*std::log(1.0+std::exp(tmp1));
    }else if (dist==2){
      ll = Rf_dlnorm(t, -th1, std::exp(-th2), true);
    }else{
      ll = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - std::exp(tmp1);
    }
  }
  if(ll>LOGSYSMIN) {
    return(ll);
  } else {
    return(LOGSYSMIN);
  }
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AFT_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                     bool BP, int dist, double xibeta){
  double ll = xibeta + logf0BP(exp(xibeta)*t, th1, th2, w, BP, dist);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double AFT_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                      bool BP, int dist, double xibeta){
  double ll = std::log(S0BP(exp(xibeta)*t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double AFT_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta){
  double ll = std::log(1.0-S0BP(exp(xibeta)*t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double AFT_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                          bool BP, int dist, double xibeta){
  double St1 = S0BP(exp(xibeta)*t1, th1, th2, w, BP, dist);
  double St2 = S0BP(exp(xibeta)*t2, th1, th2, w, BP, dist);
  double ll = std::log(std::abs(St1 - St2));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void AFT_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                   bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += AFT_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += AFT_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += AFT_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else{
      ll += AFT_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -AFT_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void AFT_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                         const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                         bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                         int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += AFT_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += AFT_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += AFT_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else{
      ll += AFT_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -AFT_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec AFT_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-AFT_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-AFT_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-AFT_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-AFT_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(AFT_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]));
  } 
  return(res);
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PH_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                    bool BP, int dist, double xibeta){
  double ll = xibeta + logf0BP(t, th1, th2, w, BP, dist);
  ll += (exp(xibeta)-1.0)*log(S0BP(t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double PH_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta){
  double ll = exp(xibeta)*log(S0BP(t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double PH_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                    bool BP, int dist, double xibeta){
  double ll = log( 1.0-exp( exp(xibeta)*log(S0BP(t, th1, th2, w, BP, dist)) ) );
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double PH_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta){
  double St1 = exp( exp(xibeta)*log(S0BP(t1, th1, th2, w, BP, dist)) );
  double St2 = exp( exp(xibeta)*log(S0BP(t2, th1, th2, w, BP, dist)) );
  double ll = std::log(std::abs(St1 - St2));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void PH_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                  const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                  bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += PH_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += PH_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += PH_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else{
      ll += PH_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -PH_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void PH_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                        int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += PH_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += PH_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += PH_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else{
      ll += PH_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -PH_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec PH_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-PH_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-PH_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-PH_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-PH_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(PH_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]));
  } 
  return(res);
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PO_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                    bool BP, int dist, double xibeta){
  double ll = logf0BP(t, th1, th2, w, BP, dist)-xibeta;
  ll += -2.0*log(1.0+(exp(-xibeta)-1.0)*S0BP(t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double PO_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta){
  double S0t = S0BP(t, th1, th2, w, BP, dist);
  double ll = log(S0t)-xibeta-log(1.0+(exp(-xibeta)-1.0)*S0t);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double PO_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                    bool BP, int dist, double xibeta){
  double S0t = S0BP(t, th1, th2, w, BP, dist);
  double ll = log(1.0-S0t)-log(1.0+(exp(-xibeta)-1.0)*S0t);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double PO_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta){
  double S0t1 = S0BP(t1, th1, th2, w, BP, dist);
  double S0t2 = S0BP(t2, th1, th2, w, BP, dist);
  double exibeta = exp(-xibeta);
  double ll = log(std::abs(exibeta*S0t1/(1.0+(exibeta-1.0)*S0t1) - exibeta*S0t2/(1.0+(exibeta-1.0)*S0t2)));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void PO_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                  const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                  bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += PO_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += PO_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += PO_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else{
      ll += PO_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -PO_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void PO_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                        int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += PO_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += PO_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += PO_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else{
      ll += PO_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -PO_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec PO_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-PO_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-PO_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-PO_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-PO_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(PO_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]));
  } 
  return(res);
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AH_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                    bool BP, int dist, double xibeta){
  double ll = logf0BP(exp(xibeta)*t, th1, th2, w, BP, dist);
  ll += (exp(-xibeta)-1.0)*log(S0BP(exp(xibeta)*t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double AH_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta){
  double ll = exp(-xibeta)*log(S0BP(exp(xibeta)*t, th1, th2, w, BP, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double AH_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                    bool BP, int dist, double xibeta){
  double ll = std::log( 1.0-exp(exp(-xibeta)*log(S0BP(exp(xibeta)*t, th1, th2, w, BP, dist))) );
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double AH_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta){
  double St1 = exp(exp(-xibeta)*log(S0BP(exp(xibeta)*t1, th1, th2, w, BP, dist)));
  double St2 = exp(exp(-xibeta)*log(S0BP(exp(xibeta)*t2, th1, th2, w, BP, dist)));
  double ll = std::log(std::abs(St1 - St2));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void AH_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                  const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                  bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += AH_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += AH_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += AH_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }else{
      ll += AH_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -AH_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void AH_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                        int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += AH_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += AH_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += AH_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }else{
      ll += AH_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -AH_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec AH_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-AH_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-AH_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-AH_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-AH_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(AH_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta[i]));
  } 
  return(res);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////// Super model: PH_PO_AFT //////////////////////////
/////////////////////////////////////////////////////////////////////////
// log density of t given xi
double PHPOAFT_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                         bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q){
  double S0t = S0BP(exp(xibeta_q)*t, th1, th2, w, BP, dist);
  double logh0t = logf0BP(exp(xibeta_q)*t, th1, th2, w, BP, dist)-log(S0t);
  double den = exp(xibeta_o+xibeta_q)*(1.0-S0t) + exp(xibeta_h)*S0t;
  double ll = xibeta_o+xibeta_h+xibeta_q+logh0t-log(den);
  double tmp2 = 1+exp(xibeta_o-xibeta_h+xibeta_q)*(1.0/S0t-1.0);
  if(tmp2>SYSMAX) tmp2=SYSMAX;
  ll+= -exp(xibeta_h-xibeta_q)*log(tmp2);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double PHPOAFT_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                          bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q){
  double tmp1 = 1.0/S0BP(exp(xibeta_q)*t, th1, th2, w, BP, dist)-1.0;
  double tmp2 = 1+exp(xibeta_o-xibeta_h+xibeta_q)*tmp1;
  if(tmp2>SYSMAX) tmp2=SYSMAX;
  double ll=-exp(xibeta_h-xibeta_q)*log(tmp2);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double PHPOAFT_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q){
  double tmp1 = 1.0/S0BP(exp(xibeta_q)*t, th1, th2, w, BP, dist)-1.0;
  double tmp2 = 1+exp(xibeta_o-xibeta_h+xibeta_q)*tmp1;
  if(tmp2>SYSMAX) tmp2=SYSMAX;
  double ll=-exp(xibeta_h-xibeta_q)*log(tmp2);
  ll = log(1.0-exp(ll));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi)) 
double PHPOAFT_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                              bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q){
  double tmp1 = 1.0/S0BP(exp(xibeta_q)*t1, th1, th2, w, BP, dist)-1.0;
  double tmp2 = 1+exp(xibeta_o-xibeta_h+xibeta_q)*tmp1;
  double ll1=-exp(xibeta_h-xibeta_q)*log(tmp2);
  tmp1 = 1.0/S0BP(exp(xibeta_q)*t2, th1, th2, w, BP, dist)-1.0;
  tmp2 = 1+exp(xibeta_o-xibeta_h+xibeta_q)*tmp1;
  double ll2=-exp(xibeta_h-xibeta_q)*log(tmp2);
  double ll = log(std::abs(exp(ll1)-exp(ll2)));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data and parameters 
void PHPOAFT_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,  
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta_h, const Rcpp::NumericVector& Xbeta_o,
                       const Rcpp::NumericVector& Xbeta_q, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += PHPOAFT_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]);
    }else if(type[i]==1){
      ll += PHPOAFT_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]);
    }else if(type[i]==2){
      ll += PHPOAFT_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]);
    }else{
      ll += PHPOAFT_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]);
    }
    if(ltr[i]>0) ll += -PHPOAFT_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec PHPOAFT_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                            const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w, 
                            bool BP, int dist, const Rcpp::NumericVector& Xbeta_h, const Rcpp::NumericVector& Xbeta_o, 
                            const Rcpp::NumericVector& Xbeta_q){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(-PHPOAFT_BP_logsurv(t1[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]));
    }else if(type[i]==1){
      res[i] = exp(-PHPOAFT_BP_logpdf(t1[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]));
    }else if(type[i]==2){
      res[i] = exp(-PHPOAFT_BP_logcdf(t2[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]));
    }else{
      res[i] = exp(-PHPOAFT_BP_logsurvdiff(t1[i], t2[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]));
    }
    if(ltr[i]>0) res[i] *= exp(PHPOAFT_BP_logsurv(ltr[i], th1, th2, w, BP, dist, Xbeta_h[i], Xbeta_o[i], Xbeta_q[i]));
  }
  return(res);
}
