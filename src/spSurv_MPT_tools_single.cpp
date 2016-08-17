#include "spSurv_MPT_tools_single.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PHlogpdf(double t, double th1, double th2, Rcpp::NumericVector probs,
                int maxL, bool MPT, int dist, double xibeta){
  double ll = xibeta + logf0MPT(t, th1, th2, probs, maxL, MPT, dist);
  ll += (exp(xibeta)-1.0)*log(S0MPT(t, th1, th2, probs, maxL, MPT, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double PHlogsurv(double t, double th1, double th2, Rcpp::NumericVector probs, 
                 int maxL, bool MPT, int dist, double xibeta){
  double ll = exp(xibeta)*log(S0MPT(t, th1, th2, probs, maxL, MPT, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double PHlogcdf(double t, double th1, double th2, Rcpp::NumericVector probs, 
                int maxL, bool MPT, int dist, double xibeta){
  double ll = log( 1.0-exp( exp(xibeta)*log(S0MPT(t, th1, th2, probs, maxL, MPT, dist)) ) );
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double PHlogsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector probs, 
                     int maxL, bool MPT, int dist, double xibeta){
  double St1 = exp( exp(xibeta)*log(S0MPT(t1, th1, th2, probs, maxL, MPT, dist)) );
  double St2 = exp( exp(xibeta)*log(S0MPT(t2, th1, th2, probs, maxL, MPT, dist)) );
  double ll = log(std::abs(St1 - St2));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void PHloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
              const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
              int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += PHlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += PHlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += PHlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else{
      ll += PHlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -PHlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void PHloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                    const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
                    int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll, 
                    int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += PHlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += PHlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += PHlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else{
      ll += PHlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -PHlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec PHinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs, 
                   int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-PHlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-PHlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-PHlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-PHlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(PHlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
  }
  return(res);
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double POlogpdf(double t, double th1, double th2, Rcpp::NumericVector probs,
                int maxL, bool MPT, int dist, double xibeta){
  double ll = logf0MPT(t, th1, th2, probs, maxL, MPT, dist)-xibeta;
  ll += -2.0*log(1.0+(exp(-xibeta)-1.0)*S0MPT(t, th1, th2, probs, maxL, MPT, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double POlogsurv(double t, double th1, double th2, Rcpp::NumericVector probs, 
                 int maxL, bool MPT, int dist, double xibeta){
  double S0t = S0MPT(t, th1, th2, probs, maxL, MPT, dist);
  double ll = log(S0t)-xibeta-log(1.0+(exp(-xibeta)-1.0)*S0t);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double POlogcdf(double t, double th1, double th2, Rcpp::NumericVector probs, 
                int maxL, bool MPT, int dist, double xibeta){
  double S0t = S0MPT(t, th1, th2, probs, maxL, MPT, dist);
  double ll = log(1.0-S0t)-log(1.0+(exp(-xibeta)-1.0)*S0t);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double POlogsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector probs, 
                     int maxL, bool MPT, int dist, double xibeta){
  double S0t1 = S0MPT(t1, th1, th2, probs, maxL, MPT, dist);
  double S0t2 = S0MPT(t2, th1, th2, probs, maxL, MPT, dist);
  double exibeta = exp(-xibeta);
  double ll = log(std::abs(exibeta*S0t1/(1.0+(exibeta-1.0)*S0t1) - exibeta*S0t2/(1.0+(exibeta-1.0)*S0t2)));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void POloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
              const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
              int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += POlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += POlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += POlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else{
      ll += POlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -POlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void POloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                    const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
                    int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll, 
                    int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += POlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += POlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += POlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else{
      ll += POlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -POlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec POinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs, 
                   int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-POlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-POlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-POlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-POlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(POlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
  }
  return(res);
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AFTlogpdf(double t, double th1, double th2, Rcpp::NumericVector probs,
                int maxL, bool MPT, int dist, double xibeta){
  double ll = xibeta + logf0MPT(exp(xibeta)*t, th1, th2, probs, maxL, MPT, dist);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
double AFTlogsurv(double t, double th1, double th2, Rcpp::NumericVector probs, 
                 int maxL, bool MPT, int dist, double xibeta){
  double ll = log(S0MPT(exp(xibeta)*t, th1, th2, probs, maxL, MPT, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
double AFTlogcdf(double t, double th1, double th2, Rcpp::NumericVector probs, 
                int maxL, bool MPT, int dist, double xibeta){
  double ll = log(1.0-S0MPT(exp(xibeta)*t, th1, th2, probs, maxL, MPT, dist));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
double AFTlogsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector probs, 
                     int maxL, bool MPT, int dist, double xibeta){
  double St1 = S0MPT(exp(xibeta)*t1, th1, th2, probs, maxL, MPT, dist);
  double St2 = S0MPT(exp(xibeta)*t2, th1, th2, probs, maxL, MPT, dist);
  double ll = log(std::abs(St1 - St2));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log likelihood given data, frailties and parameters 
void AFTloglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
              const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
              int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll){
  ll=0;
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      ll += AFTlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else if(type[i]==1){
      ll += AFTlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else if(type[i]==2){
      ll += AFTlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }else{
      ll += AFTlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
    }
    if(ltr[i]>0) ll += -AFTlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]);
  }
}
// log likelihood given frailties, parameters and data of block i
void AFTloglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                    const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs,
                    int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta, double& ll, 
                    int ind1, int ind2, double vi){
  ll=0;
  for(int i=ind1; i<=ind2; ++i){
    if(type[i]==0){
      ll += AFTlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else if(type[i]==1){
      ll += AFTlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else if(type[i]==2){
      ll += AFTlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }else{
      ll += AFTlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
    }
    if(ltr[i]>0) ll += -AFTlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]+vi);
  }
}
// Calculate 1.0/likelihood for CPO
arma::vec AFTinvLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector probs, 
                   int maxL, bool MPT, int dist, const Rcpp::NumericVector& Xbeta){
  arma::vec res(type.size());
  for(int i=0; i<type.size(); ++i){
    if(type[i]==0){
      res[i] = exp(0.0-AFTlogsurv(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else if(type[i]==1){
      res[i] = exp(0.0-AFTlogpdf(t1[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else if(type[i]==2){
      res[i] = exp(0.0-AFTlogcdf(t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }else{
      res[i] = exp(0.0-AFTlogsurvdiff(t1[i], t2[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
    }
    if(ltr[i]>0) res[i] *= exp(AFTlogsurv(ltr[i], th1, th2, probs, maxL, MPT, dist, Xbeta[i]));
  } 
  return(res);
}
