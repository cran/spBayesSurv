#include "spSurv_Coxph_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

//////////////////////////////////////////////////////////////////////
// Cox PH common
/////////////////////////////////////////////////////////////////////
// Get cutpoints from Exp(hcen)
Rcpp::NumericVector Cutpoints(double hcen, int M1){
  Rcpp::NumericVector d(M1); d[0] = 0; d[M1-1] = R_PosInf; 
  for(int k=1; k<(M1-1); ++k) d[k] = -std::log(1.0-(double)k/(M1-1.0))/hcen;
  return(d);
}

// Calculate picewise constant baseline cumulative hazard function
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
double Lambda0t(double t, Rcpp::NumericVector h, Rcpp::NumericVector d){
  if (t<=0) return(0);
  double Lamb = 0;
  int k = 1;
  while ( (t>d[k]) ){
    Lamb += (d[k]-d[k-1])*h[k];
    ++k;
  }
  Lamb += (t-d[k-1])*h[k];
  return(Lamb);
}

// Calculate picewise constant baseline cumulative hazard function at t=(t_1, ..., t_n)
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
arma::vec Lambda0tvec(Rcpp::NumericVector t, Rcpp::NumericVector h, Rcpp::NumericVector d){
  int n = t.size();
  arma::vec res(n);
  for (int i=0; i<n; ++i){
    res[i] = Lambda0t(t[i],h,d);
  }
  return(res);
}

// Claculate cdf of survival time given xi
double Foft(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta){
  double res = 1- std::exp( -Lambda0t(t, h, d)*std::exp(xibeta) );
  if(res<ESMALL) res=ESMALL;
  if(res>(1.0-ESMALL)) res=1.0-ESMALL;
  return(res);
}

// Claculate pdf of survival time given xi
double foft(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta){
  if (t<0) return(0);
  int k = 1;
  while ( (t>d[k]) ){
    ++k;
  }
  double res = (1.0-Foft(t,h,d,xibeta))*h[k]*std::exp(xibeta);
  return(res);
}

// Claculate hazard function of survival time given xi
double lambdat(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta){
  if (t<0) return(0);
  int k = 1;
  while ( (t>d[k]) ){
    ++k;
  }
  double res = h[k]*std::exp(xibeta);
  return(res);
}

//calculate Fxi^{-1}(u)
double Finvofu(double u, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta, double lower, double upper){
  double err = 1e-8;
  double tl = lower;
  double Fl = Foft(tl, h, d, xibeta) - u;
  double tr = upper;
  double Fr = Foft(tr, h, d, xibeta) - u;
  if (Fl>0) return(lower);
  if (Fr<0) return(upper);
  double tm = (tl+tr)*0.5;
  double Fm = Foft(tm, h, d, xibeta) - u;
  while( std::abs(Fm)>err ){
    //R_CheckUserInterrupt();
    if (Fl*Fm>0){
      tl = tm;
      Fl = Fm;
    } else {
      tr = tm;
    }
    tm = (tl + tr)*0.5;
    Fm = Foft(tm, h, d, xibeta) - u;
  }
  return(tm);
}

// Calculate  M(ti), i=1, ..., n;
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
void GetMt(Rcpp::IntegerVector& Mt, const Rcpp::NumericVector& t, Rcpp::NumericVector d){
  int n = t.size();
  for (int i=0; i<n; ++i){
    int k = 1;
    while ( (t[i]>d[k]) ){
      ++k;
    }
    Mt[i] = k;
  }
}

// Calculate mk = sum_i I(M(ti)=k), k=1, ..., M with m0=0;
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
void Getmk(Rcpp::IntegerVector& mk, const Rcpp::IntegerVector& Mt){
  int n = Mt.size();
  std::fill(mk.begin(), mk.end(), 0);
  for (int i=0; i<n; ++i){
    int k = Mt[i];
    mk[k] +=1;
  }
}

// Calculate lk, k=1, ..., M with m0=0;
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
void Getlk(Rcpp::NumericVector& lk, const Rcpp::IntegerVector& Mt, int M1, Rcpp::NumericVector d, 
           const Rcpp::NumericVector& t, const Rcpp::NumericVector& Xbeta){
  int n = Mt.size();
  std::fill(lk.begin(), lk.end(), 0);
  for (int k=1; k<M1; ++k){
    for (int i=0; i<n; ++i){
      if(Mt[i]>=k) lk[k] += (std::min(d[k],t[i])-d[k-1])*std::exp(Xbeta[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Independent Cox PH
/////////////////////////////////////////////////////////////////////
// Calculate CPO for Independent Cox PH
arma::vec LinvIndeptCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, Rcpp::NumericVector Xbeta,
                        Rcpp::NumericVector h, Rcpp::NumericVector d){
  int n = delta.size();
  arma::vec res(n);
  for(int i=0; i<n; ++i){
    if(delta[i]==0){
      res[i] = 1.0/(1.0-Foft(tobs[i], h, d, Xbeta[i]));
    }else{
      res[i] = 1.0/foft(tobs[i], h, d, Xbeta[i]);
    }
  }
  return(res);
}

//////////////////////////////////////////////////////////////////////
// spatial Copula Cox PH
/////////////////////////////////////////////////////////////////////
// Claculate transformed survival time vector z
void Getz(arma::vec& z, const Rcpp::NumericVector& t, Rcpp::NumericVector h, Rcpp::NumericVector d, 
          const Rcpp::NumericVector& Xbeta, int n){
  for(int i=0; i<n; ++i){
    double Fti = Foft(t[i], h, d, Xbeta[i]);
    z[i] = Rf_qnorm5(Fti, 0, 1, true, false);
  }
}

// Calculate CPO for spatial Copula Cox PH
arma::vec LinvSpCopulaCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, Rcpp::NumericVector Xbeta, 
                          Rcpp::NumericVector h, Rcpp::NumericVector d, arma::mat Cinv, arma::vec z){
  int n = delta.size();
  arma::vec res(n);
  for(int i=0; i<n; ++i){
    double Cii = Cinv(i,i);
    double s2i = 1.0/Cii;
    double mui = -s2i*( arma::dot(Cinv.col(i), z) - Cii*z[i] );
    double si = std::sqrt(s2i);
    double Fi = Foft(tobs[i], h, d, Xbeta[i]);
    double PinvFyi = Rf_qnorm5(Fi, 0, 1, true, false);
    double newzi = (PinvFyi-mui)/si;
    if(delta[i]==0){
      double St = 1.0-Rf_pnorm5(newzi, 0, 1, true, false);
      //Rprintf( "St=%f\n", St );
      res(i) = 1.0/St;
    } else {
      double fi = foft(tobs[i], h, d, Xbeta[i]);
      double diff = -0.5*std::pow(newzi, 2) + 0.5*std::pow(PinvFyi, 2);
      double ft = (1.0/si*std::exp(diff)*fi);
      res(i) = 1.0/ft;
      //Rprintf( "ft=%f\n", ft );
    }
  }
  return(res);
}
