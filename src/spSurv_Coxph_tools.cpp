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
  Rcpp::NumericVector d(M1); d[0] = 0; d[M1-1] = ELARGE; 
  for(int k=1; k<(M1-1); ++k) d[k] = -std::log(1.0-(double)k/(M1-1.0))/hcen;
  return(d);
}

// Calculate picewise constant baseline cumulative hazard function
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
double Lambda0t(double t, Rcpp::NumericVector h, Rcpp::NumericVector d){
  if (t<=0) return(0);
  double Lamb = 0;
  int i = 1;
  while ( (t>d[i]) ){
    Lamb += (d[i]-d[i-1])*h[i];
    ++i;
  }
  Lamb += (t-d[i-1])*h[i];
  return(Lamb);
}

// Calculate picewise constant baseline cumulative hazard function at t=(t_1, ..., t_n)
// where h=(h0, h1, ..., hM) with h0=0 and d=(d0, d1, ..., dM) with d0=0, dM=R_PosInf
arma::vec Lambda0tvec(Rcpp::NumericVector t, Rcpp::NumericVector h, Rcpp::NumericVector d){
  int n = t.size();
  arma::vec res(n);
  for (int i=0; i<n; ++i){
    if (t[i]<=0) res[i] = 0;
    double Lamb = 0;
    int k = 1;
    while ( (t[i]>d[k]) ){
      Lamb += (d[k]-d[k-1])*h[k];
      ++k;
    }
    Lamb += (t[i]-d[k-1])*h[k];
    res[i] = Lamb;
  }
  return(res);
}

// Claculate cdf of survival time given xi
double Foft(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta){
  double res = 1- std::exp( -Lambda0t(t, h, d)*std::exp(xibeta) );
  return(res);
}

// Claculate pdf of survival time given xi
double foft(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta){
  if (t<=0) return(0);
  int k = 1;
  while ( (t>d[k]) ){
    ++k;
  }
  double res = std::exp( -Lambda0t(t, h, d)*std::exp(xibeta) )*h[k]*std::exp(xibeta);
  return(res);
}

// Claculate hazard function of survival time given xi
double lambdat(double t, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta){
  if (t<=0) return(0);
  int k = 1;
  while ( (t>d[k]) ){
    ++k;
  }
  double res = h[k]*std::exp(xibeta);
  return(res);
}

//calculate Fxi^{-1}(u)
double Finvofu(double u, Rcpp::NumericVector h, Rcpp::NumericVector d, double xibeta, double lower, double upper){
  double err = 10e-8;
  double tl = lower;
  double Fl = Foft(tl, h, d, xibeta) - u;
  double tstep = 10;
  double tr, Fr;
  if ( (Fl>0) | (std::abs(Fl)<err) ){
    return(lower);
  } else {
    tr=tl+10; Fr = Foft(tr, h, d, xibeta) - u;
    while( Fr<=0 ){
      if((std::abs(Fl)<err) | (tr>upper) ) return(tr);
      tl=tr; Fl=Fr;
      tr += tstep;
      Fr = Foft(tr, h, d, xibeta) - u;
      R_CheckUserInterrupt();
    }
  }
  double tm = (tl+tr)*0.5;
  while( (tr-tl)>err ){
    R_CheckUserInterrupt();
    tm = (tl + tr)*0.5;
    double Fm = Foft(tm, h, d, xibeta) - u;
    //if (std::abs(Fm)<err) break;
    if (Fl*Fm>0){
      tl = tm;
      Fl = Fm;
    } else {
      tr = tm;
    }
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
           const Rcpp::NumericVector& t, const arma::vec& Xbeta){
  int n = Mt.size();
  std::fill(lk.begin(), lk.end(), 0);
  for (int k=1; k<M1; ++k){
    for (int i=0; i<n; ++i){
      if(Mt[i]>=k) lk[k] += (std::min(d[k],t[i])-d[k-1])*std::exp(Xbeta[i]);
    }
  }
}

//Sample hcen using adaptive M-H for Cox PH with random cutpoints;
void sample_hcen(double& hcen, int& rejhcen, double& hSnew, double& hbarnew, Rcpp::NumericVector h, double r0, double h0, double V0,
                 double hl0, double hs0, double hadapter, int iscan){
  int M1 = h.size();
  double M = M1-1;
  double sumlogh=0;
  for(int k=1; k<M1; ++k) sumlogh += std::log(h[k]);
  double tempold = M*r0*hcen*std::log(r0) - M*Rf_lgammafn(r0*hcen) + r0*hcen*sumlogh - 0.5*std::pow((hcen-h0),2)/V0;
  double hcenold = hcen;
  if(iscan>hl0){
      hcen = Rf_rnorm(hcenold, std::sqrt(hSnew));
  }else{
      hcen = Rf_rnorm(hcenold, hs0);
  }
  if (hcen<ESMALL) {
    hcen = hcenold; ++rejhcen;
  } else {
    double tempnew = M*r0*hcen*std::log(r0) - M*Rf_lgammafn(r0*hcen) + r0*hcen*sumlogh - 0.5*std::pow((hcen-h0),2)/V0;
    double ratio = std::exp( tempnew-tempold );
    double uu = unif_rand();
    if (uu>ratio) {
      hcen = hcenold; ++rejhcen;
    }
  }
  double nn = iscan+1;
  double hbarold = hbarnew;
  hbarnew = (nn)/(nn+1.0)*hbarold + hcen/(nn+1.0);
  double hSold = hSnew;
  hSnew = (nn-1.0)/nn*hSold + hadapter/nn*(nn*hbarold*hbarold - (nn+1.0)*hbarnew*hbarnew + hcen*hcen + 0.01 );
}

//////////////////////////////////////////////////////////////////////
// Independent Cox PH
/////////////////////////////////////////////////////////////////////
//Sample beta using adaptive M-H for independent Cox PH;
void indept_sample_beta(arma::vec& beta, int& rejbeta, arma::mat& Snew, arma::vec& betabarnew, const arma::mat& X, 
                 const arma::vec& Lamb0, arma::vec mu0, arma::mat Sig0, int p, double l0, arma::mat S0, 
                 double adapter, int iscan){
  const arma::mat Ip = eye(p,p);
  arma::vec Xbeta = X.t()*beta;
  double tempold = arma::sum( -Lamb0%arma::exp(Xbeta) + Xbeta ) - 0.5*arma::dot(beta-mu0, arma::solve(Sig0, beta-mu0) );
  arma::vec betaold = beta;
  if(iscan>l0){
      beta = mvrnorm(betaold, Snew);
  }else{
      beta = mvrnorm(betaold, S0);
  }
  Xbeta = X.t()*beta;  
  double tempnew = arma::sum( -Lamb0%arma::exp(Xbeta) + Xbeta ) - 0.5*arma::dot(beta-mu0, arma::solve(Sig0, beta-mu0) );
  double ratio = std::exp( tempnew-tempold );
  // Rprintf( "%f\n", ratio );
  double uu = unif_rand();
  if (uu>ratio) {
    beta = betaold; ++rejbeta;
  }
  double nn = iscan+1;
  arma::vec betabarold = betabarnew;
  betabarnew = (nn)/(nn+1.0)*betabarold + beta/(nn+1.0);
  arma::mat Sold = Snew;
  Snew = (nn-1.0)/nn*Sold + adapter/nn*(nn*betabarold*betabarold.t() - (nn+1.0)*betabarnew*betabarnew.t() + beta*beta.t() + 0.05*Ip );
}

// Calculate CPO for Independent Cox PH
arma::vec LinvIndeptCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, arma::vec Xbeta, Rcpp::NumericVector h, Rcpp::NumericVector d){
  int n = delta.size();
  arma::vec res(n);
  for(int i=0; i<n; ++i){
    if(delta[i]==0){
      double St = 1.0-Foft(tobs[i], h, d, Xbeta[i]);
      res(i) = 1.0/St;
    } else {
      double ft = foft(tobs[i], h, d, Xbeta[i]);
      res(i) = 1/ft;
    }
  }
  return(res);
}

//////////////////////////////////////////////////////////////////////
// spatial Copula Cox PH
/////////////////////////////////////////////////////////////////////
// Claculate transformed survival time vector z
void Getz(arma::vec& z, const Rcpp::NumericVector& t, Rcpp::NumericVector h, Rcpp::NumericVector d, const arma::vec& Xbeta, int n){
  for(int i=0; i<n; ++i){
    double Fti = 1- std::exp( -Lambda0t(t[i], h, d)*std::exp(Xbeta[i]) );
    z[i] = Rf_qnorm5(Fti, 0, 1, true, false);
  }
}

//Sample beta using adaptive M-H for spatial Copula Cox PH;
void spCopula_sample_beta(arma::vec& beta, int& rejbeta, arma::mat& Snew, arma::vec& betabarnew, const arma::mat& X, 
                 const arma::vec& Lamb0, arma::vec mu0, arma::mat Sig0, int p, double l0, arma::mat S0, 
                 double adapter, int iscan, arma::vec& z, const arma::mat& Cinv, 
                 const Rcpp::NumericVector& t, Rcpp::NumericVector h, Rcpp::NumericVector d, int n){
  arma::vec Xbeta = X.t()*beta;
  double tempold = arma::sum( -Lamb0%arma::exp(Xbeta) + Xbeta ) - 0.5*arma::dot(beta-mu0, arma::solve(Sig0, beta-mu0) ) 
                   -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  arma::vec betaold = beta;
  if(iscan>l0){
    beta = mvrnorm(betaold, Snew);
  }else{
    beta = mvrnorm(betaold, S0);
  }
  Xbeta = X.t()*beta;
  Getz(z, t, h, d, Xbeta, n);
  double tempnew = arma::sum( -Lamb0%arma::exp(Xbeta) + Xbeta ) - 0.5*arma::dot(beta-mu0, arma::solve(Sig0, beta-mu0) )
                   -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  double ratio = std::exp( tempnew-tempold );
  // Rprintf( "%f\n", ratio );
  double uu = unif_rand();
  if (uu>ratio) {
    beta = betaold; ++rejbeta;
  }
  double nn = iscan+1;
  arma::vec betabarold = betabarnew;
  betabarnew = (nn)/(nn+1.0)*betabarold + beta/(nn+1.0);
  arma::mat Sold = Snew;
  Snew = (nn-1.0)/nn*Sold + adapter/nn*(nn*betabarold*betabarold.t() - (nn+1.0)*betabarnew*betabarnew.t() + beta*beta.t() + 0.05*eye(p,p) );
}

// Calculate CPO for spatial Copula Cox PH
arma::vec LinvSpCopulaCox(Rcpp::NumericVector tobs, Rcpp::IntegerVector delta, arma::vec Xbeta, Rcpp::NumericVector h, Rcpp::NumericVector d, 
                          arma::mat Cinv, arma::vec z){
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

// Get density or survival Plots for Cox PH
SEXP CoxPHplots(SEXP xpred_, SEXP tgrid_, SEXP beta_, SEXP h_, SEXP d_, SEXP probs_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  arma::mat xpred = as<mat>(xpred_); 
  arma::vec tgrid = as<vec>(tgrid_);
  Rcpp::NumericMatrix h(h_);
  Rcpp::NumericMatrix d(d_);
  arma::mat beta = as<mat>(beta_);
  NumericVector probs(probs_);
  int nsave = h.ncol();
  int ngrid = tgrid.size();
  int npred = xpred.n_rows;
  int low = nsave*probs[0]-1;
  int up = nsave*probs[1]-1;
  
  // Temp variables
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
  
  for(int i=0; i<nsave; ++i){
    arma::vec xibeta = xpred*beta.col(i);
    Rcpp::NumericVector hi = h(_,i);
    Rcpp::NumericVector di = d(_,i);
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        estf(k, i, j) = foft(tgrid[k], hi, di, xibeta[j]);
        estS(k, i, j) = 1.0 - Foft(tgrid[k], hi, di, xibeta[j]);
        esth(k, i, j) = lambdat(tgrid[k], hi, di, xibeta[j]);
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
