#include "spSurv_DDP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

//////////////////////////////////////////////////////////////////////
// DDP common
/////////////////////////////////////////////////////////////////////
// Caculate CDF at points ygrid
arma::vec Fmix(arma::vec ygrid, arma::vec mu, arma::vec sig, arma::vec w){
  int ngrid = ygrid.size();
  int N = sig.size();
  arma::vec res(ngrid);
  for (int i=0; i<ngrid; ++i){
  arma::vec Norm(N);
    for (int k=0; k<N; ++k){
      Norm[k] = Rf_pnorm5(ygrid[i], mu[k], sig[k], true, false);
    }
  res[i] = arma::dot(w, Norm);
  }
  return(res);
}

// Caculate pdf at points ygrid
arma::vec fmix(arma::vec ygrid, arma::vec mu, arma::vec sig, arma::vec w){
  int ngrid = ygrid.size();
  int N = sig.size();
  arma::vec res(ngrid);
  for (int i=0; i<ngrid; ++i){
  arma::vec Norm(N);
    for (int k=0; k<N; ++k){
      Norm[k] = Rf_dnorm4(ygrid[i], mu[k], sig[k], false);
    }
	res[i] = arma::dot(w, Norm);
  }
  return(res);
}

// calculate marginal CDF of y: a mixture of normals
double Fofy(double y, arma::vec w, arma::vec mu, arma::vec sig){
  int N = sig.size();
  arma::vec Phi(N);
  for (int k=0; k<N; ++k){
    Phi[k] = Rf_pnorm5(y, mu[k], sig[k], true, false);
  }
  double res = arma::dot(w, Phi);
  return(res);
}

// calculate marginal pdf of y: a mixture of normals
double fofy(double y, arma::vec w, arma::vec mu, arma::vec sig){
  int N = sig.size();
  arma::vec Phi(N);
  for (int k=0; k<N; ++k){
    Phi[k] = Rf_dnorm4(y, mu[k], sig[k], false);
  }
  double res = arma::dot(w, Phi);
  return(res);
}

// calculate F^{-1}(u)
double DDP_Finvofu(double u, arma::vec wma, arma::vec mu, arma::vec sig, double lower, double upper){
  double err = 1e-8;
  double yl = lower;
  double Fl = Fofy(yl, wma, mu, sig) - u;
  double yr = upper;
  double Fr = Fofy(yr, wma, mu, sig) - u;
  if (Fl>0) return(lower);
  if (Fr<0) return(upper);
  double ym = (yl+yr)*0.5;
  double Fm = Fofy(ym, wma, mu, sig) - u;
  while( std::abs(Fm)>err ){
    R_CheckUserInterrupt();
    if (Fl*Fm>0){
      yl = ym;
      Fl = Fm;
    } else {
      yr = ym;
    }
    ym = (yl + yr)*0.5;
    Fm = Fofy(ym, wma, mu, sig) - u;
  }
  return(ym);
}

// From V to w
void DDP_Vtow(arma::vec& w, Rcpp::NumericVector V, int N){
  double temp=0.0;
  w[0] = V[0];
  for (int k=1; k<N; ++k){
    temp += std::log(1.0-V[k-1]);
    w[k] = std::max(std::exp(temp+std::log(V[k])), 1e-320);
    //w[k] = std::exp(temp+std::log(V[k]));
  }
}

// sample(1:N, prob=w), where w.size()=N
int DDP_sample(arma::vec w){
  int k = 1;
  double u = unif_rand();;
  double cdf = w[0];
  while(u>cdf){
    cdf += w[k];
    ++k;
  }
  return (k);
}

//Sample K;
void DDP_sample_K(Rcpp::IntegerVector& K, const Rcpp::NumericVector& y, const arma::mat& Xbeta, 
                  arma::vec w, Rcpp::NumericVector tau, int n, int N){
  for (int i=0; i<n; ++i){
      arma::vec phi(N);
      for (int k=0; k<N; ++k){
        phi[k] = std::max( Rf_dnorm4(y[i], Xbeta(i,k), 1.0/tau[k], false), ESMALL);
      }
      arma::vec wi = w%phi/(arma::sum(w%phi));
      K[i] = DDP_sample(wi);
  }
}

// calculate pnorm((y-mu)/sig) for vectors of mu and sig
arma::vec Phivec(double y, arma::vec mu, NumericVector sig){
  int N = sig.size();
  arma::vec res(N);
  for (int k=0; k<N; ++k){
    res[k] = Rf_pnorm5((y-mu[k])/sig[k], 0, 1.0, true, false);
  }
  return(res);
}

//////////////////////////////////////////////////////////////////////
// ANOVA DDP
/////////////////////////////////////////////////////////////////////

// Sample beta;
void anovaDDP_sample_beta(arma::mat& beta, const Rcpp::NumericVector& y, const arma::mat& X, 
                     const Rcpp::NumericVector& tau2, const Rcpp::IntegerVector& nK, 
                     const Rcpp::IntegerMatrix& Kind, arma::vec mu, arma::mat Sig, 
                     arma::mat invSig, int N, int p){
  for (int k=0; k<N; ++k){
    if (nK[k]>0){
      arma::mat xxnk(p,p); xxnk.fill(0.0);
      arma::vec xynk(p); xynk.fill(0.0);
      for (int ii=0; ii<nK[k]; ++ii){
        int i = Kind(ii, k);
        arma::vec Xi = X.col(i);
        xxnk += Xi*(arma::trans(Xi));
        xynk += Xi*y[i];
      }
      arma::mat Sigk = arma::inv_sympd(invSig + xxnk*tau2[k]);
      arma::vec muk = Sigk*(invSig*mu + xynk*tau2[k]);
      beta.col(k) = mvrnorm(muk, Sigk);
      }
    else{
      beta.col(k) = mvrnorm(mu, Sig);
    }
  }
}

//Sample simga2;
void anovaDDP_sample_sigma2(Rcpp::NumericVector& tau2, const Rcpp::NumericVector& y, const arma::mat& Xbeta, 
                     const Rcpp::IntegerVector& nK, const Rcpp::IntegerMatrix& Kind, double nua, double nub, int N){
  for (int k=0; k<N; ++k){
    if (nK[k]>0){
      double sumtemp=0;
      for (int ii=0; ii<nK[k]; ++ii){
        int i = Kind(ii, k);
        sumtemp  += pow( (y[i]-Xbeta(i,k)), 2);
      }
      double nuak = nua + nK[k]*0.5;
      double nubk = nub + 0.5*sumtemp;
      tau2[k] = Rf_rgamma(nuak, 1.0/nubk);
      }
    else{
      tau2[k] = Rf_rgamma(nua, 1.0/nub);
    }
  }
}

// Calculate CPO for anovaDDP
arma::vec anovaDDP_Linv(Rcpp::NumericVector yobs, Rcpp::IntegerVector delta, arma::mat X, arma::mat beta, arma::vec sig, arma::vec wei ){
  int n = delta.size();
  arma::vec res(n);
  arma::mat xbeta = arma::trans( beta )*X;
  for(int i=0; i<n; ++i){
    if(delta[i]==0){
      double Sy = 1.0-Fofy(yobs(i), wei, xbeta.col(i), sig);
      res(i) = 1.0/Sy;
    } else {
      double fy = fofy(yobs(i), wei, xbeta.col(i), sig);
      res(i) = std::exp(yobs(i))/fy;
    }
  }
  return(res);
}

//////////////////////////////////////////////////////////////////////
// spatial Copula DDP
/////////////////////////////////////////////////////////////////////
// Sample y_i when delta_i=0
void spCopula_sample_y(Rcpp::NumericVector& y, Rcpp::NumericVector& rejy, arma::mat& zPhi, arma::vec& z, arma::vec w, 
      const Rcpp::NumericVector& yobs, const Rcpp::IntegerVector& delta, const arma::mat& Xbeta, Rcpp::NumericVector tau, 
      Rcpp::IntegerVector K, const arma::mat& Cinv, int n, int N, int iscan, int nburn){
  for (int i=0; i<n; ++i){
    if(delta[i]==0){
      double yold = y[i];
      double zold = z[i];
      arma::rowvec zPhiold = zPhi.row(i);
      double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double mui = Xbeta(i, K[i]-1);
      y[i] = trun_rnorm(mui, 1.0/tau[K[i]-1], yobs[i], R_PosInf);
      for (int k=0; k<N;++k){
        zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*tau[k], 0, 1, true, false);
      }
      z[i] = Rf_qnorm5(arma::as_scalar(zPhi.row(i)*w), 0, 1, true, false);
      double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tempnew-tempold);
      double uu = unif_rand();
      if (uu>ratio) {
        y[i] = yold; zPhi.row(i)=zPhiold; z[i] = zold;
        if(iscan>=nburn) rejy[i]+=1.0;
      }
    }
  }
}

// Sample beta;
void spCopula_sample_beta(arma::mat& beta, Rcpp::NumericVector& rejbeta, arma::mat& zPhi, arma::vec& z, arma::vec w, 
                          const Rcpp::NumericVector& y, const arma::mat& X, Rcpp::NumericVector tau2, 
                          const Rcpp::IntegerVector& nK, const Rcpp::IntegerMatrix& Kind, arma::vec mu, 
                          arma::mat Sig, arma::mat invSig, const arma::mat& Cinv, int n, int N, int p, int iscan, int nburn){
        
  Rcpp::NumericVector tau = Rcpp::sqrt(tau2);
  for (int k=0; k<N; ++k){
    if (nK[k]>0){
      arma::mat xxnk(p,p); xxnk.fill(0.0);
      arma::vec xynk(p); xynk.fill(0.0);
      for (int ii=0; ii<nK[k]; ++ii){
        int i = Kind(ii, k);
        arma::vec Xi = X.col(i);
        xxnk += Xi*(trans(Xi));
        xynk += Xi*y[i];
      }
      arma::mat Sigk = arma::inv_sympd(invSig + xxnk*tau2[k]);
      arma::vec muk = Sigk*(invSig*mu + xynk*tau2[k]);
      arma::vec betakold = beta.col(k);
      arma::vec zold = z;
      arma::vec zPhikold = zPhi.col(k);
      beta.col(k) = mvrnorm(muk, Sigk);
      double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      for (int i=0; i<n; ++i){
        zPhi(i,k) = Rf_pnorm5((y[i]-arma::dot(X.col(i),beta.col(k)))*tau[k], 0, 1.0, true, false);
      }
      z = qnormvec( zPhi*w );
      double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tempnew-tempold);
      double uu = unif_rand(); 
      if (uu>ratio){
        arma::vec betaknew2 = mvrnorm(betakold, Sigk);
        double tmpnew2 = -0.5*arma::dot( (betaknew2-muk), solve(Sigk, (betaknew2-muk)) );
        double tmpold0 = -0.5*arma::dot( (betakold-muk), solve(Sigk, (betakold-muk)) );
        beta.col(k) = betaknew2;
        for (int i=0; i<n; ++i){
          zPhi(i,k) = Rf_pnorm5((y[i]-arma::dot(X.col(i),betaknew2))*tau[k], 0, 1.0, true, false);
        }
        z = qnormvec( zPhi*w );
        double tempnew2 = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
        double ratio2 = std::exp( tmpnew2 - tmpold0 )*(std::exp(tempnew2-tempnew)-1)/(std::exp(tempold-tempnew)-1);
        double uu2 = unif_rand();
        if (uu2>ratio2) {
          beta.col(k) = betakold; zPhi.col(k)=zPhikold; z=zold;
          if(iscan>=nburn) rejbeta[k]+=1.0;
        } 
      }
    } else{
      arma::vec betakold = beta.col(k);
      arma::vec zold = z;
      arma::vec zPhikold = zPhi.col(k);
      beta.col(k) = mvrnorm(mu, Sig);
      double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      for (int i=0; i<n; ++i){
        zPhi(i,k) = Rf_pnorm5((y[i]-arma::dot(X.col(i),beta.col(k)))*tau[k], 0, 1.0, true, false);
      }
      z = qnormvec( zPhi*w );
      double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tempnew-tempold);
      double uu = unif_rand(); 
      if (uu>ratio) {
        beta.col(k) = betakold; zPhi.col(k)=zPhikold; z=zold;
        if(iscan>=nburn) rejbeta[k]+=1.0;
      }
    }
  }
}

// Sample beta blockwise;
void spCopula_sample_beta_block(arma::mat& beta, Rcpp::NumericVector& rejbeta, arma::mat& zPhi, arma::vec& z, arma::vec w, 
      const Rcpp::NumericVector& y, const arma::mat& X, Rcpp::NumericVector tau2, const Rcpp::IntegerVector& nK, 
      const Rcpp::IntegerMatrix& Kind, arma::vec mu, arma::mat Sig, arma::mat invSig, const arma::mat& Cinv, int n, int N, int p){
  
  arma::mat betaold = beta;
  arma::vec zold = z;
  arma::mat zPhiold = zPhi;
  for (int k=0; k<N; ++k){
    if (nK[k]>0){
      arma::mat xxnk(p,p); xxnk.fill(0.0);
      arma::vec xynk(p); xynk.fill(0.0);
      for (int ii=0; ii<nK[k]; ++ii){
        int i = Kind(ii, k);
        arma::vec Xi = X.col(i);
        xxnk += Xi*(arma::trans(Xi));
        xynk += Xi*y[i];
      }
      arma::mat Sigk = arma::inv_sympd(invSig + xxnk*tau2[k]);
      arma::vec muk = Sigk*(invSig*mu + xynk*tau2[k]);
      beta.col(k) = mvrnorm(muk, Sigk);
      }
    else{
      beta.col(k) = mvrnorm(mu, Sig);
    }
  }
  double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  arma::mat Xbeta = X.t()*beta;
  for (int i=0; i<n; ++i){
    for (int k=0; k<N; ++k){
      zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*std::sqrt(tau2[k]), 0, 1.0, true, false);
    }
  }
  z = qnormvec( zPhi*w );
  double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  double ratio = std::exp(tempnew-tempold);
  double uu = unif_rand();
  double phinew = 0;
  double phiold = 0;
  if (uu>ratio){
    for (int k=0; k<N; ++k){
      arma::vec betakold = betaold.col(k);
      if (nK[k]>0){
        arma::mat xxnk(p,p); xxnk.fill(0.0);
        arma::vec xynk(p); xynk.fill(0.0);
        for (int ii=0; ii<nK[k]; ++ii){
          int i = Kind(ii, k);
          arma::vec Xi = X.col(i);
          xxnk += Xi*(arma::trans(Xi));
          xynk += Xi*y[i];
        }
        arma::mat Sigk = arma::inv_sympd(invSig + xxnk*tau2[k]);
        arma::vec muk = Sigk*(invSig*mu + xynk*tau2[k]);
        phiold += -0.5*arma::dot( (betakold-muk), solve(Sigk, (betakold-muk)) );
        beta.col(k) = mvrnorm(betakold, Sigk);
        phinew += -0.5*arma::dot( (beta.col(k)-muk), solve(Sigk, (beta.col(k)-muk)) );
      }
      else{
        //phiold += -0.5*arma::dot( (betakold-mu), solve(Sig, (betakold-mu)) );
        beta.col(k) = mvrnorm(betakold, Sig);
        //phinew += -0.5*arma::dot( (beta.col(k)-mu), solve(Sig, (beta.col(k)-mu)) );
      }
    }
    arma::mat Xbeta = X.t()*beta;
    for (int i=0; i<n; ++i){
      for (int k=0; k<N; ++k){
        zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*std::sqrt(tau2[k]), 0, 1.0, true, false);
      }
    }
    z = qnormvec( zPhi*w );
    double tempnew2 = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
    double ratio2 = std::exp( phinew - phiold )*(std::exp(tempnew2-tempnew)-1)/(std::exp(tempold-tempnew)-1);
    double uu2 = unif_rand();
    if (uu2>ratio2) {beta = betaold; rejbeta=rejbeta+1.0; zPhi=zPhiold; z=zold;}
  }
}

//Sample simga2;
void spCopula_sample_sigma2(Rcpp::NumericVector& tau2, Rcpp::NumericVector& rejsigma, arma::mat& zPhi, arma::vec& z, 
                            arma::vec w, const Rcpp::NumericVector& y, const arma::mat& Xbeta, const Rcpp::IntegerVector& nK, 
                            const Rcpp::IntegerMatrix& Kind, double nua, double nub, const arma::mat& Cinv, int n, int N, 
                            int iscan, int nburn){
  for (int k=0; k<N; ++k){
    if (nK[k]>0){
      double sumtemp=0;
      for (int ii=0; ii<nK[k]; ++ii){
        int i = Kind(ii, k);
        sumtemp  += pow( (y[i]-Xbeta(i,k)), 2);
      }
      double nuak = nua + nK[k]*0.5;
      double nubk = nub + 0.5*sumtemp;
      double tau2kold = tau2[k];
      arma::vec zold = z;
      arma::vec zPhikold = zPhi.col(k);
      tau2[k] = Rf_rgamma(nuak, 1.0/nubk);
      double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      for (int i=0; i<n; ++i){
        zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*std::sqrt(tau2[k]), 0, 1.0, true, false);
      }
      z = qnormvec( zPhi*w );
      double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tempnew-tempold);
      double uu = unif_rand();
      if (uu>ratio) {
        tau2[k] = tau2kold; zPhi.col(k)=zPhikold; z=zold; 
        if(iscan>=nburn) rejsigma[k]+=1.0;
      }
    } else{
      double tau2kold = tau2[k];
      arma::vec zold = z;
      arma::vec zPhikold = zPhi.col(k);
      tau2[k] = Rf_rgamma(nua, 1/nub);
      double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      for (int i=0; i<n; ++i){
        zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*std::sqrt(tau2[k]), 0, 1.0, true, false);
      }
      z = qnormvec( zPhi*w );
      double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tempnew-tempold);
      double uu = unif_rand();
      if (uu>ratio) {
        tau2[k] = tau2kold; zPhi.col(k)=zPhikold; z=zold; 
        if(iscan>=nburn) rejsigma[k]+=1.0;
      }
    }
  }
}

//Sample simga2 blockwise;
void spCopula_sample_sigma2_block(Rcpp::NumericVector& tau2, Rcpp::NumericVector& rejsigma, arma::mat& zPhi, arma::vec& z, 
      arma::vec w, const Rcpp::NumericVector& y, const arma::mat& Xbeta, const Rcpp::IntegerVector& nK, 
      const Rcpp::IntegerMatrix& Kind, double nua, double nub, const arma::mat& Cinv, int n, int N){
  Rcpp::NumericVector tau2old = tau2;
  arma::vec zold = z;
  arma::mat zPhiold = zPhi;
  for (int k=0; k<N; ++k){
    if (nK[k]>0){
      double sumtemp=0;
      for (int ii=0; ii<nK[k]; ++ii){
        int i = Kind(ii, k);
        sumtemp  += pow( (y[i]-Xbeta(i,k)), 2);
      }
      double nuak = nua + nK[k]*0.5;
      double nubk = nub + 0.5*sumtemp;
      tau2[k] = Rf_rgamma(nuak, 1.0/nubk);
      }
    else{
      tau2[k] = Rf_rgamma(nua, 1.0/nub);
    }
  }
  double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  for (int i=0; i<n; ++i){
    for (int k=0; k<N; ++k){
      zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*std::sqrt(tau2[k]), 0, 1.0, true, false);
    }
  }
  z = qnormvec( zPhi*w );
  double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  double ratio = std::exp(tempnew-tempold);
  double uu = unif_rand();
  if (uu>ratio) {tau2 = tau2old; rejsigma=rejsigma+1.0; zPhi=zPhiold; z=zold;}
}

//Sample V;
void spCopula_sample_V(Rcpp::NumericVector& V, Rcpp::NumericVector& rejV, arma::mat& zPhi, arma::vec& z, arma::vec& w, 
                       const Rcpp::IntegerVector& nK, double alpha, const arma::mat& Cinv, int N, int iscan, int nburn){
  arma::vec nkk = as<vec>(nK);
  for (int k=0; k<(N-1); ++k){
    double alphak = alpha + arma::sum(nkk.subvec(k+1, N-1))+ESMALL;
    double aa = nK[k] + 1.0;
    double Vkold = V[k];
    arma::vec zold = z;
    V[k] = Rf_rbeta(aa, alphak);
    double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
    DDP_Vtow(w, V, N); // From V to w
    z = qnormvec( zPhi*w );
    double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
    double ratio = std::exp(tempnew-tempold);
    double uu = unif_rand();
    if (uu>ratio) {
      V[k] = Vkold; z=zold; 
      if(iscan>=nburn) rejV[k]+=1.0;
    }    
  }
}

//Sample V blockwise;
void spCopula_sample_V_block(Rcpp::NumericVector& V, Rcpp::NumericVector& rejV, arma::mat& zPhi, arma::vec& z, arma::vec& w, 
      const Rcpp::IntegerVector& nK, double alpha, const arma::mat& Cinv, int N){
  arma::vec nkk = as<vec>(nK);
  Rcpp::NumericVector Vkold=V;
  arma::vec zold = z;
  for (int k=0; k<(N-1); ++k){
    double alphak = alpha + arma::sum(nkk.subvec(k+1, N-1))+ESMALL;
    double aa = nK[k] + 1.0;
    V[k] = Rf_rbeta(aa, alphak);
  }
  double tempold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  DDP_Vtow(w, V, N); // From V to w
  z = qnormvec( zPhi*w );
  double tempnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
  double ratio = std::exp(tempnew-tempold);
  double uu = unif_rand();
  if (uu>ratio) {V = Vkold; rejV=rejV+1.0; z=zold; }
}

// Calculate CPO for spatial Copula DDP
arma::vec spCopula_Linv(Rcpp::NumericVector yobs, Rcpp::IntegerVector delta, arma::mat X, arma::mat beta, arma::vec sig, arma::vec wei, 
                    arma::mat Cinv, arma::vec zj){
  int n = delta.size();
  arma::vec res(n);
  arma::mat xbeta = arma::trans( beta )*X;
  for(int i=0; i<n; ++i){
    double Cii = Cinv(i,i);
    double s2i = 1.0/Cii;
    double si = std::sqrt(s2i);
    double mui = -s2i*( arma::dot(Cinv.col(i), zj) - Cii*zj(i) );
    double Fyi = Fofy(yobs(i), wei, xbeta.col(i), sig);
    double PinvFyi = Rf_qnorm5(Fyi, 0, 1, true, false);
    double newzi = (PinvFyi-mui)/si;
    if(delta[i]==0){
      double Sy = 1.0 - Rf_pnorm5(newzi, 0, 1, true, false);
      res(i) = 1.0/Sy;
    } else {
      double fyi = fofy(yobs(i), wei, xbeta.col(i), sig);
      double diff = -0.5*std::pow(newzi, 2) + 0.5*std::pow(PinvFyi, 2);
      double fy = 1.0/si*std::exp(diff)*fyi;
      res(i) = std::exp(yobs(i))/fy;
    }
  }
  return(res);
}

// using anovaDDP to get initial chain
void spCopulaInitial(Rcpp::IntegerVector& K, Rcpp::NumericVector& y, arma::mat& beta, Rcpp::NumericVector& tau2, Rcpp::NumericVector& V,
                  arma::vec& w, double& alpha, arma::vec& mu, arma::mat& Sig, arma::mat& invSig, const Rcpp::NumericVector& yobs,
                  const Rcpp::IntegerVector& delta, const arma::mat& X, arma::vec m0, arma::mat S0, arma::mat Sig0, int k0, 
                  double a0, double b0, double nua, double nub, arma::mat invS0, arma::vec invS0m0 ){
  int N = tau2.size();
  int n = X.n_cols;
  int p = X.n_rows;
  IntegerVector nK(N);
  IntegerMatrix Kind(n, N);
  Rcpp::NumericVector tau = Rcpp::sqrt(tau2);
  arma::mat Xbeta = X.t()*beta;
  for (int iscan=0; iscan<50000; ++iscan){
	  R_CheckUserInterrupt();
    
    //Sample K;
    DDP_sample_K(K, y, Xbeta, w, tau, n, N);
  
    // Calculate nK and Kind
    std::fill(Kind.begin(), Kind.end(), 0);
    std::fill(nK.begin(), nK.end(), 0);
    for(int i=0; i<n; ++i){
      int k = K[i]-1;
      nK[k] += 1;
      Kind(nK[k]-1, k) = i;
    }
    
    // Sample y_i when delta_i=0
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        y[i] = trun_rnorm(Xbeta(i, K[i]-1), 1.0/tau[K[i]-1], yobs[i], R_PosInf);
      }
    }
    
    //Sample V;
    vec nkk = as<vec>(nK);
    for (int k=0; k<(N-1); ++k){
      double alphak = alpha + arma::sum(nkk.subvec(k+1, N-1));
      double aa = nK[k] + 1.0;
      V[k] = Rf_rbeta(aa, alphak);
    }
    // From V to w
    DDP_Vtow(w, V, N);

    // Sample beta;
    anovaDDP_sample_beta(beta, y, X, tau2, nK, Kind, mu, Sig, invSig, N, p);
    Xbeta = X.t()*beta;

    //Sample simga2;
    anovaDDP_sample_sigma2(tau2, y, Xbeta, nK, Kind, nua, nub, N);
    tau = Rcpp::sqrt(tau2);

    //Sample alpha;
    double a0star = a0+N-1;
    double b0star = b0-log(w[N-1]);
    if(b0star>(b0+740.0)){
      // Rprintf( "b0star = %f\n", b0star );
      b0star = b0+(N-1.0)/alpha;
    }
    alpha = Rf_rgamma(a0star, 1.0/b0star);
  
    //Sample mu;
    arma::mat S0star = inv_sympd( invS0 + (double)N*invSig );
    arma::vec m0star = S0star*( invS0m0 + invSig*arma::sum(beta, 1) );
    mu = mvrnorm(m0star, S0star);
  
    //Sample Sig^{-1}
    arma::mat mu_beta(p,p); mu_beta.fill(0.0);
    for (int k=0; k<N; ++k){
      mu_beta += (mu-beta.col(k))*(mu-beta.col(k)).t();
    }
    invSig = rwish( inv_sympd((double)k0*Sig0 + mu_beta), k0+N );
    Sig = inv_sympd(invSig);
  }
}

// Get density or survival Plots for DDP
SEXP DDPplots(SEXP xpred_, SEXP tgrid_, SEXP beta_, SEXP sigma2_, SEXP w_, SEXP CI_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  arma::mat xpred = as<mat>(xpred_); 
  arma::vec tgrid = as<vec>(tgrid_);
  arma::mat sigma2 = as<mat>(sigma2_);
  arma::mat w = as<mat>(w_);
  NumericVector vecArray(beta_);
  IntegerVector arrayDims = vecArray.attr("dim");
  arma::cube beta(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  double CI = as<double>(CI_);
  int nsave = w.n_cols;
  int ngrid = tgrid.size();
  int npred = xpred.n_rows;
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
  
  for(int i=0; i<nsave; ++i){
    arma::vec wei = w.col(i);
    arma::mat xbeta = arma::trans( xpred*beta.slice(i) );
    arma::vec sig = arma::sqrt(sigma2.col(i));
    for(int j=0; j<npred; ++j){
      (estf.slice(j)).col(i) = fmix(ygrid, xbeta.col(j), sig, wei)/tgrid;
      (estS.slice(j)).col(i) = 1.0 - Fmix(ygrid, xbeta.col(j), sig, wei);
      (esth.slice(j)).col(i) = (estf.slice(j)).col(i) / (estS.slice(j)).col(i);
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

// Spatial Maps for transformed spatial process using LDDPM-spatial model
SEXP PredMapsZ(SEXP ds0n_, SEXP dnn_, SEXP beta_, SEXP sigma2_, SEXP w_, SEXP theta1_, SEXP theta2_, SEXP z_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  arma::mat ds0n = as<mat>(ds0n_); // n by npred
  arma::mat dnn = as<mat>(dnn_);
  NumericVector vecArray(beta_);
  IntegerVector arrayDims = vecArray.attr("dim");
  arma::cube beta(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  arma::mat sigma2 = as<mat>(sigma2_);
  arma::mat w = as<mat>(w_);
  arma::vec theta1 = as<vec>(theta1_);
  arma::vec theta2 = as<vec>(theta2_);
  arma::mat z = as<mat>(z_);
  int nsave = w.n_cols;
  int n = dnn.n_rows;
  int npred = ds0n.n_cols;
  
  // temp variables;
  arma::mat Cinv=arma::eye(n,n);;
  double logdetC=0;
  // things to save;
  NumericMatrix Zpred(npred, nsave);
  
  RNGScope scope;
  arma::mat clustindx=arma::eye(n,n);
  for(int i=0; i<nsave; ++i){
    R_CheckUserInterrupt();
    arma::mat Cnn = arma::exp(-theta2(i)*dnn); 
    arma::mat Cnm = arma::exp(-theta2(i)*dnn);
    arma::mat Cmm = arma::exp(-theta2(i)*dnn);
    inv_FSA(theta1[i], Cnn, Cnm, Cmm, clustindx, Cinv, logdetC);
    arma::vec wei = w.col(i);
    arma::vec sig = arma::sqrt(sigma2.col(i));
    for(int j=0; j<npred; ++j){
      arma::vec hs0 = theta1(i)*arma::exp(-theta2(i)*ds0n.col(j));
      double mus0 = arma::dot( hs0, Cinv*z.col(i) );
      double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
      Zpred(j, i) = Rf_rnorm(mus0, sigs0);
    }
  }
  return List::create(Named("Zpred")=Zpred);
  END_RCPP
}
