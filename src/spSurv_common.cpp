#include "spSurv_common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma ;
using namespace Rcpp ;
using namespace std;

//Truncated normal N(y;mu,s)I(a<y<b); where b can be R_PosInf
// Rejection algorithm with a truncated expoential proposal for N(0,1)I(a<x<b) when a is very large: |a| < |b|
double rtexp(double a, double b){
  int stop = false;
  double twoasp = 2*std::pow(a,2);
  double expab = std::exp(-a*(b-a)) - 1;
  double z, e;
  while(!stop){
    R_CheckUserInterrupt();
    z = std::log(1 + unif_rand()*expab);
    e = -std::log(unif_rand());
    stop = (twoasp*e > std::pow(z,2));
  }
  return (a - z/a);
}
double trun_rnorm(const double mu, const double s, double a, double b){
  double xmin = -2.00443204036;                 // Left bound
  double xmax =  3.48672170399;                 // Right bound
  int stop = false;
  double r;
  //if( mu+ELARGE<0 ) return(a);
  //scalling
  if(mu!=0 || s!=1){
    a=(a-mu)/s;
    b=(b-mu)/s;
  }
  // Check if a < b
  if(a>=b){
    Rprintf( "*** B must be greater than A ! ***" ); return(NA_REAL);
  }
  else if(std::abs(a)>std::abs(b)) r = -trun_rnorm(0, 1, -b, -a);
  // If a in the right tail (a > xmax), use rejection algorithm with a truncated exponential proposal  
  else if(a>xmax) r = rtexp(a,b);
  // If a in the left tail (a < xmin), use rejection algorithm with a Gaussian proposal
  else if(a<xmin){
    while(!stop){
      R_CheckUserInterrupt();
      r = norm_rand();
      stop = (r>=a) && (r<=b);
    }
  }  
  // In other cases (xmin < a < xmax)
  else{
    double CDFa = Rf_pnorm5(a, 0, 1.0, true, false);
    double CDFb = Rf_pnorm5(b, 0, 1.0, true, false);
    double u = unif_rand();
    double CDFxi = CDFa + u*(CDFb - CDFa);
    r = Rf_qnorm5(CDFxi, 0, 1, true, false);
  }
  // Scaling
  if(mu!=0 || s!=1)
  r = r*s + mu;
  return r;
}

// generate multivariate normal (mu, sigma)
arma::vec mvrnorm(arma::vec mu, arma::mat sigma) {
  int ncols = mu.size();
  arma::vec Y(ncols);
  for (int i=0; i<ncols; i++){
    Y[i] = norm_rand();
  }
  arma::mat temp = ((arma::chol(sigma)).t())*Y;
  arma::vec res = mu + temp.col(0);
  return( res );
}

// density of multivariate normal (mu, sigma)
double mvdnorm(arma::vec x, arma::vec mu, arma::mat sigma, bool logd=true) { 
  int xdim = x.size();
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(arma::log(rooti.diag()));
  double constants = -(double)(xdim)*0.5*std::log(2.0*M_PI);
  arma::vec z = rooti*( x - mu) ;    
  double res = constants - 0.5*arma::sum(z%z) + rootisum;     
  if (logd == false) {
    res = std::exp(res);
  }
  return(res);
}

// generate Wishart random matrices
arma::mat rwish(arma::mat Sig, int n) {
  int ncols = Sig.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::mat X = Y * arma::chol(Sig);
  return( X.t()*X );
}

// sample(Nseq, prob=w), where Nseq is n-dim vector, and w.siz()=n
int sample(Rcpp::IntegerVector Nseq, Rcpp::NumericVector w){
  int k = 0;
  double u = unif_rand();;
  double cdf = w[0];
  while(u>cdf){
    cdf += w[k];
    ++k;
  }
  return (Nseq[k]);
}

// calculate qnorm(x) for a vector of x
arma::vec qnormvec(arma::vec x){
  int n = x.size();
  arma::vec res(n);
  for (int i=0; i<n; ++i){
    res[i] = std::min( Rf_qnorm5(x[i], 0, 1, true, false), 8.209536);
  }
  return (res);
}

// calculate pnorm(x) for a vector of x
arma::vec pnormvec(arma::vec x){
  int n = x.size();
  arma::vec res(n);
  for (int i=0; i<n; ++i){
    res[i] = Rf_pnorm5(x[i], 0, 1, true, false);
  }
  return (res);
}

// Gaussian correlation function
double rho_Gau(double dis, double phi){
  return ( std::exp(-std::pow(phi*dis, 2)) );
}

// Exponential correlation function
double rho_Exp(double dis, double phi){
  return ( std::exp(-std::abs(phi*dis)) );
}

// Powered Exponential
double pow_exp(double dis, double phi, double nu){
  return ( std::exp(-std::pow(std::abs(phi*dis), nu)) );
}

// Preprocess R^{-1} to get Rinv using FSA
void inv_FSA(double sill, const arma::mat& Cnn, const arma::mat& Cnm, const arma::mat& Cmm,
              const arma::mat& clustindx, arma::mat& Cinv, double& logdetC){
  int n=Cnn.n_cols;
  int nclust = clustindx.n_cols;
  arma::mat Cs = sill*(Cnn-Cnm*arma::solve(Cmm, Cnm.t())) + (1.0-sill)*arma::eye(n,n);
  arma::mat invCs = arma::diagmat(1.0/Cs.diag());
  double logdet3 = arma::sum(arma::log(Cs.diag()));
  if(nclust<n){
    logdet3=0;
    for (int i=0; i<nclust; i++){
      arma::uvec vindx = arma::find(clustindx.col(i));
      invCs.submat(vindx, vindx) = arma::inv_sympd(Cs.submat(vindx, vindx));
      double val0, sign0;
      arma::log_det(val0, sign0, Cs.submat(vindx, vindx));
      logdet3 += val0;
    }
  }
  arma::mat invCsCnm = invCs*Cnm;
  arma::mat More = Cmm + sill*Cnm.t()*invCsCnm;
  Cinv = invCs - sill*invCsCnm*arma::solve(More,invCsCnm.t());
  double val0, sign0;
  arma::log_det(val0, sign0, More);
  double logdet1 = val0;
  arma::log_det(val0, sign0, Cmm);
  double logdet2 = val0;
  logdetC = logdet1 - logdet2 + logdet3;
}
