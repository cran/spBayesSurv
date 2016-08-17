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

// Matern correlation function
double rho_Matern(double dis, double nu, double phi){
  if(nu==0.5){
    return (std::exp(-std::abs(dis/phi)));
  }else if (nu==1.5){
    return ( std::exp(-std::abs(dis/phi))*(1+std::abs(dis/phi)) );
  }else{
    return ( std::exp(-0.5*std::pow(dis/phi, 2)) );
  }
}

// process convolution bivariate Gaussian kernel
double kernel_G(double dis, double phi){
  return ( std::exp(-0.5*std::pow(dis/phi, 2)) );
}

// Powered Exponential
double pow_exp(double dis, double phi, double nu){
  return ( std::exp(-0.5*std::pow(std::abs(dis/phi), nu)) );
}

