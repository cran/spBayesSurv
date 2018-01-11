#include "spSurv_spatialtools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma ;
using namespace Rcpp ;
using namespace std;

// inverse of (ZK^{-1}Z' + In/lambda2) and its determinant using low-rank kriging
// note: for log determinant, we remove a constant (m-n)log(lambda2);
void Inverse_kriging(const arma::mat& Z, const arma::mat& K, double lambda2, const arma::mat& In,
                     arma::mat& Cinv, double& logdetC){
  arma::mat tempM = K/lambda2+Z.t()*Z;
  Cinv = lambda2*(In - Z*arma::solve(tempM, Z.t()));
  double val0, sign0;
  arma::log_det(val0, sign0, tempM);
  logdetC = val0; 
  arma::log_det(val0, sign0, K);
  logdetC -= val0;
}

// Get distance matrix for s1 s2 ... sn and s1 s2 ... sm, where s1 is column vector with coordinates
Rcpp::NumericMatrix Dist(Rcpp::NumericMatrix si, Rcpp::NumericMatrix sj){
  int ni = si.ncol();
  int nj = sj.ncol();
  NumericMatrix res(ni, nj);
  for(int i=0; i<ni; ++i){
    for(int j=0; j<nj; ++j){
      res(i,j) = std::sqrt(Rcpp::sum(Rcpp::pow(si(_,i)-sj(_,j),2)));
    }
  }
  return(res);
}

// Get distance matrix for s1 s2 ... sn and s1 s2 ... sm, where s1 is column vector with coordinates
SEXP DistMat(SEXP si_, SEXP sj_){
  Rcpp::NumericMatrix si(si_);
  Rcpp::NumericMatrix sj(sj_);
  int ni = si.ncol();
  int nj = sj.ncol();
  NumericMatrix res(ni, nj);
  for(int i=0; i<ni; ++i){
    for(int j=0; j<nj; ++j){
      res(i,j) = std::sqrt(Rcpp::sum(Rcpp::pow(si(_,i)-sj(_,j),2)));
    }
  }
  return(res);
}

//////////////////////////////////////////////////////////////////////
// spatial density things
/////////////////////////////////////////////////////////////////////
// M-distance function
double Mdist(arma::vec x1, arma::vec x2, const arma::mat& Sinv, double phi){
  return( std::exp( -phi*std::sqrt(arma::dot(x1-x2, Sinv*(x1-x2))) ) );
}

// log density for new x0 given the data (X,y);
void logf_spatdens(double y0, arma::vec x0, const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                   double th1, double exp_th2, double phi, const arma::mat& Sinv, 
                   Rcpp::IntegerMatrix& kyindex, double& logf){
  int n = y.size();
  Rcpp::IntegerVector kyindex0(J);
  int kJ = (int)(std::pow(2, J)*Rf_pnorm5(y0, th1, exp_th2, true, false));
  for(int j=0; j<J; ++j){
    kyindex0[J-j-1] += kJ; 
    kJ = (int)(kJ/2.0);
  }
  NumericVector distXx0(n);
  for(int i=0; i<n; ++i){
    distXx0[i] = Mdist(X.col(i), x0, Sinv, phi);
  }
  NumericVector ycount(J, 0.0);
  for(int j=0; j<J; ++j){
    for(int i=0; i<n; ++i){
      if(kyindex(i,j)==kyindex0[j]) ycount[j] += distXx0[i];
    }
  }
  logf = Rf_dnorm4(y0, th1, exp_th2, true);
  for(int j=1; j<J; ++j){
    logf += log(cpar*pow(j+1,2) + ycount[j]) - log(cpar*pow(j+1,2) + 0.5*ycount[j-1]);
  }
  logf += log(cpar+ycount[0]) - log(cpar+0.5*Rcpp::sum(distXx0));
}

// log likelihood given the data (X,y);
void loglik_spatdens(const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                     double phi, const arma::mat& Sinv, Rcpp::NumericVector& lognormy, 
                     Rcpp::IntegerMatrix& kyindex, double& logf){
  int n = y.size();
  logf = lognormy[0];
  for(int i=1; i<n; ++i){
    logf += lognormy[i];
    NumericVector distXx0(i);
    for(int i2=0; i2<i; ++i2){
      distXx0[i2] = Mdist(X.col(i2), X.col(i), Sinv, phi);
    }
    NumericVector ycount(J, 0.0);
    for(int j=0; j<J; ++j){
      for(int i2=0; i2<i; ++i2){
        if(kyindex(i2,j)==kyindex(i,j)) ycount[j] += distXx0[i2];
      }
    }
    for(int j=1; j<J; ++j){
      logf += log(cpar*pow(j+1,2) + ycount[j]) - log(cpar*pow(j+1,2) + 0.5*ycount[j-1]);
    }
    logf += log(cpar+ycount[0]) - log(cpar+0.5*Rcpp::sum(distXx0));
  }
}

// log likelihood given the data (X,y) with the parametric parts removed;
void loglik_spatdens_q(const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                       double phi, const arma::mat& Sinv, Rcpp::IntegerMatrix& kyindex, double& logf){
  int n = y.size();
  logf = 0;
  for(int i=1; i<n; ++i){
    NumericVector distXx0(i);
    for(int i2=0; i2<i; ++i2){
      distXx0[i2] = Mdist(X.col(i2), X.col(i), Sinv, phi);
    }
    NumericVector ycount(J, 0.0);
    for(int j=0; j<J; ++j){
      for(int i2=0; i2<i; ++i2){
        if(kyindex(i2,j)==kyindex(i,j)) ycount[j] += distXx0[i2];
      }
    }
    for(int j=1; j<J; ++j){
      logf += log(cpar*pow(j+1,2) + ycount[j]) - log(cpar*pow(j+1,2) + 0.5*ycount[j-1]);
    }
    logf += log(cpar+ycount[0]) - log(cpar+0.5*Rcpp::sum(distXx0));
  }
}

// log posterior of y_i with the parametric part removed; used for updating censored outcomes
void logq_yi_spatdens(double y0, arma::vec x0, int indx, const Rcpp::NumericVector& y, const arma::mat& X, 
                      int J, double cpar, double th1, double exp_th2, double phi, const arma::mat& Sinv, 
                      Rcpp::IntegerMatrix& kyindex, double& logf){
  int n = y.size();
  Rcpp::IntegerVector kyindex0(J);
  int kJ = (int)(std::pow(2, J)*Rf_pnorm5(y0, th1, exp_th2, true, false));
  for(int j=0; j<J; ++j){
    kyindex0[J-j-1] += kJ; 
    kJ = (int)(kJ/2.0);
  }
  NumericVector distXx0(n, 0.0);
  for(int i=0; (i<n)&&(i!=indx); ++i){
    distXx0[i] = Mdist(X.col(i), x0, Sinv, phi);
  }
  NumericVector ycount(J, 0.0);
  for(int j=0; j<J; ++j){
    for(int i=0; i<n; ++i){
      if((kyindex(i,j)==kyindex0[j])&&(i!=indx)) ycount[j] += distXx0[i];
    }
  }
  logf = 0;
  for(int j=1; j<J; ++j){
    logf += log(cpar*pow(j+1,2) + ycount[j]) - log(cpar*pow(j+1,2) + 0.5*ycount[j-1]);
  }
  logf += log(cpar+ycount[0]) - log(cpar+0.5*Rcpp::sum(distXx0));
}

// posterior prob of phi=0; used for updating phi
void prob_phi0_spatdens(const Rcpp::NumericVector& y, const arma::mat& X, int J, double cpar,
                        const Rcpp::NumericVector& phiseq, const arma::mat& Sinv,
                        Rcpp::IntegerMatrix& kyindex, double q0phi, double a0phi, double b0phi,
                        double& ratio){
  double integralpart = 0;
  double logprob_temp = 0;
  for(int k=1; k<phiseq.size(); ++k){
    loglik_spatdens_q(y, X, J, cpar, phiseq[k], Sinv, kyindex, logprob_temp);
    integralpart += exp(logprob_temp+Rf_dgamma(phiseq[k], a0phi, 1.0/b0phi, true)-log(phiseq[k]-phiseq[k-1]));
  }
  loglik_spatdens_q(y, X, J, cpar, 0.0, Sinv, kyindex, logprob_temp);
  double prob0 = exp( logprob_temp );
  ratio = q0phi*prob0/(q0phi*prob0+(1.0-q0phi)*integralpart);
}

//////////////////////////////////////////////////////////////////////
// spatial Copula things
/////////////////////////////////////////////////////////////////////
// make transformation on theta: log(theta1/(1.0-theta1)); std::log(theta2);
arma::vec trans_theta(arma::vec theta){
  arma::vec trans(2);
  trans[0] = std::log(theta[0]/(1.0-theta[0]));
  trans[1] = std::log(theta[1]);
  return(trans);
}

// transform back to get theta;
arma::vec trans_theta_inv(arma::vec trans){
  arma::vec theta(2);
  theta[0] = 1.0/(1.0+std::exp(-trans[0]));
  theta[1] = std::exp(trans[1]);
  return(theta);
}

//Sample theta using adaptive M-H for spatial Copula Model using FSA;
void spCopula_sample_theta_FSA(arma::vec& theta, double& rejtheta, arma::mat& spSnew, arma::vec& thetabarnew, 
                               arma::mat& Cinv, double& logdetC, double theta1a, double theta1b, 
                               double theta2a, double theta2b, double spl0, arma::mat spShat, const arma::mat& dnn, 
                               double spadapter, int iscan, int nburn, const arma::vec& z, int n, const arma::mat& dnm, 
                               const arma::mat& dmm, const arma::mat clustindx){
  if(theta1a<0){
    double Snew = spSnew(1,1);
    double barnew = thetabarnew[1];
    double theta2 = theta[1];
    double trans2 = std::log(theta2);
    double tempold = -0.5*logdetC - 0.5*arma::dot(z, Cinv*z) + theta2a*std::log(theta2) - theta2b*theta2;
    double thetaold = theta2;
    double transold = trans2;
    arma::mat Cinvold = Cinv;
    double logdetCold = logdetC;
    if(iscan>spl0){
        trans2 = Rf_rnorm(transold, Snew);
    }else{
        trans2 = Rf_rnorm(transold, spShat(1,1));
    }
    theta2 = std::exp(trans2);
    arma::mat Cnn = arma::exp(-theta2*dnn); 
    arma::mat Cnm = arma::exp(-theta2*dnm);
    arma::mat Cmm = arma::exp(-theta2*dmm);
    inv_FSA(theta[0], Cnn, Cnm, Cmm, clustindx, Cinv, logdetC);
    double tempnew = -0.5*logdetC - 0.5*arma::dot(z, Cinv*z) + theta2a*std::log(theta2) - theta2b*theta2;
    double ratio = std::exp( tempnew-tempold );
    // Rprintf( "%f\n", ratio );
    double uu = unif_rand();
    if (uu>ratio) {
      theta2 = thetaold; Cinv=Cinvold; logdetC=logdetCold; trans2=transold;
      if(iscan>=nburn) rejtheta+=1.0;
    }
    double nn = iscan+1;
    double barold = barnew;
    barnew = (nn)/(nn+1.0)*barold + trans2/(nn+1.0);
    spSnew(1,1) = (nn-1.0)/nn*Snew+spadapter/nn*(nn*barold*barold-(nn+1.0)*barnew*barnew+trans2*trans2 + ESMALL);
    thetabarnew[1] = barnew;
    theta[1] = theta2;
  }else{
    arma::vec trans = trans_theta(theta);
    double tempold = -0.5*logdetC - 0.5*arma::dot(z, Cinv*z) + (theta1a+1.0)*std::log(theta[0]) + (theta1b-1.0)*std::log(1.0-theta[0]) 
                     - trans[0] + theta2a*std::log(theta[1]) - theta2b*theta[1];
    arma::vec thetaold = theta;
    arma::vec transold = trans;
    arma::mat Cinvold = Cinv;
    double logdetCold = logdetC;
    if(iscan>spl0){
        trans = mvrnorm(transold, spSnew);
    }else{
        trans = mvrnorm(transold, spShat);
    }
    theta = trans_theta_inv(trans);
    arma::mat Cnn = arma::exp(-theta[1]*dnn); 
    arma::mat Cnm = arma::exp(-theta[1]*dnm);
    arma::mat Cmm = arma::exp(-theta[1]*dmm);
    inv_FSA(theta[0], Cnn, Cnm, Cmm, clustindx, Cinv, logdetC);
    double tempnew = -0.5*logdetC - 0.5*arma::dot(z, Cinv*z) + (theta1a+1.0)*std::log(theta[0]) + (theta1b-1.0)*std::log(1.0-theta[0]) 
                     - trans[0] + theta2a*std::log(theta[1]) - theta2b*theta[1];
    double ratio = std::exp( tempnew-tempold );
    // Rprintf( "%f\n", ratio );
    double uu = unif_rand();
    if (uu>ratio) {
      theta = thetaold; Cinv=Cinvold; logdetC=logdetCold; trans=transold;
      if(iscan>=nburn) rejtheta+=1.0;
    }
    double nn = iscan+1;
    arma::vec thetabarold = thetabarnew;
    thetabarnew = (nn)/(nn+1.0)*thetabarold + trans/(nn+1.0);
    spSnew = (nn-1.0)/nn*spSnew + spadapter*0.5/nn*(nn*thetabarold*thetabarold.t() - (nn+1.0)*thetabarnew*thetabarnew.t() 
            + trans*trans.t() + ESMALL*arma::eye(2,2) );
  }
}

