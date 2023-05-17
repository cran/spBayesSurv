#include "spSurv_spatialtools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

RcppExport SEXP SpatDens(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
                         SEXP y_, SEXP y1_, SEXP y2_, SEXP type_, SEXP X_, SEXP theta_, SEXP maxJ_,
                         SEXP cpar_, SEXP a0_, SEXP b0_, SEXP theta0_, 
                         SEXP V0inv_, SEXP Vhat_, SEXP l0_, SEXP adapter_, 
                         SEXP Sinv_, SEXP phi_, SEXP q0phi_, SEXP a0phi_, SEXP b0phi_, SEXP perm_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  const int ndisplay = as<int>(ndisplay_);
  const Rcpp::NumericVector y1(y1_); // n by 1;
  const Rcpp::NumericVector y2(y2_); // n by 1;
  const Rcpp::IntegerVector type(type_);// n by 1;
  const int maxJ = as<int>(maxJ_);
  const int n = y1.size();
  const int l0 = as<int>(l0_);
  const double adapter = as<double>(adapter_);
  const arma::mat Sinv = as<mat>(Sinv_); // p by p;
  arma::mat X = as<mat>(X_); // p by n;
  const int perm=as<int>(perm_);
  
  // prameters to be updated
  Rcpp::NumericVector theta(theta_); // 2 by 1
  double cpar = as<double>(cpar_);
  double phi = as<double>(phi_);
  Rcpp::NumericVector y(y_);
  
  // hyperparameters
  // prior of cpar
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  // prior of theta
  const arma::vec theta0 = as<arma::vec>(theta0_); // 2 by 1
  const arma::mat V0inv = as<arma::mat>(V0inv_);   // 2 by 2
  const arma::mat Vhat = as<arma::mat>(Vhat_);     // 2 by 2
  // prior of phi
  const double q0phi = as<double>(q0phi_);
  const double a0phi = as<double>(a0phi_);
  const double b0phi = as<double>(b0phi_);
  
  // temp variables
  const int nscan = nburn + (nskip+1)*nsave;
  int skiptally=0; 
  int isave=0;
  int distally=0;
  
  // things to save;
  Rcpp::NumericVector cpar_save(nsave);
  Rcpp::NumericMatrix theta_save(2, nsave);
  Rcpp::NumericMatrix y_save(n, nsave);
  Rcpp::NumericVector phi_save(nsave);
  double rejtheta=0;
  double rejc=0;
  double rejphi=0;
  Rcpp::NumericVector rejy(n, 0.0);
  
  // Make arma objects
  arma::vec theta_r(theta.begin(), 2, false, true);
  
  // Working temp variables
  double llold, llnew;
  double ratio, uu, nn;
  int phi_zero_count=0;
  Rcpp::IntegerMatrix kyindex(n, maxJ);
  arma::imat kyindex_r(kyindex.begin(), n, maxJ, false, true);
  arma::imat kyindexold(n, maxJ);
  Rcpp::NumericVector lognormy(n);
  Rcpp::NumericVector lognormyold(n);
  double maxJ2 = std::pow(2, maxJ);
  int kJ=0;
  int K=20;
  Rcpp::NumericVector phiseq(K+1, 0.0);
  kyindex_r.fill(0);
  for(int i=0; i<n; ++i){
    lognormy[i] = Rf_dnorm4(y[i], theta[0], exp(theta[1]), true);
    kJ = (int)(maxJ2*Rf_pnorm5(y[i], theta[0], exp(theta[1]), true, false));
    for(int j=0; j<maxJ; ++j){
      kyindex(i, maxJ-j-1) += kJ; 
      kJ = (int)(kJ/2.0);
    }
  }
  for(int k=1; k<=K; ++k){
    //phiseq[k] = -log(1.0-k/(K+1.0))/b0phi;
    phiseq[k] = Rf_qgamma((k+0.0)/(K+1.0), a0phi, 1.0/b0phi, true, false);
  }
  
  // for theta
  arma::vec thetaold(2); 
  arma::vec thBarold(2);
  arma::vec thBarnew=arma::zeros<arma::vec>(2);
  arma::mat thSnew=arma::zeros<arma::mat>(2,2);
  arma::mat I2 = ESMALL*arma::eye(2,2);
  // for cpa
  double cbarnew=0, cbarold=0, csnew=0, ctemp=0, cshat=0.16;
  // for phi
  double phibarnew=0, phibarold=0, phisnew=0, phinew=0, phishat=0.01;
  double phi_non0=0;
  phi_non0 = phi;
  
  // temp variables for permutations
  arma::vec y_r=as<arma::vec>(y);
  arma::uvec indx(n);
  for(int i=0; i<n; ++i) indx[i]=i;
  arma::uvec indx_u(n);
  
  RNGScope scope;
  
  // Set the Armadillo seed from R's 
  // int seed = (int)Rf_runif(0.0, 10000.0);
  // std::srand(seed);
  
  ////////////////////////////////////////////////////////////////////////
  // Start MCMC
  ////////////////////////////////////////////////////////////////////////
  for (int iscan=0; iscan<nscan; iscan++){
    R_CheckUserInterrupt();
    //Rprintf( "iscan = %d\n", iscan );
    //Rprintf( "lambda = %f\n", lambda );
    //Rprintf( "phi = %f\n", phi );
    
    ///////////////////////////////////////////////
    // update censored outcomes
    //////////////////////////////////////////////
    for (int i=0; i<n; ++i){
      if(type[i]!=1){
        double yiold = y[i];
        logq_yi_spatdens(y[i], X.col(i), i, y, X, maxJ, cpar, theta[0], exp(theta[1]), phi, Sinv, kyindex, llold);
        y[i] = trun_rnorm(theta[0], exp(theta[1]), y1[i], y2[i]);
        logq_yi_spatdens(y[i], X.col(i), i, y, X, maxJ, cpar, theta[0], exp(theta[1]), phi, Sinv, kyindex, llnew);
        ratio = exp(llnew-llold);
        uu = unif_rand();
        if(uu>ratio){
          y[i]=yiold; 
          if(iscan>=nburn) rejy[i]+=1.0;
        }
        y_r(i) = y[i];
      }
    }
    
    ///////////////////////////////////////////////
    // data permutation
    //////////////////////////////////////////////
    if(perm==1){
      indx_u=arma::shuffle(indx);
      y_r = y_r(indx_u);
      for(int i=0; i<n; ++i) y[i]=y_r(i);
      X = X.cols(indx_u);
      kyindex_r.fill(0);
      for(int i=0; i<n; ++i){
        lognormy[i] = Rf_dnorm4(y[i], theta[0], exp(theta[1]), true);
        kJ = (int)(maxJ2*Rf_pnorm5(y[i], theta[0], exp(theta[1]), true, false));
        for(int j=0; j<maxJ; ++j){
          kyindex(i, maxJ-j-1) += kJ; 
          kJ = (int)(kJ/2.0);
        }
      }
    }
    
    ///////////////////////////////////////////////
    // update theta
    //////////////////////////////////////////////
    if(V0inv(0,0)<SYSMAX){
      loglik_spatdens(y, X, maxJ, cpar, phi, Sinv, lognormy, kyindex, llold);
      llold += -0.5*arma::dot( (theta_r-theta0), (V0inv*(theta_r-theta0)) );
      thetaold = theta_r;
      kyindexold = kyindex_r;
      for(int i=0; i<n; ++i) lognormyold[i]=lognormy[i];
      if(iscan>l0){
        theta_r = mvrnorm(thetaold, thSnew);
      }else{
        theta_r = mvrnorm(thetaold, Vhat);
      }
      kyindex_r.fill(0);
      for(int i=0; i<n; ++i){
        lognormy[i] = Rf_dnorm4(y[i], theta[0], exp(theta[1]), true);
        kJ = (int)(maxJ2*Rf_pnorm5(y[i], theta[0], exp(theta[1]), true, false));
        for(int j=0; j<maxJ; ++j){
          kyindex(i, maxJ-j-1) += kJ; 
          kJ = (int)(kJ/2.0);
        }
      }
      loglik_spatdens(y, X, maxJ, cpar, phi, Sinv, lognormy, kyindex, llnew);
      llnew += -0.5*arma::dot( (theta_r-theta0), (V0inv*(theta_r-theta0)) );
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if(uu>ratio){
        theta_r=thetaold; kyindex_r=kyindexold;
        for(int i=0; i<n; ++i) lognormy[i]=lognormyold[i];
        if(iscan>=nburn) rejtheta+=1.0;
      }
      nn = iscan+1;
      thBarold = thBarnew;
      thBarnew = (nn)/(nn+1.0)*thBarold + theta_r/(nn+1.0);
      thSnew = (nn-1.0)/nn*thSnew + adapter/(2.0)/nn*(nn*thBarold*thBarold.t() 
                                                        - (nn+1.0)*thBarnew*thBarnew.t() + theta_r*theta_r.t() + I2 );
    }
    
    ///////////////////////////////////////////////
    // update phi
    //////////////////////////////////////////////
    if(q0phi<0.9999999999){
      if(q0phi<ESMALL){
        ratio=0.0;
      }else{
        //Rprintf( "phi = %f\n", phi );
        prob_phi0_spatdens(y, X, maxJ, cpar, phiseq, Sinv, kyindex, q0phi, a0phi, b0phi, ratio);
      }
      uu = unif_rand();
      if(uu<ratio){
        phi = 0.0; phi_zero_count +=1;
      }else{
        if(a0phi>0){
          if(iscan>l0){
            phinew = Rf_rnorm(phi_non0, std::sqrt(phisnew));
          }else{
            phinew = Rf_rnorm(phi_non0, std::sqrt(phishat));
          }
          if(phinew<SYSMIN){
            if(iscan>=nburn) rejphi+=1.0;
          }else{
            loglik_spatdens_q(y, X, maxJ, cpar, phi_non0, Sinv, kyindex, llold);
            llold += (a0phi-1.0)*log(phi_non0) - b0phi*phi_non0;
            loglik_spatdens_q(y, X, maxJ, cpar, phinew, Sinv, kyindex, llnew);
            llnew += (a0phi-1.0)*log(phinew) - b0phi*phinew;
            ratio = exp(llnew-llold);
            uu = unif_rand();
            if(uu<ratio){
              phi=phinew; phi_non0=phinew;
            }else{
              if(iscan>=nburn) rejphi+=1.0;
            }
          }
          nn = iscan-phi_zero_count+1;
          phibarold = phibarnew;
          phibarnew = (nn)/(nn+1.0)*phibarold + phi_non0/(nn+1.0);
          phisnew = (nn-1.0)/nn*phisnew + adapter/nn*(nn*pow(phibarold,2) - (nn+1.0)*pow(phibarnew,2)
                                                        + pow(phi_non0,2) + ESMALL );
        }
      }
    }else{
      phi=0.0;
    }
    
    ///////////////////////////////////////////////
    // cpar
    //////////////////////////////////////////////
    if(a0>0){
      if(iscan>l0){
        ctemp = Rf_rnorm(cpar, std::sqrt(csnew));
      }else{
        ctemp = Rf_rnorm(cpar, std::sqrt(cshat));
      }
      if(ctemp<ESMALL){
        if(iscan>=nburn) rejc += 1;
      }else{
        loglik_spatdens_q(y, X, maxJ, cpar, phi, Sinv, kyindex, llold);
        llold  += (a0-1.0)*log(cpar)-b0*cpar; 
        loglik_spatdens_q(y, X, maxJ, ctemp, phi, Sinv, kyindex, llnew);
        llnew += (a0-1.0)*log(ctemp)-b0*ctemp;
        ratio = exp(llnew-llold);
        uu = unif_rand();
        if(uu<ratio){
          cpar = ctemp; 
        }else{
          if(iscan>=nburn) rejc += 1;
        }
      }
      nn = iscan+1;
      cbarold = cbarnew;
      cbarnew = (nn)/(nn+1.0)*cbarold + cpar/(nn+1.0);
      csnew = (nn-1.0)/nn*csnew + adapter/nn*(nn*pow(cbarold,2) - (nn+1.0)*pow(cbarnew,2) + pow(cpar,2) + ESMALL );
    }
    
    // permutate back to original order;
    if(perm==1){
      indx_u = arma::sort_index(indx_u);
      y_r = y_r(indx_u);
      for(int i=0; i<n; ++i) y[i]=y_r(i);
      X = X.cols(indx_u);
      kyindex_r.fill(0);
      for(int i=0; i<n; ++i){
        //lognormy[i] = Rf_dnorm4(y[i], theta[0], exp(theta[1]), true);
        kJ = (int)(maxJ2*Rf_pnorm5(y[i], theta[0], exp(theta[1]), true, false));
        for(int j=0; j<maxJ; ++j){
          kyindex(i, maxJ-j-1) += kJ; 
          kJ = (int)(kJ/2.0);
        }
      }
    }
    
    ///////////////////////////////////////////////
    // Save the sample
    //////////////////////////////////////////////
    if(iscan>=nburn){
      ++skiptally;
      if(skiptally>nskip){
        // save data
        y_save(_,isave) = y;
        theta_save(_,isave) = theta;
        cpar_save[isave] = cpar;
        phi_save[isave] = phi;
        
        ++isave;
        ++distally;
        if(distally>=ndisplay){
          Rprintf( "scan = %d\n", isave );
          distally = 0;
        }
        skiptally=0;
      }
    }
  }
  
  // get acceptance rate
  double totscans = nscan-nburn+0.0;
  double ratetheta = 1.0 - rejtheta/totscans;
  double ratec = 1.0 - rejc/totscans;
  double ratephi = 1.0 - rejphi/(totscans-phi_zero_count);
  Rcpp::NumericVector ratey = 1.0 - rejy/totscans;
  
  return List::create(Named("y")=y_save,
                      Named("theta")=theta_save,
                      Named("cpar")=cpar_save,
                      Named("phi")=phi_save,
                      Named("ratetheta")=ratetheta,
                      Named("ratec")=ratec,
                      Named("ratey")=ratey,
                      Named("ratephi")=ratephi);
  END_RCPP
}

// Get density Plots 
RcppExport SEXP SpatDens_plots(SEXP ygrid_, SEXP xpred_, SEXP theta_, SEXP cpar_, SEXP phi_,
                               SEXP maxJ_, SEXP y_, SEXP X_, SEXP Sinv_, SEXP CI_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  Rcpp::NumericVector ygrid(ygrid_);
  arma::mat xpred = as<arma::mat>(xpred_); // p by npred;
  Rcpp::NumericMatrix theta(theta_); // 2 by nsave;
  Rcpp::NumericVector cpar(cpar_); // nsave by 1;
  Rcpp::NumericVector phi(phi_); // nsave by 1;
  const arma::mat X = as<mat>(X_); // p by n;
  const Rcpp::NumericMatrix y(y_); // n by nsave;
  const int maxJ = as<int>(maxJ_);
  const arma::mat Sinv = as<mat>(Sinv_); // p by p;
  double CI = as<double>(CI_);
  int nsave = cpar.size();
  int ngrid = ygrid.size();
  int npred = xpred.n_cols;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  int n = X.n_cols;
  
  // Temp variables
  Rcpp::NumericVector estfArray(nsave*ngrid*npred);
  arma::cube estf(estfArray.begin(), ngrid, nsave, npred, false, true);
  Rcpp::IntegerMatrix kyindex(n, maxJ);
  arma::imat kyindex_r(kyindex.begin(), n, maxJ, false, true);
  double maxJ2 = std::pow(2, maxJ);
  double logf=0;
  
  // things to save;
  arma::mat fhat(ngrid, npred);
  arma::mat fhatup(ngrid, npred);
  arma::mat fhatlow(ngrid, npred);
  
  for(int i=0; i<nsave; ++i){
    Rcpp::NumericVector yi = y(_,i);
    double th1 = theta(0,i);
    double exp_th2 = exp(theta(1,i));
    kyindex_r.fill(0);
    int kJ=0;
    for(int ii=0; ii<n; ++ii){
      kJ = (int)(maxJ2*Rf_pnorm5(yi[ii], th1, exp_th2, true, false));
      for(int j=0; j<maxJ; ++j){
        kyindex(ii, maxJ-j-1) += kJ; 
        kJ = (int)(kJ/2.0);
      }
    }
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        logf_spatdens(ygrid[k], xpred.col(j), yi, X, maxJ, cpar[i],
                      th1, exp_th2, phi[i], Sinv, kyindex, logf);
        estf(k, i, j) = std::exp(logf);
      }
    }
  }
  for(int j=0; j<npred; ++j){
    fhat.col(j) = arma::mean(estf.slice(j), 1);
    arma::mat temp = arma::sort(estf.slice(j),"ascend", 1);
    fhatlow.col(j) = temp.col(low);
    fhatup.col(j) = temp.col(up);
  }
  return List::create(Named("fhat")=fhat,
                      Named("fhatlow")=fhatlow,
                      Named("fhatup")=fhatup);
  END_RCPP
}

// Get empirical BF and permutation p-value for uncensored data
RcppExport SEXP SpatDens_BF(SEXP y_, SEXP X_, SEXP Sinv_, SEXP theta_, SEXP maxJ_, SEXP cpar_, 
                            SEXP a0_, SEXP b0_, SEXP phi_, SEXP q0phi_, SEXP a0phi_, SEXP b0phi_,
                            SEXP nperm_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  const Rcpp::NumericVector y(y_); // n by 1;
  const int maxJ = as<int>(maxJ_);
  const int n = y.size();
  const arma::mat X = as<mat>(X_); // p by n;
  const arma::mat Sinv = as<mat>(Sinv_); // p by p;
  const Rcpp::NumericVector theta(theta_); // 2 by 1
  const double mu = theta[0];
  const double sig = exp(theta[1]);
  const Rcpp::NumericVector cpar(cpar_); // m1 by 1;
  const Rcpp::NumericVector phi(phi_); // m2 by 1;
  const int m1 = cpar.size();
  const int m2 = phi.size();
  const int nperm = as<int>(nperm_);
  
  // hyperparameters
  // prior of cpar
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  // prior of phi
  const double q0phi = as<double>(q0phi_);
  const double a0phi = as<double>(a0phi_);
  const double b0phi = as<double>(b0phi_);
  
  // Working temp variables
  Rcpp::IntegerMatrix kyindex(n, maxJ);
  arma::imat kyindex_r(kyindex.begin(), n, maxJ, false, true);
  //Rcpp::NumericVector lognormy(n);
  double maxJ2 = std::pow(2, maxJ);
  int kJ=0;
  kyindex_r.fill(0);
  for(int i=0; i<n; ++i){
    //lognormy[i] = Rf_dnorm4(y[i], mu, sig, true);
    kJ = (int)(maxJ2*Rf_pnorm5(y[i], mu, sig, true, false));
    for(int j=0; j<maxJ; ++j){
      kyindex(i, maxJ-j-1) += kJ; 
      kJ = (int)(kJ/2.0);
    }
  }
  
  // temp variables for permutations
  arma::vec y_r(n);
  arma::vec yperm_r(n);
  Rcpp::NumericVector yperm(n);
  for(int i=0; i<n; ++i){
    y_r(i)=y[i];
    yperm_r(i) = y[i];
    yperm[i] = y[i];
  }
  
  arma::uvec indx(n);
  for(int i=0; i<n; ++i) indx[i]=i;
  arma::uvec indx_u(n);
  
  RNGScope scope;
  
  // Calculate BF
  arma::mat res_num(m1, m2);
  for(int i=0; i<m1; ++i){
    for(int j=0; j<m2; ++j){
      double lltemp = 0;
      loglik_spatdens_q(y, X, maxJ, cpar[i], phi[j], Sinv, kyindex, lltemp);
      lltemp += (a0-1.0)*log(cpar[i])-b0*cpar[i];
      lltemp += (a0phi-1.0)*log(phi[i]) - b0phi*phi[i];
      res_num(i,j) = lltemp;
    }
  }
  arma::vec res_den(m1);
  for(int i=0; i<m1; ++i){
    double lltemp = 0; 
    loglik_spatdens_q(y, X, maxJ, cpar[i], 0.0, Sinv, kyindex, lltemp);
    lltemp += (a0-1.0)*log(cpar[i])-b0*cpar[i];
    res_den(i) = lltemp;
  }
  double BF = exp(res_num.max()-res_den.max())*q0phi/(1.0-q0phi);
  //Rprintf( "BF = %f\n", BF );
  
  // Calculate permutation p-value
  arma::vec BF_perm(nperm);
  for (int iscan=0; iscan<nperm; iscan++){
    R_CheckUserInterrupt();
    //Rprintf( "iscan = %d\n", iscan );
    indx_u=arma::shuffle(indx);
    yperm_r = y_r(indx_u);
    for(int i=0; i<n; ++i) yperm[i]=yperm_r(i);
    kyindex_r.fill(0);
    for(int i=0; i<n; ++i){
      //lognormy[i] = Rf_dnorm4(yperm[i], mu, sig, true);
      kJ = (int)(maxJ2*Rf_pnorm5(yperm[i], mu, sig, true, false));
      for(int j=0; j<maxJ; ++j){
        kyindex(i, maxJ-j-1) += kJ; 
        kJ = (int)(kJ/2.0);
      }
    }
    for(int i=0; i<m1; ++i){
      for(int j=0; j<m2; ++j){
        double lltemp = 0;
        loglik_spatdens_q(yperm, X, maxJ, cpar[i], phi[j], Sinv, kyindex, lltemp);
        lltemp += (a0-1.0)*log(cpar[i])-b0*cpar[i];
        lltemp += (a0phi-1.0)*log(phi[i]) - b0phi*phi[i];
        res_num(i,j) = lltemp;
      }
    }
    for(int i=0; i<m1; ++i){
      double lltemp = 0; 
      loglik_spatdens_q(yperm, X, maxJ, cpar[i], 0.0, Sinv, kyindex, lltemp);
      lltemp += (a0-1.0)*log(cpar[i])-b0*cpar[i];
      res_den(i) = lltemp;
    }
    BF_perm(iscan) = exp(res_num.max()-res_den.max())*q0phi/(1.0-q0phi);
  }
  
  return List::create(Named("BF")=BF,
                      Named("BFperm")=BF_perm);
  END_RCPP
}
