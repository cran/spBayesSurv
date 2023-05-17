#include "spSurv_Coxph_tools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

///////////////////////////// Random cutpoints //////////////////////////////////
RcppExport SEXP spCopulaCoxph(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
                              SEXP tobs_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, 
                              SEXP r0_, SEXP hcen_, SEXP h0_, SEXP v0_, SEXP vhat_,
                              SEXP beta_, SEXP beta0_, SEXP S0inv_, SEXP Shat_, 
                              SEXP l0_, SEXP adapter_, SEXP xpred_, SEXP ds0n_, 
                              SEXP dnn_, SEXP theta_, SEXP theta0_, SEXP spS0_,
                              SEXP dnm_, SEXP dmm_, SEXP clustindx_, SEXP ds0m_, SEXP ds0block_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  const int ndisplay = as<int>(ndisplay_);
  const Rcpp::NumericVector tobs(tobs_);
  const Rcpp::IntegerVector delta(delta_);// n by 1;
  const NumericMatrix X(X_);
  const int n = delta.size();
  const int p = X.ncol();
  const int l0 = as<int>(l0_);
  const double adapter = as<double>(adapter_);
  const arma::mat xpred = as<arma::mat>(xpred_); //npred by p;
  const int npred = xpred.n_rows;
  
  // things about spatial copula
  arma::mat ds0block = as<mat>(ds0block_); // n by npred;
  arma::mat ds0m = as<mat>(ds0m_); // m by npred
  arma::mat ds0n = as<mat>(ds0n_); // n by npred
  const arma::mat dnn = as<mat>(dnn_); // n by n
  const arma::mat dnm = as<mat>(dnm_); // n by m
  const arma::mat dmm = as<mat>(dmm_); // m by m
  const arma::mat clustindx = as<mat>(clustindx_); // n by nclust;
  arma::vec theta = as<vec>(theta_);
  const NumericVector theta0(theta0_);
  const double theta1a = theta0[0];
  const double theta1b = theta0[1];
  const double theta2a = theta0[2];
  const double theta2b = theta0[3];
  const arma::mat spS0 = as<mat>(spS0_); // 2 by 2
  
  // prameters to be updated
  Rcpp::NumericVector beta(beta_); // p by 1
  double hcen=as<double>(hcen_);  
  Rcpp::NumericVector h(h_); // M1 by 1;
  Rcpp::NumericVector d(d_); // M1 by 1;
  const int M1 = h.size(); // M1=M+1
  
  // hyperparameters
  // prior of beta
  const arma::vec beta0 = as<arma::vec>(beta0_); // p by 1
  const arma::mat S0inv = as<arma::mat>(S0inv_); // p by p
  const arma::mat Shat = as<arma::mat>(Shat_);   // p by p
  // prior of hcen
  const double h0 = as<double>(h0_);
  const double v0 = as<double>(v0_);
  const double vhat = as<double>(vhat_);
  const double r0 = as<double>(r0_);
  
  // temp variables
  const int nscan = nburn + (nskip+1)*nsave;
  Rcpp::NumericVector Xbeta(n, 0.0);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  
  // things to save;
  Rcpp::NumericVector hcen_save(nsave);
  Rcpp::NumericMatrix t_save(n, nsave);
  Rcpp::NumericMatrix h_save(M1, nsave);
  Rcpp::NumericMatrix d_save(M1, nsave);
  Rcpp::NumericMatrix beta_save(p, nsave);
  Rcpp::NumericMatrix Tpred(npred, nsave);
  Rcpp::NumericMatrix theta_save(2, nsave);
  arma::mat z_save(n, nsave);
  NumericMatrix Zpred(npred, nsave);
  double rejhcen = 0;
  double rejbeta = 0;
  // spatial related
  arma::mat spSnew(2,2); spSnew.fill(0.0);
  arma::vec thetabarnew(2); thetabarnew.fill(0.0);
  arma::vec z(n);
  arma::mat Cinv = arma::eye(n,n);
  double logdetC=0;
  Rcpp::NumericVector rejh(M1);
  double rejtheta=0;
  
  // Make arma objects
  arma::mat X_r(const_cast<NumericMatrix&>(X).begin(), n, p, false, true);
  arma::vec beta_r(beta.begin(), p, false, true);
  arma::vec Xbeta_r(Xbeta.begin(), n, false, true);
  
  // Working temp variables
  arma::mat Linv=arma::zeros<arma::mat>(n, nsave);
  Rcpp::NumericVector t(n); 
  for (int i=0; i<n; ++i){
    if(delta[i]==0){
      t[i] = tobs[i]*1.5;
    }else{
      t[i] = tobs[i];
    }
  }
  Rcpp::IntegerVector Mt(n);
  Rcpp::IntegerVector mk(M1);
  Rcpp::NumericVector lk(M1);  
  double llold, llnew;
  double ratio, uu, nn;
  // for beta
  arma::vec betaold(p);
  arma::vec beBarold(p);
  arma::vec beBarnew=arma::zeros<arma::vec>(p);
  arma::mat beSnew=arma::zeros<arma::mat>(p,p);
  arma::mat Ip = ESMALL*arma::eye(p,p);
  // for hcen
  double hbarold=0, hbarnew=0, hsnew=0, htemp=0;
  
  RNGScope scope;
  
  ////////////
  // Start MCMC
  ////////////
  Xbeta_r = X_r*(beta_r);
  //d = Cutpoints(hcen, M1);
  arma::mat Cnn = arma::exp(-theta[1]*dnn); 
  arma::mat Cnm = arma::exp(-theta[1]*dnm);
  arma::mat Cmm = arma::exp(-theta[1]*dmm);
  Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
  inv_FSA(theta[0], Cnn, Cnm, Cmm, clustindx, Cinv, logdetC);
  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    //Rprintf( "scan = %d\n", iscan );
    
    // Sample hcen 
    if(v0>ESMALL){
      if(iscan>l0){
        htemp = Rf_rnorm(hcen, std::sqrt(hsnew));
      }else{
        htemp = Rf_rnorm(hcen, std::sqrt(vhat));
      }
      if(htemp<ESMALL){
        if(iscan>=nburn) rejhcen += 1;
      }else{
        double sumlogh=0;
        double M = M1-1;
        for(int k=1; k<M1; ++k) sumlogh += std::log(h[k]);
        llold = M*r0*hcen*std::log(r0) - M*Rf_lgammafn(r0*hcen) + r0*hcen*sumlogh - 0.5*std::pow((hcen-h0),2)/v0;
        llnew = M*r0*htemp*std::log(r0) - M*Rf_lgammafn(r0*htemp) + r0*htemp*sumlogh - 0.5*std::pow((htemp-h0),2)/v0;
        ratio = exp(llnew-llold);
        uu = unif_rand();
        if(uu<ratio){
          hcen = htemp; 
        }else{
          if(iscan>=nburn) rejhcen += 1;
        }
      }
      nn = iscan+1;
      hbarold = hbarnew;
      hbarnew = (nn)/(nn+1.0)*hbarold + hcen/(nn+1.0);
      hsnew = (nn-1.0)/nn*hsnew + adapter/nn*(nn*pow(hbarold,2) - (nn+1.0)*pow(hbarnew,2) + pow(hcen,2) + ESMALL );
      d = Cutpoints(hcen, M1);
      Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    }
    
    // Sample t_i when for censored data
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        double zobsi1 = Rf_qnorm5(Foft(tobs[i], h, d, Xbeta[i]), 0, 1, true, false);
        double Cii = Cinv(i,i);
        double s2i = 1.0/Cii;
        double mui = -s2i*( arma::dot(Cinv.col(i), z) - Cii*z[i] );
        double zi = trun_rnorm(mui, std::sqrt(s2i), zobsi1, R_PosInf); 
        double ui = Rf_pnorm5(zi, 0, 1.0, true, false);
        t[i] = Finvofu(ui, h, d, Xbeta[i], tobs[i], ELARGE);
        z[i] = Rf_qnorm5( Foft(t[i], h, d, Xbeta[i]), 0, 1, true, false);
      }
    }
    //Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
    // Sample h_k for k =1, ..., M. In following, M1=M+1
    GetMt(Mt, t, d);
    Getmk(mk, Mt);
    Getlk(lk, Mt, M1, d, t, Xbeta);
    for (int k=1; k<M1; ++k){
      double ha = r0*hcen + mk[k];
      double hb = r0 + lk[k];
      double tmpold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double hkold = h[k]; 
      h[k] = Rf_rgamma(ha, 1.0/hb);
      Getz(z, t, h, d, Xbeta, n);
      double tmpnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tmpnew-tmpold);
      double uu = unif_rand();
      if(uu>ratio) {
        h[k] = hkold; 
        if(iscan>=nburn) rejh[k]+=1.0;
      }
    }
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
    // Sample beta;
    arma::vec Lamb0 = Lambda0tvec(t, h, d);
    llold = arma::sum( -Lamb0%arma::exp(Xbeta_r) + Xbeta_r ) -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) )
      -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
    betaold = beta_r;
    if(iscan>l0){
      beta_r = mvrnorm(betaold, beSnew);
    }else{
      beta_r = mvrnorm(betaold, Shat);
    }
    Xbeta_r = X_r*(beta_r);
    Getz(z, t, h, d, Xbeta, n);
    llnew = arma::sum( -Lamb0%arma::exp(Xbeta_r) + Xbeta_r ) -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) )
      -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
    ratio = exp(llnew-llold);
    uu = unif_rand();
    if(uu>ratio){
      beta_r=betaold;
      if(iscan>=nburn) rejbeta+=1.0;
    }
    nn = iscan+1;
    beBarold = beBarnew;
    beBarnew = (nn)/(nn+1.0)*beBarold + beta_r/(nn+1.0);
    beSnew = (nn-1.0)/nn*beSnew + adapter/(p+0.0)/nn*(nn*beBarold*beBarold.t() 
                                                        - (nn+1.0)*beBarnew*beBarnew.t() + beta_r*beta_r.t() + Ip );
    Xbeta_r = X_r*(beta_r);
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
    // Sample theta;
    spCopula_sample_theta_FSA(theta, rejtheta, spSnew, thetabarnew, Cinv, logdetC, 
                              theta1a, theta1b, theta2a, theta2b, l0, spS0, dnn, 
                              adapter, iscan, nburn, z, n, dnm, dmm, clustindx);
    
    // Save the sample
    if (iscan>=nburn) {
      ++skiptally;
      if (skiptally>nskip){
        // calculate Linv
        Linv.col(isave) = LinvSpCopulaCox(tobs, delta, Xbeta, h, d, Cinv, z);
        // save samples
        t_save(_,isave) = t;
        h_save(_,isave) = h;
        d_save(_,isave) = d;
        beta_save(_,isave) = beta;
        theta_save(0,isave) = theta[0];
        theta_save(1,isave) = theta[1];
        z_save.col(isave) = z;
        hcen_save[isave] = hcen;
        // prediction
        arma::vec xbeta = xpred*beta_r;
        arma::mat rho(n, npred);
        if(((int)dmm.n_cols)==n){
          rho = arma::exp(-theta[1]*ds0n);
        }else{
          arma::mat rho1 = arma::exp(-theta[1]*dnm)*arma::solve(arma::exp(-theta[1]*dmm), arma::exp(-theta[1]*ds0m));
          arma::mat rho2 = (arma::exp(-theta[1]*ds0n)-rho1)%ds0block;
          rho = rho1 + rho2;
        }
        for(int j=0; j<npred; ++j){
          arma::vec hs0 = theta[0]*rho.col(j);
          double mus0 = arma::dot( hs0, Cinv*z );
          double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
          double znew = Rf_rnorm(mus0, sigs0);
          double u = Rf_pnorm5(znew, 0, 1.0, true, false);
          Zpred(j, isave) = znew;
          Tpred(j, isave) = Finvofu(u, h, d, xbeta[j], ESMALL, ELARGE);
        }
        
        ++isave;
        ++distally;
        if (distally>=ndisplay){
          Rprintf( "scan = %d\n", isave );
          distally = 0;
        }
        skiptally=0;
      }
    }
  }
  // get acceptance rate
  double totscans = nscan-nburn+0.0;
  double ratebeta = 1.0 - rejbeta/totscans;
  double ratetheta = 1.0-rejtheta/totscans;
  double ratehcen = 1.0 - rejhcen/totscans;
  NumericVector rateh = 1.0-rejh/totscans;
  
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("h")=h_save,
                      Named("d")=d_save,
                      Named("beta")=beta_save,
                      Named("theta")=theta_save,
                      Named("t")=t_save,
                      Named("z")=z_save,
                      Named("hcen")=hcen_save,
                      Named("ratebeta")=ratebeta,
                      Named("ratetheta")=ratetheta,
                      Named("rateh")=rateh,
                      Named("ratehcen")=ratehcen,
                      Named("cpo")=cpo,
                      Named("Tpred")=Tpred,
                      Named("Zpred")=Zpred);
  END_RCPP
}

