#include "spSurv_Coxph_tools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

///////////////////////////// Random cutpoints //////////////////////////////////
RcppExport SEXP indeptCoxph(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
                             SEXP tobs_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, 
                             SEXP r0_, SEXP hcen_, SEXP h0_, SEXP v0_, SEXP vhat_,
                             SEXP beta_, SEXP beta0_, SEXP S0inv_, SEXP Shat_, 
                             SEXP l0_, SEXP adapter_, SEXP xpred_) {
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
  double rejhcen = 0;
  double rejbeta = 0;
  
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
  // d = Cutpoints(hcen, M1);
  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    // Rprintf( "scan = %d\n", iscan );
    
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
    }
    
    // Sample t_i when for censored data
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        double ui = Rf_runif( Foft(tobs[i], h, d, Xbeta[i]), 1.0-ESMALL );
        t[i] = Finvofu(ui, h, d, Xbeta[i], tobs[i], ELARGE);
      }
    }
    
    // Sample h_k for k =1, ..., M. In following, M1=M+1
    GetMt(Mt, t, d);
    Getmk(mk, Mt);
    Getlk(lk, Mt, M1, d, t, Xbeta);
    for (int k=1; k<M1; ++k){
      double ha = r0*hcen + mk[k];
      double hb = r0 + lk[k];
      h[k] = Rf_rgamma(ha, 1.0/hb);
    }
    
    // Sample beta;
    arma::vec Lamb0 = Lambda0tvec(t, h, d);
    llold = arma::sum( -Lamb0%arma::exp(Xbeta_r) + Xbeta_r ) -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
    betaold = beta_r;
    if(iscan>l0){
      beta_r = mvrnorm(betaold, beSnew);
    }else{
      beta_r = mvrnorm(betaold, Shat);
    }
    Xbeta_r = X_r*(beta_r);
    llnew = arma::sum( -Lamb0%arma::exp(Xbeta_r) + Xbeta_r ) -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
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
    
    // Save the sample
    if (iscan>=nburn) {
      ++skiptally;
      if (skiptally>nskip){
        // calculate Linv
        Linv.col(isave) = LinvIndeptCox(tobs, delta, Xbeta, h, d);
        // save samples
        t_save(_,isave) = t;
        h_save(_,isave) = h;
        d_save(_,isave) = d;
        beta_save(_,isave) = beta;
        hcen_save[isave] = hcen;
        // prediction
        arma::vec xbeta = xpred*beta_r;
        for(int j=0; j<npred; ++j){
          double u = unif_rand();
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
  double ratehcen = 1.0 - rejhcen/totscans;
  
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("h")=h_save,
                      Named("d")=d_save,
                      Named("beta")=beta_save,
                      Named("t")=t_save,
                      Named("hcen")=hcen_save,
                      Named("ratebeta")=ratebeta,
                      Named("ratehcen")=ratehcen,
                      Named("cpo")=cpo,
                      Named("Tpred")=Tpred);
	END_RCPP
}

// Get density or survival Plots for Cox PH
RcppExport SEXP CoxPHplots(SEXP xpred_, SEXP tgrid_, SEXP beta_, SEXP h_, SEXP d_, SEXP CI_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  Rcpp::NumericVector tgrid(tgrid_);
  arma::mat xpred = as<arma::mat>(xpred_); // npred by p;
  arma::mat beta = as<arma::mat>(beta_); // p by nsave;
  Rcpp::NumericMatrix h(h_);
  Rcpp::NumericMatrix d(d_);
  double CI = as<double>(CI_);
  int nsave = h.ncol();
  int ngrid = tgrid.size();
  int npred = xpred.n_rows;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  
  // Temp variables
  arma::vec xbeta(npred);
  Rcpp::NumericVector estfArray(nsave*ngrid*npred);
  arma::cube estf(estfArray.begin(), ngrid, nsave, npred, false, true);
  Rcpp::NumericVector estSArray(nsave*ngrid*npred);
  arma::cube estS(estSArray.begin(), ngrid, nsave, npred, false, true);
  Rcpp::NumericVector esthArray(nsave*ngrid*npred);
  arma::cube esth(esthArray.begin(), ngrid, nsave, npred, false, true);
  
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
    xbeta = xpred*beta.col(i);
    Rcpp::NumericVector hi = h(_,i);
    Rcpp::NumericVector di = d(_,i);
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        estf(k, i, j) = foft(tgrid[k], hi, di, xbeta[j]);
        estS(k, i, j) = 1.0 - Foft(tgrid[k], hi, di, xbeta[j]);
        esth(k, i, j) = lambdat(tgrid[k], hi, di, xbeta[j]);
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
