#include "spSurv_Coxph_tools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

///////////////////////////// Fixed cutpoints //////////////////////////////////
RcppExport SEXP indeptCoxph( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
		SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, 
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_) {
	BEGIN_RCPP
	
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  int ndisplay = as<int>(ndisplay_);
  const NumericVector tobs(t_);
  IntegerVector delta(delta_);
  const NumericMatrix Xr(X_);
  const int n = tobs.size();
  const int l0 = as<int>(l0_);
  const arma::mat S0 = as<mat>(S0_);
  const double adapter = as<double>(adapter_);
  
  // prameters to be updated
  NumericVector beta_r(beta_);
  const int p = beta_r.size();
  NumericVector h(h_);
  NumericVector d(d_);
  const int M1 = h.size(); // M1=M+1
  
  // hyperparameters
  const double h0 = as<double>(h0_);
  const double r0 = as<double>(r0_);
  const arma::vec mu0 = as<vec>(mu0_);
  const arma::mat Sig0 = as<mat>(Sig0_);
  const int nscan = nburn + (nskip+1)*nsave;
  arma::mat xpred = as<mat>(xpred_); 
  int npred = xpred.n_rows;
  
  // Temp variables
  arma::mat Snew(p,p); Snew.fill(0.0);
  arma::vec betabarnew(p); betabarnew.fill(0.0);
  NumericVector t(n); for (int i=0; i<n; ++i) t[i] = tobs[i];
  IntegerVector Mt(n);
  Rcpp::IntegerVector mk(M1);
  Rcpp::NumericVector lk(M1);  
  int skiptally=0; 
  int isave=0;
  int distally=0;
  arma::mat Linv(n, nsave);
  int rejbeta = 0;

  // Make arma objects
  arma::mat X(const_cast<NumericMatrix&>(Xr).begin(), p, n, false);
  arma::vec beta(beta_r.begin(), p, false);
  
  // things to save;
  NumericMatrix t_save(n, nsave);
  NumericMatrix h_save(M1, nsave);
  NumericMatrix d_save(M1, nsave);
  arma::mat beta_save(p, nsave);
  NumericMatrix Tpred(npred, nsave);
  
  RNGScope scope;
	
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);

  ////////////
  // Start MCMC
  ////////////

  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    // Rprintf( "scan = %d\n", iscan );
    // Sample t_i when delta_i=0
    arma::vec Xbeta = X.t()*beta;
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        double ui = Rf_runif( Foft(tobs[i], h, d, Xbeta[i]), 1.0 );
        t[i] = Finvofu(ui, h, d, Xbeta[i], tobs[i], ELARGE);
      }
    }
    
    // Sample h_k for k =1, ..., M. In following, M1=M+1
    GetMt(Mt, t, d);
    Getmk(mk, Mt);
    Getlk(lk, Mt, M1, d, t, Xbeta);
    for (int k=1; k<M1; ++k){
      double ha = r0*h0 + mk[k];
      double hb = r0 + lk[k];
      h[k] = Rf_rgamma(ha, 1.0/hb);
    }
  
    // Sample beta;
    arma::vec Lamb0 = Lambda0tvec(t, h, d);
    indept_sample_beta(beta, rejbeta, Snew, betabarnew, X, 
                Lamb0, mu0, Sig0, p, l0, S0, adapter, iscan);
  
    // Save the sample
    if (iscan>=nburn) {
      ++skiptally;
      if (skiptally>nskip){
        // calculate Linv
        Linv.col(isave) = LinvIndeptCox(tobs, delta, X.t()*beta, h, d);
        // save samples
        t_save(_,isave) = t;
        h_save(_,isave) = h;
        d_save(_,isave) = d;
        beta_save.col(isave) = beta;
        // prediction
        arma::vec xbeta = xpred*beta;
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
  double ratebeta = 1.0-(double)rejbeta/(nscan+0.0);
  
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("h")=h_save,
                      Named("d")=d_save,
                      Named("beta")=beta_save,
                      Named("t")=t_save,
                      Named("ratebeta")=ratebeta,
                      Named("cpo")=cpo,
                      Named("Tpred")=Tpred);
	END_RCPP
}

///////////////////////////// Random cutpoints //////////////////////////////////
RcppExport SEXP indeptCoxphR( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, SEXP V0_,
    SEXP hl0_, SEXP hs0_, SEXP hadapter_,
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_) {
	BEGIN_RCPP
	
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  int ndisplay = as<int>(ndisplay_);
  const NumericVector tobs(t_);
  IntegerVector delta(delta_);
  const NumericMatrix Xr(X_);
  const int n = tobs.size();
  const int l0 = as<int>(l0_);
  const arma::mat S0 = as<mat>(S0_);
  const double adapter = as<double>(adapter_);

  // prameters to be updated
  NumericVector beta_r(beta_);
  const int p = beta_r.size();
  NumericVector h(h_);
  NumericVector d(d_);
  const int M1 = h.size(); // M1=M+1
  
  // hyperparameters
  const double h0 = as<double>(h0_);
  const double r0 = as<double>(r0_);
  double hcen = h0;
  const double V0 = as<double>(V0_);
  const int hl0 = as<int>(hl0_);
  const double hs0 = as<double>(hs0_);
  const double hadapter = as<double>(hadapter_);
  const arma::vec mu0 = as<vec>(mu0_);
  const arma::mat Sig0 = as<mat>(Sig0_);
  const int nscan = nburn + (nskip+1)*nsave;
  arma::mat xpred = as<mat>(xpred_); 
  int npred = xpred.n_rows;
  
  // Temp variables
  double hSnew = 0;
  double hbarnew = 0;
  int rejhcen = 0;
  arma::mat Snew(p,p); Snew.fill(0.0);
  arma::vec betabarnew(p); betabarnew.fill(0.0);
  NumericVector t(n); for (int i=0; i<n; ++i) t[i] = tobs[i];
  IntegerVector Mt(n);
  Rcpp::IntegerVector mk(M1);
  Rcpp::NumericVector lk(M1);  
  int skiptally=0; 
  int isave=0;
  int distally=0;
  arma::mat Linv(n, nsave);
  int rejbeta = 0;

  // Make arma objects
  arma::mat X(const_cast<NumericMatrix&>(Xr).begin(), p, n, false);
  arma::vec beta(beta_r.begin(), p, false);
  
  // things to save;
  NumericVector hcen_save(nsave);
  NumericMatrix t_save(n, nsave);
  NumericMatrix h_save(M1, nsave);
  NumericMatrix d_save(M1, nsave);
  arma::mat beta_save(p, nsave);
  NumericMatrix Tpred(npred, nsave);
  
  RNGScope scope;
	
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);

  ////////////
  // Start MCMC
  ////////////

  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    // Rprintf( "scan = %d\n", iscan );
    
    // Update hcen and d
    sample_hcen(hcen, rejhcen, hSnew, hbarnew, h, r0, h0, V0, hl0, hs0, hadapter, iscan);
    d = Cutpoints(hcen, M1);
    
    // Sample t_i when delta_i=0
    arma::vec Xbeta = X.t()*beta;
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        double ui = Rf_runif( Foft(tobs[i], h, d, Xbeta[i]), 1.0 );
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
    indept_sample_beta(beta, rejbeta, Snew, betabarnew, X, 
                Lamb0, mu0, Sig0, p, l0, S0, 
                adapter, iscan);
  
    // Save the sample
    if (iscan>=nburn) {
      ++skiptally;
      if (skiptally>nskip){
        // calculate Linv
        Linv.col(isave) = LinvIndeptCox(tobs, delta, X.t()*beta, h, d);
        // save samples
        t_save(_,isave) = t;
        h_save(_,isave) = h;
        d_save(_,isave) = d;
        beta_save.col(isave) = beta;
        hcen_save[isave] = hcen;
        // prediction
        arma::vec xbeta = xpred*beta;
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
  double ratebeta = 1.0-(double)rejbeta/(nscan+0.0);
  double ratehcen = 1.0-(double)rejhcen/(nscan+0.0);
  
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

