#include "spSurv_spCopulaCoxph.h"
#include "spSurv_Coxph_tools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

///////////////////////////// Fixed cutpoints //////////////////////////////////
SEXP spCopulaCoxph( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
		SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, 
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, SEXP theta0_, 
    SEXP spl0_, SEXP spS0_, SEXP spadapter_ ) {
	BEGIN_RCPP
	
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  int ndisplay = as<int>(ndisplay_);
  const NumericVector tobs(t_);
  IntegerVector delta(delta_);
  const arma::mat X = as<mat>(X_);
  const int n = tobs.size();
  const int l0 = as<int>(l0_);
  const arma::mat S0 = as<mat>(S0_);
  const double adapter = as<double>(adapter_);
  
  // things about spatial copula
  arma::mat ds0n = as<mat>(ds0n_);
  const arma::mat dnn = as<mat>(dnn_);
  arma::vec theta = as<vec>(theta_);
  const NumericVector theta0(theta0_);
  const double theta1a = theta0[0];
  const double theta1b = theta0[1];
  const double theta2a = theta0[2];
  const double theta2b = theta0[3];
  const int spl0 = as<int>(spl0_);
  const arma::mat spS0 = as<mat>(spS0_);
  const double spadapter = as<double>(spadapter_);
  
  // prameters to be updated
  arma::vec beta = as<vec>(beta_);
  const int p = beta.size();
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
  double Maxobs = std::max(Rcpp::max(tobs)*5, 1200.0);
  arma::vec Xbeta = X.t()*beta;
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
  // spatial related
  arma::mat spSnew(2,2); spSnew.fill(0.0);
  arma::vec thetabarnew(2); thetabarnew.fill(0.0);
  arma::vec z(n);
  arma::mat Cinv = arma::eye(n,n);
  double logdetC=0;
  Rcpp::NumericVector rejh(M1);
  int rejtheta=0;
  
  // things to save;
  NumericMatrix t_save(n, nsave);
  NumericMatrix h_save(M1, nsave);
  NumericMatrix d_save(M1, nsave);
  arma::mat beta_save(p, nsave);
  NumericMatrix Tpred(npred, nsave);
  NumericVector theta1_save(nsave);
  NumericVector theta2_save(nsave);
  arma::mat z_save(n, nsave);
  NumericMatrix Zpred(npred, nsave);
  
  RNGScope scope;
  
	// Set the Armadillo seed from R's 
	// int seed = (int)Rf_runif(0.0, 10000.0);
	// std::srand(seed);

  ////////////
  // Start MCMC
  ////////////
  Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
  GetCinv(n, theta[0], theta[1], dnn, Cinv, logdetC);

  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    //Rprintf( "scan = %d\n", iscan );
    
    // Sample t_i when delta_i=0
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        double Ftobsi = Foft(tobs[i], h, d, Xbeta[i]);
        double zobsi = Rf_qnorm5(Ftobsi, 0, 1, true, false);
        double Cii = Cinv(i,i);
        double s2i = 1.0/Cii;
        double mui = -s2i*( arma::dot(Cinv.col(i), z) - Cii*z[i] );
        double zi = trun_rnorm(mui, std::sqrt(s2i), zobsi, R_PosInf); 
        double ui = Rf_pnorm5(zi, 0, 1, true, false);
        t[i] = Finvofu(ui, h, d, Xbeta[i], tobs[i], Maxobs);
        z[i] = Rf_qnorm5( Foft(t[i], h, d, Xbeta[i]), 0, 1, true, false);
      }
    }
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
    // Sample h_k for k =1, ..., M. In following, M1=M+1
    GetMt(Mt, t, d);
    Getmk(mk, Mt);
    Getlk(lk, Mt, M1, d, t, Xbeta);
    for (int k=1; k<M1; ++k){
      double ha = r0*h0 + mk[k];
      double hb = r0 + lk[k];
      double tmpold = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double hkold = h[k]; 
      h[k] = Rf_rgamma(ha, 1.0/hb); 
      Getz(z, t, h, d, Xbeta, n);
      double tmpnew = -0.5*arma::dot(z, Cinv*z) + 0.5*arma::dot(z, z);
      double ratio = std::exp(tmpnew-tmpold);
      double uu = unif_rand();
      if(uu>ratio) {h[k] = hkold; ++rejh[k];}
    }
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
  
    // Sample beta;
    arma::vec Lamb0 = Lambda0tvec(t, h, d);
    spCopula_sample_beta(beta, rejbeta, Snew, betabarnew, X, 
                Lamb0, mu0, Sig0, p, l0, S0, adapter, iscan,
                z, Cinv, t, h, d, n);
    Xbeta = X.t()*beta;
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
    // Sample theta;
    //Rprintf( "%f\n", spSnew(0,0) );
    spCopula_sample_theta(theta, rejtheta, spSnew, thetabarnew, Cinv, logdetC, 
                 theta1a, theta1b, theta2a, theta2b, spl0, spS0, dnn, spadapter, iscan, z, n);
                 
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
        beta_save.col(isave) = beta;
        theta1_save[isave] = theta[0];
        theta2_save[isave] = theta[1];
        z_save.col(isave) = z;
        // prediction
        arma::vec xbeta = xpred*beta;
        for(int j=0; j<npred; ++j){
          arma::vec hs0 = theta[0]*arma::exp(-theta[1]*ds0n.col(j));
          double mus0 = arma::dot( hs0, Cinv*z );
          double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
          double znew = Rf_rnorm(mus0, sigs0);
          double u = Rf_pnorm5(znew, 0, 1.0, true, false);
          Zpred(j, isave) = znew;
          Tpred(j, isave) = Finvofu(u, h, d, xbeta[j], ESMALL, Maxobs);
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
  double ratetheta = 1.0-(double)rejtheta/(nscan+0.0);
  NumericVector rateh = 1.0-rejh/(nscan+0.0);
  
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("h")=h_save,
                      Named("d")=d_save,
                      Named("beta")=beta_save,
                      Named("theta1")=theta1_save,
                      Named("theta2")=theta2_save,
                      Named("t")=t_save,
                      Named("z")=z_save,
                      Named("ratebeta")=ratebeta,
                      Named("ratetheta")=ratetheta,
                      Named("rateh")=rateh,
                      Named("cpo")=cpo,
                      Named("Tpred")=Tpred,
                      Named("Zpred")=Zpred);
	END_RCPP
	}

///////////////////////////// Random cutpoints //////////////////////////////////
SEXP spCopulaCoxphR( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, SEXP V0_, 
    SEXP hl0_, SEXP hs0_, SEXP hadapter_, 
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, SEXP theta0_, 
    SEXP spl0_, SEXP spS0_, SEXP spadapter_ ) {
	BEGIN_RCPP
	
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  int ndisplay = as<int>(ndisplay_);
  const NumericVector tobs(t_);
  IntegerVector delta(delta_);
  const arma::mat X = as<mat>(X_);
  const int n = tobs.size();
  const int l0 = as<int>(l0_);
  const arma::mat S0 = as<mat>(S0_);
  const double adapter = as<double>(adapter_);
  
  // things about spatial copula
  arma::mat ds0n = as<mat>(ds0n_);
  const arma::mat dnn = as<mat>(dnn_);
  arma::vec theta = as<vec>(theta_);
  const NumericVector theta0(theta0_);
  const double theta1a = theta0[0];
  const double theta1b = theta0[1];
  const double theta2a = theta0[2];
  const double theta2b = theta0[3];
  const int spl0 = as<int>(spl0_);
  const arma::mat spS0 = as<mat>(spS0_);
  const double spadapter = as<double>(spadapter_);
  
  // prameters to be updated
  arma::vec beta = as<vec>(beta_);
  const int p = beta.size();
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
  double Maxobs = std::max(Rcpp::max(tobs)*5, 1200.0);
  arma::vec Xbeta = X.t()*beta;
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
  // spatial related
  arma::mat spSnew(2,2); spSnew.fill(0.0);
  arma::vec thetabarnew(2); thetabarnew.fill(0.0);
  arma::vec z(n);
  arma::mat Cinv = arma::eye(n,n);
  double logdetC=0;
  Rcpp::NumericVector rejh(M1);
  int rejtheta=0;
  
  // things to save;
  NumericVector hcen_save(nsave);
  NumericMatrix t_save(n, nsave);
  NumericMatrix h_save(M1, nsave);
  NumericMatrix d_save(M1, nsave);
  arma::mat beta_save(p, nsave);
  NumericMatrix Tpred(npred, nsave);
  NumericVector theta1_save(nsave);
  NumericVector theta2_save(nsave);
  arma::mat z_save(n, nsave);
  NumericMatrix Zpred(npred, nsave);
  
  RNGScope scope;
  
	// Set the Armadillo seed from R's 
	// int seed = (int)Rf_runif(0.0, 10000.0);
	// std::srand(seed);

  ////////////
  // Start MCMC
  ////////////
  Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
  GetCinv(n, theta[0], theta[1], dnn, Cinv, logdetC);

  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    //Rprintf( "scan = %d\n", iscan );
    
    // Update hcen and d
    sample_hcen(hcen, rejhcen, hSnew, hbarnew, h, r0, h0, V0, hl0, hs0, hadapter, iscan);
    d = Cutpoints(hcen, M1);
    Getz(z, t, h, d, Xbeta, n);
    
    // Sample t_i when delta_i=0
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        double Ftobsi = Foft(tobs[i], h, d, Xbeta[i]);
        double zobsi = Rf_qnorm5(Ftobsi, 0, 1, true, false);
        double Cii = Cinv(i,i);
        double s2i = 1.0/Cii;
        double mui = -s2i*( arma::dot(Cinv.col(i), z) - Cii*z[i] );
        double zi = trun_rnorm(mui, std::sqrt(s2i), zobsi, R_PosInf); 
        double ui = Rf_pnorm5(zi, 0, 1.0, true, false);
        t[i] = Finvofu(ui, h, d, Xbeta[i], tobs[i], Maxobs);
        z[i] = Rf_qnorm5( Foft(t[i], h, d, Xbeta[i]), 0, 1, true, false);
      }
    }
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
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
      if(uu>ratio) {h[k] = hkold; ++rejh[k];}
    }
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
  
    // Sample beta;
    arma::vec Lamb0 = Lambda0tvec(t, h, d);
    spCopula_sample_beta(beta, rejbeta, Snew, betabarnew, X, 
                Lamb0, mu0, Sig0, p, l0, S0, adapter, iscan,
                z, Cinv, t, h, d, n);
    Xbeta = X.t()*beta;
    Getz(z, t, h, d, Xbeta, n); // Get transformed survival time z
    
    // Sample theta;
    //Rprintf( "%f\n", spSnew(0,0) );
    spCopula_sample_theta(theta, rejtheta, spSnew, thetabarnew, Cinv, logdetC, 
                 theta1a, theta1b, theta2a, theta2b, spl0, spS0, dnn, spadapter, iscan, z, n);
    //Rprintf( "%f\n", spSnew(1,1) );
    
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
        beta_save.col(isave) = beta;
        theta1_save[isave] = theta[0];
        theta2_save[isave] = theta[1];
        z_save.col(isave) = z;
        hcen_save[isave] = hcen;
        // prediction
        arma::vec xbeta = xpred*beta;
        for(int j=0; j<npred; ++j){
          arma::vec hs0 = theta[0]*arma::exp(-theta[1]*ds0n.col(j));
          double mus0 = arma::dot( hs0, Cinv*z );
          double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
          double znew = Rf_rnorm(mus0, sigs0);
          double u = Rf_pnorm5(znew, 0, 1.0, true, false);
          Zpred(j, isave) = znew;
          Tpred(j, isave) = Finvofu(u, h, d, xbeta[j], ESMALL, Maxobs);
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
  double ratetheta = 1.0-(double)rejtheta/(nscan+0.0);
  double ratehcen = 1.0-(double)rejhcen/(nscan+0.0);
  NumericVector rateh = 1.0-rejh/(nscan+0.0);
  
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("h")=h_save,
                      Named("d")=d_save,
                      Named("beta")=beta_save,
                      Named("theta1")=theta1_save,
                      Named("theta2")=theta2_save,
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

