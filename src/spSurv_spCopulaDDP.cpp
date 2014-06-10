#include "spSurv_spCopulaDDP.h"
#include "spSurv_DDP_tools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

SEXP spCopulaDDP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
    SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
  	SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
		SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
    SEXP nua_, SEXP nub_, SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, 
    SEXP theta0_, SEXP spl0_, SEXP spS0_, SEXP spadapter_) {
	BEGIN_RCPP
	
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  int ndisplay = as<int>(ndisplay_);
  const int N = as<int>(N_);
  const NumericVector yobs(y_);
  IntegerVector delta(delta_);
  const arma::mat X = as<mat>(X_); // p by n
  const int p = X.n_rows;
  const int n = X.n_cols;
  
  // things about spatial copula
  arma::mat ds0n = as<mat>(ds0n_); // n by npred
  const arma::mat dnn = as<mat>(dnn_); // n by n
  arma::vec theta = as<vec>(theta_);
  const NumericVector theta0(theta0_);
  const double theta1a = theta0[0];
  const double theta1b = theta0[1];
  const double theta2a = theta0[2];
  const double theta2b = theta0[3];
  const int spl0 = as<int>(spl0_);
  const arma::mat spS0 = as<mat>(spS0_); // 2 by 2
  const double spadapter = as<double>(spadapter_);
  
  // prameters to be updated
  arma::mat beta = as<mat>(beta_); //p by N
  arma::vec mu = as<vec>(mu_); // p by 1
  arma::mat Sig = as<mat>(Sig_); // p by p
  arma::vec w = as<vec>(w_);
  NumericVector tau2(tau2_);
  NumericVector tau = Rcpp::sqrt(tau2);
  IntegerVector K(K_);
  NumericVector V(V_);
  double alpha = as<double>(alpha_);

  // hyperparameters 
  const arma::vec m0 = as<vec>(m0_);
  const arma::mat S0 = as<mat>(S0_);
  const arma::mat Sig0 = as<mat>(Sig0_);
  const int k0 = as<int>(k0_);
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  const double nua = as<double>(nua_);
  const double nub = as<double>(nub_);
  const int nscan = nburn + (nskip+1)*nsave;
  arma::mat xpred = as<mat>(xpred_); 
  int npred = xpred.n_rows;
  
  // Temp variables
  double MinRes = Rcpp::min(yobs)-3.0;
  double MaxRes = Rcpp::max(yobs)+3.0;
  arma::mat Xbeta = X.t()*beta;
  NumericVector y(n);  for (int i=0; i<n; ++i) y[i] = yobs[i];
  IntegerVector nK(N);
  IntegerMatrix Kind(n, N);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  arma::mat Linv(n, nsave);
  // spatial related
  arma::mat spSnew(2,2); spSnew.fill(0.0);
  arma::vec thetabarnew(2); thetabarnew.fill(0.0);
  arma::vec z(n);
  arma::mat zPhi(n, N);
  arma::mat Cinv = arma::eye(n,n);
  double logdetC=0;
  NumericVector rejy(n);
  NumericVector rejbeta(N);
  NumericVector rejsigma(N);
  NumericVector rejV(N-1);
  int rejtheta=0;

  // Make arma objects
  arma::mat invSig = inv_sympd(Sig);
  const arma::mat invS0 = arma::inv_sympd(S0);
  const arma::vec invS0m0 = arma::solve(S0,m0);
  
  // things to save;
  NumericVector betaArray(nsave*N*p);
  arma::cube beta_save(betaArray.begin(), p, N, nsave, false);
  arma::mat w_save(N, nsave);
  NumericMatrix sigma2_save(N, nsave);
  NumericVector alpha_save(nsave);
  NumericMatrix y_save(n, nsave);
  NumericMatrix Ypred(npred, nsave);
  NumericVector theta1_save(nsave);
  NumericVector theta2_save(nsave);
  arma::mat z_save(n, nsave);
  NumericMatrix Zpred(npred, nsave);
  
  GetRNGstate();
	
	// Set the Armadillo seed from R's 
	// int seed = (int)Rf_runif(0.0, 10000.0);
	// std::srand(seed);

  // From V to w
  DDP_Vtow(w, V, N);

  // initial period
  spCopulaInitial( K, y, beta, tau2, V, w, alpha, mu, Sig, invSig, yobs, delta, X, m0, S0, Sig0, k0, 
                   a0, b0, nua, nub, invS0, invS0m0);
  Xbeta = X.t()*beta;
  tau = Rcpp::sqrt(tau2);
  // get transformed survival time z
  for (int i=0; i<n; ++i){
    for (int k=0; k<N; ++k){
      zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*tau[k], 0, 1.0, true, false);
    }
  }
  z = qnormvec( zPhi*w );
  GetCinv(n, theta[0], theta[1], dnn, Cinv, logdetC); // Preprocess C^{-1} to get Cinv

  //////////////////////////////////////////////////////////////
  // Start MCMC
  //////////////////////////////////////////////////////////////
  
  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    // Rprintf( "scan = %d\n", iscan );

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
    spCopula_sample_y(y, rejy, zPhi, z, w, yobs, delta, Xbeta, tau, K, Cinv, n, N);
  
    // Sample beta;
    spCopula_sample_beta(beta, rejbeta, zPhi, z, w, y, X, tau2, nK, Kind, mu, Sig, invSig, Cinv, n, N, p);
    Xbeta = X.t()*beta;
  
    //Sample simga2;
    spCopula_sample_sigma2(tau2, rejsigma, zPhi, z, w, y, Xbeta, nK, Kind, nua, nub, Cinv, n, N);
    tau = Rcpp::sqrt(tau2);
  
    //Sample V;
    spCopula_sample_V(V, rejV, zPhi, z, w, nK, alpha, Cinv, N);
    DDP_Vtow(w, V, N); // From V to w
  
    //Sample alpha;
    double a0star = a0+N-1;
    double b0star = b0-std::log(w[N-1]);
    alpha = Rf_rgamma(a0star, 1/b0star);
  
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
  
    //Sample theta1 and theta2;
    spCopula_sample_theta(theta, rejtheta, spSnew, thetabarnew, Cinv, logdetC, 
                 theta1a, theta1b, theta2a, theta2b, spl0, spS0, dnn, spadapter, iscan, z, n);
  
    // Save the sample
    if (iscan>=nburn) {
      ++skiptally;
      if (skiptally>nskip){
        // calculate Linv
        arma::vec tautmp = as<vec>(tau);
        arma::vec sig = 1.0/tautmp;
        Linv.col(isave) = spCopula_Linv(yobs, delta, X, beta, sig, w, Cinv, z);
        // save samples
        alpha_save[isave] = alpha;
        y_save(_,isave) = y;
        for (int k=0; k<N; ++k){
          (beta_save.slice(isave)).col(k) = beta.col(k);
          sigma2_save(k, isave) = 1.0/tau2[k];
          w_save(k, isave) = w[k];
        }
        theta1_save[isave] = theta[0];
        theta2_save[isave] = theta[1];
        z_save.col(isave) = z;
        // prediction
        arma::mat xbeta = arma::trans( xpred*beta );
        for(int j=0; j<npred; ++j){
          arma::vec hs0 = theta[0]*arma::exp(-theta[1]*ds0n.col(j));
          double mus0 = arma::dot( hs0, Cinv*z );
          double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
          double znew = Rf_rnorm(mus0, sigs0);
          double u = Rf_pnorm5(znew, 0, 1.0, true, false);
          Zpred(j, isave) = znew;
          Ypred(j, isave) = DDP_Finvofu(u, w, xbeta.col(j), sig, MinRes, MaxRes);
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
  NumericVector ratey = 1.0- rejy/(nscan+0.0);
  NumericVector ratebeta = 1.0-rejbeta/(nscan+0.0);
  NumericVector ratesigma = 1.0-rejsigma/(nscan+0.0);
  NumericVector rateV = 1.0-rejV/(nscan+0.0);
  double ratetheta = 1.0-(double)rejtheta/(nscan+0.0);
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("beta")=beta_save,
                      Named("sigma2")=sigma2_save,
                      Named("w")=w_save,
                      Named("alpha")=alpha_save,
                      Named("theta1")=theta1_save,
                      Named("theta2")=theta2_save,
                      Named("z") = z_save,
                      Named("y")=y_save,
                      Named("ratey") = ratey, 
                      Named("ratebeta")=ratebeta,
                      Named("ratesigma")=ratesigma,
                      Named("rateV")=rateV,
                      Named("ratetheta")=ratetheta,
                      Named("cpo")=cpo,
                      Named("Ypred")=Ypred,
                      Named("Zpred")=Zpred);
	PutRNGstate();
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////// spatial Copula DDP using FSA /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
SEXP spCopulaDDP_FSA( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
    SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
    SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
		SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
    SEXP nua_, SEXP nub_, SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, 
    SEXP theta0_, SEXP spl0_, SEXP spS0_, SEXP spadapter_,
    SEXP dnm_, SEXP dmm_, SEXP blocki_, SEXP ds0m_, SEXP ds0block_) {
	BEGIN_RCPP
	
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  int ndisplay = as<int>(ndisplay_);
  const int N = as<int>(N_);
  const NumericVector yobs(y_);
  IntegerVector delta(delta_);
  const arma::mat X = as<mat>(X_); // p by n
  const int p = X.n_rows;
  const int n = X.n_cols;
  
  // things about spatial copula
  arma::mat ds0block = as<mat>(ds0block_); // n by npred;
  arma::mat ds0m = as<mat>(ds0m_); // m by npred
  arma::mat ds0n = as<mat>(ds0n_); // n by npred
  const arma::mat dnn = as<mat>(dnn_); // n by n
  const arma::mat dnm = as<mat>(dnm_); // n by m
  const arma::mat dmm = as<mat>(dmm_); // m by m
  const IntegerVector blocki(blocki_); 
  arma::vec theta = as<vec>(theta_);
  const NumericVector theta0(theta0_);
  const double theta1a = theta0[0];
  const double theta1b = theta0[1];
  const double theta2a = theta0[2];
  const double theta2b = theta0[3];
  const int spl0 = as<int>(spl0_);
  const arma::mat spS0 = as<mat>(spS0_); // 2 by 2
  const double spadapter = as<double>(spadapter_);
  
  // prameters to be updated
  arma::mat beta = as<mat>(beta_); //p by N
  arma::vec mu = as<vec>(mu_); // p by 1
  arma::mat Sig = as<mat>(Sig_); // p by p
  arma::vec w = as<vec>(w_);
  NumericVector tau2(tau2_);
  NumericVector tau = Rcpp::sqrt(tau2);
  IntegerVector K(K_);
  NumericVector V(V_);
  double alpha = as<double>(alpha_);

  // hyperparameters 
  const arma::vec m0 = as<vec>(m0_);
  const arma::mat S0 = as<mat>(S0_);
  const arma::mat Sig0 = as<mat>(Sig0_);
  const int k0 = as<int>(k0_);
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  const double nua = as<double>(nua_);
  const double nub = as<double>(nub_);
  const int nscan = nburn + (nskip+1)*nsave;
  arma::mat xpred = as<mat>(xpred_); 
  int npred = xpred.n_rows;
  
  // Temp variables
  double MinRes = Rcpp::min(yobs)-3.0;
  double MaxRes = Rcpp::max(yobs)+3.0;
  arma::mat Xbeta = X.t()*beta;
  NumericVector y(n);  for (int i=0; i<n; ++i) y[i] = yobs[i];
  IntegerVector nK(N);
  IntegerMatrix Kind(n, N);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  arma::mat Linv(n, nsave);
  // spatial related
  arma::mat spSnew(2,2); spSnew.fill(0.0);
  arma::vec thetabarnew(2); thetabarnew.fill(0.0);
  arma::vec z(n);
  arma::mat zPhi(n, N);
  arma::mat Cinv = arma::eye(n,n);
  double logdetC=0;
  NumericVector rejy(n);
  NumericVector rejbeta(N);
  NumericVector rejsigma(N);
  NumericVector rejV(N-1);
  int rejtheta=0;

  // Make arma objects
  arma::mat invSig = inv_sympd(Sig);
  const arma::mat invS0 = arma::inv_sympd(S0);
  const arma::vec invS0m0 = arma::solve(S0,m0);
  
  // things to save;
  NumericVector betaArray(nsave*N*p);
  arma::cube beta_save(betaArray.begin(), p, N, nsave, false);
  arma::mat w_save(N, nsave);
  NumericMatrix sigma2_save(N, nsave);
  NumericVector alpha_save(nsave);
  NumericMatrix y_save(n, nsave);
  NumericMatrix Ypred(npred, nsave);
  NumericVector theta1_save(nsave);
  NumericVector theta2_save(nsave);
  arma::mat z_save(n, nsave);
  NumericMatrix Zpred(npred, nsave);
  
  GetRNGstate();
	
	// Set the Armadillo seed from R's 
	// int seed = (int)Rf_runif(0.0, 10000.0);
	// std::srand(seed);

  // From V to w
  DDP_Vtow(w, V, N);

  // initial period
  spCopulaInitial( K, y, beta, tau2, V, w, alpha, mu, Sig, invSig, yobs, delta, X, m0, S0, Sig0, k0, 
                   a0, b0, nua, nub, invS0, invS0m0);
  Xbeta = X.t()*beta;
  tau = Rcpp::sqrt(tau2);
  // get transformed survival time z
  for (int i=0; i<n; ++i){
    for (int k=0; k<N; ++k){
      zPhi(i,k) = Rf_pnorm5((y[i]-Xbeta(i,k))*tau[k], 0, 1.0, true, false);
    }
  }
  z = qnormvec( zPhi*w );
  GetCinv_FSA(n, theta[0], theta[1], dnn, dnm, dmm, blocki, Cinv, logdetC); // Preprocess C^{-1} to get Cinv

  //////////////////////////////////////////////////////////////
  // Start MCMC
  //////////////////////////////////////////////////////////////
  
  for (int iscan=0; iscan<nscan; ++iscan){
    R_CheckUserInterrupt();
    // Rprintf( "scan = %d\n", iscan );

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
    spCopula_sample_y(y, rejy, zPhi, z, w, yobs, delta, Xbeta, tau, K, Cinv, n, N);
  
    // Sample beta;
    spCopula_sample_beta(beta, rejbeta, zPhi, z, w, y, X, tau2, nK, Kind, mu, Sig, invSig, Cinv, n, N, p);
    Xbeta = X.t()*beta;
  
    //Sample simga2;
    spCopula_sample_sigma2(tau2, rejsigma, zPhi, z, w, y, Xbeta, nK, Kind, nua, nub, Cinv, n, N);
    tau = Rcpp::sqrt(tau2);
  
    //Sample V;
    spCopula_sample_V(V, rejV, zPhi, z, w, nK, alpha, Cinv, N);
    DDP_Vtow(w, V, N); // From V to w
  
    //Sample alpha;
    double a0star = a0+N-1;
    double b0star = b0-std::log(w[N-1]);
    alpha = Rf_rgamma(a0star, 1/b0star);
  
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
  
    //Sample theta1 and theta2;
    spCopula_sample_theta_FSA(theta, rejtheta, spSnew, thetabarnew, Cinv, logdetC, 
                 theta1a, theta1b, theta2a, theta2b, spl0, spS0, dnn, spadapter, iscan, z, n, dnm, dmm, blocki);
  
    // Save the sample
    if (iscan>=nburn) {
      ++skiptally;
      if (skiptally>nskip){
        // calculate Linv
        arma::vec tautmp = as<vec>(tau);
        arma::vec sig = 1.0/tautmp;
        Linv.col(isave) = spCopula_Linv(yobs, delta, X, beta, sig, w, Cinv, z);
        // save samples
        alpha_save[isave] = alpha;
        y_save(_,isave) = y;
        for (int k=0; k<N; ++k){
          (beta_save.slice(isave)).col(k) = beta.col(k);
          sigma2_save(k, isave) = 1.0/tau2[k];
          w_save(k, isave) = w[k];
        }
        theta1_save[isave] = theta[0];
        theta2_save[isave] = theta[1];
        z_save.col(isave) = z;
        // prediction
        arma::mat xbeta = arma::trans( xpred*beta );
        arma::mat rho1 = arma::exp(-theta[1]*dnm)*arma::solve(arma::exp(-theta[1]*dmm), arma::exp(-theta[1]*ds0m));
        arma::mat rho2 = (arma::exp(-theta[1]*ds0n)-rho1)%ds0block;
        arma::mat rho = rho1 + rho2;
        for(int j=0; j<npred; ++j){
          // arma::vec hs0 = theta[0]*arma::exp(-theta[1]*ds0n.col(j));
          arma::vec hs0 = theta[0]*rho.col(j);
          double mus0 = arma::dot( hs0, Cinv*z );
          double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
          double znew = Rf_rnorm(mus0, sigs0);
          double u = Rf_pnorm5(znew, 0, 1.0, true, false);
          Zpred(j, isave) = znew;
          Ypred(j, isave) = DDP_Finvofu(u, w, xbeta.col(j), sig, MinRes, MaxRes);
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
  NumericVector ratey = 1.0- rejy/(nscan+0.0);
  NumericVector ratebeta = 1.0-rejbeta/(nscan+0.0);
  NumericVector ratesigma = 1.0-rejsigma/(nscan+0.0);
  NumericVector rateV = 1.0-rejV/(nscan+0.0);
  double ratetheta = 1.0-(double)rejtheta/(nscan+0.0);
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("beta")=beta_save,
                      Named("sigma2")=sigma2_save,
                      Named("w")=w_save,
                      Named("alpha")=alpha_save,
                      Named("theta1")=theta1_save,
                      Named("theta2")=theta2_save,
                      Named("z") = z_save,
                      Named("y")=y_save,
                      Named("ratey") = ratey, 
                      Named("ratebeta")=ratebeta,
                      Named("ratesigma")=ratesigma,
                      Named("rateV")=rateV,
                      Named("ratetheta")=ratetheta,
                      Named("cpo")=cpo,
                      Named("Ypred")=Ypred,
                      Named("Zpred")=Zpred);
	PutRNGstate();
	END_RCPP
}

