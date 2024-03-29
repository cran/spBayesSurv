#include "spSurv_DDP_tools.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////
///////////////// spatial Copula DDP using FSA /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
RcppExport SEXP spCopulaDDP(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
                            SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
                            SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
                            SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
                            SEXP nua_, SEXP nub_, SEXP xpred_, SEXP ds0n_, SEXP dnn_, 
                            SEXP theta_, SEXP theta0_, SEXP spS0_, SEXP l0_, SEXP adapter_, 
                            SEXP dnm_, SEXP dmm_, SEXP clustindx_, SEXP ds0m_, SEXP ds0block_){
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
  const int l0 = as<int>(l0_);
  const double adapter = as<double>(adapter_);
  
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
  double rejtheta=0;

  // Make arma objects
  arma::mat invSig = inv_sympd(Sig);
  const arma::mat invS0 = arma::inv_sympd(S0);
  const arma::vec invS0m0 = arma::solve(S0,m0);
  
  // things to save;
  NumericVector betaArray(nsave*N*p);
  arma::cube beta_save(betaArray.begin(), p, N, nsave, false, true);
  arma::mat w_save(N, nsave);
  NumericMatrix sigma2_save(N, nsave);
  NumericVector alpha_save(nsave);
  NumericMatrix y_save(n, nsave);
  NumericMatrix Ypred(npred, nsave);
  Rcpp::NumericMatrix theta_save(2, nsave);
  arma::mat z_save(n, nsave);
  NumericMatrix Zpred(npred, nsave);
  arma::mat V_save(N, nsave);
  IntegerMatrix K_save(n, nsave);
  
  RNGScope scope;
	
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
  arma::mat Cnn = arma::exp(-theta[1]*dnn); 
  arma::mat Cnm = arma::exp(-theta[1]*dnm);
  arma::mat Cmm = arma::exp(-theta[1]*dmm);
  inv_FSA(theta[0], Cnn, Cnm, Cmm, clustindx, Cinv, logdetC); // Preprocess C^{-1} to get Cinv

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
    
    //Sample theta1 and theta2;
    spCopula_sample_theta_FSA(theta, rejtheta, spSnew, thetabarnew, Cinv, logdetC, 
                              theta1a, theta1b, theta2a, theta2b, l0, spS0, dnn, 
                              adapter, iscan, nburn, z, n, dnm, dmm, clustindx);

    // Sample y_i when delta_i=0
    spCopula_sample_y(y, rejy, zPhi, z, w, yobs, delta, Xbeta, tau, K, Cinv, n, N, iscan, nburn);
    
    //Sample V;
    spCopula_sample_V(V, rejV, zPhi, z, w, nK, alpha, Cinv, N, iscan, nburn);
    //spCopula_sample_V_block(V, rejV, zPhi, z, w, nK, alpha, Cinv, N);
    DDP_Vtow(w, V, N); // From V to w
  
    // Sample beta;
    spCopula_sample_beta(beta, rejbeta, zPhi, z, w, y, X, tau2, nK, Kind, mu, Sig, invSig, Cinv, n, N, p, iscan, nburn);
    Xbeta = X.t()*beta;
  
    //Sample simga2;
    spCopula_sample_sigma2(tau2, rejsigma, zPhi, z, w, y, Xbeta, nK, Kind, nua, nub, Cinv, n, N, iscan, nburn);
    tau = Rcpp::sqrt(tau2);

    //Sample alpha;
    double a0star = a0+N-1;
    double b0star = b0-log(w[N-1]);
    if(b0star>(b0+740.0)){
      // Rprintf( "b0star = %f\n", b0star );
      b0star = b0+1e5; // b0star = b0+(N-1.0)/alpha;
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
        K_save(_,isave) = K;
        for (int k=0; k<N; ++k){
          (beta_save.slice(isave)).col(k) = beta.col(k);
          sigma2_save(k, isave) = 1.0/tau2[k];
          w_save(k, isave) = w[k];
          V_save(k,isave) = V[k];
        }
        theta_save(0,isave) = theta[0];
        theta_save(1,isave) = theta[1];
        z_save.col(isave) = z;
        // prediction
        arma::mat xbeta = arma::trans( xpred*beta );
        arma::mat rho(n, npred);
        if(((int)dmm.n_cols)==n){
          rho = arma::exp(-theta[1]*ds0n);
        }else{
          arma::mat rho1 = arma::exp(-theta[1]*dnm)*arma::solve(arma::exp(-theta[1]*dmm), arma::exp(-theta[1]*ds0m));
          arma::mat rho2 = (arma::exp(-theta[1]*ds0n)-rho1)%ds0block;
          rho = rho1 + rho2;
        }
        for(int j=0; j<npred; ++j){
          // arma::vec hs0 = theta[0]*arma::exp(-theta[1]*ds0n.col(j));
          arma::vec hs0 = theta[0]*rho.col(j);
          double mus0 = arma::dot( hs0, Cinv*z );
          double sigs0 = std::sqrt( 1.0 - arma::dot( hs0, Cinv*hs0) );
          double znew = Rf_rnorm(mus0, sigs0);
          double u = Rf_pnorm5(znew, 0, 1.0, true, false);
          Zpred(j, isave) = znew;
          Ypred(j, isave) = DDP_Finvofu(u, w, xbeta.col(j), sig, log(ESMALL), log(ELARGE));
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
  NumericVector ratey = 1.0- rejy/(nscan+0.0);
  NumericVector ratebeta = 1.0-rejbeta/totscans;
  NumericVector ratesigma = 1.0-rejsigma/totscans;
  NumericVector rateV = 1.0-rejV/totscans;
  double ratetheta = 1.0-rejtheta/totscans;
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  // save current values
  List state; 
  state["K"] = K; state["y"] = y; state["V"] = V;
  state["w"] = w; state["beta"] = beta; state["sigma2"] = 1.0/tau2;
  state["alpha"] = alpha; state["mu"] = mu; state["Sig"] = Sig;
  state["theta"] = theta; state["spSnew"] = spSnew;
  
  return List::create(Named("beta")=beta_save,
                      Named("sigma2")=sigma2_save,
                      Named("w")=w_save,
                      Named("alpha")=alpha_save,
                      Named("theta")=theta_save,
                      Named("z") = z_save,
                      Named("y")=y_save,
                      Named("ratey") = ratey, 
                      Named("ratebeta")=ratebeta,
                      Named("ratesigma")=ratesigma,
                      Named("rateV")=rateV,
                      Named("ratetheta")=ratetheta,
                      Named("cpo")=cpo,
                      Named("Ypred")=Ypred,
                      Named("Zpred")=Zpred,
                      Named("V")=V_save,
                      Named("K")=K_save,
                      Named("state")=state);
	END_RCPP
}

