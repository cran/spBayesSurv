#include "spSurv_anovaDDP.h"
#include "spSurv_DDP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

SEXP anovaDDP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
		SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
		SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
		SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
    SEXP nua_, SEXP nub_, SEXP xpred_ ) {
	BEGIN_RCPP
  
	// Transfer R variables into C++;
	const int nburn = as<int>(nburn_);
	const int nsave = as<int>(nsave_);
	const int nskip = as<int>(nskip_);
	const int ndisplay = as<int>(ndisplay_);
	const int N = as<int>(N_);
	const NumericVector yobs(y_);
	IntegerVector delta(delta_);
  const arma::mat X = as<mat>(X_); // p by n
  const int p = X.n_rows;
  const int n = X.n_cols;
  
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
  
  // temp variables
  double MinRes = Rcpp::min(yobs)-3.0;
  double MaxRes = Rcpp::max(yobs)+3.0;
  arma::mat Xbeta = X.t()*beta;
  NumericVector y(n);	for (int i=0; i<n; ++i) y[i] = yobs[i];
	IntegerVector nK(N);
	IntegerMatrix Kind(n, N);
  int skiptally=0; 
  int isave=0;
	int distally=0;
  arma::mat Linv(n, nsave);
	
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
  arma::mat V_save(N, nsave);
  IntegerMatrix K_save(n, nsave);
	
	RNGScope scope;
  
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);
	
	// From V to w
  DDP_Vtow(w, V, N);
  
	////////////
	// Start MCMC
	////////////
	for (int iscan=0; iscan<nscan; iscan++){
	R_CheckUserInterrupt();
      
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
    for (int i=0; i<n; ++i){
      if(delta[i]==0){
        y[i] = trun_rnorm(Xbeta(i, K[i]-1), 1.0/tau[K[i]-1], yobs[i], R_PosInf);
      }
    }
    
    //Sample V;
    vec nkk = as<vec>(nK);
    for (int k=0; k<(N-1); ++k){
      double alphak = alpha + arma::sum(nkk.subvec(k+1, N-1))+ESMALL;
      double aa = nK[k] + 1.0;
      V[k] = Rf_rbeta(aa, alphak);
    }
    // From V to w
    DDP_Vtow(w, V, N);

    // Sample beta;
    anovaDDP_sample_beta(beta, y, X, tau2, nK, Kind, mu, Sig, invSig, N, p);
    Xbeta = X.t()*beta;

    //Sample simga2;
    anovaDDP_sample_sigma2(tau2, y, Xbeta, nK, Kind, nua, nub, N);
    tau = Rcpp::sqrt(tau2);

    //Sample alpha;
    double a0star = a0+N-1;
    double b0star = b0-log(w[N-1]);
    //Rprintf( "alpha = %f\n", alpha );
    //Rprintf( "b0star = %f\n", b0star );
    if(b0star>(b0+740.0)){
      Rprintf( "b0star = %f\n", b0star );
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
        Linv.col(isave) = anovaDDP_Linv(yobs, delta, X, beta, sig, w);
        // save samples
        alpha_save[isave] = alpha;
        y_save(_,isave) = y;
        K_save(_,isave) = K;
        for (int k=0; k<N; ++k){
          (beta_save.slice(isave)).col(k) = beta.col(k);
          sigma2_save(k,isave) = 1.0/tau2[k];
          w_save(k,isave) = w[k];
          V_save(k,isave) = V[k];
        }
        // prediction
        arma::mat xbeta = arma::trans( xpred*beta );
        for(int j=0; j<npred; ++j){
          double u = unif_rand();
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
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  // save current values
  List state; 
  state["K"] = K; state["y"] = y; state["V"] = V;
  state["w"] = w; state["beta"] = beta; state["sigma2"] = 1.0/tau2;
  state["alpha"] = alpha; state["mu"] = mu; state["Sig"] = Sig;
  
  return List::create(Named("beta")=beta_save,
                      Named("sigma2")=sigma2_save,
                      Named("w")=w_save,
                      Named("alpha")=alpha_save,
                      Named("y")=y_save,
                      Named("cpo")=cpo,
                      Named("Ypred")=Ypred,
                      Named("V")=V_save,
                      Named("K")=K_save,
                      Named("state")=state);
	END_RCPP
	}



