#include "spSurv_LDTFP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

RcppExport SEXP frailtyLDTFP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
		SEXP tobs_, SEXP type_, SEXP xce_, SEXP xtf_, SEXP alpha_, SEXP betace_, 
    SEXP betatf_, SEXP sigma2_, SEXP y_, SEXP v_, SEXP blocki_, SEXP tau2_, SEXP W_, SEXP D_, SEXP maxL_,
    SEXP a0_, SEXP b0_, SEXP m0_, SEXP S0inv_, SEXP gprior_, SEXP a0sig_, SEXP b0sig_, 
    SEXP a0tau_, SEXP b0tau_) {
	BEGIN_RCPP
  
	// Transfer R variables into C++;
	const int nburn = as<int>(nburn_);
	const int nsave = as<int>(nsave_);
	const int nskip = as<int>(nskip_);
	const int ndisplay = as<int>(ndisplay_);
  const NumericMatrix tobs(tobs_); // nrec by 2;
  const IntegerVector type(type_);// nrec by 1;
  const arma::mat xce = as<mat>(xce_); // pce by nrec;
  const arma::mat xtf = as<mat>(xtf_); // ptf by nrec;
  const int pce = xce.n_rows;
  const int ptf = xtf.n_rows;
  const int nrec = xce.n_cols;
  
	// prameters to be updated
  double alpha = as<double>(alpha_);
  arma::vec betace = as<vec>(betace_); //pce by 1
  arma::mat betatf = as<mat>(betatf_); // ptf by ntlr;
  double sigma2 = as<double>(sigma2_); 
  Rcpp::NumericVector y(y_); // nrec by 1;
  arma::vec v = as<vec>(v_); // m by 1;
  Rcpp::IntegerVector blocki(blocki_); // m+1 by 1;
  double tau2 = as<double>(tau2_);
  const arma::mat W = as<mat>(W_); // m by m;
  const arma::vec D = as<vec>(D_); // m by 1;
  const int m = v.size();  
  int frailty=1;
  if (arma::sum(D)==0.0) frailty=2;
  
	// hyperparameters
  const int maxL = as<int>(maxL_);
  // prior of alpha
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  // prior of betace
	const arma::vec m0 = as<vec>(m0_); // pce by 1
	const arma::mat S0inv = as<mat>(S0inv_); // pce by pce
  // prior of betatf
  const arma::mat gprior = as<mat>(gprior_); // ptf by ptf;
  // prior of sigma2
  const double a0sig = as<double>(a0sig_);
  const double b0sig = as<double>(b0sig_);
  // prior of tau2^{-1}
  const double a0tau = as<double>(a0tau_);
  const double b0tau = as<double>(b0tau_);
	
  // temp variables
  const int maxL1 = maxL+1;
	const int nscan = nburn + (nskip+1)*nsave;
  const int ntprob = std::pow(2,maxL1)-2;
  const int ntlr = ntprob/2;
  v = v - arma::mean(v);
  arma::vec vn(nrec);
  for(int i=0; i<m; ++i){
    int ind1 = blocki[i];
    int ind2 = blocki[i+1]-1;
    (vn.subvec(ind1, ind2)).fill(v[i]);
  }
  arma::vec xbetace = xce.t()*betace; // nrec by 1
  arma::mat xbetatf = betatf.t()*xtf; // ntlr by nrec
  Rcpp::IntegerVector nobsbc(ntprob);
  Rcpp::IntegerMatrix obsbc(ntprob,nrec);
  arma::mat Linv(nrec, nsave);
  int skiptally=0; 
  int isave=0;
	int distally=0;
	
	// things to save;
	NumericVector betaArray(nsave*ntlr*ptf);
	arma::cube betatf_save(betaArray.begin(), ptf, ntlr, nsave, false);
  arma::mat beta_save(pce, nsave);
	NumericVector sigma2_save(nsave);
	NumericVector alpha_save(nsave);
	NumericMatrix y_save(nrec, nsave);
  arma::mat v_save(m, nsave);
  NumericVector tau2_save(nsave);
  NumericVector rejbetatf(ntlr);
  double rejbetace=0;
	
	RNGScope scope;
  
	// Set the Armadillo seed from R's 
	//int seed = (int)Rf_runif(0.0, 10000.0);
	//std::srand(seed);
	
  // Working temp variables
  double logliko;
  double liminf, limsup, llim, rlim, gllim, grlim;
  double xx0, xx1, gxx0, gxx1;
  double logy, uwork, win;
  int i1, evali;
  int JJ, KK, mm;
  double rej;
  
	////////////////////////////////////////////////////////////////////////
  // Initial State
	////////////////////////////////////////////////////////////////////////
  loglikldtfp(y, xbetace+vn, xbetatf, sigma2, nobsbc, obsbc, logliko, maxL);
  i1=1;
  for(int i=2; i<maxL1; ++i){
    int j1 = std::pow(2,i-1);
    for(int j=1; j<=j1; ++j){
      int k1=i1+j;
      arma::mat c0=gprior/(alpha*(std::pow(i,2)+0.0));
      int ii=(k1-1)*2+1;
      int jj=(k1-1)*2+2;
      int n1=nobsbc[ii-1];
      int n2=nobsbc[jj-1];
      if((n1>0)&(n2>0)){
        startlrcoefldtfp(5, k1-1, ii-1, jj-1, n1, n2, obsbc, betatf, xtf, c0);
      }
    }
    i1 += j1;
  }
  xbetatf = betatf.t()*xtf;
  loglikldtfp(y, xbetace+vn, xbetatf, sigma2, nobsbc, obsbc, logliko, maxL);
  
	////////////////////////////////////////////////////////////////////////
	// Start MCMC
	////////////////////////////////////////////////////////////////////////
  win=1; mm=10;
	for (int iscan=0; iscan<nscan; iscan++){
  	R_CheckUserInterrupt();
    //Rprintf( "iscan = %d\n", iscan );
    
    ///////////////////////////////////////////////
    // imputing censored data
    //////////////////////////////////////////////
    for(int i=0; i<nrec; ++i){
      if(type[i]!=1) {
        if(type[i]==2) {
          liminf=-999.0;
          limsup=std::log(tobs(i,1));
        } else if(type[i]==0) {
          liminf=std::log(tobs(i,0)); 
          limsup=999.0;
        } else {
          liminf=std::log(tobs(i,0)); 
          limsup=std::log(tobs(i,1));
        }
        double xbetavi = xbetace[i] + vn[i];
        arma::vec xbetatfi = xbetatf.col(i);
        evali=1;
        xx0=y[i];
        ldensldtfp(xx0, xbetavi, xbetatfi, sigma2, gxx0, maxL);
        logy = gxx0-exp_rand();
        uwork=unif_rand();
        llim=xx0-win*uwork;
        rlim=llim+win;
        uwork=unif_rand();
        JJ = (int)(mm*uwork);
        KK = (mm-1)-JJ;
        // repeat 
        if(llim<liminf) llim=liminf;
        if(rlim>limsup) rlim=limsup;
        ++evali;
        ldensldtfp(llim, xbetavi, xbetatfi, sigma2, gllim, maxL);
        ++evali;
        ldensldtfp(rlim, xbetavi, xbetatfi, sigma2, grlim, maxL);
        while( (JJ>0)&(gllim>logy) ){
          llim -= win;
          JJ -= 1;
          if(llim<liminf){
            llim=liminf;
            gllim=logy-1.0;
          } else {
            ++evali;
            ldensldtfp(llim, xbetavi, xbetatfi, sigma2, gllim, maxL);
          }
        }
        while( (KK>0)&(grlim>logy) ){
          rlim += win;
          KK -= 1;
          if(rlim>limsup){
            rlim=limsup;
            grlim=logy-1.0;
          } else {
            ++evali;
            ldensldtfp(rlim, xbetavi, xbetatfi, sigma2, grlim, maxL);
          }
        }
        xx1 = llim +(rlim-llim)*unif_rand();
        ++evali;
        ldensldtfp(xx1, xbetavi, xbetatfi, sigma2, gxx1, maxL);
        while(gxx1<logy){
          if(xx1>xx0) rlim=xx1;
          if(xx1<xx0) llim=xx1;
          if(llim<liminf) llim=liminf;
          if(rlim>limsup) rlim=limsup;
          xx1 = llim +(rlim-llim)*unif_rand();
          ++evali;
          ldensldtfp(xx1, xbetavi, xbetatfi, sigma2, gxx1, maxL);
        }
        y[i] = xx1;
      }
    }
    
    ///////////////////////////////////////////////
    // frailties v
    //////////////////////////////////////////////
    if(frailty==1){
      int indm1 = blocki[m-1];
      int indm2 = blocki[m]-1;
      double vm = 0;
      for(int i=0; i<m-1; ++i){
        double wistar = D[i] + D[m-1] + 2.0*W(m-1, i);
        double wi1 = arma::as_scalar(W.row(i)*v)-v[m-1]*W(i,m-1);
        double wi2 = (D[m-1]+W(m-1,i))*(arma::sum(v)-v[i]-v[m-1]);
        double wi3 = arma::as_scalar(W.row(m-1)*v)-v[i]*W(m-1,i);
        double meanvi = ( wi1-wi2-wi3 )/wistar;
        double sdvi = tau2/wistar;
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        evali=1; 
        xx0 = v[i]; vm = v[m-1];
        loglikldtfpvi(xx0, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, gxx0, maxL, vm, indm1, indm2);
        logy = gxx0-exp_rand();
        uwork=unif_rand();
        llim=xx0-win*uwork;
        rlim=llim+win;
        uwork=unif_rand();
        JJ = (int)(mm*uwork);
        KK = (mm-1)-JJ;
        ++evali; vm = v[m-1] + (xx0-llim);
        loglikldtfpvi(llim, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, gllim, maxL, vm, indm1, indm2);
        ++evali; vm = v[m-1] + (xx0-rlim);
        loglikldtfpvi(rlim, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, grlim, maxL, vm, indm1, indm2);
        while((JJ>0)&(gllim>logy)){
          llim -= win;
          JJ -= 1;
          ++evali; vm = v[m-1] + (xx0-llim);
          loglikldtfpvi(llim, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, gllim, maxL, vm, indm1, indm2);
        }
        while((KK>0)&(grlim>logy)){
          rlim += win;
          KK -= 1;
          ++evali; vm = v[m-1] + (xx0-rlim);
          loglikldtfpvi(rlim, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, grlim, maxL, vm, indm1, indm2);
        }
        xx1 = llim +(rlim-llim)*unif_rand();
        ++evali; vm = v[m-1] + (xx0-xx1);
        loglikldtfpvi(xx1, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, gxx1, maxL, vm, indm1, indm2);
        while(gxx1<logy){
          if(xx1>xx0) rlim=xx1;
          if(xx1<xx0) llim=xx1;
          xx1 = llim +(rlim-llim)*unif_rand();
          ++evali; vm = v[m-1] + (xx0-xx1);
          loglikldtfpvi(xx1, meanvi, sdvi, ind1, ind2, y, xbetace, xbetatf, sigma2, gxx1, maxL, vm, indm1, indm2);
        }
        v[i]=xx1; v[m-1]=vm;
      }
      // transfter from v to vn
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        (vn.subvec(ind1, ind2)).fill(v[i]);
      }
      // tau2
      double a0taustar = a0tau+0.5*(m-1.0);
      double b0taustar = b0tau+0.5*(arma::dot(v,D%v)-arma::dot(v,W*v)); 
      tau2 = 1.0/Rf_rgamma( a0taustar, 1.0/b0taustar ); 
    }else{
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        evali=1; 
        xx0 = v[i];
        loglikldtfpvi2(xx0, ind1, ind2, y, xbetace, xbetatf, sigma2, gxx0, maxL);
        gxx0 += -0.5*std::pow(xx0, 2)/tau2;
        logy = gxx0-exp_rand();
        uwork=unif_rand();
        llim=xx0-win*uwork;
        rlim=llim+win;
        uwork=unif_rand();
        JJ = (int)(mm*uwork);
        KK = (mm-1)-JJ;
        ++evali;
        loglikldtfpvi2(llim, ind1, ind2, y, xbetace, xbetatf, sigma2, gllim, maxL);
        gllim += -0.5*std::pow(llim, 2)/tau2;
        ++evali; 
        loglikldtfpvi2(rlim, ind1, ind2, y, xbetace, xbetatf, sigma2, grlim, maxL);
        grlim += -0.5*std::pow(rlim, 2)/tau2;
        while((JJ>0)&(gllim>logy)){
          llim -= win;
          JJ -= 1;
          ++evali; 
          loglikldtfpvi2(llim, ind1, ind2, y, xbetace, xbetatf, sigma2, gllim, maxL);
          gllim += -0.5*std::pow(llim, 2)/tau2;
        }
        while((KK>0)&(grlim>logy)){
          rlim += win;
          KK -= 1;
          ++evali; 
          loglikldtfpvi2(rlim, ind1, ind2, y, xbetace, xbetatf, sigma2, grlim, maxL);
          grlim += -0.5*std::pow(rlim, 2)/tau2;
        }
        xx1 = llim +(rlim-llim)*unif_rand();
        ++evali; 
        loglikldtfpvi2(xx1, ind1, ind2, y, xbetace, xbetatf, sigma2, gxx1, maxL);
        gxx1 += -0.5*std::pow(xx1, 2)/tau2;
        while(gxx1<logy){
          if(xx1>xx0) rlim=xx1;
          if(xx1<xx0) llim=xx1;
          xx1 = llim +(rlim-llim)*unif_rand();
          ++evali; 
          loglikldtfpvi2(xx1, ind1, ind2, y, xbetace, xbetatf, sigma2, gxx1, maxL);
          gxx1 += -0.5*std::pow(xx1, 2)/tau2;
        }
        v[i]=xx1;
      }
      // transfter from v to vn
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        (vn.subvec(ind1, ind2)).fill(v[i]);
      }
      // tau2
      double a0taustar = a0tau+0.5*(m+0.0);
      double b0taustar = b0tau+0.5*(arma::dot(v,v)); 
      tau2 = 1.0/Rf_rgamma( a0taustar, 1.0/b0taustar ); 
    }
    
    ///////////////////////////////////////////////
    // baseline reg coefficients
    //////////////////////////////////////////////
    for(int i=0; i<pce; ++i){
      evali=1;
      xx0=betace[i];
      logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, gxx0, maxL);
      logy = gxx0-exp_rand();
      uwork=unif_rand();
      llim=xx0-win*uwork;
      rlim=llim+win;
      uwork=unif_rand();
      JJ = (int)(mm*uwork);
      KK = (mm-1)-JJ;
      // repeat
      ++evali;
      betace(i)=llim;
      logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, gllim, maxL);
      ++evali;
      betace(i)=rlim;
      logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, grlim, maxL);
      while( (JJ>0)&(gllim>logy) ){
        llim -= win;
        JJ -= 1;
        ++evali;
        betace(i)=llim;
        logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, gllim, maxL);
      }
      while( (KK>0)&(grlim>logy) ){
        rlim += win;
        KK -= 1;
        ++evali;
        betace(i)=rlim;
        logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, grlim, maxL);
      }
      xx1 = llim +(rlim-llim)*unif_rand();
      ++evali;
      betace(i)=xx1;
      logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, gxx1, maxL);
      while(gxx1<logy){
        if(xx1>xx0) rlim=xx1;
        if(xx1<xx0) llim=xx1;
        xx1 = llim +(rlim-llim)*unif_rand();
        ++evali;
        betace(i)=xx1;
        logposldtfp(betace, betatf, y, xce, vn, xtf, sigma2, m0, S0inv, nobsbc, obsbc, gxx1, maxL);
      }
      betace(i)=xx1;
    }
    //if(iscan>=nburn) rejbetace += rej;
    xbetace = xce.t()*betace;
    
    ///////////////////////////////////////////////
    // baseline variance sigma2
    //////////////////////////////////////////////
    evali=1;
    xx0 = sigma2;
    loglikldtfpsig(y, xbetace+vn, xbetatf, xx0, nobsbc, obsbc, gxx0, maxL, a0sig, b0sig);
    logy = gxx0-exp_rand();
    uwork=unif_rand();
    llim=xx0-win*uwork;
    rlim=llim+win;
    uwork=unif_rand();
    JJ = (int)(mm*uwork);
    KK = (mm-1)-JJ;
    // repeat
    if(llim<1e-5) llim=1e-5;
    ++evali;
    loglikldtfpsig(y, xbetace+vn, xbetatf, llim, nobsbc, obsbc, gllim, maxL, a0sig, b0sig);
    ++evali;
    loglikldtfpsig(y, xbetace+vn, xbetatf, rlim, nobsbc, obsbc, grlim, maxL, a0sig, b0sig);
    while( (JJ>0)&(gllim>logy) ){
      llim -= win;
      JJ -= 1;
      if(llim<1e-5) {
        llim=1e-5;
        gllim=logy-1.0;
      } else{
        ++evali;
        loglikldtfpsig(y, xbetace+vn, xbetatf, llim, nobsbc, obsbc, gllim, maxL, a0sig, b0sig);
      }
    }
    while( (KK>0)&(grlim>logy) ){
      rlim += win;
      KK -= 1;
      ++evali;
      loglikldtfpsig(y, xbetace+vn, xbetatf, rlim, nobsbc, obsbc, grlim, maxL, a0sig, b0sig);
    }
    xx1 = llim +(rlim-llim)*unif_rand();
    ++evali;
    loglikldtfpsig(y, xbetace+vn, xbetatf, xx1, nobsbc, obsbc, gxx1, maxL, a0sig, b0sig);
    while(gxx1<logy){
      if(xx1>xx0) rlim=xx1;
      if(xx1<xx0) llim=xx1;
      if(llim<1e-5) llim=1e-5;
      xx1 = llim +(rlim-llim)*unif_rand();
      ++evali;
      loglikldtfpsig(y, xbetace+vn, xbetatf, xx1, nobsbc, obsbc, gxx1, maxL, a0sig, b0sig);
    }
    sigma2 = xx1;
    
    ///////////////////////////////////////////////
    // tf logistic regressions
    //////////////////////////////////////////////
    i1=1;
    for(int i=2; i<=maxL; ++i){
      int j1=std::pow(2,i-1);
      for(int j=1; j<=j1; ++j){
        int k1=i1+j;
        arma::mat c0=gprior/(alpha*(std::pow(i,2)+0.0));
        arma::mat invc0 = arma::inv_sympd(c0);
        int ii=(k1-1)*2+1;
        int jj=(k1-1)*2+2;
        int n1=nobsbc[ii-1];
        int n2=nobsbc[jj-1];
        int ntot = n1+n2;
        rej=0;
        arma::vec betatfkk = betatf.col(k1-1);
        if(ntot>0){
          update_tfcoeff_iwls(ii-1, jj-1, n1, n2, obsbc, xtf, invc0, betatfkk, rej);
        }else{
          arma::vec mu0(ptf); mu0.fill(0.0);
          betatfkk = mvrnorm(mu0, c0);
        }
        betatf.col(k1-1) = betatfkk;
        if(iscan>=nburn) rejbetatf[k1-1] += rej;
      }
      i1 += j1;
    }
    xbetatf = betatf.t()*xtf;
    
    ///////////////////////////////////////////////
    // alpha
    //////////////////////////////////////////////
    if(a0>0){
      arma::mat xtx = arma::inv_sympd(gprior);
      arma::mat c0(ptf, ptf);
      i1=1;
      double tmp1=0.0;
      for(int i=2; i<=maxL; ++i){
        int j1=std::pow(2,i-1);
        for(int j=1; j<=j1; ++j){
          int k1=i1+j;
          c0 = xtx*(std::pow(i,2)+0.0);
          tmp1 += arma::dot(betatf.col(k1-1), c0*betatf.col(k1-1));
        }
        i1 += j1;
      }
      alpha = Rf_rgamma(a0+(ptf+0.0)*(ntlr-1.0)*0.5, 1.0/(b0+0.5*tmp1));
    }
    
    ///////////////////////////////////////////////
    // Save the sample
    //////////////////////////////////////////////
    if(iscan>=nburn){
  	  ++skiptally;
  	  if(skiptally>nskip){
        // calculate Linv
        Linv.col(isave) = spldtfp_Linv(tobs, type, xbetace+vn, xbetatf, sigma2, maxL);
    		// save regression coefficient
        beta_save.col(isave) = betace;
        // centering variance
        sigma2_save[isave] = sigma2;
        // precision parameter
    		alpha_save[isave] = alpha;
        // imputed log survival times
    		y_save(_,isave) = y;
        // frailties
        v_save.col(isave) = v;
        tau2_save[isave] = tau2;
        // tf parameters
        betatf_save.slice(isave) = betatf;
    				
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
  
  NumericVector ratebetatf = 1.0- rejbetatf/(nscan-nburn+0.0);
  double ratebetace = 1.0- rejbetace/(nscan-nburn+0.0);
  
  // get CPO
  arma::vec Linvmean = arma::mean(Linv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  return List::create(Named("beta")=beta_save,
                      Named("sigma2")=sigma2_save,
                      Named("alpha")=alpha_save,
                      Named("betatf")=betatf_save,
                      Named("y")=y_save,
                      Named("v")=v_save,
                      Named("tau2")=tau2_save,
                      Named("cpo")=cpo,
                      Named("ratebetatf")=ratebetatf,
                      Named("ratebetace")=ratebetace);
	END_RCPP
}



