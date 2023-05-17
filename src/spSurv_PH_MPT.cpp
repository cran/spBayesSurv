#include "spSurv_MPT_tools_single.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

RcppExport SEXP PH_MPT(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_, SEXP ltr_, SEXP subjecti_,
                       SEXP t1_, SEXP t2_, SEXP type_, SEXP X_, SEXP theta_, SEXP beta_, 
                       SEXP cpar_, SEXP Ys_, SEXP maxL_, SEXP a0_, SEXP b0_, SEXP theta0_, 
                       SEXP V0inv_,SEXP Vhat_, SEXP beta0_, SEXP S0inv_, SEXP Shat_, 
                       SEXP l0_, SEXP adapter_, SEXP gamma_, SEXP p0gamma_, SEXP selection_,
                       SEXP frailty_, SEXP v_, SEXP blocki_, SEXP W_, 
                       SEXP lambda_, SEXP a0lambda_, SEXP b0lambda_, SEXP dist_){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  const int ndisplay = as<int>(ndisplay_);
  const Rcpp::NumericVector ltr(ltr_); // n by 1;
  const Rcpp::IntegerVector subjecti(subjecti_); //nsub+1 by 1;
  const Rcpp::NumericVector t1(t1_); // n by 1;
  const Rcpp::NumericVector t2(t2_); // n by 1;
  const Rcpp::IntegerVector type(type_);// n by 1;
  const Rcpp::NumericMatrix X(X_); // n by p;
  const int nsub = subjecti.size()-1;
  const int n = type.size();
  const int maxL = as<int>(maxL_);
  const int p = X.ncol();
  const int l0 = as<int>(l0_);
  const double adapter = as<double>(adapter_);
  const int selection = as<int>(selection_);
  const int dist = as<int>(dist_);
  // frailty model related
  const int frailty=as<int>(frailty_);
  Rcpp::NumericVector v(v_); // m by 1;
  const Rcpp::IntegerVector blocki(blocki_); //m+1 by 1;
  const arma::mat W = as<mat>(W_); // m by m;
  const int m = v.size();
  const arma::vec D = arma::sum(W, 1); // m by 1;
  
  // prameters to be updated
  double cpar = as<double>(cpar_);
  const bool MPT = (cpar<R_PosInf);
  Rcpp::NumericVector theta(theta_); // 2 by 1
  Rcpp::NumericVector beta(beta_); // p by 1
  Rcpp::NumericVector Ys(Ys_); //  2**(maxL)-1 by 1
  const int nYs = Ys.size();
  Rcpp::NumericVector gamma(gamma_); // p by 1
  double lambda = as<double>(lambda_);
  
  // hyperparameters
  // prior of cpar
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  // prior of theta
  const arma::vec theta0 = as<arma::vec>(theta0_); // 2 by 1
  const arma::mat V0inv = as<arma::mat>(V0inv_);   // 2 by 2
  const arma::mat Vhat = as<arma::mat>(Vhat_);     // 2 by 2
  // prior of beta
  const arma::vec beta0 = as<arma::vec>(beta0_); // p by 1
  const arma::mat S0inv = as<arma::mat>(S0inv_); // p by p
  const arma::mat Shat = as<arma::mat>(Shat_);   // p by p
  // prior of gamma
  Rcpp::NumericVector p0gamma(p0gamma_); // p by 1
  // prior of lambda
  const double a0lambda = as<double>(a0lambda_);
  const double b0lambda = as<double>(b0lambda_);
  
  // temp variables
  const int nscan = nburn + (nskip+1)*nsave;
  const int nprob = nYs+1;
  Rcpp::NumericVector Xbeta(n, 0.0);
  Rcpp::NumericVector probs(nprob,0.0);
  Rcpp::NumericVector probsold(nprob,0.0);
  Ys_to_probs(Ys, probs, maxL);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  Rcpp::NumericVector vn(n,0.0);
  
  // things to save;
  Rcpp::NumericVector cpar_save(nsave);
  Rcpp::NumericMatrix theta_save(2, nsave);
  Rcpp::NumericMatrix beta_save(p, nsave);
  Rcpp::NumericMatrix Ys_save(nYs, nsave);
  Rcpp::NumericMatrix gamma_save(p, nsave);
  double rejtheta=0;
  double rejbeta=0;
  double rejYs=0;
  double rejc=0;
  Rcpp::NumericVector rejv(m, 0.0);
  Rcpp::NumericMatrix v_save(m, nsave);
  Rcpp::NumericVector lambda_save(nsave);
  
  // Make arma objects
  arma::mat X_r(const_cast<NumericMatrix&>(X).begin(), n, p, false, true);
  arma::vec theta_r(theta.begin(), 2, false, true);
  arma::vec beta_r(beta.begin(), p, false, true);
  arma::vec gamma_r(gamma.begin(), p, false, true);
  arma::vec Xbeta_r(Xbeta.begin(), n, false, true);
  arma::vec Ys_r(Ys.begin(), nYs, false, true);
  arma::vec v_r(v.begin(), m, false, true);
  arma::vec vn_r(vn.begin(), n, false, true);
  
  // Working temp variables
  arma::mat Linv=arma::zeros<arma::mat>(n, nsave);
  arma::vec Dvalues = arma::zeros<arma::vec>(nsave);
  Rcpp::NumericVector sumtheta(2, 0.0); // 2 by 1
  Rcpp::NumericVector sumbeta(p, 0.0); // p by 1
  Rcpp::NumericVector sumYs(nYs, 0.0); //  2**(maxL)-1 by 1
  Rcpp::NumericVector sumvn(n, 0.0); 
  // for theta
  arma::vec thetaold(2); 
  arma::vec thBarold(2);
  arma::vec thBarnew=arma::zeros<arma::vec>(2);
  arma::mat thSnew=arma::zeros<arma::mat>(2,2);
  arma::mat I2 = ESMALL*arma::eye(2,2);
  // for beta
  arma::vec betaold(p);
  arma::vec beBarold(p);
  arma::vec beBarnew=arma::zeros<arma::vec>(p);
  arma::mat beSnew=arma::zeros<arma::mat>(p,p);
  arma::mat Ip = ESMALL*arma::eye(p,p);
  // for Ys
  arma::vec YsBarold(nYs);
  arma::vec YsBarnew=arma::zeros<arma::vec>(nYs);
  arma::mat YsSnew=arma::zeros<arma::mat>(nYs,nYs);
  arma::mat InYs = ESMALL*arma::eye(nYs,nYs);
  arma::mat YsShat=arma::zeros<arma::mat>(nYs,nYs);
  for(int j=0; j<nYs; ++j) YsShat(j,j) = 0.01*std::pow((int)(log(j+1)/log(2))+1.0, 2);
  arma::vec tmpnew=arma::zeros<arma::vec>(nYs);
  arma::vec tmpold=arma::zeros<arma::vec>(nYs);
  // for cpa
  double cbarnew=0, cbarold=0, csnew=0, ctemp=0, cshat=0.16;
  double llold, llnew;
  double ratio, uu, nn, j2;
  // for v
  if(frailty==1)  {
    v_r = v_r - arma::mean(v_r);
  } else if(frailty==0){
    for(int i=0; i<m; ++i) v[i]=0;
  }
  for(int i=0; i<m; ++i){
    int ind1 = blocki[i];
    int ind2 = blocki[i+1]-1;
    (vn_r.subvec(ind1, ind2)).fill(v[i]);
  }
  
  RNGScope scope;
  
  // Set the Armadillo seed from R's 
  // int seed = (int)Rf_runif(0.0, 10000.0);
  // std::srand(seed);
  
  ////////////////////////////////////////////////////////////////////////
  // Start MCMC
  ////////////////////////////////////////////////////////////////////////
  Xbeta_r = X_r*(beta_r%gamma_r);
  for (int iscan=0; iscan<nscan; iscan++){
    R_CheckUserInterrupt();
    // Rprintf( "iscan = %d\n", iscan );
    
    ///////////////////////////////////////////////
    // update PT conditional probabilities
    //////////////////////////////////////////////
    if(MPT){
      PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llold);
      tmpold = arma::log(Ys_r/(1.0-Ys_r));
      for(int i=0; i<nprob; ++i) probsold[i] = probs[i];
      if(iscan>l0){
        tmpnew = mvrnorm(tmpold, YsSnew);
      }else{
        tmpnew = mvrnorm(tmpold, YsShat);
      }
      Ys_r = 1.0/(1.0+arma::exp(-tmpnew));
      Ys_to_probs(Ys, probs, maxL);
      PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llnew);
      for(int i=0; i<nYs; ++i){
        j2 = std::pow((int)(log(i+1)/log(2))+1, 2);
        llold += cpar*j2*tmpold[i]-2.0*cpar*j2*log(1+exp(tmpold[i]));
        llnew += cpar*j2*tmpnew[i]-2.0*cpar*j2*log(1+exp(tmpnew[i]));
      }
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if(uu>ratio){
        for(int i=0; i<nYs; ++i){
          tmpnew[i] = tmpold[i];
          Ys[i] = 1.0/(1.0+exp(-tmpold[i]));
        }
        if(iscan>=nburn) rejYs += 1.0;
        for(int i=0; i<nprob; ++i) probs[i] = probsold[i];
      }
      nn = iscan+1;
      YsBarold = YsBarnew;
      YsBarnew = (nn)/(nn+1.0)*YsBarold + tmpnew/(nn+1.0);
      YsSnew = (nn-1.0)/nn*YsSnew + adapter/(double)(nYs)/nn*(nn*YsBarold*YsBarold.t() 
                                                                - (nn+1.0)*YsBarnew*YsBarnew.t() + tmpnew*tmpnew.t() + InYs );
    }
    
    ///////////////////////////////////////////////
    // update beta
    //////////////////////////////////////////////
    PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llold);
    llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
    betaold = beta_r;
    if(iscan>l0){
      beta_r = mvrnorm(betaold, beSnew);
    }else{
      beta_r = mvrnorm(betaold, Shat);
    }
    Xbeta_r = X_r*(beta_r%gamma_r);
    PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llnew);
    llnew += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
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
    Xbeta_r = X_r*(beta_r%gamma_r);
    
    ///////////////////////////////////////////////
    // update theta
    //////////////////////////////////////////////
    if(V0inv(0,0)<SYSMAX){
      PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llold);
      llold += -0.5*arma::dot( (theta_r-theta0), (V0inv*(theta_r-theta0)) );
      thetaold = theta_r;
      if(iscan>l0){
        theta_r = mvrnorm(thetaold, thSnew);
      }else{
        theta_r = mvrnorm(thetaold, Vhat);
      }
      PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llnew);
      llnew += -0.5*arma::dot( (theta_r-theta0), (V0inv*(theta_r-theta0)) );
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if(uu>ratio){
        theta_r=thetaold;
        if(iscan>=nburn) rejtheta+=1.0;
      }
      nn = iscan+1;
      thBarold = thBarnew;
      thBarnew = (nn)/(nn+1.0)*thBarold + theta_r/(nn+1.0);
      thSnew = (nn-1.0)/nn*thSnew + adapter/(2.0)/nn*(nn*thBarold*thBarold.t() 
                                                        - (nn+1.0)*thBarnew*thBarnew.t() + theta_r*theta_r.t() + I2 );
    }
    
    ///////////////////////////////////////////////
    // update v and lambda
    //////////////////////////////////////////////
    if(frailty==1){
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        double meanvi = arma::as_scalar(W.row(i)*v_r)/D[i];
        double sdvi = std::sqrt(1.0/(D[i]*lambda));
        double viold = v[i];
        PHloglikblocki(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta, llold, ind1, ind2, v[i]);
        llold += -0.5*D[i]*lambda*std::pow(v[i]-meanvi,2);
        v[i] = Rf_rnorm(viold, sdvi);
        PHloglikblocki(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta, llnew, ind1, ind2, v[i]);
        llnew += -0.5*D[i]*lambda*std::pow(v[i]-meanvi,2);
        ratio = exp(llnew-llold);
        uu = unif_rand();
        if(uu>ratio){
          v[i]=viold;
          if(iscan>=nburn) rejv[i]+=1.0;
        }
        v_r = v_r - arma::mean(v_r);
      }
      // transfter from v to vn
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        (vn_r.subvec(ind1, ind2)).fill(v[i]);
      }
      // lambda
      double a0lambdastar = a0lambda+0.5*(m-1.0);
      double b0lambdastar = b0lambda+0.5*(arma::dot(v_r,D%v_r)-arma::dot(v_r,W*v_r)); 
      lambda = Rf_rgamma( a0lambdastar, 1.0/b0lambdastar ); 
    } else if(frailty==2){
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        double viold = v[i];
        PHloglikblocki(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta, llold, ind1, ind2, v[i]);
        llold += -0.5*lambda*std::pow(v[i],2);
        v[i] = Rf_rnorm(viold, std::sqrt(1.0/lambda));
        PHloglikblocki(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta, llnew, ind1, ind2, v[i]);
        llnew += -0.5*lambda*std::pow(v[i],2);
        ratio = exp(llnew-llold);
        uu = unif_rand();
        if(uu>ratio){
          v[i]=viold;
          if(iscan>=nburn) rejv[i]+=1.0;
        }
      }
      // transfter from v to vn
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        (vn_r.subvec(ind1, ind2)).fill(v[i]);
      }
      // lambda
      double a0lambdastar = a0lambda+0.5*(m+0.0);
      double b0lambdastar = b0lambda+0.5*(arma::dot(v_r,v_r)); 
      lambda = Rf_rgamma( a0lambdastar, 1.0/b0lambdastar );       
    }
    
    ///////////////////////////////////////////////
    // gamma
    //////////////////////////////////////////////
    if(selection==1){
      for(int j=0; j<p; ++j){
        arma::vec gamma_tmp = gamma_r;
        gamma_tmp[j] = 1.0;
        Xbeta_r = X_r*(beta_r%gamma_tmp);
        PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llold);
        gamma_tmp[j] = 0.0;
        Xbeta_r = X_r*(beta_r%gamma_tmp);
        PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, llnew);
        double pstar = p0gamma[j]/( p0gamma[j]+(1.0-p0gamma[j])*exp(llnew-llold) );
        uu = unif_rand();
        //Rprintf( "p = %f\n", pstar );
        if(uu<pstar){
          gamma_r[j]=1.0;
        }else{
          gamma_r[j]=0.0;
        }
        Xbeta_r = X_r*(beta_r%gamma_r);
      }
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
        llold=0; llnew=0;
        for(int i=0; i<nYs; ++i){
          j2 = std::pow((int)(log(i+1)/log(2))+1, 2);
          llold += Rf_lgammafn(2.0*cpar*j2) - 2.0*Rf_lgammafn(cpar*j2) + cpar*j2*(log(Ys[i])+log(1.0-Ys[i]));
          llnew += Rf_lgammafn(2.0*ctemp*j2) - 2.0*Rf_lgammafn(ctemp*j2) + ctemp*j2*(log(Ys[i])+log(1.0-Ys[i]));
        }
        llold += (a0-1.0)*log(cpar)-b0*cpar;
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
    
    ///////////////////////////////////////////////
    // Save the sample
    //////////////////////////////////////////////
    if(iscan>=nburn){
      ++skiptally;
      if(skiptally>nskip){
        // calculate 1.0/likelihood
        Linv.col(isave) = PHinvLik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn);
        // calculate -2loglikelihood
        double tempLik = 0;
        PHloglik(t1, t2, ltr, type, theta[0], theta[1], probs, maxL, MPT, dist, Xbeta+vn, tempLik);
        Dvalues(isave) = -2.0*tempLik;
        // save regression coefficient
        theta_save(_,isave) = theta;
        sumtheta += theta; 
        beta_save(_,isave) = beta;
        sumbeta += beta*gamma;
        gamma_save(_, isave) = gamma;
        // precision parameter
        cpar_save[isave] = cpar;
        // PT probs
        Ys_save(_,isave) = Ys;
        sumYs += Ys;
        // frailties
        v_save(_,isave) = v;
        sumvn += vn;
        lambda_save[isave] = lambda;
        
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
  double ratebeta = 1.0 - rejbeta/totscans;
  double rateYs = 1.0 - rejYs/totscans;
  double ratec = 1.0 - rejc/totscans;
  Rcpp::NumericVector ratev = 1.0 - rejv/totscans;
  
  // get CPO
  arma::mat finv(nsub, nsave);
  for(int i=0; i<nsub; ++i){
    int ind1 = subjecti[i];
    int ind2 = subjecti[i+1]-1;
    finv.row(i) = arma::prod(Linv.rows(ind1, ind2), 0);
  }
  arma::vec Linvmean = arma::mean(finv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  // get DIC
  double meanD = arma::mean(Dvalues);
  Xbeta_r = X_r*as<arma::vec>(sumbeta)/(nsave+0.0);
  Ys_to_probs(sumYs/(nsave+0.0), probs, maxL);
  double Dmean = 0;
  PHloglik(t1, t2, ltr, type, sumtheta[0]/(nsave+0.0), sumtheta[1]/(nsave+0.0), 
           probs, maxL, MPT, dist, Xbeta+sumvn/(nsave+0.0), Dmean);
  double pD = meanD + 2.0*Dmean;
  double DIC = meanD + pD; 
  
  return List::create(Named("theta")=theta_save,
                      Named("beta")=beta_save,
                      Named("gamma")=gamma_save,
                      Named("cpar")=cpar_save,
                      Named("Ys")=Ys_save,
                      Named("cpo")=cpo,
                      Named("pD")=pD,
                      Named("DIC")=DIC,
                      Named("ratetheta")=ratetheta,
                      Named("ratebeta")=ratebeta,
                      Named("rateYs")=rateYs,
                      Named("ratec")=ratec,
                      Named("v")=v_save,
                      Named("lambda")=lambda_save,
                      Named("ratev")=ratev);
  END_RCPP
}

// Get density or survival Plots for frailty LDTFP PH
RcppExport SEXP PH_MPT_plots(SEXP tgrid_, SEXP xpred_, SEXP theta_, SEXP beta_, SEXP Ys_, 
                             SEXP maxL_, SEXP CI_, SEXP dist_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  Rcpp::NumericVector tgrid(tgrid_);
  arma::mat xpred = as<arma::mat>(xpred_); // npred by p;
  Rcpp::NumericMatrix theta(theta_); // 2 by nsave;
  Rcpp::NumericMatrix beta(beta_); // p by nsave;
  Rcpp::NumericMatrix Ys(Ys_); // 2**maxL-1 by nsave;
  int maxL = as<int>(maxL_);
  double CI = as<double>(CI_);
  const int dist = as<int>(dist_);
  int nsave = beta.ncol();
  int ngrid = tgrid.size();
  int npred = xpred.n_rows;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  int nYs = Ys.nrow();
  int nprob = nYs+1;
  int p = beta.nrow();
  
  // Temp variables
  arma::vec xbeta(npred);
  Rcpp::NumericVector probs(nprob,0.0);
  Rcpp::NumericVector estfArray(nsave*ngrid*npred);
  arma::cube estf(estfArray.begin(), ngrid, nsave, npred, false, true);
  Rcpp::NumericVector estSArray(nsave*ngrid*npred);
  arma::cube estS(estSArray.begin(), ngrid, nsave, npred, false, true);
  Rcpp::NumericVector esthArray(nsave*ngrid*npred);
  arma::cube esth(esthArray.begin(), ngrid, nsave, npred, false, true);
  
  // Make arma objects
  arma::mat beta_r(beta.begin(), p, nsave, false, true);
  
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
    xbeta = xpred*beta_r.col(i);
    Rcpp::NumericVector Ysi = Ys(_,i);
    Ys_to_probs(Ysi, probs, maxL);
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        estf(k, i, j) = std::exp(PHlogpdf(tgrid[k], theta(0,i), theta(1,i), probs, maxL, true, dist,
                                 xbeta[j]));
        estS(k, i, j) = std::exp(PHlogsurv(tgrid[k], theta(0,i), theta(1,i), probs, maxL, true, dist,
                                 xbeta[j]));
        esth(k, i, j) = estf(k, i, j)/estS(k, i, j);
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
