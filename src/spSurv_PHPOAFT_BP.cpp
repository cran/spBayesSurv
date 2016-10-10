#include "spSurv_BP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std; 

RcppExport SEXP PHPOAFT_BP(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_, SEXP ltr_, SEXP subjecti_,
                           SEXP t1_, SEXP t2_, SEXP type_, SEXP X_, SEXP theta_, SEXP beta_, 
                           SEXP weight_, SEXP cpar_, SEXP a0_, SEXP b0_, SEXP theta0_, SEXP V0inv_, SEXP Vhat_, 
                           SEXP beta0_, SEXP S0inv_, SEXP Shat_, 
                           SEXP l0_, SEXP adapter_, SEXP dist_){
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
  const int l0 = as<int>(l0_);
  const double adapter = as<double>(adapter_);
  const int dist = as<int>(dist_);
  const int n = type.size();
  const int p = X.ncol();
  const int p3 = 3*p;
  // prameters to be updated and priors
  // prior of theta
  Rcpp::NumericVector theta(theta_); // 2 by 1
  const arma::vec theta0 = as<arma::vec>(theta0_); // 2 by 1
  const arma::mat V0inv = as<arma::mat>(V0inv_);   // 2 by 2
  const arma::mat Vhat = as<arma::mat>(Vhat_);     // 2 by 2
  // prior of beta
  Rcpp::NumericVector beta(beta_); // 3*p by 1
  const arma::vec beta0 = as<arma::vec>(beta0_); // 3*p by 1
  const arma::mat S0inv = as<arma::mat>(S0inv_); // 3*p by 3*p
  const arma::mat Shat = as<arma::mat>(Shat_);   // 3*p by 3*p
  // prior of cpar
  double cpar = as<double>(cpar_);
  const bool BP = (cpar<R_PosInf);
  Rcpp::NumericVector weight(weight_); //  maxL by 1
  const int maxL = weight.size();
  const double a0 = as<double>(a0_); 
  const double b0 = as<double>(b0_);
  const int nYs = maxL-1;
  // temp variables
  const int nscan = nburn + (nskip+1)*nsave;
  Rcpp::NumericVector beta_h(p, 0.0);
  Rcpp::NumericVector beta_o(p, 0.0);
  Rcpp::NumericVector beta_q(p, 0.0);
  Rcpp::NumericVector Xbeta_h(n, 0.0);
  Rcpp::NumericVector Xbeta_o(n, 0.0);
  Rcpp::NumericVector Xbeta_q(n, 0.0);
  Rcpp::NumericVector Ys(nYs);
  for(int i=0; i<nYs; ++i) Ys[i] = log(weight[i])-log(weight[nYs]);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  
  // things to save;
  Rcpp::NumericVector cpar_save(nsave);
  Rcpp::NumericMatrix theta_save(2, nsave);
  Rcpp::NumericMatrix beta_h_save(p, nsave);
  Rcpp::NumericMatrix beta_o_save(p, nsave);
  Rcpp::NumericMatrix beta_q_save(p, nsave);
  Rcpp::NumericMatrix weight_save(maxL, nsave);
  double rejtheta=0;
  double rejbeta=0;
  double rejYs=0;
  double rejc=0;
  
  // Make arma objects
  arma::mat X_r(const_cast<NumericMatrix&>(X).begin(), n, p, false);
  arma::vec theta_r(theta.begin(), 2, false);
  arma::vec beta_r(beta.begin(), p3, false);
  arma::vec beta_h_r(beta_h.begin(), p, false);
  arma::vec beta_o_r(beta_o.begin(), p, false);
  arma::vec beta_q_r(beta_q.begin(), p, false);
  arma::vec Xbeta_h_r(Xbeta_h.begin(), n, false);
  arma::vec Xbeta_o_r(Xbeta_o.begin(), n, false);
  arma::vec Xbeta_q_r(Xbeta_q.begin(), n, false);
  arma::vec Ys_r(Ys.begin(), nYs, false);
  
  // Working temp variables
  arma::mat Linv=arma::zeros<arma::mat>(n, nsave);
  arma::vec Dvalues = arma::zeros<arma::vec>(nsave);
  Rcpp::NumericVector sumtheta(2, 0.0); // 2 by 1
  Rcpp::NumericVector sumbeta_h(p, 0.0); // p by 1
  Rcpp::NumericVector sumbeta_o(p, 0.0); // p by 1
  Rcpp::NumericVector sumbeta_q(p, 0.0); // p by 1
  Rcpp::NumericVector sumweight(maxL, 0.0); //  maxL by 1
  // for theta
  arma::vec thetaold(2); 
  arma::vec thBarold(2);
  arma::mat I2 = ESMALL*arma::eye(2,2);
  arma::vec thBarnew=arma::zeros<arma::vec>(2);
  arma::mat thSnew=arma::zeros<arma::mat>(2,2);
  // for beta
  arma::vec betaold(p3);
  arma::vec beBarold(p3);
  arma::mat Ip3 = ESMALL*arma::eye(p3,p3);
  arma::vec beBarnew=arma::zeros<arma::vec>(p3);
  arma::mat beSnew=arma::zeros<arma::mat>(p3,p3);
  // for Ys
  arma::vec Ysold(nYs);
  arma::vec YsBarold(nYs);
  arma::vec YsBarnew=arma::zeros<arma::vec>(nYs);
  arma::mat YsSnew=arma::zeros<arma::mat>(nYs,nYs);
  arma::mat InYs = ESMALL*arma::eye(nYs,nYs);
  arma::mat YsShat=0.16*arma::eye(nYs,nYs);
  // for cpa
  double cbarnew=0, cbarold=0, csnew=0, ctemp=0, cshat=0.16;
  double llold, llnew;
  double ratio, uu, nn;
  
  RNGScope scope;
  
  // Set the Armadillo seed from R's 
  // int seed = (int)Rf_runif(0.0, 10000.0);
  // std::srand(seed);
  
  ////////////////////////////////////////////////////////////////////////
  // Start MCMC
  ////////////////////////////////////////////////////////////////////////
  beta_h_r = beta_r.subvec(0, p-1);
  beta_o_r = beta_r.subvec(p, 2*p-1);
  beta_q_r = beta_r.subvec(2*p, p3-1);
  Xbeta_h_r = X_r*beta_h_r;
  Xbeta_o_r = X_r*beta_o_r;
  Xbeta_q_r = X_r*beta_q_r;
  for (int iscan=0; iscan<nscan; iscan++){
    R_CheckUserInterrupt();
    // Rprintf( "iscan = %d\n", iscan );
    
    ///////////////////////////////////////////////
    // update PT conditional probabilities
    //////////////////////////////////////////////
    if(BP){
      PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, llold);
      llold += Rf_lgammafn(cpar*maxL) - maxL*Rf_lgammafn(cpar) + cpar*Rcpp::sum(Rcpp::log(weight)); 
      Ysold=Ys_r;
      if(iscan>l0){
        Ys_r = mvrnorm(Ysold, YsSnew);
      }else{
        Ys_r = mvrnorm(Ysold, YsShat);
      }
      Ys_to_weight(Ys, weight);
      PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, llnew);
      llnew += Rf_lgammafn(cpar*maxL) - maxL*Rf_lgammafn(cpar) + cpar*Rcpp::sum(Rcpp::log(weight));
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if(uu>ratio){
        Ys_r=Ysold;
        if(iscan>=nburn) rejYs += 1.0;
        Ys_to_weight(Ys, weight);
      }
      nn = iscan+1;
      YsBarold = YsBarnew;
      YsBarnew = (nn)/(nn+1.0)*YsBarold + Ys_r/(nn+1.0);
      YsSnew = (nn-1.0)/nn*YsSnew + adapter/(double)(nYs)/nn*(nn*YsBarold*YsBarold.t() 
                                                                - (nn+1.0)*YsBarnew*YsBarnew.t() + Ys_r*Ys_r.t() + InYs );
    }
    
    ///////////////////////////////////////////////
    // update beta
    //////////////////////////////////////////////
    PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, llold);
    llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
    betaold = beta_r;
    if(iscan>l0){
      beta_r = mvrnorm(betaold, beSnew);
    }else{
      beta_r = mvrnorm(betaold, Shat);
    }
    beta_h_r = beta_r.subvec(0, p-1);
    beta_o_r = beta_r.subvec(p, 2*p-1);
    beta_q_r = beta_r.subvec(2*p, p3-1);
    Xbeta_h_r = X_r*beta_h_r;
    Xbeta_o_r = X_r*beta_o_r;
    Xbeta_q_r = X_r*beta_q_r;
    PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, llnew);
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
    beSnew =(nn-1.0)/nn*beSnew+adapter/(p3+0.0)/nn*(nn*beBarold*beBarold.t()-(nn+1.0)*beBarnew*beBarnew.t()
                                                         +beta_r*beta_r.t()+Ip3);
    beta_h_r = beta_r.subvec(0, p-1);
    beta_o_r = beta_r.subvec(p, 2*p-1);
    beta_q_r = beta_r.subvec(2*p, p3-1);
    Xbeta_h_r = X_r*beta_h_r;
    Xbeta_o_r = X_r*beta_o_r;
    Xbeta_q_r = X_r*beta_q_r;
    
    ///////////////////////////////////////////////
    // update theta
    //////////////////////////////////////////////
    if(V0inv(0,0)<SYSMAX){
      PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, llold);
      llold += -0.5*arma::dot( (theta_r-theta0), (V0inv*(theta_r-theta0)) );
      thetaold = theta_r;
      if(iscan>l0){
        theta_r = mvrnorm(thetaold, thSnew);
      }else{
        theta_r = mvrnorm(thetaold, Vhat);
      }
      PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, llnew);
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
        double tempsum = Rcpp::sum(Rcpp::log(weight));
        llold = Rf_lgammafn(cpar*maxL) - maxL*Rf_lgammafn(cpar) + (cpar-1.0)*tempsum + (a0-1.0)*log(cpar)-b0*cpar; 
        llnew = Rf_lgammafn(ctemp*maxL) - maxL*Rf_lgammafn(ctemp) + (ctemp-1.0)*tempsum + (a0-1.0)*log(ctemp)-b0*ctemp;
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
        Linv.col(isave) = PHPOAFT_BP_invLik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q);
        // calculate -2loglikelihood
        double tempLik = 0;
        PHPOAFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, tempLik);
        Dvalues(isave) = -2.0*tempLik;
        // save regression coefficient
        theta_save(_,isave) = theta;
        sumtheta += theta; 
        beta_h_save(_,isave) = beta_h;
        sumbeta_h += beta_h;
        beta_o_save(_,isave) = beta_o;
        sumbeta_o += beta_o;
        beta_q_save(_,isave) = beta_q;
        sumbeta_q += beta_q;
        // precision parameter
        cpar_save[isave] = cpar;
        // PT probs
        weight_save(_,isave) = weight;
        sumweight += weight;
        
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
  Xbeta_h_r = X_r*as<arma::vec>(sumbeta_h)/(nsave+0.0);
  Xbeta_o_r = X_r*as<arma::vec>(sumbeta_o)/(nsave+0.0);
  Xbeta_q_r = X_r*as<arma::vec>(sumbeta_q)/(nsave+0.0);
  double Dmean = 0;
  PHPOAFT_BP_loglik(t1, t2, ltr, type, sumtheta[0]/(nsave+0.0), sumtheta[1]/(nsave+0.0), 
                    sumweight/(nsave+0.0), BP, dist, Xbeta_h, Xbeta_o, Xbeta_q, Dmean);
  double pD = meanD + 2.0*Dmean;
  double DIC = meanD + pD; 
  
  return List::create(Named("theta")=theta_save,
                      Named("beta_h")=beta_h_save,
                      Named("beta_o")=beta_o_save,
                      Named("beta_q")=beta_q_save,
                      Named("cpar")=cpar_save,
                      Named("weight")=weight_save,
                      Named("cpo")=cpo,
                      Named("pD")=pD,
                      Named("DIC")=DIC,
                      Named("ratetheta")=ratetheta,
                      Named("ratebeta")=ratebeta,
                      Named("rateYs")=rateYs,
                      Named("ratec")=ratec);
  END_RCPP
}

// Get density or survival Plots for frailty LDTFP PO
RcppExport SEXP PHPOAFT_BP_plots(SEXP tgrid_, SEXP xpred_, SEXP theta_, SEXP beta_h_, SEXP beta_o_, SEXP beta_q_, 
                                 SEXP weight_, SEXP CI_, SEXP dist_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  Rcpp::NumericVector tgrid(tgrid_);
  arma::mat xpred = as<arma::mat>(xpred_); // npred by p;
  Rcpp::NumericMatrix theta(theta_); // 2 by nsave;
  Rcpp::NumericMatrix beta_h(beta_h_); // p by nsave;
  Rcpp::NumericMatrix beta_o(beta_o_); // p by nsave;
  Rcpp::NumericMatrix beta_q(beta_q_); // p by nsave;
  Rcpp::NumericMatrix weight(weight_); // maxL by nsave;
  double CI = as<double>(CI_);
  const int dist = as<int>(dist_);
  int nsave = theta.ncol();
  int ngrid = tgrid.size();
  int npred = xpred.n_rows;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  int p = beta_h.nrow();
  
  // Temp variables
  arma::vec xbeta_h(npred);
  arma::vec xbeta_o(npred);
  arma::vec xbeta_q(npred);
  Rcpp::NumericVector estfArray(nsave*ngrid*npred);
  arma::cube estf(estfArray.begin(), ngrid, nsave, npred, false);
  Rcpp::NumericVector estSArray(nsave*ngrid*npred);
  arma::cube estS(estSArray.begin(), ngrid, nsave, npred, false);
  Rcpp::NumericVector esthArray(nsave*ngrid*npred);
  arma::cube esth(esthArray.begin(), ngrid, nsave, npred, false);
  
  // Make arma objects
  arma::mat beta_h_r(beta_h.begin(), p, nsave, false);
  arma::mat beta_o_r(beta_o.begin(), p, nsave, false);
  arma::mat beta_q_r(beta_q.begin(), p, nsave, false);
  
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
    xbeta_h = xpred*beta_h_r.col(i);
    xbeta_o = xpred*beta_o_r.col(i);
    xbeta_q = xpred*beta_q_r.col(i);
    Rcpp::NumericVector wi = weight(_,i);
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        estf(k, i, j) = std::exp(PHPOAFT_BP_logpdf(tgrid[k], theta(0,i), theta(1,i), wi, true, dist,
                                 xbeta_h[j], xbeta_o[j], xbeta_q[j]));
        estS(k, i, j) = std::exp(PHPOAFT_BP_logsurv(tgrid[k], theta(0,i), theta(1,i), wi, true, dist,
                                 xbeta_h[j], xbeta_o[j], xbeta_q[j]));
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
