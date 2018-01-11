#include "spSurv_BP_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

RcppExport SEXP AFT_BP(SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_, SEXP ltr_, SEXP subjecti_,
                      SEXP t1_, SEXP t2_, SEXP type_, SEXP X_, SEXP theta_, SEXP beta_, 
                      SEXP weight_, SEXP cpar_, SEXP a0_, SEXP b0_, SEXP theta0_, 
                      SEXP V0inv_, SEXP Vhat_, SEXP beta0_, SEXP S0inv_, SEXP Shat_, 
                      SEXP l0_, SEXP adapter_, SEXP gamma_, SEXP p0gamma_, SEXP selection_,
                      SEXP frailty_, SEXP v_, SEXP blocki_, SEXP W_, SEXP clustindx_,
                      SEXP Dmm_, SEXP Dmr_, SEXP Drr_, SEXP phi_, SEXP nu_, SEXP a0phi_, SEXP b0phi_,
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
  const arma::mat Dmm = as<mat>(Dmm_); // m by m;
  const arma::mat Dmr = as<mat>(Dmr_); // m by r; 
  const arma::mat Drr = as<mat>(Drr_); // r by r; r is # of knots
  const arma::mat clustindx = as<mat>(clustindx_); // m by nclust;
  const int r = Drr.n_cols;
  const bool FSA = (r<m);
  
  // prameters to be updated
  double cpar = as<double>(cpar_);
  const bool BP = (cpar<R_PosInf);
  Rcpp::NumericVector theta(theta_); // 2 by 1
  Rcpp::NumericVector beta(beta_); // p by 1
  Rcpp::NumericVector weight(weight_); //  maxL by 1
  const int maxL = weight.size();
  Rcpp::NumericVector gamma(gamma_); // p by 1
  double lambda = as<double>(lambda_);
  double phi = as<double>(phi_);
  double nu = as<double>(nu_);
  
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
  // prior of phi
  const double a0phi = as<double>(a0phi_);
  const double b0phi = as<double>(b0phi_);
  
  // temp variables
  const int nscan = nburn + (nskip+1)*nsave;
  const int nYs = maxL-1;
  Rcpp::NumericVector Xbeta(n, 0.0);
  Rcpp::NumericVector Ys(nYs);
  for(int i=0; i<nYs; ++i) Ys[i] = log(weight[i])-log(weight[nYs]);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  Rcpp::NumericVector vn(n,0.0);
  
  // things to save;
  Rcpp::NumericVector cpar_save(nsave);
  Rcpp::NumericMatrix theta_save(2, nsave);
  Rcpp::NumericMatrix beta_save(p, nsave);
  Rcpp::NumericMatrix weight_save(maxL, nsave);
  Rcpp::NumericMatrix gamma_save(p, nsave);
  double rejtheta=0;
  double rejbeta=0;
  double rejYs=0;
  double rejc=0;
  double rejphi=0;
  Rcpp::NumericVector rejv(m, 0.0);
  Rcpp::NumericMatrix v_save(m, nsave);
  Rcpp::NumericVector lambda_save(nsave);
  Rcpp::NumericVector phi_save(nsave);
  
  // Make arma objects
  arma::mat X_r(const_cast<NumericMatrix&>(X).begin(), n, p, false);
  arma::vec theta_r(theta.begin(), 2, false);
  arma::vec beta_r(beta.begin(), p, false);
  arma::vec gamma_r(gamma.begin(), p, false);
  arma::vec Xbeta_r(Xbeta.begin(), n, false);
  arma::vec Ys_r(Ys.begin(), nYs, false);
  arma::vec v_r(v.begin(), m, false);
  arma::vec vn_r(vn.begin(), n, false);
  
  // Working temp variables
  arma::mat logLik=arma::zeros<arma::mat>(n, nsave);
  arma::mat St1=arma::zeros<arma::mat>(nsub, nsave);
  arma::mat St2=arma::zeros<arma::mat>(nsub, nsave);
  arma::vec Dvalues = arma::zeros<arma::vec>(nsave);
  Rcpp::NumericVector sumtheta(2, 0.0); // 2 by 1
  Rcpp::NumericVector sumbeta(p, 0.0); // p by 1
  Rcpp::NumericVector sumweight(maxL, 0.0); //  maxL by 1
  Rcpp::NumericVector sumvn(n, 0.0); 
  double llold, llnew;
  double ratio, uu, nn;
  //double phi_min = std::pow(-log(0.001), 1.0/nu)/Dmm.max();
  double phi_min = ESMALL;
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
  arma::vec Ysold(nYs);
  arma::vec YsBarold(nYs);
  arma::vec YsBarnew=arma::zeros<arma::vec>(nYs);
  arma::mat YsSnew=arma::zeros<arma::mat>(nYs,nYs);
  arma::mat InYs = ESMALL*arma::eye(nYs,nYs);
  arma::mat YsShat=0.16*arma::eye(nYs,nYs);
  // for cpa
  double cbarnew=0, cbarold=0, csnew=0, ctemp=0, cshat=0.16;
  // for log(phi)
  double phibarnew=0, phibarold=0, phisnew=0, phinew=0, phishat=0.01;
  // for sill
  double sill=1.0-ESMALL;
  
  // for v
  arma::mat Im = arma::eye(m,m);
  arma::mat Rmm = arma::zeros<arma::mat>(m,m); // Correlation matrix of v
  arma::mat Rmr = arma::zeros<arma::mat>(m,r); 
  arma::mat Rrr = arma::zeros<arma::mat>(r,r); 
  arma::mat Rmminv = arma::zeros<arma::mat>(m,m); 
  double logdetR = 0;
  if(frailty==3){
    for(int i=0; i<m; ++i){
      for(int j=0; j<m; ++j){
        Rmm(i,j) = pow_exp(Dmm(i,j),phi,nu);
      }
    }
    if(FSA){
      for(int i=0; i<m; ++i){
        for(int j=0; j<r; ++j){
          Rmr(i,j) = pow_exp(Dmr(i,j),phi,nu);
        }
      }
      for(int i=0; i<r; ++i){
        for(int j=0; j<r; ++j){
          Rrr(i,j) = pow_exp(Drr(i,j),phi,nu);
        }
      }
      inv_FSA(sill, Rmm, Rmr, Rrr, clustindx, Rmminv, logdetR);
    }else{
      Rmminv = arma::inv_sympd(sill*Rmm+(1.0-sill)*Im);
      double sign0;
      arma::log_det(logdetR, sign0, sill*Rmm+(1.0-sill)*Im);
    }
  }
  if(frailty==1)  {
    v_r = v_r - arma::mean(v_r);
  }
  if((frailty==1)|(frailty==2)|(frailty==3)){
    for(int i=0; i<m; ++i){
      int ind1 = blocki[i];
      int ind2 = blocki[i+1]-1;
      (vn_r.subvec(ind1, ind2)).fill(v[i]);
    }
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
    //Rprintf( "iscan = %d\n", iscan );
    //Rprintf( "lambda = %f\n", lambda );
    //Rprintf( "phi = %f\n", phi );
    
    ///////////////////////////////////////////////
    // update BP weights
    //////////////////////////////////////////////
    if(BP){
      AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llold);
      llold += cpar*Rcpp::sum(Rcpp::log(weight)); 
      Ysold=Ys_r;
      if(iscan>l0){
        Ys_r = mvrnorm(Ysold, YsSnew);
      }else{
        Ys_r = mvrnorm(Ysold, YsShat);
      }
      Ys_to_weight(Ys, weight);
      AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llnew);
      llnew += cpar*Rcpp::sum(Rcpp::log(weight));
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
    AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llold);
    llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
    betaold = beta_r;
    if(iscan>l0){
      beta_r = mvrnorm(betaold, beSnew);
    }else{
      beta_r = mvrnorm(betaold, Shat);
    }
    Xbeta_r = X_r*(beta_r%gamma_r);
    AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llnew);
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
      AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llold);
      llold += -0.5*arma::dot( (theta_r-theta0), (V0inv*(theta_r-theta0)) );
      thetaold = theta_r;
      if(iscan>l0){
        theta_r = mvrnorm(thetaold, thSnew);
      }else{
        theta_r = mvrnorm(thetaold, Vhat);
      }
      AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llnew);
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
    if(frailty==1){// CAR frailty
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        double meanvi = arma::as_scalar(W.row(i)*v_r)/D[i];
        double sdvi = std::sqrt(1.0/(D[i]*lambda));
        double viold = v[i];
        AFT_BP_loglikblocki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta, llold, ind1, ind2, v[i]);
        llold += -0.5*D[i]*lambda*std::pow(v[i]-meanvi,2);
        v[i] = Rf_rnorm(viold, sdvi);
        AFT_BP_loglikblocki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta, llnew, ind1, ind2, v[i]);
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
    } else if(frailty==2){// IID frailty
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        double viold = v[i];
        AFT_BP_loglikblocki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta, llold, ind1, ind2, v[i]);
        llold += -0.5*lambda*std::pow(v[i],2);
        v[i] = Rf_rnorm(viold, std::sqrt(1.0/lambda));
        AFT_BP_loglikblocki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta, llnew, ind1, ind2, v[i]);
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
    } else if (frailty==3){// Gaussian random field frailty for clustered point-referenced data.
      for(int i=0; i<m; ++i){
        int ind1 = blocki[i];
        int ind2 = blocki[i+1]-1;
        double Rii = Rmminv(i,i);
        double meanvi = -(arma::as_scalar(Rmminv.row(i)*v_r)-Rii*v[i])/Rii;
        double sdvi = std::sqrt(1.0/(Rii*lambda));
        double viold = v[i];
        AFT_BP_loglikblocki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta, llold, ind1, ind2, v[i]);
        llold += -0.5*Rii*lambda*std::pow(v[i]-meanvi,2);
        v[i] = Rf_rnorm(viold, sdvi);
        AFT_BP_loglikblocki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta, llnew, ind1, ind2, v[i]);
        llnew += -0.5*Rii*lambda*std::pow(v[i]-meanvi,2);
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
      double b0lambdastar = b0lambda+0.5*(arma::dot(v_r,Rmminv*v_r)); 
      lambda = Rf_rgamma( a0lambdastar, 1.0/b0lambdastar );
      // phi
      if(a0phi>0){
        if(iscan>l0){
          phinew = Rf_rnorm(phi, std::sqrt(phisnew));
        }else{
          phinew = Rf_rnorm(phi, std::sqrt(phishat));
        }
        if(phinew<phi_min){
          if(iscan>=nburn) rejphi+=1.0;
        }else{
          arma::mat Rmminvnew(m,m);
          double logdetRnew=0;
          llold = -0.5*logdetR - 0.5*lambda*arma::dot(v_r, Rmminv*v_r);
          llold += (a0phi-1.0)*log(phi) - b0phi*phi;
          for(int i=0; i<m; ++i){
            for(int j=0; j<m; ++j){
              Rmm(i,j) = pow_exp(Dmm(i,j),phinew,nu);
            }
          }
          if(FSA){
            for(int i=0; i<m; ++i){
              for(int j=0; j<r; ++j){
                Rmr(i,j) = pow_exp(Dmr(i,j),phinew,nu);
              }
            }
            for(int i=0; i<r; ++i){
              for(int j=0; j<r; ++j){
                Rrr(i,j) = pow_exp(Drr(i,j),phinew,nu);
              }
            }
            inv_FSA(sill, Rmm, Rmr, Rrr, clustindx, Rmminvnew, logdetRnew);
          }else{
            Rmminvnew = arma::inv_sympd(sill*Rmm+(1.0-sill)*Im);
            double sign0;
            arma::log_det(logdetRnew, sign0, sill*Rmm+(1.0-sill)*Im);
          }
          llnew = -0.5*logdetRnew - 0.5*lambda*arma::dot(v_r, Rmminvnew*v_r);
          llnew += (a0phi-1.0)*log(phinew) - b0phi*phinew;
          ratio = exp(llnew-llold);
          uu = unif_rand();
          if(uu<ratio){
            phi=phinew; Rmminv = Rmminvnew; logdetR=logdetRnew;
          }else{
            if(iscan>=nburn) rejphi+=1.0;
          }
        }
        nn = iscan+1;
        phibarold = phibarnew;
        phibarnew = (nn)/(nn+1.0)*phibarold + phi/(nn+1.0);
        phisnew = (nn-1.0)/nn*phisnew + adapter/nn*(nn*pow(phibarold,2) - (nn+1.0)*pow(phibarnew,2)
                                                      + pow(phi,2) + ESMALL );
      }
    }
    
    ///////////////////////////////////////////////
    // gamma
    //////////////////////////////////////////////
    if(selection==1){
      for(int j=0; j<p; ++j){
        arma::vec gamma_tmp = gamma_r;
        gamma_tmp[j] = 1.0;
        Xbeta_r = X_r*(beta_r%gamma_tmp);
        AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llold);
        gamma_tmp[j] = 0.0;
        Xbeta_r = X_r*(beta_r%gamma_tmp);
        AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, llnew);
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
        // calculate loglikelihood
        logLik.col(isave) = AFT_BP_logliki(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn);
        // calculate -2loglikelihood
        double tempLik = 0;
        AFT_BP_loglik(t1, t2, ltr, type, theta[0], theta[1], weight, BP, dist, Xbeta+vn, tempLik);
        Dvalues(isave) = -2.0*tempLik;
        // Cox-Snell residuals
        for(int i=0; i<nsub; ++i){
          int ind1 = subjecti[i];
          int ind2 = subjecti[i+1]-1;
          double logSt1tmp=0;
          double logSt2tmp=0;
          for(int j=ind1; j<=ind2; ++j){
            double logSt0tmp=0;
            if(ltr[j]>0){
              logSt0tmp = -AFT_BP_logsurv(ltr[j], theta[0], theta[1], weight, BP, dist, Xbeta[j]+vn[j]);
            }
            if((type[j]==0)|(type[j]==1)){
              double logtmp = AFT_BP_logsurv(t1[j], theta[0], theta[1], weight, BP, dist, Xbeta[j]+vn[j]);
              logSt1tmp += logtmp; logSt1tmp += logSt0tmp;
              logSt2tmp += logtmp; logSt2tmp += logSt0tmp;
            }else if(type[j]==2){
              double logtmp = AFT_BP_logsurv(t2[j], theta[0], theta[1], weight, BP, dist, Xbeta[j]+vn[j]);
              logSt1tmp += logtmp; logSt1tmp += logSt0tmp;
              logSt2tmp += logtmp; logSt2tmp += logSt0tmp;
            }else{
              double logtmp1 = AFT_BP_logsurv(t1[j], theta[0], theta[1], weight, BP, dist, Xbeta[j]+vn[j]);
              double logtmp2 = AFT_BP_logsurv(t2[j], theta[0], theta[1], weight, BP, dist, Xbeta[j]+vn[j]);
              logSt1tmp += logtmp1; logSt1tmp += logSt0tmp;
              logSt2tmp += logtmp2; logSt2tmp += logSt0tmp;
            }
          }
          St1(i,isave) = exp(logSt1tmp);
          St2(i,isave) = exp(logSt2tmp);
        }
        // save regression coefficient
        theta_save(_,isave) = theta;
        sumtheta += theta; 
        beta_save(_,isave) = beta;
        sumbeta += beta*gamma;
        gamma_save(_, isave) = gamma;
        // precision parameter
        cpar_save[isave] = cpar;
        // PT probs
        weight_save(_,isave) = weight;
        sumweight += weight;
        // frailties
        v_save(_,isave) = v;
        sumvn += vn;
        lambda_save[isave] = lambda;
        // phi
        phi_save[isave] = phi;
        
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
  double ratephi = 1.0 - rejphi/totscans;
  
  // get cpo
  arma::mat logLik_subject(nsub, nsave);
  for(int i=0; i<nsub; ++i){
    int ind1 = subjecti[i];
    int ind2 = subjecti[i+1]-1;
    logLik_subject.row(i) = arma::sum(logLik.rows(ind1, ind2), 0);
  }
  arma::mat fpostinv = arma::exp(-logLik_subject);
  arma::vec Linvmean = arma::mean(fpostinv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  // get DIC
  double meanD = arma::mean(Dvalues);
  Xbeta_r = X_r*as<arma::vec>(sumbeta)/(nsave+0.0);
  double Dmean = 0;
  AFT_BP_loglik(t1, t2, ltr, type, sumtheta[0]/(nsave+0.0), sumtheta[1]/(nsave+0.0), 
               sumweight/(nsave+0.0), BP, dist, Xbeta+sumvn/(nsave+0.0), Dmean);
  double pD = meanD + 2.0*Dmean;
  double DIC = meanD + pD; 
  arma::vec DIC_pD(2); DIC_pD[0]=DIC; DIC_pD[1]=pD;
  
  //get stabilized cpo
  arma::mat fpost = arma::exp(logLik_subject);
  arma::mat fpostinv_bd = arma::zeros<mat>(nsub, nsave);
  fpostinv_bd.each_col() = std::sqrt(nsave+0.0)*Linvmean;
  fpostinv_bd = arma::min(fpostinv, fpostinv_bd);
  arma::vec cpo_stab = arma::mean(fpostinv_bd%fpost, 1)/arma::mean(fpostinv_bd, 1);
  
  // get WAIC
  arma::vec fpostmean=arma::mean(fpost, 1);
  double lpd = arma::sum(arma::log(fpostmean));
  arma::vec logfpostvar=arma::var(logLik_subject, 0, 1);
  double pwaic = arma::sum(logfpostvar);
  double WAIC = -2.0*lpd + 2.0*pwaic;
  arma::vec WAIC_pwaic(2); WAIC_pwaic[0]=WAIC; WAIC_pwaic[1]=pwaic;
  
  // get Cox-Snell residuals
  arma::vec r1 = -arma::log(arma::mean(St1, 1));
  arma::vec r2 = -arma::log(arma::mean(St2, 1));
  
  return List::create(Named("theta")=theta_save,
                      Named("beta")=beta_save,
                      Named("gamma")=gamma_save,
                      Named("cpar")=cpar_save,
                      Named("weight")=weight_save,
                      Named("cpo")=cpo,
                      Named("cpo_stab")=cpo_stab,
                      Named("DIC_pD")=DIC_pD,
                      Named("WAIC_pwaic")=WAIC_pwaic,
                      Named("resid1")=r1,
                      Named("resid2")=r2,
                      Named("ratetheta")=ratetheta,
                      Named("ratebeta")=ratebeta,
                      Named("rateYs")=rateYs,
                      Named("ratec")=ratec,
                      Named("v")=v_save,
                      Named("lambda")=lambda_save,
                      Named("ratev")=ratev,
                      Named("phi")=phi_save,
                      Named("ratephi")=ratephi);
  END_RCPP
}

// Get LOO and WAIC
RcppExport SEXP AFT_BP_loo_waic(SEXP ltr_, SEXP subjecti_, SEXP t1_, SEXP t2_, SEXP type_, 
                               SEXP X_, SEXP theta_, SEXP beta_, SEXP vn_, SEXP weight_, SEXP dist_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  const Rcpp::NumericVector ltr(ltr_); // n by 1;
  const Rcpp::IntegerVector subjecti(subjecti_); //nsub+1 by 1;
  const Rcpp::NumericVector t1(t1_); // n by 1;
  const Rcpp::NumericVector t2(t2_); // n by 1;
  const Rcpp::IntegerVector type(type_);// n by 1;
  const arma::mat X=as<arma::mat>(X_); // n by p;
  const int nsub = subjecti.size()-1;
  const arma::mat theta=as<arma::mat>(theta_); // 2 by nsave;
  const arma::mat beta=as<arma::mat>(beta_); // p by nsave;
  const arma::mat vn=as<arma::mat>(vn_); // n by nsave;
  const Rcpp::NumericMatrix weight(weight_); // maxL by nsave;
  const int dist = as<int>(dist_);
  int nsave = beta.n_cols;
  int n = X.n_rows;
  
  // Temp variables
  arma::mat logLik=arma::zeros<arma::mat>(n, nsave);
  Rcpp::NumericVector xbeta(n, 0.0);
  arma::vec xbeta_r(xbeta.begin(), n, false);
  
  // get LOO and WAIC
  for(int isave=0; isave<nsave; ++isave){
    double th1 = theta(0,isave);
    double th2 = theta(1,isave);
    Rcpp::NumericVector wi = weight(_,isave);
    xbeta_r = X*beta.col(isave)+vn.col(isave);
    // calculate loglikelihood
    logLik.col(isave) = AFT_BP_logliki(t1, t2, ltr, type, th1, th2, wi, true, dist, xbeta);
  }
  
  // get cpo
  arma::mat logLik_subject(nsub, nsave);
  for(int i=0; i<nsub; ++i){
    int ind1 = subjecti[i];
    int ind2 = subjecti[i+1]-1;
    logLik_subject.row(i) = arma::sum(logLik.rows(ind1, ind2), 0);
  }
  arma::mat fpostinv = arma::exp(-logLik_subject);
  arma::vec Linvmean = arma::mean(fpostinv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  //get stabilized cpo
  arma::mat fpost = arma::exp(logLik_subject);
  arma::mat fpostinv_bd = arma::zeros<mat>(nsub, nsave);
  fpostinv_bd.each_col() = std::sqrt(nsave+0.0)*Linvmean;
  fpostinv_bd = arma::min(fpostinv, fpostinv_bd);
  arma::vec cpo_stab = arma::mean(fpostinv_bd%fpost, 1)/arma::mean(fpostinv_bd, 1);
  
  // get WAIC
  arma::vec fpostmean=arma::mean(fpost, 1);
  double lpd = arma::sum(arma::log(fpostmean));
  arma::vec logfpostvar=arma::var(logLik_subject, 0, 1);
  double pwaic = arma::sum(logfpostvar);
  double WAIC = -2.0*lpd + 2.0*pwaic;
  arma::vec WAIC_pwaic(2); WAIC_pwaic[0]=WAIC; WAIC_pwaic[1]=pwaic;
  
  return List::create(Named("cpo")=cpo,
                      Named("cpo_stab")=cpo_stab,
                      Named("WAIC_pwaic")=WAIC_pwaic,
                      Named("fpost")=fpostmean);
  END_RCPP
}

// Get Cox-Snell residuals
RcppExport SEXP AFT_BP_cox_snell(SEXP ltr_, SEXP subjecti_, SEXP t1_, SEXP t2_, SEXP type_, 
                                SEXP X_, SEXP theta_, SEXP beta_, SEXP vn_, SEXP weight_, SEXP dist_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  const Rcpp::NumericVector ltr(ltr_); // n by 1;
  const Rcpp::IntegerVector subjecti(subjecti_); //nsub+1 by 1;
  const Rcpp::NumericVector t1(t1_); // n by 1;
  const Rcpp::NumericVector t2(t2_); // n by 1;
  const Rcpp::IntegerVector type(type_);// n by 1;
  const arma::mat X=as<arma::mat>(X_); // n by p;
  const int nsub = subjecti.size()-1;
  const arma::mat theta=as<arma::mat>(theta_); // 2 by nsave;
  const arma::mat beta=as<arma::mat>(beta_); // p by nsave;
  const arma::mat vn=as<arma::mat>(vn_); // n by nsave;
  const Rcpp::NumericMatrix weight(weight_); // maxL by nsave;
  const int dist = as<int>(dist_);
  int nsave = beta.n_cols;
  
  // Temp variables
  arma::mat St1=arma::zeros<arma::mat>(nsub, nsave);
  arma::mat St2=arma::zeros<arma::mat>(nsub, nsave);
  
  // get residuals
  for(int isave=0; isave<nsave; ++isave){
    double th1 = theta(0,isave);
    double th2 = theta(1,isave);
    Rcpp::NumericVector wi = weight(_,isave);
    arma::vec xbeta = X*beta.col(isave)+vn.col(isave);
    for(int i=0; i<nsub; ++i){
      int ind1 = subjecti[i];
      int ind2 = subjecti[i+1]-1;
      double logSt1tmp=0;
      double logSt2tmp=0;
      for(int j=ind1; j<=ind2; ++j){
        double logSt0tmp=0;
        if(ltr[j]>0){
          logSt0tmp = -AFT_BP_logsurv(ltr[j], th1, th2, wi, true, dist, xbeta[j]);
        }
        if((type[j]==0)|(type[j]==1)){
          double logtmp = AFT_BP_logsurv(t1[j], th1, th2, wi, true, dist, xbeta[j]);
          logSt1tmp += logtmp; logSt1tmp += logSt0tmp;
          logSt2tmp += logtmp; logSt2tmp += logSt0tmp;
        }else if(type[j]==2){
          double logtmp = AFT_BP_logsurv(t2[j], th1, th2, wi, true, dist, xbeta[j]);
          logSt1tmp += logtmp; logSt1tmp += logSt0tmp;
          logSt2tmp += logtmp; logSt2tmp += logSt0tmp;
        }else{
          double logtmp1 = AFT_BP_logsurv(t1[j], th1, th2, wi, true, dist, xbeta[j]);
          double logtmp2 = AFT_BP_logsurv(t2[j], th1, th2, wi, true, dist, xbeta[j]);
          logSt1tmp += logtmp1; logSt1tmp += logSt0tmp;
          logSt2tmp += logtmp2; logSt2tmp += logSt0tmp;
        }
      }
      St1(i,isave) = exp(logSt1tmp);
      St2(i,isave) = exp(logSt2tmp);
    }
  }
  // get Cox-Snell residuals
  //arma::vec r1 = -arma::log(arma::mean(St1, 1));
  //arma::vec r2 = -arma::log(arma::mean(St2, 1));
  
  return List::create(Named("St1")=St1,
                      Named("St2")=St2);
  END_RCPP
}


// Get density or survival Plots for frailty LDTFP AFT
RcppExport SEXP AFT_BP_plots(SEXP tgrid_, SEXP xpred_, SEXP theta_, SEXP beta_, SEXP v_,
                            SEXP weight_, SEXP CI_, SEXP dist_){
  BEGIN_RCPP
  // Transfer R variables into C++;
  Rcpp::NumericVector tgrid(tgrid_);
  arma::mat xpred = as<arma::mat>(xpred_); // npred by p;
  arma::mat beta = as<arma::mat>(beta_); // p by nsave;
  arma::mat v = as<arma::mat>(v_); // npred by nsave;
  Rcpp::NumericMatrix theta(theta_); // 2 by nsave;
  Rcpp::NumericMatrix weight(weight_); // maxL by nsave;
  double CI = as<double>(CI_);
  const int dist = as<int>(dist_);
  int nsave = beta.n_cols;
  int ngrid = tgrid.size();
  int npred = xpred.n_rows;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  
  // Temp variables
  arma::vec xbeta(npred);
  Rcpp::NumericVector estfArray(nsave*ngrid*npred);
  arma::cube estf(estfArray.begin(), ngrid, nsave, npred, false);
  Rcpp::NumericVector estSArray(nsave*ngrid*npred);
  arma::cube estS(estSArray.begin(), ngrid, nsave, npred, false);
  Rcpp::NumericVector esthArray(nsave*ngrid*npred);
  arma::cube esth(esthArray.begin(), ngrid, nsave, npred, false);
  
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
    xbeta = xpred*beta.col(i)+v.col(i);
    Rcpp::NumericVector wi = weight(_,i);
    for(int j=0; j<npred; ++j){
      for(int k=0; k<ngrid; ++k){
        estf(k, i, j) = std::exp(AFT_BP_logpdf(tgrid[k], theta(0,i), theta(1,i), wi, true, dist,
                                 xbeta[j]));
        estS(k, i, j) = std::exp(AFT_BP_logsurv(tgrid[k], theta(0,i), theta(1,i), wi, true, dist,
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
