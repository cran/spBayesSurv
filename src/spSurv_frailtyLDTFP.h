#ifndef _spSurv_frailtyLDTFP_h
#define _spSurv_frailtyLDTFP_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

RcppExport SEXP frailtyLDTFP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP tobs_, SEXP type_, SEXP xce_, SEXP xtf_, SEXP alpha_, SEXP betace_, 
    SEXP betatf_, SEXP sigma2_, SEXP y_, SEXP v_, SEXP blocki_, SEXP tau2_, SEXP W_, SEXP D_, SEXP maxL_,
    SEXP a0_, SEXP b0_, SEXP m0_, SEXP S0inv_, SEXP gprior_, SEXP a0sig_, SEXP b0sig_, 
    SEXP a0tau_, SEXP b0tau_);


#endif
