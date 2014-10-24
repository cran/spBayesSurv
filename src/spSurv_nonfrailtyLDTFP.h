#ifndef _spSurv_nonfrailtyLDTFP_h
#define _spSurv_nonfrailtyLDTFP_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

RcppExport SEXP nonfrailtyLDTFP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP tobs_, SEXP type_, SEXP xce_, SEXP xtf_, SEXP alpha_, SEXP betace_, 
    SEXP betatf_, SEXP sigma2_, SEXP y_, SEXP maxL_,
    SEXP a0_, SEXP b0_, SEXP m0_, SEXP S0inv_, SEXP gprior_, SEXP a0sig_, SEXP b0sig_);


#endif
