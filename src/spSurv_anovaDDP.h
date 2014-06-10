#ifndef _spSurv_anovaDDP_h
#define _spSurv_anovaDDP_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

RcppExport SEXP anovaDDP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
		SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
		SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
    SEXP nua_, SEXP nub_, SEXP xpred_ );
#endif
