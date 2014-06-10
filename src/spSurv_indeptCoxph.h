#ifndef _spSurv_indeptCoxph_h
#define _spSurv_indeptCoxph_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

///////////////////////////// Fixed cutpoints //////////////////////////////////
RcppExport SEXP indeptCoxph( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, 
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_);

///////////////////////////// Random cutpoints //////////////////////////////////
RcppExport SEXP indeptCoxphR( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
    SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, SEXP V0_,
    SEXP hl0_, SEXP hs0_, SEXP hadapter_,
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_) ;

#endif
