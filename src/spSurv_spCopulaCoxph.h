#ifndef _spSurv_spCopulaCoxph_h
#define _spSurv_spCopulaCoxph_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

///////////////////////////// Fixed cutpoints //////////////////////////////////
RcppExport SEXP spCopulaCoxph( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
  	SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, 
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, SEXP theta0_, 
    SEXP spl0_, SEXP spS0_, SEXP spadapter_ );

///////////////////////////// Random cutpoints //////////////////////////////////
RcppExport SEXP spCopulaCoxphR( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
    SEXP t_, SEXP delta_, SEXP X_, SEXP d_, SEXP h_, SEXP r0_, SEXP h0_, SEXP V0_, 
    SEXP hl0_, SEXP hs0_, SEXP hadapter_, 
    SEXP beta_, SEXP mu0_, SEXP Sig0_, SEXP l0_, SEXP S0_, SEXP adapter_,
		SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, SEXP theta0_, 
    SEXP spl0_, SEXP spS0_, SEXP spadapter_ );
#endif
