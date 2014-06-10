#ifndef _spSurv_spCopulaDDP_h
#define _spSurv_spCopulaDDP_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"

RcppExport SEXP spCopulaDDP( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
    SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
    SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
		SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
    SEXP nua_, SEXP nub_, SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, 
    SEXP theta0_, SEXP spl0_, SEXP spS0_, SEXP spadapter_);

///////////////// spatial Copula DDP using FSA /////////////////////////////////////
RcppExport SEXP spCopulaDDP_FSA( SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_,
    SEXP y_, SEXP delta_, SEXP X_, SEXP N_, SEXP beta_, SEXP tau2_,
    SEXP K_, SEXP V_, SEXP w_, SEXP alpha_, SEXP mu_, SEXP Sig_,
  	SEXP m0_, SEXP S0_, SEXP Sig0_, SEXP k0_, SEXP a0_, SEXP b0_, 
    SEXP nua_, SEXP nub_, SEXP xpred_, SEXP ds0n_, SEXP dnn_, SEXP theta_, 
    SEXP theta0_, SEXP spl0_, SEXP spS0_, SEXP spadapter_,
    SEXP dnm_, SEXP dmm_, SEXP blocki_, SEXP ds0m_, SEXP ds0block_);

#endif
