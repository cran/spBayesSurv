"GetCurves" <- function(fit, xpred, ygrid, CI=c(0.05, 0.95)){
  if((class(fit)=="anovaDDP")|(class(fit)=="spCopulaDDP")){
    xpred1 = cbind(1.0, xpred);
    output = .DDPplots( xpred1, ygrid, fit$beta, fit$sigma2, fit$w, CI );
  }
  if((class(fit)=="indeptCoxph")|(class(fit)=="spCopulaCoxph")){
    tgrid = exp(ygrid);
    output = .CoxPHplots( xpred, tgrid, fit$beta, fit$h, fit$d, CI );
  }
  class(output) <- c("GetCurves")
  output
}