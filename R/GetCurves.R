"GetCurves" <- function(fit, xpred, ygrid, CI=c(0.05, 0.95), xtfpred=NULL, IDpred=NULL, AFTbaseline=FALSE){
  if(!is.null(xpred)){
    if(is.vector(xpred)) xpred = matrix(xpred, 1);
    npred = nrow(xpred);
  }else{
    npred = 1;
  }
  if(!is.null(xtfpred)) {
    if(is.vector(xtfpred)) xtfpred = matrix(xtfpred, 1);
    nxtfpred = nrow(xtfpred);
  }else{
    nxtfpred =1;
  }
  if((class(fit)=="anovaDDP")|(class(fit)=="spCopulaDDP")){
    xpred1 = cbind(rep(1.0,npred), xpred);
    output = .DDPplots( xpred1, ygrid, fit$beta, fit$sigma2, fit$w, CI );
  }else if((class(fit)=="indeptCoxph")|(class(fit)=="spCopulaCoxph")){
    tgrid = exp(ygrid);
    if(is.null(xpred)) stop("please specify xpred!");
    output = .CoxPHplots( xpred, tgrid, fit$beta, fit$h, fit$d, CI );
  }else if((class(fit)=="frailtyGAFT")){
    if(AFTbaseline){
      xpred1 = matrix(0, nxtfpred, dim(fit$beta)[1]);
      xtfpred1 = cbind(rep(1.0,nxtfpred), xtfpred);
      vpred = matrix(0, nxtfpred, length(fit$sigma2));
      output = .frailtyGAFTplots(ygrid, xpred1, xtfpred1, fit$beta, fit$betatf, vpred, fit$sigma2, fit$maxL, CI);
    }else{
      xpred1 = cbind(rep(1.0,npred), xpred);
      xtfpred1 = cbind(rep(1.0,npred), xtfpred);
      if(is.null(fit$frailty)){
        vpred = matrix(0, npred, length(fit$sigma2));
      }else{
        if(is.null(IDpred)) stop("please specify the ID names for each row of xpred!");
        if(!(length(IDpred)==npred)) stop("xpred and IDpred should have same nrows!");
        if(!is.character(IDpred)) stop("IDpred should be character!");
        IDind = rep(0, npred);
        for(i in 1:npred){
          IDind[i] = which(fit$IDnames==IDpred[i]);
        }
        vpred = matrix(fit$v[IDind,],npred);
      }
      output = .frailtyGAFTplots(ygrid, xpred1, xtfpred1, fit$beta, fit$betatf, vpred, fit$sigma2, fit$maxL, CI);
    }
  }else{
    output = NULL;
  }
  class(output) <- c("GetCurves")
  output
}