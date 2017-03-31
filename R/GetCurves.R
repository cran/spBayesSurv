"GetCurves" <- function (x, xpred=NULL, tgrid=NULL, CI=0.95, PLOT=FALSE, ...) {
  if((is(x,"anovaDDP"))|(is(x,"spCopulaDDP"))){
    if(is.null(tgrid)) tgrid = seq(0.01, max(x$Surv[,1], na.rm=T), length.out=200);
    if(is.null(xpred)) {
      if(x$p==0){
        xpred=matrix(1,nrow=1,ncol=1);
      }else{
        stop("please specify xpred")
      }
    }else{
      if(is.vector(xpred)) xpred=matrix(xpred, nrow=1);
      if(ncol(xpred)!=x$p) stop("please make sure the number of columns matches!");
    }
    if(x$p>0){
      X.center = attributes(x$X.scaled)$`scaled:center`;
      X.scale = attributes(x$X.scaled)$`scaled:scale`;
      xpred = cbind(xpred);
      npred = nrow(xpred);
      for(i in 1:npred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
      xpred <- cbind(rep(1,npred), xpred);
    }else{
      npred = 1;
    }
    estimates <- .Call("DDPplots", xpred, tgrid, x$beta, x$sigma2, x$w, CI,
                       PACKAGE = "spBayesSurv");
    if(PLOT){
      par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
      plot(tgrid, estimates$Shat[,1], "l", lwd=3, xlab="time", ylab="survival", main=paste(i));
      for(i in 1:npred){
        polygon(x=c(rev(tgrid),tgrid),
                y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
                border=NA,col="lightgray");
      }
      for(i in 1:npred){
        lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
      }
    }
  }else if ((is(x,"indeptCoxph"))|(is(x,"spCopulaCoxph"))){
    if(is.null(tgrid)) tgrid = seq(0.01, max(x$Surv[,1], na.rm=T), length.out=200);
    if(is.null(xpred)) {
      stop("please specify xpred")
    }else{
      if(is.vector(xpred)) xpred=matrix(xpred, nrow=1);
      if(ncol(xpred)!=x$p) stop("please make sure the number of columns matches!");
    }
    X.center = attributes(x$X.scaled)$`scaled:center`;
    X.scale = attributes(x$X.scaled)$`scaled:scale`;
    xpred = cbind(xpred);
    npred = nrow(xpred);
    for(i in 1:npred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
    betafitted = x$beta.scaled;
    estimates <- .Call("CoxPHplots", xpred, tgrid, betafitted, x$h.scaled, x$d.scaled, CI,
                       PACKAGE = "spBayesSurv");
    if(PLOT){
      par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
      plot(tgrid, estimates$Shat[,1], "l", lwd=3, xlab="time", ylab="survival", main=paste(i));
      for(i in 1:npred){
        polygon(x=c(rev(tgrid),tgrid),
                y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
                border=NA,col="lightgray");
      }
      for(i in 1:npred){
        lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
      }
    }
  }
  estimates$tgrid=tgrid;
  invisible(estimates)
}