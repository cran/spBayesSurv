"anovaDDP" <- function(formula, data, na.action, prediction=NULL, 
                       mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                       prior=NULL, state=NULL, scale.designX=TRUE) {
  #########################################################################################
  # call parameters
  #########################################################################################
  Call <- match.call(); # save a copy of the call 
  indx <- match(c("formula", "data", "na.action", "truncation_time", "subject.num"),
                names(Call), nomatch=0) 
  if (indx[1] ==0) stop("A formula argument is required");
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1L]] <- quote(stats::model.frame)
  
  special <- c("baseline", "frailtyprior", "truncation_time", "subject.num", "bspline")
  temp$formula <- if (missing(data)) 
    terms(formula, special)
  else terms(formula, special, data = data)
  
  m <- eval(temp, parent.frame())
  Terms <- attr(m, 'terms')
  
  if(any(names(m)=="(truncation_time)")){
    truncation_time = m[,"(truncation_time)"]
  }else{
    truncation_time = NULL
  }
  
  if(any(names(m)=="(subject.num)")){
    subject.num = m[,"(subject.num)"]
  }else{
    subject.num = NULL
  }
  
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  baseline0 <- attr(Terms, "specials")$baseline
  frailtyprior0<- attr(Terms, "specials")$frailtyprior
  bspline0<- attr(Terms, "specials")$bspline
  
  if (length(frailtyprior0)) {
    temp <- survival::untangle.specials(Terms, 'frailtyprior', 1)
    dropfrail <- c(temp$terms)
    frail.terms <- m[[temp$vars]]
  }else{
    dropfrail <- NULL
    frail.terms <- NULL;
  }
  if (length(baseline0)) {
    temp <- survival::untangle.specials(Terms, 'baseline', 1)
    dropXtf <- c(temp$terms)
    Xtf <- m[[temp$vars]]
  }else{
    dropXtf <- NULL
    Xtf <- NULL
  }
  if (length(bspline0)) {
    temp <- survival::untangle.specials(Terms, 'bspline', 1)
    #dropx <- c(dropx, temp$terms);
    X.bs = NULL;
    n.bs = rep(0, length(temp$vars));
    for(ii in 1:length(temp$vars)){
      X.bs = cbind(X.bs, m[[temp$vars[ii]]]);
      n.bs[ii] = ncol(m[[temp$vars[ii]]]); 
    }
  }else{
    X.bs <- NULL;
    n.bs <- NULL;
  }
  
  dropx <- c(dropfrail, dropXtf)
  if (length(dropx)) {
    newTerms <- Terms[-dropx]
    # R (version 2.7.1) adds intercept=T anytime you drop something
    attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
  } else  newTerms <- Terms
  
  X <- model.matrix(newTerms, m);
  assign <- lapply(survival::attrassign(X, newTerms)[-1], function(x) x-1)
  xlevels <- .getXlevels(newTerms, m)
  contr.save <- attr(X, 'contrasts')
  
  # drop the intercept after the fact, and also drop baseline if necessary
  adrop <- 0  #levels of "assign" to be dropped; 0= intercept
  Xatt <- attributes(X) 
  xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  
  n <- nrow(X)
  p <- ncol(X)
  if(p==0){
    #stop("covariate is required")
    X.scaled <- NULL;
    X1 = cbind(rep(1,n), X.scaled); 
  }else{
    if(scale.designX){
      X.scaled <- scale(X);
    }else{
      X.scaled <- scale(X, center=rep(0,p), scale=rep(1,p));
    }
    X.center = attributes(X.scaled)$`scaled:center`;
    X.scale = attributes(X.scaled)$`scaled:scale`;
    X1 = cbind(rep(1,n), X.scaled);
  }
  #########################################################################################
  # data structure
  #########################################################################################
  t1 = Y[,1]; t2 = Y[,1];
  type <- attr(Y, "type")
  exactsurv <- Y[,ncol(Y)] ==1
  if (any(exactsurv)) {
    t1[exactsurv]=Y[exactsurv,1];
    t2[exactsurv]=Y[exactsurv,1];
  }
  if (type== 'counting') stop ("Invalid survival type")
  if (type=='interval') {
    intsurv <- Y[,3]==3;
    if (any(intsurv)){
      t1[intsurv]=Y[intsurv,1];
      t2[intsurv]=Y[intsurv,2];
    }
  } 
  delta = Y[,ncol(Y)];
  if (!all(is.finite(Y))) {
    stop("Invalid survival times for this distribution")
  } else {
    if (type=='left') delta <- 2- delta;
  }
  
  ##############################################
  ### Currently it only supports right censored data 
  ##############################################
  model.name <- "Bayesian Nonparametric Survival Model"
  if(sum(delta%in%c(0,1))!=n) stop("This function currently only supports right-censored data.")
  
  #########################################################################################
  # prediction
  #########################################################################################
  if(p==0){
    nxpred = 1;
    xpred = cbind(rep(1,nxpred));
  }else{
    xpred <- prediction$xpred;
    if(is.null(xpred)){
      xpred = X;
    }
    if(is.vector(xpred)) xpred=matrix(xpred, nrow=1);
    if(ncol(xpred)!=p) stop("Please make sure the number of columns for xpred equals the number of covariates!");
    xpred = cbind(xpred);
    nxpred = nrow(xpred);
    for(i in 1:nxpred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
    xpred <- cbind(rep(1,nxpred), xpred);
  }
  
  #########################################################################################
  # initial analysis and priors
  #########################################################################################
  fit0 <- survival::survreg(formula = Surv(t1, delta) ~ X1-1, dist = "lognormal");
  muhat = as.vector(fit0$coefficients);
  sig2hat = fit0$scale^2
  Sighat = as.matrix(fit0$var[(1:(p+1)),(1:(p+1))]); Sigscale=30
  N <- prior$N; if(is.null(N)) N <- 10;
  m0 <- prior$m0; if(is.null(m0)) m0 <- muhat;
  S0 <- prior$S0; if(is.null(S0)) S0 <- Sighat;
  Sig0 <- prior$Sig0; if(is.null(Sig0)) Sig0 <- Sigscale*Sighat; #Sig0 <- diag(rep(1e5,p+1), nrow=p+1, ncol=p+1); 
  k0 <- prior$k0; if(is.null(k0)) k0 <- p+1+5;
  nua <-prior$nua; nub <- prior$nub;
  if(is.null(nua)) nua=2+1; #nua=2+sig2hat/4; #
  if(is.null(nub)) nub=sig2hat; #nub=sig2hat/4*(nua-1); #
  a0 <-prior$a0; b0 <- prior$b0;
  if(is.null(a0)) a0=2; if(is.null(b0)) b0=2;
  
  #########################################################################################
  # current state and mcmc specification
  #########################################################################################
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  currenty = log(t1);
  mu <- state$mu; if(is.null(mu)) mu <- muhat;
  Sig <- state$Sig; if(is.null(Sig)) Sig <- 25*Sighat;
  beta<- state$beta; if(is.null(beta)) beta = matrix(muhat, p+1, N);
  sigma2<- state$sigma2; if(is.null(sigma2)) sigma2 = rep(sig2hat/2,N);
  alpha <- state$alpha; if(is.null(alpha)) alpha <- 2;
  K = sample(1:N,n, replace=T);
  V = rbeta(N, 1, alpha); V[N] =1;
  w = V; 
  for (k in 2:N){
    w[k] = max(exp( sum(log(1-V[1:(k-1)]))+log(V[k]) ), 1e-320);
    #w[k] = max( (1 - sum(w[1:(k-1)]))*V[k], 1e-320);
  }
  
  #########################################################################################
  # calling the c++ code
  #########################################################################################
  foo <- .Call("anovaDDP", nburn_ = nburn, nsave_ = nsave, nskip_ = nskip, ndisplay_ = ndisplay,
               y_ = currenty, delta_ = delta, X_ = as.matrix(t(X1)), N_ = N, beta_ = beta, 
               tau2_ = 1.0/sigma2, K_ = K, V_ = V, w_ = w, alpha_ = alpha, mu_ = mu, Sig_ = Sig,
               m0_ = m0, S0_ = S0, Sig0_ = Sig0, k0_ = k0, a0_ = a0, b0_ = b0, nua_ = nua, 
               nub_ = nub, xpred_ = as.matrix(xpred), PACKAGE = "spBayesSurv")
  
  #########################################################################################
  # output
  #########################################################################################
  output <- list(modelname=model.name,
                 terms=m,
                 call=Call,
                 n=n,
                 p=p,
                 Surv=survival::Surv(t1, delta),
                 X.scaled=X.scaled,
                 X = X,
                 beta = foo$beta,
                 sigma2 = foo$sigma2,
                 w = foo$w,
                 alpha = foo$alpha,
                 #survt = exp(foo$y),
                 cpo = foo$cpo,
                 Tpred = exp(foo$Ypred),
                 V = foo$V,
                 K = foo$K,
                 state = foo$state);
  class(output) <- c("anovaDDP")
  output
}

"plot.anovaDDP" <- function (x, xnewdata, tgrid=NULL, CI=0.95, PLOT=TRUE, ...) {
  if(is.null(tgrid)) tgrid = seq(0.01, max(x$Surv[,1], na.rm=T), length.out=200);
  if(x$p==0){
    xpred = cbind(1); nxpred=1;
    rnames = "baseline"
  }else{
    if(missing(xnewdata)) {
      stop("please specify xnewdata")
    }else{
      rnames = row.names(xnewdata)
      m = x$terms
      Terms = attr(m, 'terms')
      baseline0 <- attr(Terms, "specials")$baseline
      frailtyprior0<- attr(Terms, "specials")$frailtyprior
      dropx <- NULL
      if (length(frailtyprior0)) {
        temp <- survival::untangle.specials(Terms, 'frailtyprior', 1)
        dropx <- c(dropx, temp$terms)
        frail.terms <- m[[temp$vars]]
      }else{
        frail.terms <- NULL;
      }
      if (length(baseline0)) {
        temp <- survival::untangle.specials(Terms, 'baseline', 1)
        dropx <- c(dropx, temp$terms)
        Xtf <- m[[temp$vars]]
      }else{
        Xtf <- NULL;
      }
      if (length(dropx)) {
        newTerms <- Terms[-dropx]
        # R (version 2.7.1) adds intercept=T anytime you drop something
        attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
      } else  newTerms <- Terms
      newTerms <- delete.response(newTerms)
      mnew <- model.frame(newTerms, xnewdata, na.action = na.omit, xlev = .getXlevels(newTerms, m))
      Xnew <- model.matrix(newTerms, mnew);
      assign <- lapply(survival::attrassign(Xnew, newTerms)[-1], function(x) x-1)
      xlevels <- .getXlevels(newTerms, mnew)
      contr.save <- attr(Xnew, 'contrasts')
      # drop the intercept after the fact, and also drop baseline if necessary
      adrop <- 0  #levels of "assign" to be dropped; 0= intercept
      Xatt <- attributes(Xnew) 
      xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
      Xnew <- Xnew[, !xdrop, drop=FALSE]
      attr(Xnew, "assign") <- Xatt$assign[!xdrop]
      xpred = Xnew
      if(ncol(xpred)!=x$p) stop("please make sure the number of columns matches!");
      X.center = attributes(x$X.scaled)$`scaled:center`;
      X.scale = attributes(x$X.scaled)$`scaled:scale`;
      xpred = cbind(xpred);
      nxpred = nrow(xpred);
      for(i in 1:nxpred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
      xpred <- cbind(rep(1,nxpred), xpred);
    }
  }
  estimates <- .Call("DDPplots", xpred, tgrid, x$beta, x$sigma2, x$w, CI,
                     PACKAGE = "spBayesSurv");
  if(PLOT){
    par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
    plot(tgrid, estimates$Shat[,1], "l", lwd=3, xlab="time", ylab="survival", 
         xlim=c(0, max(tgrid)), ylim=c(0,1));
    for(i in 1:nxpred){
      polygon(x=c(rev(tgrid),tgrid),
              y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
              border=NA,col="lightgray");
    }
    for(i in 1:nxpred){
      lines(tgrid, estimates$Shat[,i], lty=i, lwd=3, col=i);
    }
    legend("topright", rnames, col = 1:nxpred, lty=1:nxpred, ...)
  }
  estimates$tgrid=tgrid;
  invisible(estimates)
}
