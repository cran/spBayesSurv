"SpatDensReg" <- function (formula, data, na.action, prior=NULL, state=NULL,
                           mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500),
                           permutation=TRUE, fix.theta=TRUE) {
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
  
  if (is.R()) 
    m <- eval(temp, parent.frame())
  else m <- eval(temp, sys.parent())
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
    if (is.R()) attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
  } else  newTerms <- Terms
  
  X <- model.matrix(newTerms, m);
  if (is.R()) {
    assign <- lapply(survival::attrassign(X, newTerms)[-1], function(x) x-1)
    xlevels <- .getXlevels(newTerms, m)
    contr.save <- attr(X, 'contrasts')
  }else {
    assign <- lapply(attr(X, 'assign')[-1], function(x) x -1)
    xvars <- as.character(attr(newTerms, 'variables'))
    xvars <- xvars[-attr(newTerms, 'response')]
    if (length(xvars) >0) {
      xlevels <- lapply(m[xvars], levels)
      xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
      if(length(xlevels) == 0)
        xlevels <- NULL
    } else xlevels <- NULL
    contr.save <- attr(X, 'contrasts')
  }
  
  # drop the intercept after the fact, and also drop baseline if necessary
  adrop <- 0  #levels of "assign" to be dropped; 0= intercept
  Xatt <- attributes(X) 
  xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  X <- X[, !xdrop, drop=FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  
  n <- nrow(X)
  p <- ncol(X)
  if(p==0){
    stop("covariate is required; you could creat a working covariate and set phi=0")
  }
  Sinv = solve(var(X));
  
  # find the maximumu M distance between X[i,] and colMeans(X)
  distseq = rep(0, n);
  Xbar = colMeans(X);
  for(i in 1:n) distseq[i] = sqrt(as.vector((X[i,]-Xbar)%*%Sinv%*%(X[i,]-Xbar)))
  maxdist = max(distseq)
  phi0 = (-log(0.001))/maxdist;
  
  #########################################################################################
  # data structure
  #########################################################################################
  y1 = Y[,1]; y2 = Y[,1];
  type <- attr(Y, "type")
  exactsurv <- Y[,ncol(Y)] ==1
  if (any(exactsurv)) {
    y1[exactsurv]=Y[exactsurv,1];
    y2[exactsurv]=Y[exactsurv,1];
  }
  if (type== 'counting') stop ("Invalid survival type")
  if (type=='interval') {
    intsurv <- Y[,3]==3;
    if (any(intsurv)){
      y1[intsurv]=Y[intsurv,1];
      y2[intsurv]=Y[intsurv,2];
    }
  } 
  delta = Y[,ncol(Y)];
  if (!all(is.finite(Y))) {
    stop("Invalid survival times for this distribution")
  } else {
    if (type=='left') delta <- 2- delta;
  }
  
  #########################################################################################
  # initial MLE analysis and mcmc parameters
  #########################################################################################
  fit0 <- survival::survreg(formula = survival::Surv(y1, y2, type="interval2")~1, dist="gaussian");
  theta1 = fit0$coefficients[1];
  theta2 = log(fit0$scale);
  theta0 = c(theta1, theta2); theta_prior = c(theta1, theta2);
  Vhat0 = as.matrix(fit0$var[c(1,2),c(1,2)]);
  
  #########################################################################################
  # priors and initial values
  #########################################################################################
  alpha=state$alpha; if(is.null(alpha)) alpha=1;
  theta=state$theta; if(is.null(theta)) theta=c(theta1, theta2);
  phi = state$phi; if(is.null(phi)) phi=phi0;
  y <- state$y; 
  if(is.null(y)){
    y <- rep(0, n);
    for(i in 1:n){
      if(delta[i]==0) y[i] = y1[i]+sd(y);  
      if(delta[i]==1) y[i] = y1[i]; 
      if(delta[i]==2) y[i] = y2[i]-sd(y);
      if(delta[i]==3) y[i] = mean(c(y1[i], y2[i]));
    }
  }
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  maxL <- prior$maxL; if(is.null(maxL)) maxL<-5;
  a0=prior$a0; if(is.null(a0)) a0=5;
  b0=prior$b0; if(is.null(b0)) b0=1;
  if(fix.theta){
    V0_prior = diag(0, 2);
  }else{
    V0_prior = 10*Vhat0;
  }
  theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- theta_prior;
  V0 <- prior$V0; if(is.null(V0)) V0 <- V0_prior;
  if(sum(abs(V0))==0){
    V0inv <- diag(c(Inf,Inf));
  }else {
    V0inv <- solve(V0);
  }
  Vhat <- prior$Vhat; if(is.null(Vhat)) Vhat <- Vhat0;
  phiq0 = prior$phiq0; if(is.null(phiq0)) phiq0=0.5;
  phia0 = prior$phia0; if(is.null(phia0)) phia0=2;
  phib0 = prior$phib0; if(is.null(phib0)) phib0=1/phi0;
  ## save to output list
  mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
  theta_initial=theta+0;
  alpha_initial=alpha+0;
  phi_initial = phi+0;
  initial.values = list(alpha=alpha_initial, theta=theta_initial, phi=phi_initial);
  prior = list(maxL=maxL, a0=a0, b0=b0, theta0=theta0, V0=V0, Vhat=Vhat, 
               phiq0=phiq0, phia0=phia0, phib0=phib0);
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  y1new=y1; y2new=y2;
  for(i in 1:n){
    if(delta[i]==0) y2new[i] = Inf;  
    if(delta[i]==2) y1new[i] = -Inf;
  }
  model.name <- "Spatially Smoothed Polya Tree Density Estimation:";
  foo <- .Call("SpatDens", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, 
               y_=y, y1_=y1new, y2_=y2new, type_=delta, X_=t(X), theta_=theta, maxJ_=maxL, 
               cpar_=alpha, a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, 
               l0_=min(5000,nsave/2), adapter_=2.38^2, Sinv_=Sinv, phi_=phi, q0phi_=phiq0, 
               a0phi_=phia0, b0phi_=phib0, perm_=permutation+0, PACKAGE = "spBayesSurv");
  
  ## Bayes Factor for the spatial model vs. the exchangeable model
  q.bar = mean(foo$phi==0);
  BF = (phiq0*(1-q.bar))/((1-phiq0)*q.bar);
  #########################################################################################
  # save to a list
  #########################################################################################
  output <- list(modelname=model.name,
                 terms=m,
                 call=Call,
                 prior=prior,
                 mcmc=mcmc,
                 n=n,
                 p=p,
                 Surv=survival::Surv(y1, y2, type="interval2"),
                 X = X,
                 alpha = foo$cpar,
                 theta = foo$theta,
                 phi = foo$phi,
                 y = foo$y,
                 maxL = maxL,
                 ratec = foo$ratec,
                 ratetheta = foo$ratetheta,
                 ratephi = foo$ratephi,
                 ratey = foo$ratey,
                 initial.values=initial.values,
                 BF = BF);
  class(output) <- c("SpatDensReg")
  output
}

#### empirial BF and p-value for the spatial model vs. the exchangeable model
"BF.SpatDensReg" <- function (y, X, prior=NULL, nperm=100, c_seq=NULL, phi_seq=NULL) {
  n = length(y);
  X = cbind(X);
  Sinv = solve(var(X));
  # find the maximumu M distance between X[i,] and colMeans(X)
  distseq = rep(0, n);
  Xbar = colMeans(X);
  for(i in 1:n) distseq[i] = sqrt(as.vector((X[i,]-Xbar)%*%Sinv%*%(X[i,]-Xbar)))
  maxdist = max(distseq)
  phi0 = (-log(0.001))/maxdist;
  #########################################################################################
  # initial MLE analysis
  #########################################################################################
  fit0 <- survival::survreg(formula = survival::Surv(y)~1, dist="gaussian");
  theta1 = fit0$coefficients[1];
  theta2 = log(fit0$scale);
  theta = c(theta1, theta2); 
  
  #########################################################################################
  # priors and initial values
  #########################################################################################
  maxL <- prior$maxL; if(is.null(maxL)) maxL<-5;
  a0=prior$a0; if(is.null(a0)) a0=5;
  b0=prior$b0; if(is.null(b0)) b0=1;
  phiq0 = prior$phiq0; if(is.null(phiq0)) phiq0=0.5;
  phia0 = prior$phia0; if(is.null(phia0)) phia0=2;
  phib0 = prior$phib0; if(is.null(phib0)) phib0=1/phi0;
  if(is.null(c_seq)) c_seq=c(0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 1000);
  if(is.null(phi_seq)) phi_seq = qgamma((1:10)/11, phia0, phib0)
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  foo <- .Call("SpatDens_BF", y_=y, X_=t(X), Sinv_=Sinv, theta_=theta, maxJ_=maxL, 
               cpar_=c_seq, a0_=a0, b0_=b0, phi_=phi_seq, q0phi_=phiq0, 
               a0phi_=phia0, b0phi_=phib0, nperm_=nperm, PACKAGE = "spBayesSurv");
  
  ## Bayes Factor for the spatial model vs. the exchangeable model
  BF = foo$BF;
  pvalue = sum(foo$BFperm>foo$BF)/nperm;
  output <- list(BF = BF,
                 pvalue = pvalue);
  output
}

#### print, summary, plot
"print.SpatDensReg" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  cat(paste("\nBayes Factor for the spatial model vs. the exchangeable model:", sep=""), x$BF);
  cat("\nn=",x$n, "\n", sep="")
  invisible(x)
}

"plot.SpatDensReg" <- function (x, xnewdata, ygrid=NULL, CI=0.95, PLOT=TRUE, ...) {
  if(is.null(ygrid)) ygrid = seq(min(x$Surv[,1], na.rm=T)-sd(x$Surv[,1], na.rm=T), 
                                 max(x$Surv[,2], na.rm=T)+sd(x$Surv[,2], na.rm=T), length.out=200);
  if(missing(xnewdata)){
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
      if (is.R()) attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
    } else  newTerms <- Terms
    newTerms <- delete.response(newTerms)
    mnew <- model.frame(newTerms, xnewdata, na.action = na.omit, xlev = .getXlevels(newTerms, m))
    Xnew <- model.matrix(newTerms, mnew);
    if (is.R()) {
      assign <- lapply(survival::attrassign(Xnew, newTerms)[-1], function(x) x-1)
      xlevels <- .getXlevels(newTerms, mnew)
      contr.save <- attr(Xnew, 'contrasts')
    }else {
      assign <- lapply(attr(Xnew, 'assign')[-1], function(x) x -1)
      xvars <- as.character(attr(newTerms, 'variables'))
      xvars <- xvars[-attr(newTerms, 'response')]
      if (length(xvars) >0) {
        xlevels <- lapply(mnew[xvars], levels)
        xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
        if(length(xlevels) == 0)
          xlevels <- NULL
      } else xlevels <- NULL
      contr.save <- attr(Xnew, 'contrasts')
    }
    # drop the intercept after the fact, and also drop baseline if necessary
    adrop <- 0  #levels of "assign" to be dropped; 0= intercept
    Xatt <- attributes(Xnew) 
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    Xnew <- Xnew[, !xdrop, drop=FALSE]
    attr(Xnew, "assign") <- Xatt$assign[!xdrop]
    xpred = Xnew
    if(ncol(xpred)!=x$p) stop("please make sure the number of columns matches!");
  }
  xpred = cbind(xpred);
  nxpred = nrow(xpred);
  Sinv = solve(var(x$X));
  estimates <- .Call("SpatDens_plots", ygrid, t(xpred), x$theta, x$alpha, x$phi, x$maxL,
                     x$y, t(x$X), Sinv, CI, PACKAGE = "spBayesSurv");
  if(PLOT){
    par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
    plot(ygrid, estimates$fhat[,1], "l", lwd=3, xlab="log time", ylab="density", 
         xlim=c(min(ygrid), max(ygrid)), ylim=c(0,max(estimates$fhatup)));
    for(i in 1:nxpred){
      polygon(x=c(rev(ygrid),ygrid),
              y=c(rev(estimates$fhatlow[,i]),estimates$fhatup[,i]),
              border=NA,col="lightgray");
    }
    for(i in 1:nxpred){
      lines(ygrid, estimates$fhat[,i], lty=i, lwd=3, col=i);
    }
    legend("topright", rnames, col = 1:nxpred, lty=1:nxpred, ...)
  }
  estimates$ygrid=ygrid;
  invisible(estimates)
}

"summary.SpatDensReg" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### Baseline Information
  mat <- as.matrix(object$theta)
  coef.p <- apply(mat, 1, mean); names(coef.p)=c("location", "log(scale)");
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$theta.var <- coef.table
  
  ### Precision parameter
  if(object$prior$a0<=0){
    ans$alpha.var <- NULL
  }else{
    mat <- object$alpha
    coef.p <- mean(mat); names(coef.p)="alpha";    
    coef.m <- median(mat)    
    coef.sd <- sd(mat)
    limm <- as.vector(quantile(mat, probs=c((1-CI.level)/2, 1-(1-CI.level)/2)))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%CI-Low", sep=""),
                                                  paste(CI.level*100, "%CI-Upp", sep="")))
    ans$alpha.var <- coef.table
  }
  
  ### phi parameter
  mat <- object$phi
  coef.p <- mean(mat); names(coef.p)="range";
  coef.m <- median(mat)    
  coef.sd <- sd(mat)
  limm <- as.vector(quantile(mat, probs=c((1-CI.level)/2, 1-(1-CI.level)/2)))
  coef.l <- limm[1]
  coef.u <- limm[2]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$phi.var <- coef.table;
  ans$BF <- object$BF
  ans$n <- object$n
  ans$p <- object$p
  ans$prior <- object$prior
  
  ### acceptance rates
  ans$ratetheta = object$ratetheta;
  ans$ratephi = object$ratephi;
  ans$ratey = object$ratey;
  ans$ratec = object$ratec;
  
  class(ans) <- "summary.SpatDensReg"
  return(ans)
}

"print.summary.SpatDensReg"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  if(x$theta.var[1,3]==0){
    cat("\nCentering distribution parameters are fixed at:\n")
    cat("location=", x$theta.var[1,1], ", log(scale)=", x$theta.var[2,1], "\n", sep="")
  }else{
    cat("\nPosterior inference of centering distribution parameters\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratetheta, "):\n", sep="")
    print.default(format(x$theta.var, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  
  if (!is.null(x$alpha.var)) {
    cat("\nPosterior inference of precision parameter\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratec, "):\n", sep="")
    print.default(format(x$alpha.var, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nPosterior inference of distance function range phi\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratephi, "):\n", sep="")
  print.default(format(x$phi.var, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat(paste("\nBayes Factor for the spatial model vs. the exchangeable model:", sep=""), x$BF)   
  
  cat("\nNumber of subjects: n=", x$n, "\n", sep="")

  invisible(x)
}


