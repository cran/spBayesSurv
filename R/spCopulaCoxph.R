"spCopulaCoxph" <- function(formula, data, na.action, prediction=NULL, 
                            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                            prior=NULL, state=NULL, scale.designX=TRUE,
                            Coordinates, DIST=NULL, Knots=NULL) {
  #########################################################################################
  # call parameters
  #########################################################################################
  Call <- match.call(); # save a copy of the call 
  indx <- match(c("formula", "data", "na.action", "truncation_time", "subject.num"),
                names(Call), nomatch=0) 
  if (indx[1] ==0) stop("A formula argument is required");
  temp <- Call[c(1,indx)]  # only keep the arguments we wanted
  temp[[1]] <- as.name('model.frame')  # change the function called
  
  special <- c("baseline", "frailtyprior", "truncation_time", "subject.num")
  temp$formula <- if(missing(data)) terms(formula, special)
  else              terms(formula, special, data=data)
  if (is.R()) m <- eval(temp, parent.frame())
  else        m <- eval(temp, sys.parent())
  
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
  
  Terms <- attr(m, 'terms')
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
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
    stop("covariate is required")
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
  #### Distance type:
  if(is.null(DIST)){
    DIST <- function(x, y) fields::rdist(x, y)
  }
  
  ##############################################
  ### Currently it only supports PH for right censored data 
  ##############################################
  model.name <- "Spatial Copula Cox PH model with piecewise constant baseline hazards"
  if(sum(delta%in%c(0,1))!=n) stop("This function currently only supports right-censored data.")
  
  #########################################################################################
  # prediction
  #########################################################################################
  if(p==0){
    s0 <- prediction$spred;
    if(ncol(s0)!=2) stop("Make sure that prediction$spred is a matrix with two columns.")
    if(is.null(s0)) s0=Coordinates
    npred = nrow(s0);
    xpred = cbind(rep(1,npred));
  }else{
    xpred <- prediction$xpred;
    s0 <- prediction$spred;
    if(is.null(xpred)){
      xpred = X;
      s0 = Coordinates;
    }
    if(is.vector(xpred)) xpred=matrix(xpred, nrow=1);
    if(ncol(xpred)!=p) stop("Please make sure the number of columns for xpred equals the number of covariates!");
    xpred = cbind(xpred);
    npred = nrow(xpred);
    if(nrow(s0)!=npred) stop("Error: nrow(xpred) is not equal to nrow(spred) in prediction");
    for(i in 1:npred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
  }
  
  #########################################################################################
  # general setup based on copula priors
  #########################################################################################
  s = Coordinates;
  if(is.null(s)) stop("please specify Coordinates for each subject");
  if(nrow(s)!=n) stop("the number of coordinates should be equal to the sample size")
  dnn = DIST(s, s);
  if(min(dnn[row(dnn)!=col(dnn)])<=0) stop("each subject should have different Coordinates");
  if(is.null(Knots)){
    nknots = prior$nknots; if(is.null(nknots)) nknots=n;
  }else{
    nknots = nrow(Knots);
  }
  if(is.null(Knots)){
    if(nknots<n){
      ss = as.matrix(fields::cover.design(s, nd=nknots, DIST=DIST)$design);
    }else{
      ss = s;
    }
  }else{
    ss = Knots;
  }
  dnm = DIST(s, ss); 
  dmm = DIST(ss, ss);
  nblock=prior$nblock; if(is.null(nblock)) nblock=n;
  if(nblock==n){
    s0tmp = s;
    Dtmp = DIST(s, s0tmp);
    idtmp = apply(Dtmp, 1, which.min);
    clustindx=diag(1,n, n);
  }else{
    s0tmp = as.matrix(fields::cover.design(s, nd=nblock, DIST=DIST)$design);
    Dtmp = DIST(s, s0tmp);
    idtmp = apply(Dtmp, 1, which.min);
    nblock=length(table(idtmp));
    idnames = as.numeric(names(table(idtmp)))
    clustindx=matrix(0, n, nblock);
    for(jj in 1:nblock){
      clustindx[which(idtmp==idnames[jj]),jj] = 1;
    }
  }
  ds0n <- DIST(s, s0);
  ds0m <- DIST(ss, s0);
  Ds0tmps0 <- DIST(s0, s0tmp);
  idpred <- apply(Ds0tmps0, 1, which.min);
  ds0block <- matrix(0, n, npred);
  for(i in 1:n){
    for(j in 1:npred){
      ds0block[i,j] = (idtmp[i]==idpred[j])+0
    }
  }
  
  #########################################################################################
  # initial analysis and mcmc parameters
  #########################################################################################
  tbase1 = t1; tbase2 = t2; deltabase = delta;
  Xbase.scaled = X.scaled;
  for(i in 1:n){
    if(deltabase[i]==0) tbase2[i]=NA;
    if(deltabase[i]==2) tbase1[i]=NA;
  }
  fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, 
                            dist="exponential");
  
  #########################################################################################
  # priors
  #########################################################################################
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  r0 <- prior$r0; if(is.null(r0)) r0 = 1;
  h0 <- prior$h0; if(is.null(h0)) h0 = as.vector( exp( -fit0$coefficients[1] ) );
  v0 <- prior$v0; if(is.null(v0)) v0=0;#v0 = 10*as.vector( exp( -2*fit0$coefficients[1] )*fit0$var[1,1] );
  vhat <- prior$vhat; if(is.null(vhat)) vhat <- as.vector( exp( -2*fit0$coefficients[1] )*fit0$var[1,1] );
  beta0 <- prior$beta0; if(is.null(beta0)) beta0 <- rep(0,p);
  S0 <- prior$S0; if(is.null(S0)) S0=diag(1e10, p);
  S0inv <- solve(S0);
  Shat <- prior$Shat; if(is.null(Shat)) Shat <- as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))])/(fit0$scale)^2;
  M <- prior$M; if(is.null(M)) M <- 20;
  M1<- M+1;
  d <- prior$cutpoints; 
  if(is.null(d)){
    d = as.vector(quantile(t1, probs=seq(0,1,length=M1)));
    d = d[-1];
    d[M] = Inf;
  }
  d <- c(0, d);
  if(!(M1==length(d))) stop("error: M is not equal to length(cutpoints)");
  theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- c(1.0, 1.0, 1.0, 1.0);
  spS0 <- prior$spS0; if(is.null(spS0)) spS0 <- diag(c(0.5,0.1));
  
  #########################################################################################
  # current state and mcmc specification
  #########################################################################################
  h = c(0, rep(h0, M)); 
  hcen=state$hcen; if(is.null(hcen)) hcen=h0;
  beta=state$beta; if(is.null(beta)) beta=as.vector( -fit0$coefficients[-1] );
  theta = state$theta; if(is.null(theta)) theta <- c(0.98, 1);
  
  #########################################################################################
  # calling the c++ code
  #########################################################################################
  foo <- .Call("spCopulaCoxph", nburn_ = nburn, nsave_ = nsave, nskip_ = nskip, ndisplay_ = ndisplay,
               tobs_ = t1, delta_ = delta, X_=X.scaled, d_ = d, h_ = h, r0_ = r0, hcen_=hcen,
               h0_ = h0, v0_ = v0, vhat_ = vhat, beta_ = beta, beta0_ = beta0, S0inv_ = S0inv,
               Shat_=Shat, l0_ = round(min(5000,nsave/2)), adapter_ = (2.38)^2, 
               xpred_ = as.matrix(xpred), ds0n_=ds0n, dnn_=dnn, theta_=theta, theta0_=theta0,
               spS0_=spS0, dnm_=dnm, dmm_=dmm, clustindx_=clustindx, ds0m_=ds0m, 
               ds0block_=ds0block, PACKAGE = "spBayesSurv");
  #### transfer the estimates back to original scales;
  beta.scaled = matrix(foo$beta, p, nsave);
  beta.original = matrix(beta.scaled, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
  
  #### coefficients
  coeff1 <- c(apply(beta.original, 1, mean));
  coeff2 <- c(apply(foo$theta, 1, mean));
  coeff <- c(coeff1, coeff2);
  names(coeff) = c(colnames(X.scaled), "sill", "range");
  #### Save to a list
  output <- list(modelname=model.name,
                 coefficients=coeff,
                 call=Call,
                 prior=prior,
                 mcmc=mcmc,
                 n=n,
                 p=p,
                 Surv=survival::Surv(tbase1, tbase2, type="interval2"),
                 X.scaled=X.scaled,
                 X = X,
                 #survt = foo$t,
                 beta = beta.original,
                 beta.scaled = beta.scaled,
                 h.scaled = foo$h,
                 d.scaled = foo$d,
                 cutpoints = foo$d[,1],
                 #hcen.scaled = foo$hcen,
                 M=M,
                 ratebeta = foo$ratebeta,
                 #ratehcen = foo$ratehcen,
                 theta = foo$theta,
                 ratebeta = foo$ratebeta,
                 ratetheta = foo$ratetheta,
                 rateh = foo$rateh,
                 cpo = foo$cpo,
                 Coordinates = s,
                 Tpred = foo$Tpred,
                 Zpred = foo$Zpred);
  class(output) <- c("spCopulaCoxph")
  output
}

#### print, summary, plot
"print.spCopulaCoxph" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior means for regression coefficients:\n")
  if(x$p>0){
    print.default(format(x$coefficients[1:x$p], digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nLPML:", sum(log(x$cpo)))
  cat("\nn=",x$n, "\n", sep="")
  invisible(x)
}

"plot.spCopulaCoxph" <- function (x, xpred=NULL, tgrid=NULL, CI=0.95, PLOT=FALSE, ...) {
  if(is(x,"spCopulaCoxph")){
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

"summary.spCopulaCoxph" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### CPO
  ans$cpo <- object$cpo
  
  ### Median information
  mat <- as.matrix(object$beta)
  coef.p <- object$coefficients[(1:object$p)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$coeff <- coef.table
  
  ### Spatial Copula sill and range parameters
  mat <- as.matrix(object$theta)
  coef.p <- object$coefficients[-(1:object$p)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$sill.range <- coef.table
  
  ## LPML and DIC
  ans$n <- object$n
  ans$p <- object$p
  ans$LPML <- sum(log(object$cpo))
  
  ### acceptance rates
  ans$ratebeta = object$ratebeta;
  ans$ratetheta = object$ratetheta;
  class(ans) <- "summary.spCopulaCoxph"
  return(ans)
}

"print.summary.spCopulaCoxph"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  if(x$p>0){
    cat("\nPosterior inference of regression coefficients\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratebeta, "):\n", sep="")
    print.default(format(x$coeff, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nPosterior inference of spatial sill and range parameters\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratetheta, "):\n", sep="")
  print.default(format(x$sill.range, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nLog pseudo marginal likelihood: LPML=", x$LPML, sep="")
  cat("\nNumber of subjects: n=", x$n, "\n", sep="")
  invisible(x)
}

