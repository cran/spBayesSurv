"SuperSurvRegBayes" <- function (formula, data, na.action, dist="lognormal", 
                                  mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500),
                                  prior=NULL, state=NULL, truncation_time=NULL, subject.num=NULL, 
                                  InitParamMCMC=FALSE) {
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
    #X.scaled <- scale(X, center=rep(0,p), scale=rep(1,p));
    X.scaled <- scale(X);
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
  if (type== 'counting'){ #stop ("Invalid survival type")
    t1 = Y[,2]; t2 = Y[,2];
    truncation_time = as.vector(Y[,1]);
  } 
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
  if(is.null(truncation_time)) truncation_time=rep(0, n);
  if(is.null(subject.num)) subject.num=(1:n);
  distcode = switch(dist, loglogistic=1, lognormal=2, 3);
  frail.prior = colnames(frail.terms)[1];
  #### if subject.num has been ordered if time-dependent covariates are considered
  if(!(sum(order(subject.num)==(1:n))==n)){
    message("Please sort data by subject.num; otherwise, the LPML will not be adjusted for time-dependent covariates.");
  }
  subjecti = c(0, cumsum(as.vector(table(subject.num))));
  nsubject = length(subjecti)-1;
  baseindx = (subjecti+1)[-length(subjecti)];
  lastindx = subjecti[-1]
  
  #########################################################################################
  # initial MLE analysis and mcmc parameters
  #########################################################################################
  tbase1 = t1; tbase2 = t2; deltabase = delta;
  Xbase.scaled = X.scaled;
  for(i in 1:n){
    if(deltabase[i]==0) tbase2[i]=NA;
    if(deltabase[i]==2) tbase1[i]=NA;
  }
  ## initial fit for theta:
  fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, dist=dist);
  theta1 = -fit0$coefficients[1];
  theta2 = -log(fit0$scale);
  theta = c(theta1, theta2); theta_prior = c(theta1, theta2);
  Vhat0 = as.matrix(fit0$var[c(1,p+2),c(1,p+2)]);
  fit1 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, dist="weibull");
  fit2 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, dist="loglogistic");
  if(AIC(fit1)>AIC(fit2)){
    beta_h = rep(0, p); Shat0_h=as.matrix(fit1$var[c(2:(1+p)),c(2:(1+p))])/(fit1$scale)^2;
    beta_q = rep(0, p); Shat0_q=as.matrix(fit2$var[c(2:(1+p)),c(2:(1+p))]);
    beta_o = -fit2$coefficients[-1]/fit2$scale; 
    Shat0_o = as.matrix(fit2$var[c(2:(1+p)),c(2:(1+p))])/(fit2$scale)^2;
    beta = c(beta_h, beta_o, beta_q);
    Shat0 = matrix(0, 3*p, 3*p);
    Shat0[(1:p), (1:p)] = Shat0_h; 
    Shat0[(1:p)+p, (1:p)+p] = Shat0_o;
    Shat0[(1:p)+2*p, (1:p)+2*p] = Shat0_q;
  }else{
    beta_h = -fit1$coefficients[-1]; 
    Shat0_h = as.matrix(fit1$var[c(2:(1+p)),c(2:(1+p))]);
    beta_q = -fit1$coefficients[-1]; 
    Shat0_q = as.matrix(fit1$var[c(2:(1+p)),c(2:(1+p))]);
    beta_o = rep(0, p); Shat0_o=as.matrix(fit2$var[c(2:(1+p)),c(2:(1+p))])/(fit2$scale)^2;
    beta = c(beta_h, beta_o, beta_q);
    Shat0 = matrix(0, 3*p, 3*p);
    Shat0[(1:p), (1:p)] = Shat0_h; 
    Shat0[(1:p)+p, (1:p)+p] = Shat0_o;
    Shat0[(1:p)+2*p, (1:p)+2*p] = Shat0_q;
  }
  if(InitParamMCMC){
    S0param0 = 3.228191/p*n*solve(t(X.scaled)%*%X.scaled);
    S0param = matrix(0, 3*p, 3*p);
    S0param[(1:p), (1:p)] = S0param0; 
    S0param[(1:p)+p, (1:p)+p] = S0param0;
    S0param[(1:p)+2*p, (1:p)+2*p] = S0param0;
    ## initial fit for beta
    nburn0=5000; nsave0=2000;
    fit0<- .Call("PHPOAFT_BP", nburn_=nburn0, nsave_=nsave0, nskip_=1, ndisplay_=10000, 
                 ltr_=truncation_time, subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, 
                 beta_=beta, weight_=c(1), cpar_=Inf, a0_=-1, b0_=1, theta0_=theta_prior, 
                 V0inv_=solve(10*Vhat0), Vhat_=Vhat0, beta0_=rep(0,3*p), S0inv_=solve(S0param), 
                 Shat_=Shat0, l0_=3000, adapter_=2.38^2, dist_=distcode, PACKAGE = "spBayesSurv");
    beta_h = colMeans(t(matrix(fit0$beta_h, p, nsave0)));
    beta_o = colMeans(t(matrix(fit0$beta_o, p, nsave0)));
    beta_q = colMeans(t(matrix(fit0$beta_q, p, nsave0)));
    beta = c(beta_h, beta_o, beta_q);
    beta_prior = c(beta_h, beta_o, beta_q);
    Shat0 = cov(t(rbind(fit0$beta_h, fit0$beta_o, fit0$beta_q)));
    theta = colMeans(t(matrix(fit0$theta, 2, nsave0)));
    theta_prior = colMeans(t(matrix(fit0$theta, 2, nsave0)));
    Vhat0 = cov(t(matrix(fit0$theta, 2, nsave0)));
  }
  
  #########################################################################################
  # priors
  # note the priors should be based on scaled data.
  #########################################################################################
  cpar=state$cpar; if(is.null(cpar)) cpar=1;
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  maxL <- prior$maxL; if(is.null(maxL)) maxL<-15;
  weight = rep(1/maxL, maxL);
  a0=prior$a0; if(is.null(a0)) a0=1;
  b0=prior$b0; if(is.null(b0)) b0=1;
  beta0 <- prior$beta0; 
  M=prior$M; if(is.null(M)) M=10;
  q=prior$q; if(is.null(q)) q=0.9;
  if(is.null(beta0)){
    beta0 <- rep(0,3*p);
  }
  S0 <- prior$S0; 
  if(is.null(S0)) {
    gg = (log(M)/qnorm(q))^2/p;
    S0_0 = gg*n*solve(t(X.scaled)%*%X.scaled);
    S0 = matrix(0, 3*p, 3*p);
    S0[(1:p), (1:p)] = S0_0; 
    S0[(1:p)+p, (1:p)+p] = S0_0;
    S0[(1:p)+2*p, (1:p)+2*p] = S0_0;
  }
  S0inv <- solve(S0);
  theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- theta_prior;
  V0 <- prior$V0; if(is.null(V0)) V0 <- 10*Vhat0;
  if(sum(abs(V0))==0){
    V0inv <- diag(c(Inf,Inf));
  }else {
    V0inv <- solve(V0);
  }
  Vhat <- prior$Vhat; if(is.null(Vhat)) Vhat <- Vhat0;
  Shat <- prior$Shat; if(is.null(Shat)) Shat <- Shat0;
  mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
  prior = list(maxL=maxL, a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0,
               Shat=Shat, Vhat=Vhat, M=M, q=q);
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  model.name <- "Super Survival Model:";
  foo <- .Call("PHPOAFT_BP", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, 
               ltr_=truncation_time, subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, 
               beta_=beta, weight_=weight, cpar_=cpar, a0_=a0, b0_=b0, theta0_=theta0, 
               V0inv_=V0inv, Vhat_=Vhat, beta0_=beta0, S0inv_=S0inv, Shat_=Shat,
               l0_=min(5000,nsave/2), adapter_=2.38^2, dist_=distcode, PACKAGE = "spBayesSurv");
  
  #########################################################################################
  # save state
  #########################################################################################
  #### transfer the estimates back to original scales;
  theta.scaled = foo$theta;
  theta.original = foo$theta;
  beta_h.scaled = matrix(foo$beta_h, p, nsave);
  beta_o.scaled = matrix(foo$beta_o, p, nsave);
  beta_q.scaled = matrix(foo$beta_q, p, nsave);
  beta_h.original = matrix(beta_h.scaled, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
  beta_o.original = matrix(beta_o.scaled, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
  beta_q.original = matrix(beta_q.scaled, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
  
  #### coefficients
  coeff1h <- c(apply(beta_h.original, 1, mean));
  coeff1o <- c(apply(beta_o.original, 1, mean));
  coeff1q <- c(apply(beta_q.original, 1, mean));
  coeff2 <- c(apply(theta.scaled, 1, mean));
  coeff <- c(coeff1h, coeff1o, coeff1q, coeff2);
  names(coeff) = c(colnames(X.scaled), colnames(X.scaled), colnames(X.scaled), "theta1", "theta2");
  
  #### Save to a list
  output <- list(modelname=model.name,
                 dist = dist,
                 coefficients=coeff,
                 call=Call,
                 prior=prior,
                 mcmc=mcmc,
                 n=n,
                 p=p,
                 nsubject=nsubject,
                 subject.num=subject.num,
                 truncation_time=truncation_time,
                 Surv=survival::Surv(tbase1, tbase2, type="interval2"),
                 X.scaled=X.scaled,
                 X = X,
                 beta_h = beta_h.original,
                 beta_o = beta_o.original,
                 beta_q = beta_q.original,
                 beta_h.scaled = beta_h.scaled,
                 beta_o.scaled = beta_o.scaled,
                 beta_q.scaled = beta_q.scaled,
                 theta.scaled = theta.scaled,
                 cpar = foo$cpar,
                 maxL = maxL,
                 weight = foo$weight,
                 cpo = foo$cpo,
                 pD = foo$pD, 
                 DIC = foo$DIC,
                 ratetheta = foo$ratetheta,
                 ratebeta_h = foo$ratebeta,
                 ratebeta_o = foo$ratebeta,
                 ratebeta_q = foo$ratebeta,
                 rateYs = foo$rateYs,
                 ratec = foo$ratec);
  #########################################################################################
  # Calculate Bayes Factors;
  #########################################################################################
  ## PH:
  HH = rbind(foo$beta_q, foo$beta_o-foo$beta_h); 
  meanH = apply(HH, 1, mean);
  varH = var(t(HH));
  numH = (2*pi)^(-p/2)*(det(S0))^(-1/2)*(2*pi)^(-p/2)*(det(2*S0))^(-1/2);
  denH = as.vector((2*pi)^(-p)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BF_h = denH/numH; 
  ## AFT:
  HH = rbind(foo$beta_o, foo$beta_h-foo$beta_q); 
  meanH = apply(HH, 1, mean);
  varH = var(t(HH));
  numH = (2*pi)^(-p/2)*(det(S0))^(-1/2)*(2*pi)^(-p/2)*(det(2*S0))^(-1/2);
  denH = as.vector((2*pi)^(-p)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BF_q = denH/numH; 
  ## PO:
  HH = rbind(foo$beta_q, foo$beta_h); 
  meanH = apply(HH, 1, mean);
  varH = var(t(HH));
  numH = ((2*pi)^(-p/2)*(det(S0))^(-1/2))^2;
  denH = as.vector((2*pi)^(-p)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BF_o = denH/numH; 
  ## AH:
  HH = rbind(foo$beta_h, foo$beta_q+foo$beta_o); 
  meanH = apply(HH, 1, mean);
  varH = var(t(HH));
  numH = (2*pi)^(-p/2)*(det(S0))^(-1/2)*(2*pi)^(-p/2)*(det(2*S0))^(-1/2);
  denH = as.vector((2*pi)^(-p)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BF_a = denH/numH; 
  ## Yang and Prentice:
  HH = rbind(foo$beta_q); 
  meanH = apply(HH, 1, mean);
  varH = as.matrix(var(t(HH)), p, p);
  numH = ((2*pi)^(-p/2)*(det(S0))^(-1/2));
  denH = as.vector((2*pi)^(-p/2)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BF_y = denH/numH; 
  ## EH:
  HH = rbind(foo$beta_h-foo$beta_q-foo$beta_o); 
  meanH = apply(HH, 1, mean);
  varH = as.matrix(var(t(HH)), p, p);
  numH = ((2*pi)^(-p/2)*(det(3*S0))^(-1/2));
  denH = as.vector((2*pi)^(-p/2)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BF_e = denH/numH; 
  ## save:
  BayesFactors = c(BF_q, BF_h, BF_o, BF_a, BF_e, BF_y); 
  names(BayesFactors) = c( "AFT", "PH", "PO", "AH", "EH", "YP");
  output$BF = BayesFactors; 
  class(output) <- c("SuperSurvRegBayes")
  output
}

#### print, summary, plot
"print.SuperSurvRegBayes" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior means for regression coefficients:\n")
  if(x$p>0){
    print.default(format(x$coefficients[1:(3*(x$p))], digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nBayes factors for testing AFT, PH and PO:\n")       
  print.default(format(x$BF, digits = digits), print.gap = 2, 
                quote = FALSE) 
  
  cat("\nLPML:", sum(log(x$cpo)))
  cat("\nDIC:", x$DIC)
  cat("\nn=",x$n, "\n", sep="")
  invisible(x)
}

"plot.SuperSurvRegBayes" <- function (x, xpred=NULL, tgrid=NULL, CI=0.95, PLOT=FALSE, ...) {
  if(is(x,"SuperSurvRegBayes")){
    if(is.null(tgrid)) tgrid = seq(0.01, max(x$Surv[,2], na.rm=T), length.out=200);
    if(is.null(xpred)) {
      stop("please specify xpred")
    }else{
      if(is.vector(xpred)) xpred=matrix(xpred, nrow=1);
      if(ncol(xpred)!=x$p) stop("please make sure the number of columns matches!");
    };
    X.center = attributes(x$X.scaled)$`scaled:center`;
    X.scale = attributes(x$X.scaled)$`scaled:scale`;
    xpred = cbind(xpred);
    nxpred = nrow(xpred);
    for(i in 1:nxpred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
    distcode = switch(x$dist, loglogistic=1, lognormal=2, 3);
    estimates <- .Call("PHPOAFT_BP_plots", tgrid, xpred, x$theta.scaled, 
                       x$beta_h.scaled, x$beta_o.scaled, x$beta_q.scaled, x$weight, CI, 
                       dist_=distcode, PACKAGE = "spBayesSurv");
    if(PLOT){
      par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
      plot(tgrid, estimates$Shat[,1], "l", lwd=3, xlab="time", ylab="survival", main=paste(i));
      for(i in 1:nxpred){
        polygon(x=c(rev(tgrid),tgrid),
                y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
                border=NA,col="lightgray");
      }
      for(i in 1:nxpred){
        lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
      }
    }
  }
  estimates$tgrid=tgrid;
  invisible(estimates)
}

"summary.SuperSurvRegBayes" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### CPO
  ans$cpo <- object$cpo
  
  ### Median information
  mat <- as.matrix(object$beta_h)
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
  ans$coeff_h <- coef.table
  mat <- as.matrix(object$beta_o)
  coef.p <- object$coefficients[(object$p+1):(2*object$p)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$coeff_o <- coef.table
  mat <- as.matrix(object$beta_q)
  coef.p <- object$coefficients[(1:object$p)+2*object$p];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$coeff_q <- coef.table
  
  ### Baseline Information
  mat <- as.matrix(object$theta.scaled)
  coef.p <- object$coefficients[-(1:(3*object$p))];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$basepar <- coef.table
  
  ### Precision parameter
  if(object$prior$a0<=0){
    ans$prec <- NULL
  }else{
    mat <- object$cpar
    coef.p <- mean(mat)    
    coef.m <- median(mat)    
    coef.sd <- sd(mat)
    limm <- as.vector(quantile(mat, probs=c((1-CI.level)/2, 1-(1-CI.level)/2)))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%CI-Low", sep=""),
                                                  paste(CI.level*100, "%CI-Upp", sep="")))
    ans$prec <- coef.table
  }
  
  ## LPML and DIC
  ans$n <- object$n
  ans$p <- object$p
  ans$LPML <- sum(log(object$cpo))
  ans$DIC <- object$DIC
  
  ### acceptance rates
  ans$ratetheta = object$ratetheta;
  ans$ratebeta_h = object$ratebeta_h;
  ans$ratebeta_o = object$ratebeta_o;
  ans$ratebeta_q = object$ratebeta_q;
  ans$rateYs = object$rateYs;
  ans$ratec = object$ratec;
  ans$cpar = object$cpar;
  ans$BF <- object$BF;
  
  class(ans) <- "summary.SuperSurvRegBayes"
  return(ans)
}


"print.summary.SuperSurvRegBayes"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  if(x$p>0){
    cat("\nPosterior inference of beta_h \n")
    cat("(Adaptive M-H acceptance rate: ", x$ratebeta_h, "):\n", sep="")
    print.default(format(x$coeff_h, digits = digits), print.gap = 2, 
                  quote = FALSE)
    cat("\nPosterior inference of beta_o \n")
    cat("(Adaptive M-H acceptance rate: ", x$ratebeta_o, "):\n", sep="")
    print.default(format(x$coeff_o, digits = digits), print.gap = 2, 
                  quote = FALSE)
    cat("\nPosterior inference of beta_q \n")
    cat("(Adaptive M-H acceptance rate: ", x$ratebeta_q, "):\n", sep="")
    print.default(format(x$coeff_q, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  if(x$cpar[1]==Inf){
    cat("\nPosterior inference of baseline parameters\n")
    message("Note: the baseline estimates are based on centered covariates")
    cat("(Adaptive M-H acceptance rate: ", x$ratetheta, "):\n", sep="")
    print.default(format(x$basepar, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  #cat("(Adaptive M-H acceptance rate for conditional probabilities: ", x$rateYs, ")\n", sep="")
  
  if (!is.null(x$prec)) {
    cat("\nPosterior inference of precision parameter\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratec, "):\n", sep="")
    print.default(format(x$prec, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nBayes factors relative to the super model:\n")       
  print.default(format(x$BF, digits = digits), print.gap = 2, 
                quote = FALSE) 
  
  cat("\nLog pseudo marginal likelihood: LPML=", x$LPML, sep="")
  cat("\nDeviance Information Criterion: DIC=", x$DIC, sep="")
  cat("\nNumber of subjects: n=", x$n, "\n", sep="")

  invisible(x)
}


