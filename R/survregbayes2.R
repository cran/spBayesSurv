"survregbayes2" <- function (formula, data, na.action, survmodel="PH", dist="loglogistic", 
                            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                            prior=NULL, state=NULL, selection=FALSE, Proximity=NULL, 
                            truncation_time=NULL, subject.num=NULL, InitParamMCMC=TRUE,
                            scale.designX=TRUE) {
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
  # check frailty
  #########################################################################################
  if(!is.null(frail.prior)){
    if(frail.prior=="car"){
      frailtyCode = 1;
      ID = frail.terms[,1];
      orderindex = order(ID);
      if(!(sum(orderindex==(1:n))==n)) stop("please sort the data by ID");
      blocki = c(0, cumsum(as.vector(table(ID))));
      if(is.null(Proximity)) stop("please specify prxoimity matrix");
      W = Proximity;
      D = rowSums(W);
      if (any(D==0)) stop("it seems that some region does not have any neighbers, which is not allowed, pelase check");
    }else if (frail.prior=="iid") {
      frailtyCode = 2;
      ID = frail.terms[,1];
      orderindex = order(ID);
      if(!(sum(orderindex==(1:n))==n)) stop("please sort the data by ID");
      blocki = c(0, cumsum(as.vector(table(ID))));
      W = matrix(0, length(blocki)-1, length(blocki)-1);
      D = rowSums(W);
    }else {
      stop("This function only supports non-frailty, CAR frailty and IID frailty models.")
    }
  } else {
    ID = NULL;
    frailtyCode = 0;
    blocki = c(0, n);
    W = matrix(1, length(blocki)-1, length(blocki)-1);
    D = rowSums(W);
  }
  if(is.null(state$frail)) {
    v <- rep(0, length(blocki)-1);
  } else {
    v <- state$frail; if(length(v)!=(length(blocki)-1)) stop("check the length of frail");
  }
  
  #########################################################################################
  # initial MLE analysis and mcmc parameters
  #########################################################################################
  tbase1 = t1; tbase2 = t2; deltabase = delta;
  Xbase.scaled = X.scaled;
  for(i in 1:n){
    if(deltabase[i]==0) tbase2[i]=NA;
    if(deltabase[i]==2) tbase1[i]=NA;
  }
  ## initial MCMC
  if(InitParamMCMC){
    ## initial fit for theta:
    fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, dist=dist);
    theta1 = -fit0$coefficients[1];
    theta2 = -log(fit0$scale);
    theta0 = c(theta1, theta2); theta_prior = c(theta1, theta2);
    theta = c(theta1, theta2);
    Vhat0 = as.matrix(fit0$var[c(1,p+2),c(1,p+2)]);
    ## initial fit for beta;
    beta0 = -fit0$coefficients[-1]; beta_prior = -fit0$coefficients[-1];
    beta = -fit0$coefficients[-1];
    Shat0 = as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))]);
    if(survmodel=="PH"){
      fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, dist="weibull");
      beta0 = -fit0$coefficients[-1]/fit0$scale; beta_prior = -fit0$coefficients[-1]/fit0$scale;
      beta = -fit0$coefficients[-1]/fit0$scale;
      Shat0 = as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))])/(fit0$scale)^2;
    }else if (survmodel=="PO"){
      fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, dist="loglogistic");
      beta0 = -fit0$coefficients[-1]/fit0$scale; beta_prior = -fit0$coefficients[-1]/fit0$scale;
      beta = -fit0$coefficients[-1]/fit0$scale;
      Shat0 = as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))])/(fit0$scale)^2;
    }
    message("Starting initial MCMC based on parametric model:")
    nburn0=5000; nsave0=5000;
    scaleV = 10; scaleS = 100000;
    if(sum(as.numeric(names(table(delta))))==2){ ## if current status data
      scaleV= 1; scaleS = 100000;
    }
    if(survmodel=="AFT"){
      fit0<- .Call("AFT_MPT", nburn_=nburn0, nsave_=nsave0, nskip_=0, ndisplay_=1000, ltr_=truncation_time,
                   subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta0, beta_=beta0, cpar_=Inf, 
                   Ys_=rep(0.5, 2^4-1), maxL_=4, a0_=-1, b0_=1, theta0_=theta_prior, V0inv_=solve(scaleV*Vhat0), 
                   Vhat_=Vhat0, beta0_=beta_prior, S0inv_=solve(scaleS*Shat0),  Shat_=Shat0, l0_=3000, adapter_=2.38^2, 
                   gamma_=rep(1, p), p0gamma_=rep(0.5, p), selection_=0, frailty_=frailtyCode, v_=v, 
                   blocki_=blocki, W_=W, lambda_=1, a0lambda_=1, b0lambda_=1, 
                   dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(survmodel=="PH"){
      fit0<- .Call("PH_MPT", nburn_=nburn0, nsave_=nsave0, nskip_=0, ndisplay_=1000, ltr_=truncation_time,
                   subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta0, beta_=beta0, cpar_=Inf, 
                   Ys_=rep(0.5, 2^4-1), maxL_=4, a0_=-1, b0_=1, theta0_=theta_prior, V0inv_=solve(scaleV*Vhat0), 
                   Vhat_=Vhat0, beta0_=beta_prior, S0inv_=solve(scaleS*Shat0),  Shat_=Shat0, l0_=3000, adapter_=2.38^2, 
                   gamma_=rep(1, p), p0gamma_=rep(0.5, p), selection_=0, frailty_=frailtyCode, v_=v, 
                   blocki_=blocki, W_=W, lambda_=1, a0lambda_=1, b0lambda_=1, 
                   dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(survmodel=="PO"){
      fit0<- .Call("PO_MPT", nburn_=nburn0, nsave_=nsave0, nskip_=0, ndisplay_=1000, ltr_=truncation_time,
                   subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta0, beta_=beta0, cpar_=Inf, 
                   Ys_=rep(0.5, 2^4-1), maxL_=4, a0_=-1, b0_=1, theta0_=theta_prior, V0inv_=solve(scaleV*Vhat0), 
                   Vhat_=Vhat0, beta0_=beta_prior, S0inv_=solve(scaleS*Shat0),  Shat_=Shat0, l0_=3000, adapter_=2.38^2, 
                   gamma_=rep(1, p), p0gamma_=rep(0.5, p), selection_=0, frailty_=frailtyCode, v_=v, 
                   blocki_=blocki, W_=W, lambda_=1, a0lambda_=1, b0lambda_=1, 
                   dist_=distcode, PACKAGE = "spBayesSurv");
    }else{
      stop("This function only supports PH, PO or AFT");
    }
    beta = colMeans(t(matrix(fit0$beta, p, nsave0)));
    beta_prior = colMeans(t(matrix(fit0$beta, p, nsave0)));
    Shat0 = cov(t(matrix(fit0$beta, p, nsave0)));
    theta = colMeans(t(matrix(fit0$theta, 2, nsave0)));
    theta_prior = colMeans(t(matrix(fit0$theta, 2, nsave0)));
    Vhat0 = cov(t(matrix(fit0$theta, 2, nsave0)));
    lambda = mean(fit0$lambda);
    if(!is.null(frail.prior)){
      if((frail.prior=="car")|(frail.prior=="iid")){
        if(is.null(state$frail)) v = rowMeans(fit0$v);
      }
    }
    message("Starting the MCMC for the semiparametric model:")
  }else{
    ## initial fit for theta:
    frailty.ID = ID;
    if(!is.null(frail.prior)){
      fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled+survival::frailty.gaussian(frailty.ID), 
                                dist=dist);
      v_initial = fit0$frail
    }else{
      fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, 
                                dist=dist);
      v_initial = rep(0, length(blocki)-1);
    }
    theta1 = -fit0$coefficients[1];
    theta2 = -log(fit0$scale);
    theta0 = c(theta1, theta2); theta_prior = c(theta1, theta2);
    theta = c(theta1, theta2);
    Vhat0 = as.matrix(fit0$var[c(1,p+2),c(1,p+2)]);
    ## initial fit for beta;
    beta0 = -fit0$coefficients[-1]; beta_prior = -fit0$coefficients[-1];
    beta = -fit0$coefficients[-1];
    Shat0 = as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))]);
    if((survmodel=="PH")|(survmodel=="PO")){
      if(survmodel=="PH") dist0 = "weibull";
      if(survmodel=="PO") dist0 = "loglogistic";
      frailty.ID = ID;
      if(!is.null(frail.prior)){
        fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled+survival::frailty.gaussian(frailty.ID), 
                                  dist=dist0);
        v_initial = fit0$frail
      }else{
        fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, 
                                  dist=dist0);
        v_initial = rep(0, length(blocki)-1);
      }
      v_initial = v_initial/fit0$scale;
      beta0 = -fit0$coefficients[-1]/fit0$scale; beta_prior = -fit0$coefficients[-1]/fit0$scale;
      beta = -fit0$coefficients[-1]/fit0$scale;
      Shat0 = as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))])/(fit0$scale)^2;
    }
    lambda=1;
    ## initial frailties
    if(!is.null(frail.prior)){
      lambda = 1/var(v_initial);
      if(is.null(state$frail)) v = v_initial;
    }
  }
  
  #########################################################################################
  # priors
  # note the priors should be based on scaled data.
  #########################################################################################
  alpha=state$alpha; if(is.null(alpha)) alpha=1;
  tau2 = state$tau2; if(is.null(tau2)) tau2=1/lambda; lambda=1/tau2; 
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  maxL <- prior$maxL; if(is.null(maxL)) maxL<-4;
  Ys = rep(0.5, 2^maxL-1);
  a0=prior$a0; if(is.null(a0)) a0=1;
  b0=prior$b0; if(is.null(b0)) b0=1;
  V0_prior = 10*Vhat0;
  S0_prior = diag(1e10, p);
  if(sum(as.numeric(names(table(delta))))==2){ ## if current status data
    V0_prior = 1*Vhat0;
  }
  beta0 <- prior$beta0; 
  if(is.null(beta0)){
    if(selection){
      beta0 <- rep(0,p);
    }else{
      beta0 <- rep(0,p);
    }
  }
  S0 <- prior$S0; 
  if(is.null(S0)) {
    if(selection){
      #M = quantile(as.vector(exp(X.scaled%*%beta)), probs=0.95);
      M = 10;
      gg = (log(M)/qnorm(0.9))^2/p;
      S0 = gg*n*solve(t(X.scaled)%*%X.scaled)
    }else{
      S0 <- S0_prior;
    }
  }
  S0inv <- solve(S0);
  theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- theta_prior;
  V0 <- prior$V0; if(is.null(V0)) V0 <- V0_prior;
  if(sum(abs(V0))==0){
    V0inv <- diag(c(Inf,Inf));
  }else {
    V0inv <- solve(V0);
  }
  Vhat <- prior$Vhat; if(is.null(Vhat)) Vhat <- Vhat0;
  Shat <- prior$Shat; if(is.null(Shat)) Shat <- Shat0;
  p0gamma = prior$p0gamma; if(is.null(p0gamma)) p0gamma=rep(0.5, p);
  gamma0 = rep(1.0, p);
  taua0 = prior$taua0; if(is.null(taua0)) taua0=1;
  taub0 = prior$taub0; if(is.null(taub0)) taub0=1;
  mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
  theta_initial=theta+0;
  beta_initial=beta+0;
  tau2_initial=1/lambda+0;
  frail_initial=v+0;
  if(!is.null(frail.prior)){
    initial.values = list(beta=beta_initial, theta=theta_initial, tau2=tau2_initial, frail=frail_initial);
    prior = list(maxL=maxL, a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0,
                 Shat=Shat, Vhat=Vhat, taua0=taua0, taub0=taub0);
  }else{
    initial.values = list(beta=beta_initial, theta=theta_initial);
    prior = list(maxL=maxL, a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0,
                 Shat=Shat, Vhat=Vhat);
  }
  if(selection){
    prior$p0gamma=p0gamma;
  }
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  if(survmodel=="AFT"){
    model.name <- "Accelerated failure time model:";
    foo <- .Call("AFT_MPT", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, ltr_=truncation_time, subjecti_=subjecti, 
                 t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=alpha, Ys_=Ys, maxL_=maxL,
                 a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, beta0_=beta0, S0inv_=S0inv, 
                 Shat_=Shat, l0_=min(5000,nsave/2), adapter_=2.38^2, gamma_=gamma0, p0gamma_=p0gamma, selection_=selection+0,
                 frailty_=frailtyCode, v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, 
                 dist_=distcode, PACKAGE = "spBayesSurv");
  }else if(survmodel=="PO"){
    model.name <- "Proportional Odds model:";
    foo <- .Call("PO_MPT", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, ltr_=truncation_time, subjecti_=subjecti, 
                 t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=alpha, Ys_=Ys, maxL_=maxL,
                 a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, beta0_=beta0, S0inv_=S0inv, 
                 Shat_=Shat, l0_=min(5000,nsave/2), adapter_=2.38^2, gamma_=gamma0, p0gamma_=p0gamma, selection_=selection+0,
                 frailty_=frailtyCode, v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, 
                 dist_=distcode, PACKAGE = "spBayesSurv");
  }else if(survmodel=="PH"){
    model.name <- "Proportional hazards model:";
    foo <- .Call("PH_MPT", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, ltr_=truncation_time, subjecti_=subjecti, 
                 t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, cpar_=alpha, Ys_=Ys, maxL_=maxL,
                 a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, beta0_=beta0, S0inv_=S0inv, 
                 Shat_=Shat, l0_=min(5000,nsave/2), adapter_=2.38^2, gamma_=gamma0, p0gamma_=p0gamma, selection_=selection+0,
                 frailty_=frailtyCode, v_=v, blocki_=blocki, W_=W, lambda_=lambda, a0lambda_=taua0, b0lambda_=taub0, 
                 dist_=distcode, PACKAGE = "spBayesSurv");
  }else{
    stop("This function only supports PH, PO or AFT");
  }
  #########################################################################################
  # save state
  #########################################################################################
  #### transfer the estimates back to original scales;
  theta.scaled = foo$theta;
  theta.original = foo$theta;
  beta.scaled = matrix(foo$beta, p, nsave);
  beta.original = matrix(beta.scaled, p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
  
  #### coefficients
  coeff1 <- c(apply(beta.original, 1, mean));
  coeff2 <- c(apply(theta.original, 1, mean));
  coeff <- c(coeff1, coeff2);
  names(coeff) = c(colnames(X.scaled),"theta1", "theta2");
  
  #### Save to a list
  output <- list(modelname=model.name,
                 dist = dist,
                 survmodel = survmodel,
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
                 beta = beta.original,
                 theta.scaled = theta.scaled,
                 beta.scaled = beta.scaled,
                 alpha = foo$cpar,
                 maxL = maxL,
                 Ys = foo$Ys,
                 cpo = foo$cpo,
                 pD = foo$pD, 
                 DIC = foo$DIC,
                 ratetheta = foo$ratetheta,
                 ratebeta = foo$ratebeta,
                 rateYs = foo$rateYs,
                 ratec = foo$ratec,
                 frail.prior=frail.prior,
                 selection = selection,
                 initial.values=initial.values);
  if(!is.null(frail.prior)){
    output$v = foo$v;
    output$ratev = foo$ratev;
    output$tau2 = 1/foo$lambda;
    output$ID = ID;
  }
  if(selection){
    output$gamma = foo$gamma;
  }
  class(output) <- c("survregbayes2")
  output
}

#### print, summary, plot
"print.survregbayes2" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior means for regression coefficients:\n")
  if(x$p>0){
    print.default(format(x$coefficients[1:x$p], digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nLPML:", sum(log(x$cpo)))
  cat("\nDIC:", x$DIC)
  cat("\nn=",x$n, "\n", sep="")
  
  if(x$selection){
    models=apply(x$gamma, 2, function(x) paste(x, sep="", collapse=",") )
    selected = table(models)/length(models);
    selected.max = selected[which.max(selected)];
    selected.ind = scan(text = names(selected)[which.max(selected)], what = 0L, sep=",", quiet = TRUE);
    if(sum(selected.ind)==0){
      cat("\nNo covariates are selected.\n")
    }else{
      selected.name = names(x$coefficients[1:x$p])[selected.ind==1];
      cat("\nThe selected covariates are :", selected.name, "\n")
    }
  }
  
  invisible(x)
}

"plot.survregbayes2" <- function (x, xpred=NULL, tgrid=NULL, CI=0.95, PLOT=FALSE, ...) {
  if(is(x,"survregbayes2")){
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
    betafitted = x$beta.scaled;
    if(x$selection){
      betafitted = x$beta.scaled*x$gamma;
    }
    if(x$survmodel=="AFT"){
      estimates <- .Call("AFT_MPT_plots", tgrid, xpred, x$theta.scaled, betafitted, x$Ys, x$maxL, CI, 
                         dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PO"){
      estimates <- .Call("PO_MPT_plots", tgrid, xpred, x$theta.scaled, betafitted, x$Ys, x$maxL, CI, 
                         dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PH"){
      estimates <- .Call("PH_MPT_plots", tgrid, xpred, x$theta.scaled, betafitted, x$Ys, x$maxL, CI, 
                         dist_=distcode, PACKAGE = "spBayesSurv");
    }else{
      stop("This function only supports PH, PO or AFT");
    }
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

"summary.survregbayes2" <- function(object, CI.level=0.95, ...) {
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
  
  ### Baseline Information
  mat <- as.matrix(object$theta.scaled)
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
  ans$basepar <- coef.table
  
  ### Precision parameter
  if(object$prior$a0<=0){
    ans$prec <- NULL
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
    ans$prec <- coef.table
  }
  
  ### frailty variance parameter
  ans$frail.prior=object$frail.prior;
  if(is.null(object$frail.prior)){
    ans$frailvar <- NULL
  }else{
    mat <- object$tau2
    coef.p <- mean(mat); names(coef.p)="variance";
    coef.m <- median(mat)    
    coef.sd <- sd(mat)
    limm <- as.vector(quantile(mat, probs=c((1-CI.level)/2, 1-(1-CI.level)/2)))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%CI-Low", sep=""),
                                                  paste(CI.level*100, "%CI-Upp", sep="")))
    ans$frailvar <- coef.table;
  }
  
  ## LPML and DIC
  ans$n <- object$n
  ans$p <- object$p
  ans$LPML <- sum(log(object$cpo))
  ans$DIC <- object$DIC
  
  ### acceptance rates
  ans$ratetheta = object$ratetheta;
  ans$ratebeta = object$ratebeta;
  ans$rateYs = object$rateYs;
  ans$ratec = object$ratec;
  ans$selection = object$selection;
  if(object$selection){
    ans$gamma = object$gamma;
  }
  ans$alpha = object$alpha;
  
  class(ans) <- "summary.survregbayes2"
  return(ans)
}


"print.summary.survregbayes2"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  if(x$p>0){
    cat("\nPosterior inference of regression coefficients\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratebeta, "):\n", sep="")
    print.default(format(x$coeff, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  if(x$alpha[1]==Inf){
    cat("\nPosterior inference of baseline parameters\n")
    message("Note: the baseline estimates are based on scaled covariates")
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
  
  if (!is.null(x$frailvar)) {
    if(x$frail.prior=="car"){
      cat("\nPosterior inference of conditional CAR frailty variance\n")
      print.default(format(x$frailvar, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }else if(x$frail.prior=="iid"){
      cat("\nPosterior inference of frailty variance\n")
      print.default(format(x$frailvar, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
  }
  
  if(x$selection){
    models=apply(x$gamma, 2, function(x) paste(x, sep="", collapse=",") )
    selected = table(models)/length(models);
    for(i in 1:length(selected)){
      selected.ind = scan(text = names(selected)[i], what = 0L, sep=",", quiet = TRUE);
      if(sum(selected.ind)==0){
        selected.name = "None";
      }else{
        selected.name = rownames(x$coeff)[selected.ind==1];
      }
      names(selected)[i] = paste(selected.name, collapse =",");
    }
    selected = sort(selected, decreasing = TRUE);
    selected.mat = matrix( selected, nrow=1 );
    colnames(selected.mat) = names(selected);
    rownames(selected.mat) = "prop.";
    cat("\nVariable selection:\n")
    print.default(format(selected.mat, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nLog pseudo marginal likelihood: LPML=", x$LPML, sep="")
  cat("\nDeviance Information Criterion: DIC=", x$DIC, sep="")
  cat("\nNumber of subjects: n=", x$n, "\n", sep="")

  invisible(x)
}


