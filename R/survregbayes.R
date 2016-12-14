"survregbayes" <- function (formula, data, na.action, survmodel="PH", dist="loglogistic", 
                            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                            prior=NULL, state=NULL, selection=FALSE, Proximity=NULL, 
                            truncation_time=NULL, subject.num=NULL, Knots=NULL, 
                            Coordinates=NULL, DIST=NULL, InitParamMCMC=TRUE, 
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
  #### Distance type:
  if(is.null(DIST)){
    DIST <- function(x, y) fields::rdist(x, y)
  }
  
  #########################################################################################
  # general setup based on frailty priors
  #########################################################################################
  if(!is.null(frail.prior)){
    ID = frail.terms[,1];
    orderindex = order(ID);
    if(!(sum(orderindex==(1:n))==n)) stop("please sort the data by ID");
    blocki = c(0, cumsum(as.vector(table(ID))));
    nID = length(blocki)-1;
    if(is.null(Knots)){
      nknots = prior$nknots; if(is.null(nknots)) nknots=nID;
    }else{
      nknots = nrow(Knots);
    }
    if(nknots>nID) stop("the number of knots needs to be smaller than the number of region IDs");
    if(frail.prior=="car") {
      frailtyCode = 1;
      if(is.null(Proximity)) stop("please specify prxoimity matrix");
      W = Proximity;
      D = rowSums(W);
      if (any(D==0)) stop("it seems that some region does not have any neighbers, which is not allowed, pelase check");
      Dmm = matrix(1000, nrow(W), nrow(W)); diag(Dmm)=0;
      Dmr = matrix(1000, nrow(W), 1);
      Drr = matrix(1000, 1,1);
    } else if (frail.prior=="iid"){
      frailtyCode = 2;
      W = matrix(0, nID, nID);
      Dmm = matrix(1000, nrow(W), nrow(W)); diag(Dmm)=0;
      Dmr = matrix(1000, nrow(W), 1);
      Drr = matrix(1000, 1,1);
    } else if (frail.prior=="grf") {
      frailtyCode = 3;
      if(is.null(Coordinates)) stop("please specify Coordinates for each ID");
      if(nrow(Coordinates)!=nID) stop("the number of coordinates should be equal to the number of ID")
      W = matrix(0, nID, nID);
      if(is.null(Knots)){
        if(nknots<nID){
          s0 = as.matrix(fields::cover.design(Coordinates, nd=nknots, DIST=DIST)$design);
        }else{
          s0 = Coordinates;
        }
      }else{
        s0 = Knots;
      }
      Dmm = DIST(Coordinates, Coordinates);
      if(min(Dmm[row(Dmm)!=col(Dmm)])<=0) stop("each ID should have different Coordinates");
      Dmr = DIST(Coordinates, s0); 
      Drr = DIST(s0, s0);
    } else {
      stop("This function only supports non-frailty, CAR frailty, IID, and GRF frailty models.")
    }
  } else {
    ID = NULL;
    frailtyCode = 0;
    blocki = c(0, n);
    W = matrix(1, length(blocki)-1, length(blocki)-1);
    Dmm = matrix(1000, nrow(W), nrow(W)); diag(Dmm)=0;
    Dmr = matrix(1000, nrow(W), 1);
    Drr = matrix(1000, 1,1);
  }
  phi = state$phi; if(is.null(phi)) phi=1;
  phib0_prior = 1;
  nu = prior$nu; if(is.null(nu)) nu=1;
  if(!is.null(frail.prior)){
    if((frail.prior=="grf")){
      maxdis = max(Dmm);
      phi_min = (-log(0.001))^(1/nu)/maxdis
      phib0_prior = -log(.95)/phi_min;
      phi = 1/phib0_prior;
      if(!is.null(state$phi)){
        phi = state$phi;
      }
      if (phi<=0) stop("phi in state arguement should be greater than 0.")
    }
  }
  # frailty initials
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
    frailtyCode0 = frailtyCode;
    blocki0=blocki; W0=W; v0=v; subjecti0=subjecti;
    Dmm0=Dmm; Dmr0=Dmr; Drr0=Drr; clustindx0=matrix(1,length(blocki)-1,1);
    truncation_time0 = truncation_time;
    t10 = t1; t20=t2; delta0=delta; X.scaled0=X.scaled;
    if(!is.null(frail.prior)){
      if(frail.prior=="grf") {
        frailtyCode0=2;
      }
    }
    scaleV = 10; scaleS = 100000;
    if(sum(as.numeric(names(table(delta))))==2){ ## if current status data
      scaleV= 1; scaleS = 100000;
    }
    if(survmodel=="AFT"){
      fit0<- .Call("AFT_BP", nburn_=nburn0, nsave_=nsave0, nskip_=0, ndisplay_=1000, ltr_=truncation_time0, subjecti_=subjecti0,
                   t1_=t10, t2_=t20, type_=delta0, X_=X.scaled0, theta_=theta0, beta_=beta0, weight_=c(1), 
                   cpar_=Inf, a0_=-1, b0_=1, theta0_=theta_prior, V0inv_=solve(scaleV*Vhat0), Vhat_=Vhat0, 
                   beta0_=beta_prior, S0inv_=solve(scaleS*Shat0), Shat_=Shat0, l0_=3000, adapter_=2.38^2, 
                   gamma_=rep(1.0, p), p0gamma_=rep(0.5, p), selection_=0, frailty_=frailtyCode0, v_=v0, blocki_=blocki0, 
                   W_=W0, clustindx_=clustindx0, Dmm_=Dmm0, Dmr_=Dmr0, Drr_=Drr0, phi_=phi, nu_=nu, a0phi_=1, b0phi_=1, 
                   lambda_=1, a0lambda_=1, b0lambda_=1, dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(survmodel=="PO"){
      fit0<- .Call("PO_BP", nburn_=nburn0, nsave_=nsave0, nskip_=0, ndisplay_=1000, ltr_=truncation_time0, subjecti_=subjecti0,
                   t1_=t10, t2_=t20, type_=delta0, X_=X.scaled0, theta_=theta0, beta_=beta0, weight_=c(1), 
                   cpar_=Inf, a0_=-1, b0_=1, theta0_=theta_prior, V0inv_=solve(scaleV*Vhat0), Vhat_=Vhat0, 
                   beta0_=beta_prior, S0inv_=solve(scaleS*Shat0), Shat_=Shat0, l0_=3000, adapter_=2.38^2, 
                   gamma_=rep(1.0, p), p0gamma_=rep(0.5, p), selection_=0, frailty_=frailtyCode0, v_=v0, blocki_=blocki0, 
                   W_=W0, clustindx_=clustindx0, Dmm_=Dmm0, Dmr_=Dmr0, Drr_=Drr0, phi_=phi, nu_=nu, a0phi_=1, b0phi_=1, 
                   lambda_=1, a0lambda_=1, b0lambda_=1, dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(survmodel=="PH"){
      fit0<- .Call("PH_BP", nburn_=nburn0, nsave_=nsave0, nskip_=0, ndisplay_=1000, ltr_=truncation_time0, subjecti_=subjecti0,
                   t1_=t10, t2_=t20, type_=delta0, X_=X.scaled0, theta_=theta0, beta_=beta0, weight_=c(1), 
                   cpar_=Inf, a0_=-1, b0_=1, theta0_=theta_prior, V0inv_=solve(scaleV*Vhat0), Vhat_=Vhat0, 
                   beta0_=beta_prior, S0inv_=solve(scaleS*Shat0), Shat_=Shat0, l0_=3000, adapter_=2.38^2, 
                   gamma_=rep(1.0, p), p0gamma_=rep(0.5, p), selection_=0, frailty_=frailtyCode0, v_=v0, blocki_=blocki0, 
                   W_=W0, clustindx_=clustindx0, Dmm_=Dmm0, Dmr_=Dmr0, Drr_=Drr0, phi_=phi, nu_=nu, a0phi_=1, b0phi_=1, 
                   lambda_=1, a0lambda_=1, b0lambda_=1, dist_=distcode, PACKAGE = "spBayesSurv");
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
      if(is.null(state$frail)) v = rowMeans(fit0$v);
    }
    message("Starting the MCMC for the semiparametric model:")
  }else{
    ## initial fit for theta:
    fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, 
                              dist=dist);
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
      fit0 <- survival::survreg(formula = survival::Surv(tbase1, tbase2, type="interval2")~Xbase.scaled, 
                                dist=dist0);
      beta0 = -fit0$coefficients[-1]/fit0$scale; beta_prior = -fit0$coefficients[-1]/fit0$scale;
      beta = -fit0$coefficients[-1]/fit0$scale;
      Shat0 = as.matrix(fit0$var[c(2:(1+p)),c(2:(1+p))])/(fit0$scale)^2;
    }
    ## initial frailties
    lambda=1;
    if(!is.null(frail.prior)){
      if(is.null(state$frail)) v=rep(0, length(blocki)-1);
    }
  }
  
  #########################################################################################
  # priors
  # note the priors should be based on scaled data.
  #########################################################################################
  cpar=state$cpar; if(is.null(cpar)) cpar=1;
  tau2 = state$tau2; if(is.null(tau2)) tau2=1/lambda; lambda=1/tau2;
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  maxL <- prior$maxL; if(is.null(maxL)) maxL<-15;
  #if (cpar==Inf) maxL=1;
  weight = rep(1/maxL, maxL);
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
  M=prior$M; if(is.null(M)) M=10;
  q=prior$q; if(is.null(q)) q=0.9;
  if(is.null(S0)) {
    if(selection){
      gg = (log(M)/qnorm(q))^2/p;
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
  taua0 = prior$taua0; if(is.null(taua0)) taua0=.001;
  taub0 = prior$taub0; if(is.null(taub0)) taub0=.001;
  phia0 = prior$phia0; if(is.null(phia0)) phia0=1;
  phib0 = prior$phib0; if(is.null(phib0)) phib0=phib0_prior;
  mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
  theta_initial=theta+0;
  beta_initial=beta+0;
  tau2_initial=1/lambda+0;
  frail_initial=v+0;
  clustindx=matrix(1,length(blocki)-1, 1);
  if(!is.null(frail.prior)){
    initial.values = list(beta=beta_initial, theta=theta_initial, tau2=tau2_initial, frail=frail_initial);
    prior = list(maxL=maxL, a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0,
                 Shat=Shat, Vhat=Vhat, taua0=taua0, taub0=taub0);
    if((frail.prior=="grf")){
      prior$nknots=nknots; initial.values$phi=phi; prior$nu=nu;
      prior$phia0=phia0; prior$phib0=phib0;
      if(is.null(prior$nblock)) nblock=nID;
      if(nblock==nID){
        clustindx=diag(1,length(blocki)-1, length(blocki)-1);
      }else{
        s0tmp = as.matrix(fields::cover.design(Coordinates, nd=nblock, DIST=DIST)$design);
        Dtmp = DIST(Coordinates, s0tmp);
        idtmp = apply(Dtmp, 1, which.min);
        nblock=length(table(idtmp));
        idnames = as.numeric(names(table(idtmp)))
        clustindx=matrix(0, length(blocki)-1, nblock);
        for(jj in 1:nblock){
          clustindx[which(idtmp==idnames[jj]),jj] = 1;
        }
      }
      prior$nblock=nblock; 
    }
  }else{
    initial.values = list(beta=beta_initial, theta=theta_initial);
    prior = list(maxL=maxL, a0=a0, b0=b0, theta0=theta0, V0=V0, beta0=beta0, S0=S0,
                 Shat=Shat, Vhat=Vhat);
  }
  if(selection){
    prior$p0gamma=p0gamma; prior$M=M; prior$q=q;
  }
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  if(survmodel=="AFT"){
    model.name <- "Accelerated failure time model:";
    foo <- .Call("AFT_BP", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, ltr_=truncation_time,
                 subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, weight_=weight, 
                 cpar_=cpar, a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, 
                 beta0_=beta0, S0inv_=S0inv, Shat_=Shat, l0_=min(5000,nsave/2), adapter_=2.38^2, 
                 gamma_=gamma0, p0gamma_=p0gamma, selection_=selection+0,
                 frailty_=frailtyCode, v_=v, blocki_=blocki, W_=W, clustindx_=clustindx, 
                 Dmm_=Dmm, Dmr_=Dmr, Drr_=Drr, phi_=phi, nu_=nu, a0phi_=phia0, b0phi_=phib0, lambda_=lambda, 
                 a0lambda_=taua0, b0lambda_=taub0, dist_=distcode, PACKAGE = "spBayesSurv");
  }else if(survmodel=="PO"){
    model.name <- "Proportional Odds model:";
    foo <- .Call("PO_BP", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, ltr_=truncation_time,
                 subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, weight_=weight, 
                 cpar_=cpar, a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, 
                 beta0_=beta0, S0inv_=S0inv, Shat_=Shat, l0_=min(5000,nsave/2), adapter_=2.38^2, 
                 gamma_=gamma0, p0gamma_=p0gamma, selection_=selection+0,
                 frailty_=frailtyCode, v_=v, blocki_=blocki, W_=W, clustindx_=clustindx, 
                 Dmm_=Dmm, Dmr_=Dmr, Drr_=Drr, phi_=phi, nu_=nu, a0phi_=phia0, b0phi_=phib0, lambda_=lambda, 
                 a0lambda_=taua0, b0lambda_=taub0, dist_=distcode, PACKAGE = "spBayesSurv");
  }else if(survmodel=="PH"){
    model.name <- "Proportional hazards model:";
    foo <- .Call("PH_BP", nburn_=nburn, nsave_=nsave, nskip_=nskip, ndisplay_=ndisplay, ltr_=truncation_time,
                 subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, X_=X.scaled, theta_=theta, beta_=beta, weight_=weight, 
                 cpar_=cpar, a0_=a0, b0_=b0, theta0_=theta0, V0inv_=V0inv, Vhat_=Vhat, 
                 beta0_=beta0, S0inv_=S0inv, Shat_=Shat, l0_=min(5000,nsave/2), adapter_=2.38^2, 
                 gamma_=gamma0, p0gamma_=p0gamma, selection_=selection+0,
                 frailty_=frailtyCode, v_=v, blocki_=blocki, W_=W, clustindx_=clustindx, 
                 Dmm_=Dmm, Dmr_=Dmr, Drr_=Drr, phi_=phi, nu_=nu, a0phi_=phia0, b0phi_=phib0, lambda_=lambda, 
                 a0lambda_=taua0, b0lambda_=taub0, dist_=distcode, PACKAGE = "spBayesSurv");
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
  
  #### Cox-Snell residuals
  resid1 = foo$resid1; resid2 = foo$resid2;
  for(i in 1:nsubject){
    if(delta[lastindx[i]]==0){
      resid2[i]=NA;
    }else if (delta[lastindx[i]]==1){
      resid2[i]=resid1[i];
    }else if (delta[lastindx[i]]==2){
      resid1[i]=NA;
    }
  }
  Surv.cox.snell = survival::Surv(resid1, resid2, type="interval2");
  
  ### Calculate Bayes Factors for log(weight[-maxL])-log(weight[maxL]);
  HH = apply(foo$weight, 2, function(x) log(x[-maxL])-log(x[maxL]) ); 
  meanH = apply(HH, 1, mean);
  varH = var(t(HH));
  meancpar = mean(foo$cpar);
  numH = exp( lgamma(meancpar*maxL)-maxL*lgamma(meancpar)-meancpar*maxL*log(maxL) );
  denH = as.vector((2*pi)^(-(maxL-1)/2)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
  BayesFactor = numH/denH; 
  
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
                 cpar = foo$cpar,
                 maxL = maxL,
                 weight = foo$weight,
                 cpo = foo$cpo,
                 pD = foo$pD, 
                 DIC = foo$DIC,
                 Surv.cox.snell = Surv.cox.snell,
                 ratetheta = foo$ratetheta,
                 ratebeta = foo$ratebeta,
                 rateYs = foo$rateYs,
                 ratec = foo$ratec,
                 frail.prior=frail.prior,
                 selection = selection,
                 initial.values=initial.values,
                 BF = BayesFactor);
  if(!is.null(frail.prior)){
    output$v = foo$v;
    output$ratev = foo$ratev;
    output$tau2 = 1/foo$lambda;
    output$ID = ID;
    if(frail.prior=="grf"){
      output$Coordinates = Coordinates;
      output$ratephi = foo$ratephi;
      output$phi = foo$phi;
      output$Knots = s0;
    }
  }
  if(selection){
    output$gamma = foo$gamma;
  }
  class(output) <- c("survregbayes")
  output
}

#### BF for testing the centering distribution
"BF.survregbayes" <- function (x) {
  if(is(x,"survregbayes")){
    maxL = x$prior$maxL;
    HH = apply(x$weight, 2, function(y) log(y[-maxL])-log(y[maxL]) ); 
    meanH = apply(HH, 1, mean);
    varH = var(t(HH));
    meancpar = mean(x$cpar);
    numH = exp( lgamma(meancpar*maxL)-maxL*lgamma(meancpar)-meancpar*maxL*log(maxL) );
    denH = as.vector((2*pi)^(-(maxL-1)/2)*(det(varH))^(-1/2)*exp(-1/2*t(meanH)%*%solve(varH)%*%meanH));
    BayesFactor = numH/denH;
  }
  BayesFactor
}

#### print, summary, plot
"print.survregbayes" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior means for regression coefficients:\n")
  if(x$p>0){
    print.default(format(x$coefficients[1:x$p], digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat(paste("\nBayes Factor for ", x$dist, " baseline vs. Bernstein poly:", sep=""), x$BF)       
  
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

"cox.snell.survregbayes" <- function (x, ncurves=0) {
  Resids = list();
  if(is(x,"survregbayes")){
    Y = x$Surv;
    t1 = Y[,1]; t2 = Y[,1];
    type <- attr(Y, "type")
    exactsurv <- Y[,ncol(Y)] ==1
    if (any(exactsurv)) {
      t1[exactsurv]=Y[exactsurv,1];
      t2[exactsurv]=Y[exactsurv,1];
    }
    if (type== 'counting'){ #stop ("Invalid survival type")
      t1 = Y[,2]; t2 = Y[,2];
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
    truncation_time = x$truncation_time;
    subject.num = x$subject.num;
    subjecti = c(0, cumsum(as.vector(table(subject.num))));
    nsubject = length(subjecti)-1;
    baseindx = (subjecti+1)[-length(subjecti)];
    lastindx = subjecti[-1];
    if(is.null(x$ID)){
      frailn = matrix(0, x$n, x$mcmc$nsave);
    }else{
      ID = x$ID;
      freq.frail = table(ID);
      frailn = apply(x$v, 2, function(x) rep(x, freq.frail) );
    }
    distcode = switch(x$dist, loglogistic=1, lognormal=2, 3);
    if(x$selection){
      betafitted = x$beta.scaled*x$gamma;
    }else{
      betafitted = x$beta.scaled;
    }
    if(x$survmodel=="AFT"){
      foo <- .Call("AFT_BP_cox_snell", ltr_=truncation_time, subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, 
                   X_=x$X.scaled, theta_=x$theta.scaled, beta_=betafitted, vn_=frailn, weight_=x$weight, 
                   dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PO"){
      foo <- .Call("PO_BP_cox_snell", ltr_=truncation_time, subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, 
                   X_=x$X.scaled, theta_=x$theta.scaled, beta_=betafitted, vn_=frailn, weight_=x$weight, 
                   dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PH"){
      foo <- .Call("PH_BP_cox_snell", ltr_=truncation_time, subjecti_=subjecti, t1_=t1, t2_=t2, type_=delta, 
                   X_=x$X.scaled, theta_=x$theta.scaled, beta_=betafitted, vn_=frailn, weight_=x$weight, 
                   dist_=distcode, PACKAGE = "spBayesSurv");
    }else{
      stop("This function only supports PH, PO or AFT");
    }
    resid1 = -log(apply(foo$St1, 1, mean)); 
    resid2 = -log(apply(foo$St2, 1, mean)); 
    for(i in 1:nsubject){
      if(delta[lastindx[i]]==0){
        resid2[i]=NA;
      }else if (delta[lastindx[i]]==1){
        resid2[i]=resid1[i];
      }else if (delta[lastindx[i]]==2){
        resid1[i]=NA;
      }
    }
    Resids$resid=survival::Surv(resid1, resid2, type="interval2");
    if(ncurves>=1){
      res.indx = sample(x$mcmc$nsave, ncurves);
      for(k in 1:ncurves){
        resid1 = -log(foo$St1[,res.indx[k]]);
        resid2 = -log(foo$St2[,res.indx[k]]);
        for(i in 1:nsubject){
          if(delta[lastindx[i]]==0){
            resid2[i]=NA;
          }else if (delta[lastindx[i]]==1){
            resid2[i]=resid1[i];
          }else if (delta[lastindx[i]]==2){
            resid1[i]=NA;
          }
        }
        Resids[[k+1]] = survival::Surv(resid1, resid2, type="interval2");
        names(Resids)[k+1] = paste("resid", k, sep="");
      }
    }
  }
  Resids$St1 = foo$St1; Resids$St2 = foo$St2; 
  Resids
}

"plot.survregbayes" <- function (x, xpred=NULL, tgrid=NULL, frail=NULL, CI=0.95, PLOT=FALSE, ...) {
  if(is(x,"survregbayes")){
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
    nxpred = nrow(xpred);
    if(is.null(frail)){
      frail=matrix(0, nrow=nxpred, ncol=x$mcmc$nsave);
    }else{
      if(is.vector(frail)) frail=matrix(frail, nrow=1);
      if((nrow(frail)!=nxpred)|(ncol(frail)!=x$mcmc$nsave)) stop("The dim of frail should be nrow(xpred) by nsave.")
    }
    for(i in 1:nxpred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
    distcode = switch(x$dist, loglogistic=1, lognormal=2, 3);
    betafitted = x$beta.scaled;
    if(x$selection){
      betafitted = x$beta.scaled*x$gamma;
    }
    if(x$survmodel=="AFT"){
      estimates <- .Call("AFT_BP_plots", tgrid, xpred, x$theta.scaled, betafitted, frail, x$weight, CI, 
                         dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PO"){
      estimates <- .Call("PO_BP_plots", tgrid, xpred, x$theta.scaled, betafitted, frail, x$weight, CI, 
                         dist_=distcode, PACKAGE = "spBayesSurv");
    }else if(x$survmodel=="PH"){
      estimates <- .Call("PH_BP_plots", tgrid, xpred, x$theta.scaled, betafitted, frail, x$weight, CI, 
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

"summary.survregbayes" <- function(object, CI.level=0.95, ...) {
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
    mat <- object$cpar
    coef.p <- mean(mat); names(coef.p)="";    
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
    if((object$frail.prior=="grf")){
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
      ans$phivar <- coef.table;
    }
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
  ans$cpar = object$cpar;
  ans$BF = object$BF;
  
  class(ans) <- "summary.survregbayes"
  return(ans)
}


"print.summary.survregbayes"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  if(x$p>0){
    cat("\nPosterior inference of regression coefficients\n")
    cat("(Adaptive M-H acceptance rate: ", x$ratebeta, "):\n", sep="")
    print.default(format(x$coeff, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  if(x$cpar[1]==Inf){
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
    } else if((x$frail.prior=="grf")){
      cat("\nPosterior inference of frailty variance\n")
      print.default(format(x$frailvar, digits = digits), print.gap = 2, 
                    quote = FALSE)
      
      cat("\nPosterior inference of correlation function range phi\n")
      print.default(format(x$phivar, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
  }
  
  cat(paste("\nBayes Factor for ", x$dist, " baseline vs. Bernstein poly:", sep=""), x$BF)  
  
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


