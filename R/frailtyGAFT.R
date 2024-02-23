"frailtyGAFT" <- function (formula, data, na.action, 
                           mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                           prior=NULL, state=NULL, Proximity=NULL, Coordinates=NULL, 
                           DIST=NULL, scale.designX=TRUE){
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
  pce = p+1;
  if(p==0){
    X.scaled <- NULL;
    X1 = cbind(rep(1,n), X.scaled); colnames(X1)[1]="intercept";
  }else{
    if(scale.designX){
      X.scaled <- scale(X);
    }else{
      X.scaled <- scale(X, center=rep(0,p), scale=rep(1,p));
    }
    X.center = attributes(X.scaled)$`scaled:center`;
    X.scale = attributes(X.scaled)$`scaled:scale`;
    X1 = cbind(rep(1,n), X.scaled); colnames(X1)[1]="intercept";
  }
  if(is.null(Xtf)){
    Xtf.scaled <- NULL;
    Xtf1 = cbind(rep(1,n), Xtf.scaled); colnames(Xtf1)[1]="intercept";
  }else{
    if(scale.designX){
      Xtf.scaled <- scale(Xtf);
    }else{
      Xtf.scaled <- scale(Xtf, center=rep(0,ncol(Xtf)), scale=rep(1,ncol(Xtf)));
    }
    Xtf.center = attributes(Xtf.scaled)$`scaled:center`;
    Xtf.scale = attributes(Xtf.scaled)$`scaled:scale`;
    Xtf1 = cbind(rep(1,n), Xtf.scaled); colnames(Xtf1)[1]="intercept";
  }
  ptf = ncol(Xtf1);
  
  #########################################################################################
  # data structure
  #########################################################################################
  t1 = Y[,1]; t2 = Y[,1];
  type <- attr(Y, "type")
  if (type== 'counting') stop ("Invalid survival type")
  exactsurv <- Y[,ncol(Y)] ==1
  if (any(exactsurv)) {
    t1[exactsurv]=Y[exactsurv,1];
    t2[exactsurv]=Y[exactsurv,1];
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
  frail.prior = colnames(frail.terms)[1];
  ID = frail.terms[,1];
  #### Distance type:
  if(is.null(DIST)){
    DIST <- function(x, y) fields::rdist(x, y)
  }
  
  ##############################################
  ### Currently it only supports AFT ###########
  ##############################################
  #########################################################################################
  # initial MLE analysis and mcmc parameters
  #########################################################################################
  ## initial fit
  fit0 <- survival::survreg(formula = Y~X1-1, dist="lognormal");
  betace =  fit0$coefficients; 
  betaShat0 = 100*fit0$var[1:pce,1:pce];
  sigma2 = (fit0$scale)^2;
  sig2hat <- (fit0$scale)^2;
  sig2var <- 100*fit0$var[pce+1, pce+1]*4*sig2hat^2;
  sig2a0 = 2+sig2hat^2/sig2var;
  sig2b0 = sig2hat*(sig2a0-1);
  y <- state$logt; 
  if(is.null(y)){
    y <- rep(0, n);
    for(i in 1:n){
      if(delta[i]==0) y[i] = log(t1[i]+1);  
      if(delta[i]==1) y[i] = log(t1[i]); 
      if(delta[i]==2) y[i] = log(t2[i]/2);
      if(delta[i]==3) y[i] = log(mean(c(t1[i], t2[i])));
    }
  }  
  
  #########################################################################################
  # check frailty
  #########################################################################################
  if(!is.null(frail.prior)) {
    if(is.null(ID)) stop("please specify ID");
    orderindex = order(ID); 
    if(!(sum(orderindex==(1:n))==n)) stop("please sort the data by ID");
    blocki = c(0, cumsum(as.vector(table(ID))));
    if(frail.prior=="car") {
      if(is.null(Proximity)) stop("please specify prxoimity matrix");
      W = Proximity;
      D = rowSums(W);
      if (any(D==0)) stop("it seems that some region does not have any neighbers, which is not allowed, pelase check")
    }else if(frail.prior=="iid"){
      W = matrix(0, length(blocki)-1, length(blocki)-1);
      D = rep(0, length(blocki)-1);
    }else if(frail.prior=="grf"){
      if(is.null(Coordinates)) stop("please specify Coordinates for each ID");
      if(nrow(Coordinates)!=(length(blocki)-1)) stop("the number of coordinates should be equal to the number of ID")
      Dmm = DIST(Coordinates, Coordinates);
      if(min(Dmm[row(Dmm)!=col(Dmm)])<=0) stop("each ID should have different Coordinates");
    }else{
      stop("This function only supports non-frailty, car frailty, iid and grf frailty models.")
    }
  }else{
    ID = NULL;
    blocki = c(0, n);
    W = matrix(1, length(blocki)-1, length(blocki)-1);
    D = rep(1, length(blocki)-1);
    v <- rep(0, length(blocki)-1);
  }
  phi = state$phi; if(is.null(phi)) phi=1;
  phia0_prior = 2;
  nu = prior$nu; if(is.null(nu)) nu=1;
  if(!is.null(frail.prior)){
    if(frail.prior=="grf"){
      maxdis = max(Dmm);
      #phi_min = (-log(0.001))^(1/nu)/maxdis; phib0_prior = -log(.95)/phi_min; phi = 1/phib0_prior;
      #cc = sqrt((-log(0.001))^(1/nu)/maxdis); phi = maxdis*cc;
      phi = (-log(0.001))^(1/nu)/maxdis; 
      phia0_prior = 2;
      if(!is.null(state$phi)){
        phi = state$phi;
      }
      if (phi<=0) stop("phi in state arguement should be greater than 0.")
    }
  }
  phi0=phi;
  if(is.null(state$frail)) {
    v <- rep(0, length(blocki)-1);
  } else {
    v <- state$frail; if(length(v)!=(length(blocki)-1)) stop("check the length of frail");
  }
  #########################################################################################
  # priors
  #########################################################################################
  alpha=state$alpha; if(is.null(alpha)) alpha=5;
  tau2 = state$tau2; if(is.null(tau2)) tau2=0.5;
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  maxL <- prior$maxL; if(is.null(maxL)) maxL<-4;
  ntprob <- 2^(maxL+1)-2; 
  ntlr <- 2^maxL-1;
  betatf <- matrix(0,nrow=ptf,ncol=ntlr);
  a0=prior$a0; if(is.null(a0)) a0=5;
  b0=prior$b0; if(is.null(b0)) b0=1;
  if(a0<=0){
    a0=-1;alpha=state$alpha; 
    if(is.null(alpha)) stop("please specify state$alpha if prior$a0 is not positive");
  }
  m0 <- prior$m0; if(is.null(m0)) m0 <- rep(0, pce);
  S0 <- prior$S0; if(is.null(S0)) S0 <- diag(1e5, pce);
  S0inv <- solve(S0);
  siga0=prior$siga0; if(is.null(siga0)) siga0=sig2a0;
  sigb0=prior$sigb0; if(is.null(sigb0)) sigb0=sig2b0;
  gprior <- prior$gprior; if(is.null(gprior)) gprior <- 2*n*solve(t(Xtf1)%*%Xtf1);
  taua0 = prior$taua0; if(is.null(taua0)) taua0=1;
  taub0 = prior$taub0; if(is.null(taub0)) taub0=1;
  phia0 = prior$phia0; if(is.null(phia0)) phia0=phia0_prior;
  phib0 = prior$phib0; if(is.null(phib0)) phib0=(phia0-1)/phi0;
  mm = prior$mm; if (is.null(mm)) mm = 10;
  win = prior$win; if (is.null(win)) win = 1;
  mcmc = list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=ndisplay)
  if(!is.null(frail.prior)){
    prior = list(maxL=maxL, a0=a0, b0=b0, siga0=siga0, sigb0=sigb0, m0=m0, S0=S0, mm=mm, win=win);
    if((frail.prior=="iid")|(frail.prior=="car")){
      prior$taua0=taua0; prior$taub0=taub0;
    }else if (frail.prior=="grf"){
      prior$nu=nu;
      prior$taua0=taua0; prior$taub0=taub0; #prior$silla0=silla0; prior$sillb0=sillb0;
      prior$phia0=phia0; prior$phib0=phib0;
    }
  }else{
    prior = list(maxL=maxL, a0=a0, b0=b0, siga0=siga0, sigb0=sigb0, m0=m0, S0=S0, mm=mm, win=win)
  }
  
  #########################################################################################
  # calling the c++ code and # output 
  #########################################################################################
  if(!is.null(frail.prior)){
    model.name <- "Generalized accelerated failure time frailty model:";
    if((frail.prior=="iid")|(frail.prior=="car")){
      foo <- .Call("frailtyLDTFP", nburn_ = nburn, nsave_ = nsave, nskip_ = nskip, ndisplay_ = ndisplay,
                   tobs_ = cbind(t1, t2), type_ = delta, xce_ = t(X1), xtf_ = t(Xtf1), alpha_ = alpha, 
                   betace_ = betace, betatf_ = betatf,  sigma2_ = sigma2, y_ = y, v_ = v, blocki_ = blocki,
                   tau2_ = tau2, W_ = W, D_ = D, maxL_ = maxL, a0_ = a0, b0_ = b0,  m0_ = m0,  S0inv_ = S0inv, 
                   gprior_ = gprior, a0sig_ = siga0, b0sig_ = sigb0, a0tau_ = taua0, b0tau_ = taub0, 
                   win_ = win, mm_ = mm, PACKAGE = "spBayesSurv");
    }else if (frail.prior=="grf"){
      foo <- .Call("frailty_GRF_LDTFP", nburn_ = nburn, nsave_ = nsave, nskip_ = nskip, ndisplay_ = ndisplay,
                   tobs_ = cbind(t1, t2), type_ = delta, xce_ = t(X1), xtf_ = t(Xtf1), alpha_ = alpha, 
                   betace_ = betace, betatf_ = betatf,  sigma2_ = sigma2, y_ = y, v_ = v, blocki_ = blocki,
                   tau2_ = tau2, Dmm_ = Dmm, maxL_ = maxL, a0_ = a0, b0_ = b0,  
                   m0_ = m0,  S0inv_ = S0inv, gprior_ = gprior, a0sig_ = siga0, b0sig_ = sigb0, 
                   a0tau_ = taua0, b0tau_ = taub0, nu_ = nu, phi_ = phi, 
                   a0phi_ = phia0, b0phi_ = phib0, win_ = win, mm_ = mm, PACKAGE = "spBayesSurv");
    }
  }else{
    model.name <- "Generalized accelerated failure time model:";
    foo <- .Call("nonfrailtyLDTFP", nburn_ = nburn, nsave_ = nsave, nskip_ = nskip, ndisplay_ = ndisplay,
                 tobs_ = cbind(t1, t2), type_ = delta, xce_ = t(X1), xtf_ = t(Xtf1), alpha_ = alpha, 
                 betace_ = betace, betatf_ = betatf,  sigma2_ = sigma2, y_ = y, maxL_ = maxL, a0_ = a0, b0_ = b0,  
                 m0_ = m0,  S0inv_ = S0inv, gprior_ = gprior, a0sig_ = siga0, b0sig_ = sigb0, 
                 win_ = win, mm_ = mm, PACKAGE = "spBayesSurv");
  }
  #########################################################################################
  # save state
  #########################################################################################
  #### transfer the estimates back to original scales;
  beta.scaled = matrix(foo$beta, pce, nsave);  beta.original = matrix(foo$beta, pce, nsave);
  if(p>0){
    beta.original[2:pce, ] = matrix(beta.scaled[2:pce, ], p, nsave)/matrix(rep(X.scale, nsave), p, nsave);
    mat.tmp1 = matrix(beta.scaled[2:pce, ], p, nsave)
    mat.tmp2 = matrix(rep(X.center, nsave), p, nsave)
    mat.tmp3 = matrix(rep(X.scale, nsave), p, nsave)
    beta.original[1, ] = beta.scaled[1,] - colSums(mat.tmp1*mat.tmp2/mat.tmp3)
  }
  
  #### coefficients
  coeff <- c( apply(beta.original, 1, mean), mean(sqrt(foo$sigma2)), mean(foo$alpha) );
  names(coeff) = c(colnames(X1), "scale", "precision");
  
  #### Save to a list
  output <- list(modelname=model.name,
                 terms = m,
                 coefficients=coeff,
                 call=Call,
                 prior=prior,
                 mcmc=mcmc,
                 n=n,
                 pce=pce,
                 ptf=ptf,
                 Surv=Y,
                 X.scaled=X.scaled,
                 X=X,
                 Xtf.scaled=Xtf.scaled,
                 Xtf=Xtf,
                 sigma2 = foo$sigma2,
                 beta = beta.original,
                 beta.scaled = beta.scaled,
                 alpha = foo$alpha,
                 maxL = maxL,
                 betatf = foo$betatf,
                 logt = foo$y,
                 cpo = foo$cpo,
                 accept_beta = foo$ratebetace,
                 accept_betatf = foo$ratebetatf,
                 frail.prior=frail.prior);
  ## check frailties
  if(!is.null(frail.prior)){
    output$v = foo$v;
    output$ratev = foo$ratev;
    output$tau2 = foo$tau2;
    output$ID = ID;
    if (frail.prior=="grf"){
      output$Coordinates = Coordinates
      output$phi = foo$phi;
      output$ratephi = foo$ratephi;
    }
  }
  
  ### Calculate Bayes Factors for betatf
  gprior = solve(t(Xtf1)%*%Xtf1)*(2*n);
  betatfmat = matrix(as.vector(foo$betatf[,-1,]), (2^maxL-2)*ptf, nsave);
  BFs = .Call("BayesFactor", betatf_=betatfmat, maxL_=maxL, gprior_=gprior, alpha_=mean(foo$alpha), PACKAGE = "spBayesSurv");
  BayesFactors = c(as.vector( BFs$BFindividual ), BFs$BFoverallLDTFP, BFs$BFoverallParam); 
  names(BayesFactors) = c( colnames(Xtf1), "overall", "normality");
  output$BF = BayesFactors; 
  class(output) <- c("frailtyGAFT")
  output
}

#### print, summary, plot
"print.frailtyGAFT" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior means for regression coefficients:\n")
  print.default(format(x$coefficients[1:x$pce], digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nBayes factors for LDTFP covariate effects:\n")       
  print.default(format(x$BF, digits = digits), print.gap = 2, 
                quote = FALSE) 
  
  cat("\nLPML:", sum(log(x$cpo)))
  cat("\nn =",x$n, "\n")
  invisible(x)
}

"plot.frailtyGAFT" <- function (x, xnewdata, xtfnewdata, tgrid = NULL, 
                                frail = NULL, CI = 0.95, PLOT = TRUE, ...) {
  if(is.null(tgrid)) tgrid = seq(0.01, max(x$Surv[,1], na.rm=T), length.out=200);
  if(x$pce==1){
    if(is.null(frail)){
      if(missing(xtfnewdata)){
        xpred = cbind(1); nxpred=1;
        rnames = "baseline"
      }else{
        nxpred = nrow(xtfnewdata);
        xpred = matrix(1, nrow = nxpred);
        rnames = row.names(xtfnewdata)
      }
    }else{
      if(is.vector(frail)) frail=matrix(frail, nrow=1);
      nxpred = nrow(frail);
      xpred = matrix(1, nrow = nxpred);
      rnames = row.names(as.data.frame(frail))
      if(!missing(xtfnewdata)){
        if(nrow(xtfnewdata)!=nxpred) stop("xtfnewdata and frail should have the same numbers of rows")
      }
    }
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
      nxpred = nrow(xpred);
      X.center = attributes(x$X.scaled)$`scaled:center`;
      X.scale = attributes(x$X.scaled)$`scaled:scale`;
      for(i in 1:nxpred) xpred[i,] = (xpred[i,]-X.center)/X.scale;
      xpred = cbind(rep(1,nxpred),xpred);
      if(ncol(xpred)!=x$pce) stop("please make sure the number of columns matches!");
    }
  }
  
  if(x$ptf==1){
    xtfpred = matrix(1, nrow = nxpred);
  }else{
    if(missing(xtfnewdata)) {
      stop("please specify xtfnewdata")
    }else{
      newTerms = attr(x$Xtf, 'terms')
      mnew <- model.frame(newTerms, xtfnewdata, na.action = na.omit, xlev = attr(x$Xtf, "levels"))
      Xnew <- model.matrix(newTerms, mnew);
      xtfpred = cbind(Xnew);
      Xtf.center = attributes(x$Xtf.scaled)$`scaled:center`;
      Xtf.scale = attributes(x$Xtf.scaled)$`scaled:scale`;
      for(i in 1:nxpred) xtfpred[i, 2:x$ptf] = (xtfpred[i,2:x$ptf]-Xtf.center)/Xtf.scale;
      if(ncol(xtfpred)!=x$ptf) stop("please make sure the number of columns matches!");
    }
  }
  if(nrow(xpred)!=nrow(xtfpred)) stop("please make sure xnewdata and xtfnewdata have the same numbers of rows!");
  if(is.null(frail)){
    frail=matrix(0, nrow=nxpred, ncol=x$mcmc$nsave);
  }else{
    if(is.vector(frail)) frail=matrix(frail, nrow=1);
    if((nrow(frail)!=nxpred)|(ncol(frail)!=x$mcmc$nsave)) stop("The dim of frail should be nrow(xpred) by nsave.")
  }
  estimates <- .Call("frailtyGAFTplots", tgrid, xpred, xtfpred, x$beta.scaled, x$betatf, frail, 
                     x$sigma2, x$maxL, CI, PACKAGE = "spBayesSurv");
  
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

"summary.frailtyGAFT" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### CPO
  ans$cpo <- object$cpo
  
  ### Median information
  mat <- as.matrix(object$beta)
  coef.p <- object$coefficients[(1:object$pce)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(coda::HPDinterval(coda::as.mcmc(x), prob=CI.level)))
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%HPD-Low", sep=""),
                                                paste(CI.level*100, "%HPD-Upp", sep="")))
  ans$coeff <- coef.table
  
  ### Scale parameter
  mat <- c(sqrt( object$sigma2)); 
  coef.p <- object$coefficients[c(object$pce+1)];
  coef.m <- median(mat) 
  coef.sd <-sd(mat)
  limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
  coef.l <- limm[1]
  coef.u <- limm[2]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%HPD-Low", sep=""),
                                                paste(CI.level*100, "%HPD-Upp", sep="")))
  ans$scale<- coef.table
  
  ### Precision parameter
  if(object$prior$a0<=0){
    ans$prec <- NULL
  }else{
    mat <- object$alpha
    coef.p <- mean(mat); names(coef.p)="alpha";
    coef.m <- median(mat)    
    coef.sd <- sd(mat)
    limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%HPD-Low", sep=""),
                                                  paste(CI.level*100, "%HPD-Upp", sep="")))
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
    limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%HPD-Low", sep=""),
                                                  paste(CI.level*100, "%HPD-Upp", sep="")))
    ans$frailvar <- coef.table;
    if(object$frail.prior=="grf"){
      mat <- object$phi
      coef.p <- mean(mat); names(coef.p)="range";
      coef.m <- median(mat)    
      coef.sd <- sd(mat)
      limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
      coef.l <- limm[1]
      coef.u <- limm[2]
      
      coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
      dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                    paste(CI.level*100, "%HPD-Low", sep=""),
                                                    paste(CI.level*100, "%HPD-Upp", sep="")))
      ans$range <- coef.table;
    }
  }
  
  ### summaries
  ans$n <- object$n
  ans$pce <- object$pce
  ans$LPML <- sum(log(object$cpo))
  ans$BF <- object$BF;
  
  ### acceptance rates
  ans$ratephi = object$ratephi;
  
  class(ans) <- "summary.frailtyGAFT"
  return(ans)
}


"print.summary.frailtyGAFT"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  if(!is.null(x$coeff)){
    cat("\nPosterior inference of regression coefficients\n")
    print.default(format(x$coeff, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nPosterior inference of scale parameter\n")
  print.default(format(x$scale, digits = digits), print.gap = 2, quote = FALSE)
  
  if (!is.null(x$prec)) {
    cat("\nPosterior inference of precision parameter of LDTFP\n")
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
    } else if(x$frail.prior=="grf"){
      cat("\nPosterior inference of frailty variance\n")
      print.default(format(x$frailvar, digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("\nPosterior inference of correlation function range phi \n")
      print.default(format(x$range, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
  }
  
  cat("\nBayes factors for LDTFP covariate effects:\n")       
  print.default(format(x$BF, digits = digits), print.gap = 2, 
                quote = FALSE) 
  
  cat("\nLog pseudo marginal likelihood: LPML", x$LPML, sep="=")
  #cat("\nDeviance information criterion: DIC", x$DIC, sep="=")
  cat("\nNumber of subjects:", x$n, sep="=")     
  invisible(x)
}
