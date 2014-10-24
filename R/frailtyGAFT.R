"frailtyGAFT" <- function(tobs, x=NULL, xtf=NULL, prior, mcmc, state, frailty=NULL, Proximity=NULL,
                            data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
UseMethod("frailtyGAFT")

"frailtyGAFT.default" <- 
  function (tobs, 
            x=NULL, 
            xtf=NULL, 
            prior, 
            mcmc,
            state,
            frailty=NULL,
            Proximity=NULL,
            data=sys.frame(sys.parent()),
            na.action=na.fail, 
            work.dir=NULL) {
    #########################################################################################
    # call parameters
    #########################################################################################
    m <- mcall <- cl <- match.call()
    
    #########################################################################################
    # data structure
    #########################################################################################
    tobs <- as.matrix(tobs);
    nrec <- length(tobs[,1]);
    xce <- as.matrix(cbind(rep(1,nrec), x));   pce <- ncol(xce); colnames(xce)[1] = "intercept";
    xtf <- as.matrix(cbind(rep(1,nrec), xtf)); ptf <- ncol(xtf); colnames(xtf)[1] = "intercept";
    if(!is.null(frailty)) {
      ID = frailty;
      orderindex = order(ID); 
      if(!(sum(orderindex==(1:nrec))==nrec)) stop("please sort the data by ID");
      blocki = c(0, cumsum(as.vector(table(ID))));
      if(is.null(Proximity)){
        W = matrix(0, length(blocki)-1, length(blocki)-1);
        D = rep(1, length(blocki)-1);
      }else{
        W = Proximity;
        D = rowSums(W);
      }
    }
    type <- rep(2, nrec);
    type[tobs[,1]==0] <- 1;
    type[tobs[,2]==Inf] <- 3;
    type[tobs[,1]==tobs[,2]] <- 4;
    nevents = length(which(type==4));
    
    #########################################################################################
    # change working directory (if requested..)
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n Changing working directory to ",work.dir,"\n")
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    
    #########################################################################################
    # prediction
    #########################################################################################
    
    
    #########################################################################################
    # initial MLE analysis and mcmc parameters
    #########################################################################################
    tmat <- tobs; 
    tmat[which(tobs[,1]==0),1]=0;
    tmat[which(tobs[,2]==Inf),2]=1e3;
    ll <- tmat[,1]; rr <- tmat[,2];
    fit0 <- survival::survreg(formula = survival::Surv(time=ll, time2=rr, type="interval2") ~ xce-1, dist = "lognormal");
    nburn <- mcmc$nburn;
    nsave <- mcmc$nsave;
    nskip <- mcmc$nskip;
    ndisplay <- mcmc$ndisplay;
    
    #########################################################################################
    # priors
    #########################################################################################
    maxL <- prior$maxL; if(is.null(maxL)) maxL<-5;
    ntprob <- 2^(maxL+1)-2; 
    ntlr <- 2^maxL-1;
    if(is.null(prior$a0)){
      a0=-1; b0=-1;
      alpha=state$alpha; if(is.null(alpha)) stop("please specify state$alpha if prior$a0 is missed");
    } else {
      a0=prior$a0; b0=prior$b0; if(is.null(b0)) b0=1;
      alpha=state$alpha; if(is.null(alpha)) alpha=5;
    }
    m0 <- prior$m0; if(is.null(m0)) m0 <- rep(0,pce);
    S0 <- prior$S0; if(is.null(S0)) S0 <- diag(rep(1e5,pce), nrow=pce, ncol=pce);
    S0inv <- solve(S0);
    if(is.null(prior$siga0)){
      siga0=-1; sigb0=-1;
    } else {
      siga0=prior$siga0; sigb0=prior$sigb0; if(is.null(sigb0)) sigb0=1;
    }
    gprior <- prior$gprior; if(is.null(gprior)) gprior <- 2*nrec*solve(t(xtf)%*%xtf);
    
    #########################################################################################
    # current state
    #########################################################################################
    betace <- state$betace; if(is.null(betace)) betace <- coefficients(fit0);
    sigma2 <- state$sigma2; if(is.null(sigma2)) sigma2 <- fit0$scale;
    betatf <- state$betatf; if(is.null(betatf)) betatf <- matrix(0,nrow=ptf,ncol=ntlr); 
    y <- state$y; 
    if(is.null(y)){
      y <- rep(0, nrec);
      for(i in 1:nrec){
        if(type[i]==1) y[i] = log(tobs[i,2]/2);
        if(type[i]==2) y[i] = log(mean(tobs[i,]));
        if(type[i]==3) y[i] = log(tobs[i,1]+1);
        if(type[i]==4) y[i] = log(tobs[i,1]);
      }
    }
    if(!is.null(frailty)){
      if(is.null(state$frail)) {
        v <- rep(0, length(blocki)-1);
      } else {
        v <- state$frail; if(length(v)!=(length(blocki)-1)) stop("check the length of frail");
      }
      taua0 = prior$taua0; if(is.null(taua0)) taua0=0.1;
      taub0 = prior$taub0; if(is.null(taub0)) taub0=0.1;
      tau2 = state$tau2; if(is.null(tau2)) tau2=1;
    }
        
    #########################################################################################
    # calling the c++ code and # output
    #########################################################################################
    if(is.null(frailty)){
      foo <- .Call("nonfrailtyLDTFP", 
                   nburn_ = nburn, 
                   nsave_ = nsave, 
                   nskip_ = nskip, 
                   ndisplay_ = ndisplay,
                   tobs_ = tobs, 
                   type_ = type, 
                   xce_ = t( xce ),
                   xtf_ = t( xtf ),
                   alpha_ = alpha, 
                   betace_ = betace, 
                   betatf_ = betatf, 
                   sigma2_ = sigma2,
                   y_ = y, 
                   maxL_ = maxL,
                   a0_ = a0, 
                   b0_ = b0, 
                   m0_ = m0, 
                   S0inv_ = S0inv, 
                   gprior_ = gprior,
                   a0sig_ = siga0, 
                   b0sig_ = sigb0, 
                   PACKAGE = "spBayesSurv")
    }else{
      foo <- .Call("frailtyLDTFP", 
                   nburn_ = nburn, 
                   nsave_ = nsave, 
                   nskip_ = nskip, 
                   ndisplay_ = ndisplay,
                   tobs_ = tobs, 
                   type_ = type, 
                   xce_ = t( xce ),
                   xtf_ = t( xtf ),
                   alpha_ = alpha, 
                   betace_ = betace, 
                   betatf_ = betatf, 
                   sigma2_ = sigma2,
                   y_ = y, 
                   v_ = v, 
                   blocki_ = blocki,
                   tau2_ = tau2,
                   W_ = W, 
                   D_ = D,
                   maxL_ = maxL,
                   a0_ = a0, 
                   b0_ = b0, 
                   m0_ = m0, 
                   S0inv_ = S0inv, 
                   gprior_ = gprior,
                   a0sig_ = siga0, 
                   b0sig_ = sigb0, 
                   a0tau_ = taua0, 
                   b0tau_ = taub0, 
                   PACKAGE = "spBayesSurv")
    }
    
    #########################################################################################
    # save state
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n\n Changing working directory back to ",old.dir,"\n")
      setwd(old.dir)
    }
    model.name <- "Generalized accelerated failure time model for time-to-event data";
    
    ### Calculate Bayes Factors for betatf
    BFs = .BayesFactor(foo$betatf, maxL, a0, b0, gprior, mean(foo$alpha));
    BayesFactors = c(as.vector( BFs$BFindividual ), BFs$BFoverallLDTFP, BFs$BFoverallParam); 
    names(BayesFactors) = c( colnames(xtf), "Covariate Dependency", "Normality");
    
    #### coefficients
    if(is.null(frailty)){
      coeff <- c(apply(foo$beta, 1, mean), mean(foo$sigma2), mean(foo$alpha));
      names(coeff) = c(colnames(xce), "sigma2", "alpha");
    }else{
      coeff <- c(apply(foo$beta, 1, mean), mean(foo$sigma2), mean(foo$alpha), mean(foo$tau2));
      names(coeff) = c(colnames(xce), "sigma2", "alpha", "tau2");
    }
    
    #### Save to a list
    output <- list(modelname=model.name,
                   coefficients=coeff,
                   BF = BayesFactors,
                   call=cl,
                   prior=prior,
                   mcmc=mcmc,
                   nrec=nrec,
                   nevents=nevents,
                   pce=pce,
                   ptf=ptf,
                   tobs=tobs,
                   X=xce,
                   Xtf=xtf,
                   maxL = maxL,
                   frailty = frailty,
                   beta = foo$beta,
                   sigma2 = foo$sigma2,
                   alpha = foo$alpha,
                   betatf = foo$betatf,
                   y = foo$y,
                   cpo = foo$cpo,
				   accept_beta = foo$ratebetace,
				   accept_betatf = foo$ratebetatf);
    if(!is.null(frailty)){
      output$v = foo$v;
      output$tau2 = foo$tau2;
      output$IDnames = names(table(ID));
    }
    cat("\n\n")
    class(output) <- c("frailtyGAFT")
    output
  }

#### print, summary, plot
"print.frailtyGAFT" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior Inference of Parameters:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nBayes Factors for each LDFTP covariate effect:\n")       
  print.default(format(x$BF, digits = digits), print.gap = 2, 
                quote = FALSE) 
  
  cat("\nNumber of Observations:",x$nrec)
  cat("\nNumber of Events:",x$nevents,"\n") 
  invisible(x)
}

"summary.frailtyGAFT" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### CPO
  ans$cpo <- object$cpo
  
  ### Median information
  mat <- as.matrix(object$beta)
  dimen1 <- object$pce
  coef.p <- object$coefficients[1:dimen1];
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
  
  ### Baseline Information
  mat <- object$sigma2
  coef.p <- object$coefficients[(dimen1+1)]    
  coef.m <- median(mat)    
  coef.sd <- sd(mat)
  limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
  coef.l <- limm[1]
  coef.u <- limm[2]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%HPD-Low", sep=""),
                                                paste(CI.level*100, "%HPD-Upp", sep="")))
  ans$basevariance <- coef.table
  
  ### Precision parameter
  if(is.null(object$prior$a0)){
    ans$prec <- NULL
  }else{
    mat <- object$alpha
    coef.p <- object$coefficients[(dimen1+2)]    
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
  if(is.null(object$frailty)){
    ans$frailvar <- NULL
  }else{
    mat <- object$tau2
    coef.p <- object$coefficients[(dimen1+3)]    
    coef.m <- median(mat)    
    coef.sd <- sd(mat)
    limm <- as.vector(coda::HPDinterval(coda::mcmc(mat), prob=CI.level))
    coef.l <- limm[1]
    coef.u <- limm[2]
    
    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                  paste(CI.level*100, "%HPD-Low", sep=""),
                                                  paste(CI.level*100, "%HPD-Upp", sep="")))
    ans$frailvar <- coef.table
  }
  ans$nrec <- object$nrec
  ans$nevents <- object$nevents
  ans$pce <- object$pce
  ans$ptf <- object$ptf
  ans$BF <- object$BF
  
  class(ans) <- "summaryfrailtyGAFT"
  return(ans)
}


"print.summaryfrailtyGAFT"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior Inference of Regression Parameters:\n")
  print.default(format(x$coeff, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nPosterior Inference of Baseline Variance:\n")
  print.default(format(x$basevariance, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  if (!is.null(x$prec)) {
    cat("\nPrecision parameter:\n")
    print.default(format(x$prec, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$frailvar)){
    cat("\nVariance parameter of frailties:\n")
    print.default(format(x$frailvar, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  
  cat("\nBayes Factors for each LDFTP covariate effect:\n")       
  print.default(format(x$BF, digits = digits), print.gap = 2, 
                quote = FALSE) 
  
  cat("\nNumber of Observations:",x$nrec)
  cat("\nNumber of Events:",x$nevents) 
  cat("\nNumber of Predictors for the Median Regression:",x$pce)    
  cat("\nNumber of Predictors for the Tailfree Probabilities:",x$pce,"\n")    
  invisible(x)
}