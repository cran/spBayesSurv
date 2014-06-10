"anovaDDP" <- function(y, delta, x=NULL, prediction, prior, mcmc, state, 
                    data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
  UseMethod("anovaDDP")

"anovaDDP.default" <- 
  function (y, 
            delta, 
            x=NULL,  
            prediction, 
            prior, 
            mcmc,
            state,
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
    y <- as.vector(y);
    n <- length(y);
    X <- t(cbind(rep(1,n), x));
    p <- nrow(X);
    
    #########################################################################################
    # change working directory (if requested..)
    #########################################################################################
    if(!is.null(work.dir))
    {
      cat("\n Changing working directory to ",work.dir,"\n")
      old.dir <- getwd()  # by default work in current working directory
      setwd(work.dir)
    }
    model.name <- "ANOVA DDP model for point-referenced time-to-event data"
    
    #########################################################################################
    # prediction
    #########################################################################################
    xnew <- prediction$xpred;
    if(is.null(xnew)) return("please specify xpred")
    if(is.vector(xnew)) {
      npred = length(xnew)
    } else npred = nrow(xnew);
    xpred <- cbind(rep(1,npred), xnew);
    if(!(ncol(xpred)==p)) { return ("error: ncol(xpred) is not equal to ncol(x)");} 

    #########################################################################################
    # initial analysis and priors
    #########################################################################################
    fit0 <- survival::survreg(formula = Surv(exp(y), delta) ~ t(X)-1, dist = "lognormal")
    N <- prior$N; if(is.null(N)) N <- 10;
    m0 <- prior$m0; if(is.null(m0)) m0 <- rep(0,p);
    S0 <- prior$S0; if(is.null(S0)) S0 <- diag(rep(1e5,p), nrow=p, ncol=p);
    Sig0 <- prior$Sig0; if(is.null(Sig0)) Sig0 <- diag(rep(1e5,p), nrow=p, ncol=p);
    k0 <- prior$k0; if(is.null(k0)) k0 <- 7;
    nua <-prior$nua; nub <- prior$nub;
    if(is.null(nua)) nua=2; if(is.null(nub)) nub=1;
    a0 <-prior$a0; b0 <- prior$b0;
    if(is.null(a0)) a0=1; if(is.null(b0)) b0=1;
    
    #########################################################################################
    # current state and mcmc specification
    #########################################################################################
    mu <- state$mu; if(is.null(mu)) mu <- rep(0,p);
    Sig <- state$Sig; if(is.null(Sig)) Sig <- diag(rep(1e5,p), nrow=p, ncol=p); 
    beta<- state$beta; if(is.null(beta)) beta = matrix(coefficients(fit0), p, N);
    sigma2<- state$sigma2; if(is.null(sigma2)) sigma2 = rep(fit0$scale,N);
    alpha <- state$alpha; if(is.null(alpha)) alpha <- 1;
    K = sample(1:N,n, replace=T);
    V = rbeta(N, 1, alpha); V[N] =1;
    w = V; 
    for (k in 2:N){
      w[k] = max( (1 - sum(w[1:(k-1)]))*V[k], 1e-20);
    }
    nburn <- mcmc$nburn;
    nsave <- mcmc$nsave;
    nskip <- mcmc$nskip;
    ndisplay <- mcmc$ndisplay;
    
    #########################################################################################
    # calling the c++ code
    #########################################################################################
    foo <- .Call("anovaDDP", 
                 nburn_ = nburn, 
                 nsave_ = nsave, 
                 nskip_ = nskip, 
                 ndisplay_ = ndisplay,
                 y_ = y,
                 delta_ = delta, 
                 X_ = as.matrix(X), 
                 N_ = N,
                 beta_ = beta, 
                 tau2_ = 1.0/sigma2,
                 K_ = K, 
                 V_ = V, 
                 w_ = w,
                 alpha_ = alpha, 
                 mu_ = mu, 
                 Sig_ = Sig,
                 m0_ = m0, 
                 S0_ = S0, 
                 Sig0_ = Sig0, 
                 k0_ = k0,
                 a0_ = a0, 
                 b0_ = b0, 
                 nua_ = nua, 
                 nub_ = nub,
                 xpred_ = as.matrix(xpred),  
                 PACKAGE = "spBayesSurv")
    
    #########################################################################################
    # output
    #########################################################################################
    output <- list(modelname=model.name,
                   beta = foo$beta,
                   sigma2 = foo$sigma2,
                   w = foo$w,
                   alpha = foo$alpha,
                   y = foo$y,
                   cpo = foo$cpo,
                   Ypred = foo$Ypred);
    
    cat("\n\n")
    class(output) <- c("anovaDDP")
    output
  }
