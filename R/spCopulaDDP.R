"spCopulaDDP" <- function(y, delta, x=NULL, s, prediction, prior, mcmc, state, FSA = TRUE, knots,
                    data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
UseMethod("spCopulaDDP")

"spCopulaDDP.default" <- 
function (y, 
          delta, 
          x=NULL, 
          s, 
          prediction, 
          prior, 
          mcmc,
          state,
          FSA = TRUE,
          knots,
          data=sys.frame(sys.parent()),
          na.action=na.fail, 
          work.dir=NULL) {
  #########################################################################################
  # call parameters
  #########################################################################################
  m <- mcall <- cl <- match.call()
  
  #########################################################################################
  # data structure and check if FSA is needed
  #########################################################################################
  if(FSA) {
    ss <- knots$ss; if(is.null(ss)) stop("please specify the knots vector ss if FSA is used");
    ss <- t(as.matrix(ss));
    id <- knots$blockid; if(is.null(id)) stop("please specify the bolock id that each observation is located in if FSA is used ");
    orderindex = order(id); id = id[orderindex];
    y <- as.vector(y); y = y[orderindex]; n = length(y);
    delta <- as.vector(delta); delta = delta[orderindex];
    X <- t(cbind(rep(1,n), x)); X = X[,orderindex]; p <- nrow(X);
    s <- t(as.matrix(s)); s = s[,orderindex];
    dnn <- .DistMat(s, s);
    dnm <- .DistMat(s, ss);
    dmm <- .DistMat(ss, ss);
    blocki = c( 0, cumsum(as.vector(table(id))) );
    #########################################################################################
    # prediction
    #########################################################################################
    s0 <- t(prediction$spred);
    xnew <- prediction$xpred;
    if(is.null(xnew)) stop("please specify xpred");
    if(is.vector(xnew)) {
      npred = length(xnew)
    } else npred = nrow(xnew);
    xpred <- cbind(rep(1,npred), xnew);
    if(!(ncol(xpred)==p)) { stop ("error: ncol(xpred) is not equal to ncol(x)");} 
    ds0n <- .DistMat(s, s0);
    ds0m <- .DistMat(ss, s0);
    ds0block <- matrix(0, n, npred);
    predid <- prediction$predid; 
    if(is.null(predid)) stop("please specify the pred id that each observation is located in if FSA is used ");
    for(i in 1:n){
      for(j in 1:npred){
        ds0block[i,j] = (id[i]==predid[j])+0
      }
    }
  } else{
    y <- as.vector(y); n <- length(y);
    delta <- as.vector(delta);
    X <- t(cbind(rep(1,n), x)); p <- nrow(X);
    s <- t(s);
    dnn <- .DistMat(s, s);
  }

  #########################################################################################
  # change working directory (if requested..)
  #########################################################################################
  if(!is.null(work.dir))
  {
    cat("\n Changing working directory to ",work.dir,"\n")
    old.dir <- getwd()  # by default work in current working directory
    setwd(work.dir)
  }
  model.name <- "Spatial copula model for point-referenced time-to-event data"
  
  #########################################################################################
  # prediction
  #########################################################################################
  s0 <- t(prediction$spred);
  xnew <- prediction$xpred;
  if(is.null(xnew)) stop("please specify xpred");
  if(is.vector(xnew)) {
    npred = length(xnew)
  } else npred = nrow(xnew);
  xpred <- cbind(rep(1,npred), xnew);
  if(!(ncol(xpred)==p)) { stop ("error: ncol(xpred) is not equal to ncol(x)");} 
  ds0n <- .DistMat(s, s0);
  
  #########################################################################################
  # initial analysis and mcmc parameters
  #########################################################################################
  fit0 <- survival::survreg(formula = Surv(exp(y), delta) ~ t(X)-1, dist = "lognormal");
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  
  #########################################################################################
  # priors
  #########################################################################################
  N <- prior$N; if(is.null(N)) N <- 10;
  m0 <- prior$m0; if(is.null(m0)) m0 <- rep(0,p);
  S0 <- prior$S0; if(is.null(S0)) S0 <- diag(rep(1e5,p), nrow=p, ncol=p);
  Sig0 <- prior$Sig0; if(is.null(Sig0)) Sig0 <- diag(rep(1e5,p), nrow=p, ncol=p);
  k0 <- prior$k0; if(is.null(k0)) k0 <- 7;
  nua <-prior$nua; nub <- prior$nub;
  if(is.null(nua)) nua=2; if(is.null(nub)) nub=1;
  a0 <-prior$a0; b0 <- prior$b0;
  if(is.null(a0)) a0=1; if(is.null(b0)) b0=1;
  theta0 <- prior$theta0; if(is.null(theta0)) theta0 <- c(1.0, 1.0, 0.5, 0.2);
  spl0 <- prior$spl0; if(is.null(spl0)) spl0 <- round(nburn/2);
  spS0 <- prior$spS0; if(is.null(spS0)) spS0 <- diag(c(4,4));
  spadapter <- prior$spadapter; if(is.null(spadapter)) spadapter <- (2.4)^2/2;

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
  theta = state$theta; if(is.null(theta)) theta <- c(0.95, 1);

  #########################################################################################
  # calling the c++ code
  #########################################################################################
  if(FSA){
    foo <- .Call("spCopulaDDP_FSA", 
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
                 ds0n_ = ds0n,
                 dnn_ = dnn,
                 theta_ = theta,
                 theta0_ = theta0,
                 spl0_ = spl0,
                 spS0_ = spS0, 
                 spadapter_ = spadapter,
                 dnm_ = dnm, dmm_=dmm, blocki_=blocki,
                 ds0m_= ds0m, ds0block_=ds0block,
                 PACKAGE = "spBayesSurv")
  } else{
    foo <- .Call("spCopulaDDP", 
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
                 ds0n_ = ds0n,
                 dnn_ = dnn,
                 theta_ = theta,
                 theta0_ = theta0,
                 spl0_ = spl0,
                 spS0_ = spS0, 
                 spadapter_ = spadapter,
                 PACKAGE = "spBayesSurv")
  }
  #########################################################################################
  # output
  #########################################################################################
  output <- list(modelname=model.name,
                 beta = foo$beta,
                 sigma2 = foo$sigma2,
                 w = foo$w,
                 alpha = foo$alpha,
                 theta1 = foo$theta1,
                 theta2 = foo$theta2,
                 z = foo$z,
                 y = foo$y,
                 ratey = foo$ratey,
                 ratebeta = foo$ratebeta,
                 ratesigma = foo$ratesigma,
                 rateV = foo$rateV,
                 ratetheta = foo$ratetheta,
                 cpo = foo$cpo,
                 Ypred = foo$Ypred,
                 Zpred = foo$Zpred);
  
  cat("\n\n")
  class(output) <- c("spCopulaDDP")
  output
}