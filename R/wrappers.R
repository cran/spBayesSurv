
.DDPplots <- function ( xpred_, ygrid_, beta_, sigma2_, w_, probs_ ) {
  .Call("DDPplots", xpred_, ygrid_, beta_, sigma2_, w_, probs_, PACKAGE = "spBayesSurv")
}

.DistMat <- function ( si_, sj_ ) {
  .Call("DistMat", si_, sj_, PACKAGE = "spBayesSurv")
}

.CoxPHplots <- function ( xpred_, tgrid_, beta_, h_, d_, probs_ ) {
  .Call("CoxPHplots", xpred_, tgrid_, beta_, h_, d_, probs_, PACKAGE = "spBayesSurv")
}

.rho_Gau <- function(dis, phi){
  exp(-(phi*dis)^2);
}

.rho_Exp <- function(dis, phi){
  exp(-abs(phi*dis));
}

.rho_pow <- function(dis, phi, nu){
  exp(-abs(phi*dis)^nu);
}

baseline = function (...) 
{
  words <- as.character((match.call())[-1]);
  allf <- cbind(...);
  nterms <- ncol(allf);
  if (is.null(names(allf))){
    argname <- words[1:nterms]
  }else{
    argname <- ifelse(names(allf)=="", words[1:nterms], names(allf));
  }
  allf
}

frailtyprior = function (prior="car", ...) {
  x = cbind(...);
  colnames(x)=rep(prior,ncol(x));
  x
}