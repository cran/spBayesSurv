
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
  allf <- data.frame(...);
  mtf <- model.frame(as.formula(paste("~", paste(words, collapse = "+"))), data = allf)
  Terms <- attr(mtf, 'terms')
  Xtf = model.matrix(Terms, mtf);
  adrop <- 0  #levels of "assign" to be dropped; 0= intercept
  Xatt <- attributes(Xtf) 
  xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
  Xtf <- Xtf[, !xdrop, drop=FALSE]
  attr(Xtf, "assign") <- Xatt$assign[!xdrop]
  attr(Xtf, "terms") <- Terms
  levelslist = list()
  k = 0
  for(i in 1:ncol(allf)){
    if(is.factor(allf[,i])){
      k = k+1
      levelslist[[k]] = levels(factor(allf[,i]));
      names(levelslist)[k] = words[i];
    }
  }
  attr(Xtf, "levels") = levelslist
  Xtf
}

frailtyprior = function (prior="car", ...) {
  x = cbind(...);
  colnames(x)=rep(prior,ncol(x));
  x
}

