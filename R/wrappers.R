
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

"b_spline" <- function(x, df=NULL, knots=NULL, Boundary.knots = range(x)){
  if(is.null(df)){
    if(is.null(knots)){
      df = 3;
    }else{
      df = length(knots)+2;
    }
  }
  if(df<2) stop("df has to be greater than 1 for cubic splines")
  df0 = df+2; 
  nknots0 = df0-3+1;
  #if(is.null(knots)){
  #  knots0 = seq(min(x), max(x), length.out=nknots0)[-c(1,nknots0)]
  #}else{
  #  knots0 = knots
  # if(length(knots0)!=(df-2)) stop("the number of internal knots has to be equal to df-2.")
  #}
  x.bs = splines::bs(x, df=df0, knots=knots, intercept=TRUE, Boundary.knots=Boundary.knots);
  x.attr = attributes(x.bs);
  x.attr$dim[2] = x.attr$dim[2]-2
  x.attr$dimnames[[2]] = (x.attr$dimnames[[2]])[-c(df0-1,df0)];
  x.bs2 = x.bs[,-c(1,df0)];
  attributes(x.bs2)=x.attr;
  class(x.bs2) <- c("b_spline")
  x.bs2
}

"predict.b_spline" <- function(object, newx, ...){
  if(is(object,"b_spline")){
    res=b_spline(newx, knots=attr(object,"knots"), 
             Boundary.knots=attr(object,"Boundary.knots"));
  }else{
    return(NULL)
  }
  res
}
