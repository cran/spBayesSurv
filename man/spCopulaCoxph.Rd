\name{spCopulaCoxph}
\alias{spCopulaCoxph}
\alias{spCopulaCoxph.default}

\title{A Spatial Copula Approach to Cox Proportional Hazards Model}
\description{
    This function fits a spatial copula Cox PH model for survival analysis of 
    point-referenced right censored time-to-event data.
}
   
\usage{
spCopulaCoxph(y, delta, x=NULL, s, prediction, prior, mcmc, state, RandomIntervals=F,
      data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
}

\arguments{
   \item{y}{         an n by 1 vector giving the survival times.}
   
   \item{delta}{	 an n by 1 vector indicating whether it is right censored (=0) or not (=1).}
   
   \item{x}{         an n by p matrix of covariates without intercept. The default is NULL, 
					 indicating no covariates included.} 
					 
   \item{s}{         an n by d matrix of UMT coordinates, where d is the dimension of space.}

   \item{prediction}{a list giving the information used to obtain conditional
                     inferences. The list includes the following
                     elements: \code{spred} and \code{xpred} giving the n by 2 new locations 
					 and corresponding n by p covariates matrix, respectively, used for prediction .}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{M} an integer giving the total number of cut points for 
                     baseline hazard. See Zhou, Hanson and Knapp (2014+) 
                     for more detailed hyperprior specifications.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following elements: \code{nburn} an integer giving the number of burn-in 
                     scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer 
                     giving the total number of scans to be saved, \code{ndisplay} an integer giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.} 
    \item{RandomIntervals}{ indicate if random cut ponts for baseline hazard is need. The default is TRUE.}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{spCopulaCoxph} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This generic function fits a spatial copula Cox PH model (Zhou, Hanson and Knapp, 2014+),
  for (potentially) point-referenced right-censored data. 
}

\value{
  The results include the MCMC chains for the parameters discussed in Zhou, Hanson and Knapp (2014+). 
  Use \code{names} to find out what they are.

}

\seealso{
\code{\link{spCopulaDDP}}
}

\references{
Zhou, H., Hanson, T. and Knapp, R. (2014+). A Spatial Copula Approach to Fully Bayesian Nonparametric Survival Analysis. 
In preparation for publication.
}


\author{
Haiming Zhou \email{<zhouh@email.sc.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{Spatial copula Cox PH}
\keyword{Bayesian}

\examples{
\dontrun{

###############################################################
# A simulated data: spatial Copula Cox PH
###############################################################
rm(list=ls())
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(survival)
library(spBayesSurv)
## True parameters 
betaT = c(-1); 
theta1 = 0.98; theta2 = 0.1;

## generate coordinates: 
## npred is the # of locations for prediction
n = 300; npred = 30; ntot = n + npred;
ldist = 100; wdist = 40;
s1 = runif(ntot, 0, wdist); s2 = runif(ntot, 0, ldist);
s = rbind(s1,s2); #plot(s[1,], s[2,]);

## Covariance matrix
corT = matrix(1, ntot, ntot);
for (i in 1:(ntot-1)){
  for (j in (i+1):ntot){
    dij = sqrt(sum( (s[,i]-s[,j])^2 ));
    corT[i,j] = theta1*exp(-theta2*dij);
    corT[j,i] = theta1*exp(-theta2*dij);
  }
}

## Generate x 
x = runif(ntot,-1.5,1.5);
## Generate transformed log of survival times
z = mvrnorm(1, rep(0, ntot), corT);
## The CDF of Ti: Lambda(t) = t^2;
Fi = function(t, xi){
  res = 1-exp(-t^3*exp(sum(xi*betaT)));
  res[which(t<0)] = 0;
  res
}
## The pdf of Ti:
fi = function(t, xi){
  res=(1-Fi(t,xi))*3*t^2*exp(sum(xi*betaT));
  res[which(t<0)] = 0;
  res
}
#integrate(function(x) fi(x, 0), -Inf, Inf) 
## true plot
xx = seq(0, 10, 0.1)
plot(xx, fi(xx, -1), "l", lwd=2, col=2)
lines(xx, fi(xx, 1), "l", lwd=2, col=3)

## The inverse for CDF of Ti
Finvsingle = function(u, xi) {
  res = uniroot(function (x) Fi(x, xi)-u, lower=0, upper=5000);
  res$root
}
Finv = function(u, xi) {sapply(u, Finvsingle, xi)};
## Generate survival times t
u = pnorm(z);
t = rep(0, ntot);
for (i in 1:ntot){
  t[i] = Finv(u[i], x[i]);
}
tTrue = t; #plot(x,t);

## Censoring scheme
Centime = runif(ntot, 1, 5);
#Centime = 10000;
delta = (t<=Centime) +0 ;
sum(delta)/ntot;
cen = which(delta==0);
t[cen] = Centime[cen];

## make a data frame
dtotal = data.frame(s1=s1, s2=s2, t=t, logt=log(t), x=x, delta=delta, 
                    tTrue=tTrue, logtTrue=log(tTrue));
## Hold out npred=30 for prediction purpose
predindex = sample(1:ntot, npred);
dpred = dtotal[predindex,];
dtrain = dtotal[-predindex,];

# Prediction settings 
xpred = dpred$x;
s0 = cbind( dpred$s1, dpred$s2 );
prediction = list(spred=s0, xpred=xpred);

###############################################################
# spatial Copula Cox PH
###############################################################
# rename the variables 
d = dtrain;n=nrow(d); n;
s = cbind(d$s1, d$s2);
t = d$t;
x = d$x;
X = cbind(x);
p = ncol(X); # number of covariates + 1
delta =d$delta;

# Initial Exp Cox PH
fit0 = survreg( Surv(t, delta) ~ x, dist="exponential", data=d );
h0 = as.vector( exp( -fit0$coefficients[1] ) );
beta = as.vector( -fit0$coefficients[-1] );
S0 = as.matrix( fit0$var[-1,-1] );

# Prior information
prior = list(M = 10,
             h0 = h0,
             r0 = 1,
             V0 = 2*as.vector( exp( -2*fit0$coefficients[1] )*fit0$var[1,1] ),
             mu0 = as.vector( -fit0$coefficients[-1] ),
             Sig0 = diag(rep(1e5,p), nrow=p, ncol=p),
             S0 = as.matrix( fit0$var[-1,-1] ),
             adapter = 2,
             hs0 = sqrt( as.vector( exp( -2*fit0$coefficients[1] )*fit0$var[1,1] ) ),
             hadapter = (2.4)^2,
             spS0 = diag(c(2,4)),
             spadapter = (2.4)^2/2, 
             theta0=c(1.0, 1.0, 0.5, 0.2));

# current state values
state <- list(theta = c(0.9, 1));

# MCMC parameters
nburn <- 1000
nsave <- 500
nskip <- 4
ndisplay <- 100
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)

# Fit the spatial Cox PH model
res2 = spCopulaCoxph( y = t,
                    delta =delta, 
                    s = s,
                    x = X, 
                    RandomIntervals=T,
                    prediction=prediction, 
                    prior=prior, 
                    mcmc=mcmc,
                    state=state);
par(mfrow = c(2,1))
h.save = res2$h; rowMeans(h.save);
traceplot(mcmc(h.save[2,]), main="h")
beta.save = res2$beta; rowMeans(beta.save);
traceplot(mcmc(beta.save[1,]), main="beta")
res2$ratebeta;
res2$ratetheta;
res2$rateh;
## LPML
LPML2 = sum(log(res2$cpo)); LPML2;
## MSPE
mean((dpred$tTrue-apply(res2$Tpred, 1, median))^2); 
## plots
par(mfrow = c(2,1))
xnew = c(-1, 1)
xpred = cbind(xnew); 
nxpred = nrow(xpred);
ygrid = seq(-5, 1.1, 0.05);tgrid = exp(ygrid);
ngrid = length(ygrid);
estimates = GetCurves(res2, xpred, ygrid, CI=c(0.05, 0.95));
fhat = estimates$fhat; 
Shat = estimates$Shat;
## density in t
plot(tgrid, fi(tgrid, xnew[1]), "l", lwd=2,  xlim=c(0,6), main="density in y")
for(i in 1:nxpred){
  lines(tgrid, fi(tgrid, xnew[i]), lwd=2)
  lines(tgrid, fhat[,i], lty=2, lwd=2, col=4);
}
## survival in t
plot(tgrid, 1-Fi(tgrid, xnew[1]), "l", lwd=2, xlim=c(0,6), main="survival in y")
for(i in 1:nxpred){
  lines(tgrid, 1-Fi(tgrid, xnew[i]), lwd=2)
  lines(tgrid, Shat[,i], lty=2, lwd=2, col=4);
}

}
}

