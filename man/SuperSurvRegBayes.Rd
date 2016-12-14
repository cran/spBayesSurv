\name{SuperSurvRegBayes}
\alias{SuperSurvRegBayes}
\alias{print.SuperSurvRegBayes}
\alias{plot.SuperSurvRegBayes}
\alias{summary.SuperSurvRegBayes}
\alias{print.summary.SuperSurvRegBayes}

\title{Bayesian Semiparametric Super Survival Model}
\description{
This function fits a super survival model. It can fit both Case I and Case II interval censored data, as well as standard right-censored, uncensored, and mixtures of these. The Bernstein Polynomial Prior is used for fitting the baseline survival function. 
}
   
\usage{
SuperSurvRegBayes(formula, data, na.action, dist="lognormal", 
                  mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500),
                  prior=NULL, state=NULL, truncation_time=NULL, subject.num=NULL, 
                  InitParamMCMC=FALSE, scale.designX=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{dist}{centering distribution for MPT. Choices include \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The list includes: \code{maxL} an integer giving the maximum level of MPT, \code{a0} and \code{b0} gamma prior of the precision parameter of MPT, \code{Shat} the initial covariance matrix of adaptive M-H for each of the beta_h, beta_o and beta_q, \code{Vhat} the initial covariance matrix of adaptive M-H for theta, \code{beta0} and \code{S0} the prior for each of the beta_h, beta_o and beta_q, for which the default is the g-prior, \code{theta0} and \code{V0} the prior for theta.}
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{truncation_time}{a vector of left-trucation times with length n.}
   
   \item{subject.num}{a vector of suject id numbers when time dependent covariates are considered. For example, for subject 1 time dependent covariates are recorded over [0,1), [1,3), and for subject 2 covariates are recorded over [0,2), [2,3), [3,4). Suppose we only have two subjects, i.e. n=2. In this case, we save the data in the long format, set truncation_time=c(0,1,0,2,3) and subject.num=c(1,1,2,2,2).}
   
   \item{InitParamMCMC}{flag to indicate wheter an initial MCMC will be run based on the centering parametric model, where \code{TRUE} indicates yes.}
   
   \item{scale.designX}{flag to indicate wheter the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
}

\value{
  The results include the MCMC chains for the parameters; use \code{names} to find out what they are.

}

\seealso{
\code{\link{survregbayes}}
}

\references{
In preparation.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
#################################################################
# A simulated data based on PH_PO_AFT super model
#################################################################
rm(list=ls())
library(coda)
library(survival)
library(spBayesSurv)

## True coeffs
betaT_h = c(1, 1);
betaT_o = c(0, 0);
betaT_q = c(1, 1);
## Baseline Survival
f0oft = function(t) 0.5*dlnorm(t, -1, 0.5)+0.5*dlnorm(t,1,0.5);
S0oft = function(t) {
  0.5*plnorm(t, -1, 0.5, lower.tail=FALSE)+
    0.5*plnorm(t, 1, 0.5, lower.tail=FALSE)
}
h0oft = function(t) f0oft(t)/S0oft(t);
## The Survival function:
Sioft = function(t,x){
  xibeta_h = sum(x*betaT_h);
  xibeta_o = sum(x*betaT_o);
  xibeta_q = sum(x*betaT_q);
  (1+exp(xibeta_o-xibeta_h+xibeta_q)*
     (1/S0oft(exp(xibeta_q)*t)-1))^(-exp(xibeta_h-xibeta_q));
}
Fioft = function(t,x) 1-Sioft(t,x);
## The inverse for Fioft
Finv = function(u, x) uniroot(function (t) Fioft(t,x)-u, lower=1e-100, 
                              upper=1e100, extendInt ="yes", tol=1e-6)$root
### true plots
tt=seq(1e-10, 4, 0.02);
xpred1 = c(0,0);
xpred2 = c(0,1);
plot(tt, Sioft(tt, xpred1), "l", ylim=c(0,1));
lines(tt, Sioft(tt, xpred2), "l");

##-------------Generate data-------------------##
## generate x 
n = 80; 
x1 = rbinom(n, 1, 0.5); x2 = rnorm(n, 0, 1); X = cbind(x1, x2);
## generate survival times
u = runif(n);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = Finv(u[i], X[i,]);
}

### ----------- right censored -------------###
t1=tT;t2=tT;
Centime = runif(n, 2, 6);
delta = (tT<=Centime) +0 ; length(which(delta==0))/n;
rcen = which(delta==0);
t1[rcen] = Centime[rcen];
t2[rcen] = NA;
## make a data frame
d = data.frame(t1=t1, t2=t2, x1=x1, x2=x2, delta=delta, tT=tT); table(d$delta)/n;

##-------------Fit the model-------------------##
# MCMC parameters
nburn=200; nsave=500; nskip=0; niter = nburn+nsave
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(maxL=15, a0=1, b0=1, M=10, q=.9);
state <- list(cpar=1);
ptm<-proc.time()
res1 = SuperSurvRegBayes(formula = Surv(t1, t2, type="interval2")~x1+x2, data=d, 
                          prior=prior, mcmc=mcmc, state=state, dist="lognormal");
sfit=summary(res1); sfit;
systime1=proc.time()-ptm; systime1;
par(mfrow = c(3,2))
traceplot(mcmc(res1$beta_h[1,]), main="beta_h for x1");
traceplot(mcmc(res1$beta_h[2,]), main="beta_h for x2");
traceplot(mcmc(res1$beta_o[1,]), main="beta_o for x1");
traceplot(mcmc(res1$beta_o[2,]), main="beta_o for x2");
traceplot(mcmc(res1$beta_q[1,]), main="beta_q for x1");
traceplot(mcmc(res1$beta_q[2,]), main="beta_q for x2");

####################################################################
## Get curves
####################################################################
wide=0.2;
tgrid = seq(1e-10,4,wide);
ngrid = length(tgrid);
xpred = rbind(c(0,0), c(0,1)); 
estimates=plot(res1, xpred=xpred, tgrid=tgrid);
## plots
par(mfrow = c(1,1))
plot(tgrid, Sioft(tgrid, xpred[2,]), "l", lwd=3);
lines(tgrid, Sioft(tgrid, xpred[1,]), "l", lwd=3);
lines(estimates$tgrid, estimates$Shat[,1], lty=2, lwd=3)
lines(estimates$tgrid, estimates$Shat[,2], lty=2, lwd=3)
}

