\name{spCopulaDDP}
\alias{spCopulaDDP}

\title{Marginal Bayesian Nonparametric Survival Model via Spatial Copula}
\description{
This function fits a marginal Bayesian Nonparametric model (Zhou, Hanson and Knapp, 2015) for point-referenced right censored time-to-event data.
}
   
\usage{
spCopulaDDP(formula, data, na.action, prediction=NULL,
            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
            prior=NULL, state=NULL, scale.designX=TRUE,
            Coordinates, DIST=NULL, Knots=NULL)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It currently only supports right-censoring.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{prediction}{a list giving the information used to obtain conditional inferences. The list includes the following elements: \code{spred} and \code{xpred} giving the n by 2 new locations and corresponding npred by p covariates matrix, respectively, used for prediction. If \code{prediction=NULL}, \code{xpred} will be set to be the design matrix used in \code{formula}, and \code{spred} will be set to be in \code{Coordinates}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
    
    \item{prior}{a list giving the prior information. The list includes the following elements: \code{N} an integer giving the truncation of the Dirichlet process. See Zhou, Hanson and Zhang (2017) for more detailed hyperprior specifications.}

    \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.}  
    
    \item{scale.designX}{flag to indicate wheter the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
    
    \item{Coordinates}{an n by 2 coordinates matrix, where n is the sample size, 2 is the dimension of coordiantes. Note all cocordinates should be distinct. }
    
    \item{DIST}{This is a function argument, used to calculate the distance. The default is Euclidean distance (\code{fields::rdist}). This function should have two arguments (X1,X2), where X1 and X2 are matrices with coordinates as the rows. The returned value of this function should be the pairwise distance matrix. If nrow(X1)=m and nrow(X2)=n then the function should return an m by n matrix of all distances between these two sets of points.}
    
    \item{Knots}{an nknots by 2 matrix, where nknots is the number of selected knots for FSA, and 2 is the dimension of each location. If \code{Knots} is not specified, the space-filling algorithm will be used to find the knots.}
    
}

\details{
This function fits a marginal Bayesian Nonparametric model (Zhou, Hanson and Knapp, 2015) for point-referenced right censored time-to-event data. 
}

\value{
  The results include the MCMC chains for the parameters discussed in Zhou, Hanson and Zhang (2017).
  Use \code{names} to find out what they are.

}

\seealso{
\code{\link{anovaDDP}, \link{GetCurves}}
}

\references{
Zhou, H., Hanson, T., and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, 71(4): 1101-1110. 
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{Spatial copula}
\keyword{Bayesian nonparametric}


\examples{
###############################################################
# A simulated data: mixture of two normals with spatial dependence
###############################################################
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
## True parameters 
betaT = cbind(c(3.5, 0.5), c(2.5, -1)); 
wT = c(0.4, 0.6); 
sig2T = c(1^2, 0.5^2);
theta1 = 0.98; theta2 = 0.1;
n=30; npred=3; ntot = n+npred;
## The Survival function for log survival times:
fiofy = function(y, xi, w=wT){
  nw = length(w);
  ny = length(y);
  res = matrix(0, ny, nw);
  Xi = c(1,xi);
  for (k in 1:nw){
    res[,k] = w[k]*dnorm(y, sum(Xi*betaT[,k]), sqrt(sig2T[k]) )
  }
  apply(res, 1, sum)
}
fioft = function(t, xi, w=wT) fiofy(log(t), xi, w)/t;
Fiofy = function(y, xi, w=wT){
  nw = length(w);
  ny = length(y);
  res = matrix(0, ny, nw);
  Xi = c(1,xi);
  for (k in 1:nw){
    res[,k] = w[k]*pnorm(y, sum(Xi*betaT[,k]), sqrt(sig2T[k]) )
  }
  apply(res, 1, sum)
}
Fioft = function(t, xi, w=wT) Fiofy(log(t), xi, w);
## The inverse for Fioft
Finv = function(u, x) uniroot(function (y) Fiofy(y,x)-u, lower=-250, 
                              upper=250, extendInt ="yes", tol=1e-6)$root

## generate coordinates: 
## npred is the # of locations for prediction
ldist = 100; wdist = 40;
s1 = runif(ntot, 0, wdist); s2 = runif(ntot, 0, ldist);
s = cbind(s1,s2); #plot(s[,1], s[,2]);
## Covariance matrix
corT = matrix(1, ntot, ntot);
for (i in 1:(ntot-1)){
  for (j in (i+1):ntot){
    dij = sqrt(sum( (s[i,]-s[j,])^2 ));
    corT[i,j] = theta1*exp(-theta2*dij);
    corT[j,i] = theta1*exp(-theta2*dij);
  }
}

## generate x 
x1 = runif(ntot,-1.5,1.5); X = cbind(x1);
## generate transformed log of survival times
z = MASS::mvrnorm(1, rep(0, ntot), corT);
## generate survival times
u = pnorm(z);
tT = rep(0, ntot);
for (i in 1:ntot){
  tT[i] = exp(Finv(u[i], X[i,]));
}

### ----------- right-censored -------------###
t_obs=tT 
Centime = runif(ntot, 200, 500);
delta = (tT<=Centime) +0 ; 
length(which(delta==0))/ntot; # censoring rate
rcen = which(delta==0);
t_obs[rcen] = Centime[rcen]; ## observed time 
## make a data frame
dtot = data.frame(tobs=t_obs, x1=x1, delta=delta, tT=tT,
                  s1=s1, s2=s2); 
## Hold out npred for prediction purpose
predindex = sample(1:ntot, npred);
dpred = dtot[predindex,];
d = dtot[-predindex,];
# Prediction settings 
prediction = list(xpred=cbind(dpred$x1), 
                  spred=cbind(dpred$s1, dpred$s2));

###############################################################
# Independent DDP: Bayesian Nonparametric Survival Model
###############################################################
# MCMC parameters
nburn=100; nsave=100; nskip=0;
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(N=10, a0=2, b0=2, nknots=n, nblock=round(n/2));
# here nknots=n, so FSA is not used.
# If nknots<n, FSA will be used with nblock=round(n/2).
# As nknots is getting larger, the FSA is more accurate but slower
# As nblock is getting smaller, the FSA is more accurate but slower. 
# In most applications, setting nblock=n works fine, which is the
# setting by not specifying nblock. 
# If nknots is not specified or nknots=n, the exact covariance is used. 
# Fit the Cox PH model
res1 = spCopulaDDP(formula = Surv(tobs, delta)~x1, data=d, 
                     prior=prior, mcmc=mcmc, prediction=prediction,
                     Coordinates=cbind(d$s1,d$s2), Knots=NULL);
# here if prediction=NULL, prediction$xpred will be set as the design matrix
# in formula, and prediction$spred will be set as the Coordinates argument. 
# Knots=NULL is the defualt setting, for which the knots will be generated 
# using fields::cover.design() with number of knots equal to prior$nknots. 
## LPML
LPML = sum(log(res1$cpo)); LPML;
## Number of non-negligible components
quantile(colSums(res1$w>0.05))
## MSPE
mean((log(dpred$tT)-apply(log(res1$Tpred), 1, median))^2); 

## traceplot
par(mfrow = c(1,2))
traceplot(mcmc(res1$theta[1,]), main="sill")
traceplot(mcmc(res1$theta[2,]), main="range")

############################################
## Curves
############################################
ygrid = seq(0,6.0,length=100); tgrid = exp(ygrid);
ngrid = length(tgrid);
xpred = rbind(-1, 1); 
estimates=plot(res1, xpred=xpred, tgrid=tgrid);

## plot
par(mfrow = c(1,2))
plot(tgrid, 1-Fioft(tgrid, xpred[1,]), "l", lwd=3, 
     main="Survival function", xlab="time", ylab="survival");
lines(tgrid, 1-Fioft(tgrid, xpred[2,]), "l", lwd=3);
lines(estimates$tgrid, estimates$Shat[,1], lty=2, lwd=3)
lines(estimates$tgrid, estimates$Shatlow[,1], lty=3, lwd=1)
lines(estimates$tgrid, estimates$Shatup[,1], lty=3, lwd=1)
lines(estimates$tgrid, estimates$Shat[,2], lty=2, lwd=3)
lines(estimates$tgrid, estimates$Shatlow[,2], lty=3, lwd=1)
lines(estimates$tgrid, estimates$Shatup[,2], lty=3, lwd=1)
plot(log(tgrid), fiofy(log(tgrid), xpred[1,]), "l", lwd=3, 
     main="Density of log time", xlab="log time", ylab="survival");
lines(log(tgrid), fiofy(log(tgrid), xpred[2,]), "l", lwd=3);
lines(log(estimates$tgrid), estimates$fhat[,1]*tgrid, lty=2, lwd=3)
lines(log(estimates$tgrid), estimates$fhatlow[,1]*tgrid, lty=3, lwd=1)
lines(log(estimates$tgrid), estimates$fhatup[,1]*tgrid, lty=3, lwd=1)
lines(log(estimates$tgrid), estimates$fhat[,2]*tgrid, lty=2, lwd=3)
lines(log(estimates$tgrid), estimates$fhatlow[,2]*tgrid, lty=3, lwd=1)
lines(log(estimates$tgrid), estimates$fhatup[,2]*tgrid, lty=3, lwd=1)
}
