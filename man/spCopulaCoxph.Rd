\name{spCopulaCoxph}
\alias{spCopulaCoxph}
\alias{print.spCopulaCoxph}
\alias{summary.spCopulaCoxph}
\alias{print.summary.spCopulaCoxph}

\title{Marginal Bayesian Proportional Hazards Model via Spatial Copula}
\description{
This function fits a marginal Bayesian proportional hazards model (Zhou, Hanson and Zhang, 2018) for point-referenced right censored time-to-event data.
}
   
\usage{
spCopulaCoxph(formula, data, na.action, prediction=NULL, 
              mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
              prior=NULL, state=NULL, scale.designX=TRUE,
              Coordinates, DIST=NULL, Knots=NULL)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{\link[survival]{Surv}} function in the \code{survival} package. It currently only supports right-censoring.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{prediction}{a list giving the information used to obtain conditional inferences. The list includes the following elements: \code{spred} and \code{xpred} giving the n by 2 new locations and corresponding npred by p covariates matrix, respectively, used for prediction. If \code{prediction=NULL}, \code{xpred} will be set to be the design matrix used in \code{formula}, and \code{spred} will be set to be in \code{Coordinates}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
    
    \item{prior}{a list giving the prior information. See Zhou, Hanson and Zhang (2018) for more detailed hyperprior specifications.}

    \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.}  
    
    \item{scale.designX}{flag to indicate wheter the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
    
    \item{Coordinates}{an n by 2 coordinates matrix, where n is the sample size, 2 is the dimension of coordiantes. Note all cocordinates should be distinct. }
    
    \item{DIST}{This is a function argument, used to calculate the distance. The default is Euclidean distance (\code{fields::rdist}). This function should have two arguments (X1,X2), where X1 and X2 are matrices with coordinates as the rows. The returned value of this function should be the pairwise distance matrix. If nrow(X1)=m and nrow(X2)=n then the function should return an m by n matrix of all distances between these two sets of points.}
    
    \item{Knots}{an nknots by 2 matrix, where nknots is the number of selected knots for FSA, and 2 is the dimension of each location. If \code{Knots} is not specified, the space-filling algorithm will be used to find the knots.}

}

\details{
This function fits a marginal Bayesian proportional hazards model (Zhou, Hanson and Zhang, 2018) for point-referenced right censored time-to-event data.
}

\value{
The \code{spCopulaCoxph} object is a list containing at least the following components:
\item{modelname}{
the name of the fitted model
}
\item{terms}{
the \code{\link{terms}} object used
}
\item{coefficients}{
a named vector of coefficients.  
}
\item{call}{
the matched call
}
\item{prior}{
the list of hyperparameters used in all priors. 
}
\item{mcmc}{
the list of MCMC parameters used 
}
\item{n}{
the number of row observations used in fitting the model
}
\item{p}{
the number of columns in the model matrix 
}
\item{Surv}{
the \code{\link[survival]{Surv}} object used 
}
\item{X.scaled}{
the n by p scaled design matrix
}
\item{X}{
the n by p orginal design matrix
}
\item{beta}{
the p by nsave matrix of posterior samples for the coefficients in the \code{linear.predictors} 
}
\item{beta.scaled}{
the p by nsave matrix of posterior samples for the coefficients in the \code{linear.predictors}. Note that these posterior samples are based scaled design matrix.
}
\item{theta}{
the 2 by nsave matrix of posterior samples for sill and range parameters 
}
\item{ratebeta}{
the acceptance rate in the posterior sampling of beta coefficient vector
}
\item{ratetheta}{
the acceptance rate in the posterior sampling of theta 
}
\item{cpo}{
the length n vector of the stabilized estiamte of CPO; used for calculating LPML
}
\item{Coordinates}{
the \code{Coordinates} matrix used in \code{\link{survregbayes}} 
}
\item{Tpred}{
the npred by nsave predicted survival times for covariates specified in the argument \code{prediction}.
}
\item{Zpred}{
the npred by nsave predicted z values for covariates specified in the argument \code{prediction}.
}
}

\seealso{
\code{\link{spCopulaDDP}, \link{GetCurves}}
}

\references{
Zhou, H., Hanson, T., and Zhang, J. (2020). spBayesSurv: Fitting Bayesian Spatial Survival Models Using R. \emph{Journal of Statistical Software}, 92(9): 1-33.

Zhou, H., Hanson, T., and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, 71(4): 1101-1110. 
}


\author{
Haiming Zhou and Timothy Hanson
}

\keyword{Spatial copula Cox PH}
\keyword{Bayesian}

\examples{
###############################################################
# A simulated data: spatial Copula Cox PH
###############################################################
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
## True parameters 
betaT = c(1,1); 
theta1 = 0.98; theta2 = 0.1;
n=50; npred=3; ntot = n+npred;
## Baseline Survival
f0oft = function(t) 0.5*dlnorm(t, -1, 0.5)+0.5*dlnorm(t,1,0.5);
S0oft = function(t) (0.5*plnorm(t, -1, 0.5, lower.tail=FALSE)+
                       0.5*plnorm(t, 1, 0.5, lower.tail=FALSE))
## The Survival function:
Sioft = function(t,x)  exp( log(S0oft(t))*exp(sum(x*betaT)) ) ;
fioft = function(t,x) exp(sum(x*betaT))*f0oft(t)/S0oft(t)*Sioft(t,x);
Fioft = function(t,x) 1-Sioft(t,x);
## The inverse for Fioft
Finv = function(u, x) uniroot(function (t) Fioft(t,x)-u, lower=1e-100, 
                              upper=1e100, extendInt ="yes", tol=1e-6)$root

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
x1 = rbinom(ntot, 1, 0.5); x2 = rnorm(ntot, 0, 1); X = cbind(x1, x2);
## generate transformed log of survival times
z = MASS::mvrnorm(1, rep(0, ntot), corT);
## generate survival times
u = pnorm(z);
tT = rep(0, ntot);
for (i in 1:ntot){
  tT[i] = Finv(u[i], X[i,]);
}

### ----------- right-censored -------------###
t_obs=tT 
Centime = runif(ntot, 2, 6);
delta = (tT<=Centime) +0 ; 
length(which(delta==0))/ntot; # censoring rate
rcen = which(delta==0);
t_obs[rcen] = Centime[rcen]; ## observed time 
## make a data frame
dtot = data.frame(tobs=t_obs, x1=x1, x2=x2, delta=delta, tT=tT,
                  s1=s1, s2=s2); 
## Hold out npred for prediction purpose
predindex = sample(1:ntot, npred);
dpred = dtot[predindex,];
d = dtot[-predindex,];
# Prediction settings 
prediction = list(xpred=cbind(dpred$x1, dpred$x2), 
                  spred=cbind(dpred$s1, dpred$s2));

###############################################################
# Independent Cox PH
###############################################################
# MCMC parameters
nburn=500; nsave=500; nskip=0;
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(M=10, r0=1, nknots=10, nblock=n);
# here nknots=10<n, so FSA will be used with nblock=n.
# As nknots is getting larger, the FSA is more accurate but slower
# As nblock is getting smaller, the FSA is more accurate but slower. 
# In most applications, setting nblock=n works fine, which is the
# setting by not specifying nblock. 
# If nknots is not specified or nknots=n, the exact covariance is used. 
# Fit the Cox PH model
res1 = spCopulaCoxph(formula = Surv(tobs, delta)~x1+x2, data=d, 
                     prior=prior, mcmc=mcmc, prediction=prediction,
                     Coordinates=cbind(d$s1,d$s2), Knots=NULL);
# here if prediction=NULL, prediction$xpred will be set as the design matrix
# in formula, and prediction$spred will be set as the Coordinates argument. 
# Knots=NULL is the defualt setting, for which the knots will be generated 
# using fields::cover.design() with number of knots equal to prior$nknots. 
sfit1=summary(res1); sfit1;
## MSPE
mean((dpred$tT-apply(res1$Tpred, 1, median))^2); 

## traceplot
par(mfrow = c(2,2))
traceplot(mcmc(res1$beta[1,]), main="beta1")
traceplot(mcmc(res1$beta[2,]), main="beta2")
traceplot(mcmc(res1$theta[1,]), main="sill")
traceplot(mcmc(res1$theta[2,]), main="range")

############################################
## Curves
############################################
par(mfrow = c(1,1))
tgrid = seq(1e-10,4,0.1);
xpred = data.frame(x1=c(0,0), x2=c(0,1)); 
plot(res1, xnewdata=xpred, tgrid=tgrid);
}

