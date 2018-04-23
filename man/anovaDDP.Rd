\name{anovaDDP}
\alias{anovaDDP}

\title{Bayesian Nonparametric Survival Model}
\description{
This function fits a Bayesian Nonparametric model (De Iorio et al., 2009) for non-spatial right censored time-to-event data. Note that the notations are different with those presented in the original paper; see Zhou, Hanson and Zhang (2018) for new examples.
}
   
\usage{
anovaDDP(formula, data, na.action, prediction=NULL,
         mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500),
         prior=NULL, state=NULL, scale.designX=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It currently only supports right-censoring.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{prediction}{a list giving the information used to obtain conditional inferences. The list includes the following element: \code{xpred} giving the npred by p covariates matrix, used for prediction. If \code{prediction=NULL}, \code{xpred} will be set to be the design matrix used in \code{formula}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
    
    \item{prior}{a list giving the prior information. See Zhou, Hanson and Zhang (2018) for more detailed hyperprior specifications.}

    \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.}  
    
    \item{scale.designX}{flag to indicate wheter the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. All returned posterior samples fit from scaled covariates. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}

}

\details{
This function fits a Bayesian Nonparametric model (De Iorio et al., 2009) for non-spatial right censored time-to-event data. Note that the notations are different with those presented in the original paper; see Zhou, Hanson and Zhang (2018) for new examples. 
}

\value{
The \code{anovaDDP} object is a list containing at least the following components:
\item{n}{
the number of row observations used in fitting the model
}
\item{p}{
the number of columns in the model matrix
}
\item{Surv}{
the \code{\link{Surv}} object used 
}
\item{X.scaled}{
the n by p scaled design matrix
}
\item{X}{
the n by p orginal design matrix
}
\item{beta}{
the p+1 by N by nsave array of posterior samples for the coefficients
}
\item{sigma2}{
the N by nsave matrix of posterior samples for sigma2 involved in the DDP.
}
\item{w}{
the N by nsave matrix of posterior samples for weights involved in the DDP. 
}
\item{Tpred}{
the npred by nsave predicted survival times for covariates specified in the argument \code{prediction}.
}
}

\seealso{
\code{\link{spCopulaDDP}, \link{GetCurves}} 
}

\references{
Zhou, H., Hanson, T., and Zhang, J. (2018). spBayesSurv: Fitting Bayesian Spatial Survival Models Using R. \emph{Journal of Statistical Software}, accept minor.

Zhou, H., Hanson, T., and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, 71(4): 1101-1110. 

De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. \emph{Biometrics}, 65(3): 762-771. 
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{ANOVA DDP}
\keyword{Bayesian nonparametric}

\examples{
###############################################################
# A simulated data: mixture of two normals
###############################################################
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
## True parameters 
betaT = cbind(c(3.5, 0.5), c(2.5, -1)); 
wT = c(0.4, 0.6); 
sig2T = c(1^2, 0.5^2);
n=100; 
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

## generate x 
x1 = runif(n,-1.5,1.5); X = cbind(x1);
## generate survival times
u = runif(n);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = exp(Finv(u[i], X[i,]));
}

### ----------- right-censored -------------###
t_obs=tT 
Centime = runif(n, 20, 200);
delta = (tT<=Centime) +0 ; 
length(which(delta==0))/n; # censoring rate
rcen = which(delta==0);
t_obs[rcen] = Centime[rcen]; ## observed time 
## make a data frame
d = data.frame(tobs=t_obs, x1=x1, delta=delta, tT=tT); 
table(d$delta)/n;

###############################################################
# Independent DDP: Bayesian Nonparametric Survival Model
###############################################################
# MCMC parameters
nburn=500; nsave=500; nskip=0;
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(N=10, a0=2, b0=2);
# Fit the Cox PH model
res1 = anovaDDP(formula = Surv(tobs, delta)~x1, data=d, 
                prior=prior, mcmc=mcmc);
## LPML
LPML = sum(log(res1$cpo)); LPML;
## Number of non-negligible components
quantile(colSums(res1$w>0.05))

############################################
## Curves
############################################
ygrid = seq(0,6.0,length=100); tgrid = exp(ygrid);
xpred = data.frame(x1=c(-1, 1))
plot(res1, xnewdata=xpred, tgrid=tgrid);

}
