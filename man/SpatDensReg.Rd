\name{SpatDensReg}
\alias{SpatDensReg}
\alias{BF.SpatDensReg}
\alias{print.SpatDensReg}
\alias{summary.SpatDensReg}
\alias{print.summary.SpatDensReg}

\title{Bayesian Nonparametric Spatially Smoothed Density Estimation}
\description{
This function provides a Bayesian nonparametric density estimator that changes smoothly in space. The estimator is built from the predictive rule for a marginalized Polya tree (PT), modified so that observations are spatially weighted by their distance from the location of interest.
}
   
\usage{
SpatDensReg(formula, data, na.action, prior=NULL, state=NULL,
            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500),
            permutation=TRUE, fix.theta=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them. Note: for survival data, the input response should be log survival times.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{prior}{a list giving the prior information. The list includes the following parameter: \code{maxL} an integer giving the maximum level of Polya trees, (\code{a0}, \code{b0}) parameters of the gamma prior for the precision parameter alpha, (\code{theta0}, \code{V0}) mean and variance of the normal prior for the centering distribution parameter vector theta, \code{phiq0} the prior probability of phi=0, \code{phib0} the rate parameter of the exponential prior of phi. The function itself provides all default priors. See Hanson, Zhou, and de Carvalho (2018).}
   
   \item{state}{a list giving the current value of the parameters. If \code{NULL}, all values are provided based on the centering parametric model.} 
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{permutation}{flag to indicate whether a random data permutation will be implemented in the beginning of each iterate; default is \code{TRUE}.}
   
   \item{fix.theta}{flag to indicate whether the centering distribution parameters theta=(location, log(scale)) are fixed; default is \code{TRUE} indicating fixed.}
}

\value{
The \code{SpatDensReg} object is a list containing at least the following components:
\item{modelname}{
the name of the fitted model
}
\item{terms}{
the \code{\link{terms}} object used
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
the \code{\link{Surv}} object used 
}
\item{X}{
the n by p design matrix
}
\item{theta}{
the 2 by nsave matrix of posterior samples for location and log(scale) involved in the centering distribution of Polya tree prior. 
}
\item{phi}{
the vector of posterior samples for the phi parameter in the marginal PT definition. 
}
\item{alpha}{
the vector of posterior samples for the precision parameter alpha in the PT prior. 
}
\item{maxL}{
the truncation level used in the PT prior. 
}

\item{Surv.cox.snell}{
the \code{\link{Surv}} object used for Cox-Snell residual plot. This is not recommended for frailty models, for which please use the function \code{\link{cox.snell.survregbayes}}. 
}
\item{ratetheta}{
the acceptance rate in the posterior sampling of theta vector involved in the centering distribution
}
\item{ratec}{
the acceptance rate in the posterior sampling of precision parameter alpha involved in the PT prior
}
\item{ratephi}{
the acceptance rate in the posterior sampling of phi parameter 
}
\item{initial.values}{
the list of initial values used for the parameters 
}
\item{BF}{
the Bayes factor for comparing the spatial model vs. the exchangeable model
}
}


\references{
Hanson, T., Zhou, H., and de Carvalho, V. (2018). Bayesian nonparametric spatially smoothed density estimation. In \emph{New Frontiers of Biostatistics and Bioinformatics} (pp 87-105). Springer.
}


\author{
Haiming Zhou and Timothy Hanson
}

\examples{
## Simulated data
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)

## True conditional density
fiofy_x = function(y, x){
  0.5*dnorm(y, -x, 1)+0.5*dnorm(y, x, 1);
}

## Generate data
n = 200;
x = runif(n, 0, 3)
y = rep(0, n);
uu = runif(n);
for(i in 1:n){
  if(uu[i]<0.5){
    y[i] = rnorm(1, -x[i], 1);
  }else{
    y[i] = rnorm(1, x[i], 1);
  }
}

## right censored 
y1=y;y2=y;
Centime = runif(n, 2, 4); 
delta = (y<=Centime) +0 ; 
length(which(delta==0))/n; ## censoring rate
rcen = which(delta==0);
y1[rcen] = Centime[rcen];
y2[rcen] = NA;
## make a data frame
## Method 1: in the interval-censoring notation: 
## y1 is the left endpoint and y2 is the right endpoint.
## This way we could use Surv(y1, y2, type="interval2")
## Method 2: Because we have right-censored data, 
## we could use y1 as the observed survival times and delta as the indicator. 
## This way we could use Surv(y1, delta). This is the same as above. 
d = data.frame(y1=y1, y2=y2, x=x, delta=delta); 

##-------------fit the model-------------------##
# MCMC parameters
nburn=50; nsave=50; nskip=0;
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=50);
prior = list(maxL=4, phiq0=0);
# Note please set 0<phiq0<1 for a valid Bayes factor of testing 
# spatial model vs. exchangeable model.
# If the Bayes factor is not needed, setting phiq0=0 will speed up 
# the computing time about seven times. 
state = list(alpha=1);
ptm<-proc.time()
res1 = SpatDensReg(formula = Surv(y1, delta)~x, data=d, 
                   prior=prior, state=state, mcmc=mcmc, permutation = TRUE, 
                   fix.theta=FALSE);
## Or equivalently formula = Surv(y1, y2, type="interval2") can also be used.
sfit=summary(res1); sfit
systime1=proc.time()-ptm; systime1;
traceplot(mcmc(res1$theta[1,]))
traceplot(mcmc(res1$theta[2,]))
traceplot(mcmc(res1$alpha))
traceplot(mcmc(res1$phi))

## plots
ygrid = seq(-6, 6,length.out=100);
xpred = data.frame(x=c(0,1,2,3)); 
plot(res1, xnewdata=xpred, ygrid=ygrid);
}

