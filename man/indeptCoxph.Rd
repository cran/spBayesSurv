\name{indeptCoxph}
\alias{indeptCoxph}
\alias{print.indeptCoxph}
\alias{summary.indeptCoxph}
\alias{print.summary.indeptCoxph}

\title{Bayesian Proportional Hazards Model}
\description{
This function fits a Bayesian proportional hazards model (Zhou, Hanson and Zhang, 2018) for non-spatial right censored time-to-event data.
}
   
\usage{
indeptCoxph(formula, data, na.action, prediction=NULL, 
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
    
    \item{scale.designX}{flag to indicate wheter the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
}

\details{
This function fits a Bayesian proportional hazards model (Zhou, Hanson and Zhang, 2018) for non-spatial right censored time-to-event data.
}

\value{
The \code{indeptCoxph} object is a list containing at least the following components:
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
the \code{\link{Surv}} object used 
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
\item{ratebeta}{
the acceptance rate in the posterior sampling of beta coefficient vector
}
\item{cpo}{
the length n vector of the stabilized estiamte of CPO; used for calculating LPML
}
\item{Tpred}{
the npred by nsave predicted survival times for covariates specified in the argument \code{prediction}.
}
}

\seealso{
\code{\link{spCopulaCoxph}, \link{GetCurves}}
}

\references{
Zhou, H., Hanson, T., and Zhang, J. (2018). spBayesSurv: Fitting Bayesian Spatial Survival Models Using R. \emph{Journal of Statistical Software}, accept minor.

Zhou, H., Hanson, T., and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, 71(4): 1101-1110. 
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}

}

\keyword{Cox PH}
\keyword{Bayesian}

\examples{
###############################################################
# A simulated data: Cox PH
###############################################################
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
## True parameters 
betaT = c(1,1); 
n=100; 
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

## generate x 
x1 = rbinom(n, 1, 0.5); x2 = rnorm(n, 0, 1); X = cbind(x1, x2);
## generate survival times
u = runif(n);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = Finv(u[i], X[i,]);
}

### ----------- right-censored -------------###
t_obs=tT 
Centime = runif(n, 2, 6);
delta = (tT<=Centime) +0 ; 
length(which(delta==0))/n; # censoring rate
rcen = which(delta==0);
t_obs[rcen] = Centime[rcen]; ## observed time 
## make a data frame
d = data.frame(tobs=t_obs, x1=x1, x2=x2, delta=delta, tT=tT); 
table(d$delta)/n;

###############################################################
# Independent Cox PH
###############################################################
# MCMC parameters
nburn=500; nsave=300; nskip=0;
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(M=10, r0=1);
# Fit the Cox PH model
res1 = indeptCoxph(formula = Surv(tobs, delta)~x1+x2, data=d, 
                   prior=prior, mcmc=mcmc);
sfit1=summary(res1); sfit1;
## traceplot
par(mfrow = c(2,2))
traceplot(mcmc(res1$beta[1,]), main="beta1")
traceplot(mcmc(res1$beta[2,]), main="beta2")
traceplot(mcmc(res1$h.scaled[2,]), main="h")
traceplot(mcmc(res1$h.scaled[3,]), main="h")

############################################
## Curves
############################################
par(mfrow = c(1,1))
tgrid = seq(1e-10,4,0.1);
xpred = data.frame(x1=c(0,0), x2=c(0,1)); 
plot(res1, xnewdata=xpred, tgrid=tgrid);
}

