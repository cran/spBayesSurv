\name{survregbayes2}
\alias{survregbayes2}
\alias{print.survregbayes2}
\alias{plot.survregbayes2}
\alias{summary.survregbayes2}
\alias{print.summary.survregbayes2}

\title{Bayesian Semiparametric Survival Models}
\description{
This function fits mixtures of Polya trees (MPT) proportional hazards, proportional odds, and accelerated failture time models. It also allows to include exchangeable and CAR frailties for fitting clusterd survival data. The function can fit both Case I and Case II interval censored data, as well as standard right-censored, uncensored, and mixtures of these.
}
   
\usage{
survregbayes2(formula, data, na.action, survmodel="PH", dist="loglogistic", 
              mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
              prior=NULL, state=NULL, selection=FALSE, Proximity=NULL,
              truncation_time=NULL, subject.num=NULL, InitParamMCMC=TRUE, 
              scale.designX=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them. To include CAR frailties, add \code{frailtyprior("car",ID)} to the formula, where \code{ID} is an n dimensional vector of cluster ID numbers. Furthermore, use \code{frailtyprior("iid",ID)} for Gaussian IID frailties, and exclude the term \code{frailtyprior()} for non-frailty models. Note: the data need to be sorted by ID.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{survmodel}{a character string for the assumed survival model. The options include \code{"PH"} for proportional hazards, \code{"PO"} for proportional odds, and \code{"AFT"} for accelerated failture time.}
   
   \item{dist}{centering distribution for MPT. Choices include \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The list includes the following parameter: \code{maxL} an integer giving the maximum number of mixtures of beta distributions. The function itself provides all default priors. }
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{selection}{flag to indicate whether variable selection is performed, where \code{FALSE} indicates that no variable selection will be performed.}
   
   \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters/regions. If CAR frailty model is specified in the formula, \code{Proximity} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID.}

   \item{truncation_time}{a vector of left-trucation times with length n.}
   
   \item{subject.num}{a vector of suject id numbers when time dependent covariates are considered. For example, for subject 1 time dependent covariates are recorded over [0,1), [1,3), and for subject 2 covariates are recorded over [0,2), [2,3), [3,4). Suppose we only have two subjects, i.e. n=2. In this case, we save the data in the long format, set truncation_time=c(0,1,0,2,3) and subject.num=c(1,1,2,2,2).}
   
   \item{InitParamMCMC}{flag to indicate wheter an initial MCMC will be run based on the centering parametric model, where \code{TRUE} indicates yes. }
   
   \item{scale.designX}{flag to indicate wheter the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
}


\value{
  The results include the MCMC chains for the parameters; use \code{names} to find out what they are.

}

\seealso{
\code{\link{frailtyprior}}
}

\references{
Zhou, H. and Hanson, T. (2015). Bayesian spatial survival models. In Nonparametric Bayesian Inference in Biostatistics (pp 215-246). Springer International Publishing.

Zhao, L. and Hanson, T. (2011). Spatially dependent Polya tree modeling for survival data. \emph{Biometrics}, 67(2), 391-403.
}


\author{
Haiming Zhou and Timothy Hanson
}

\examples{
rm(list=ls())
library(coda)
library(survival)
library(spBayesSurv)

## True coeffs
betaT = c(1,1); 
## Baseline Survival
f0oft = function(t) 0.5*dlnorm(t, -1, 0.5)+0.5*dlnorm(t,1,0.5);
S0oft = function(t) (0.5*plnorm(t, -1, 0.5, lower.tail=FALSE)+
                       0.5*plnorm(t, 1, 0.5, lower.tail=FALSE));
## The Survival function:
Sioft = function(t,x,v=0)  exp( log(S0oft(t))*exp(sum(x*betaT)+v) ) ;
Fioft = function(t,x,v=0) 1-Sioft(t,x,v);
## The inverse for Fioft
Finv = function(u, x,v=0) uniroot(function (t) Fioft(t,x,v)-u, 
                                  lower=1e-100, upper=1e100,
                                  extendInt ="yes")$root

##-------------Generate data-------------------##
## generate x 
n = 100; 
x1 = rbinom(n, 1, 0.5); x2 = rnorm(n, 0, 1); X = cbind(x1, x2);
## generate survival times
u = runif(n);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = Finv(u[i], X[i,]);
}

### ----------- partly interval-censored -------------###
t1=rep(NA, n);t2=rep(NA, n); delta=rep(NA, n); 
n1 = floor(0.5*n); ## right-censored part
n2 = n-n1; ## interval-censored part
# right-censored part
rcen = sample(1:n, n1);
t1_r=tT[rcen];t2_r=tT[rcen];
Centime = runif(n1, 2, 6);
delta_r = (tT[rcen]<=Centime) +0 ; length(which(delta_r==0))/n1;
t1_r[which(delta_r==0)] = Centime[which(delta_r==0)];
t2_r[which(delta_r==0)] = NA;
t1[rcen]=t1_r; t2[rcen]=t2_r; delta[rcen] = delta_r;
# interval-censored part
intcen = (1:n)[-rcen];
t1_int=rep(NA, n2);t2_int=rep(NA, n2); delta_int=rep(NA, n2);
npois = rpois(n2, 2)+1;
for(i in 1:n2){
  gaptime = cumsum(rexp(npois[i], 1)); 
  pp = Fioft(gaptime, X[intcen[i],]);
  ind = sum(u[intcen[i]]>pp); 
  if (ind==0){
    delta_int[i] = 2;
    t2_int[i] = gaptime[1];
  }else if (ind==npois[i]){
    delta_int[i] = 0;
    t1_int[i] = gaptime[ind];
  }else{
    delta_int[i] = 3;
    t1_int[i] = gaptime[ind];
    t2_int[i] = gaptime[ind+1];
  }
}
t1[intcen]=t1_int; t2[intcen]=t2_int; delta[intcen] = delta_int;
## make a data frame
d = data.frame(t1=t1, t2=t2, x1=x1, x2=x2, delta=delta, tT=tT);
table(d$delta)/n;

##-------------spBayesSurv-------------------##
fit0=survreg(formula = Surv(t1, t2, type="interval2")~x1+x2, 
             data=d, dist="loglogistic");
# MCMC parameters
nburn=500; nsave=500; nskip=0; niter = nburn+nsave
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=500);
prior = list(maxL=4, a0=1, b0=1);
ptm<-proc.time()
res = survregbayes2(formula = Surv(t1, t2, type="interval2")~x1+x2, data=d, 
                     survmodel="PH", prior=prior, mcmc=mcmc, 
                    dist="loglogistic", InitParamMCMC=FALSE);
## Note InitParamMCMC=FALSE is used only speeding,
## InitParamMCMC=TRUE is recommended in general. 
sfit=summary(res); sfit;
systime=proc.time()-ptm; 

### trace plots
par(mfrow = c(2,2));
traceplot(mcmc(res$beta[1,]), main="beta1");
traceplot(mcmc(res$beta[2,]), main="beta2");

####################################################################
## Get curves
####################################################################
par(mfrow = c(1,1));
wide=0.01;
tgrid = seq(0.001,4,wide);
ngrid = length(tgrid);
xnew = c(0,1)
xpred = cbind(c(0,0), xnew); 
nxpred = nrow(xpred);
estimates=plot(res, xpred, tgrid);
## plots
## survival function when x=(0,0)
i=2
par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
plot(tgrid, Sioft(tgrid, c(0,xnew[i])), "l", lwd=3, 
     xlim=c(0,3), xlab="time", ylab="survival");
polygon(x=c(rev(tgrid),tgrid),
        y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
        border=NA,col="lightgray");
lines(tgrid, Sioft(tgrid, c(0,xnew[i])), "l", lwd=3);
lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
## survival function when x=(0,0)
i=1
par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
lines(tgrid, Sioft(tgrid, c(0,xnew[i])), "l", lwd=3, 
     xlim=c(0,3), xlab="time", ylab="survival");
polygon(x=c(rev(tgrid),tgrid),
        y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
        border=NA,col="lightgray");
lines(tgrid, Sioft(tgrid, c(0,xnew[i])), "l", lwd=3);
lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
}

