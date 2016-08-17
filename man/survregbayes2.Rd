\name{survregbayes2}
\alias{survregbayes2}
\alias{print.survregbayes2}
\alias{plot.survregbayes2}
\alias{summary.survregbayes2}
\alias{print.summary.survregbayes2}



\title{Bayesian Semiparametric Survival Models}
\description{
This function fits semiparametric proportional hazards, proportional odds, and accelerated failture time models. Both georeferenced (location observed exactly) and areally observed (location known up to a geographic unit such as a county) spatial locations can be handled. Georeferenced data are modeled a discrete process convolution whereas areal data are modeled with a Markov random field. Variable selection is also incorporated. The function can fit both Case I and Case II interval censored data, as well as standard right-censored, uncensored, and mixtures of these. The Bernstein Polynomial (BP) prior is used for fitting the baseline survival function. 
}
   
\usage{
survregbayes2(formula, data, na.action, survmodel="PH", dist="weibull", 
              mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
              prior=NULL, state=NULL, selection=FALSE, Proximity=NULL, 
              truncation_time=NULL, Knots=NULL, InitParamMCMC=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them. To include CAR frailties, add \code{frailtyprior("car",ID)} to the formula, where \code{ID} is an n dimensional vector of cluster ID numbers. Furthermore, use \code{frailtyprior("iid",ID)} for Gaussian exchangeable frailties, use \code{frailtyprior("kriging",s1,s2)} for process convolution frailties, and exclude the term \code{frailtyprior()} for non-frailty models. Note: the data need to be sorted by ID.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{survmodel}{a character string for the assumed survival model. The options include \code{"PH"} for proportional hazards, \code{"PO"} for proportional odds, and \code{"AFT"} for accelerated failture time.}
   
   \item{dist}{centering distribution for BP. Choices include \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The list includes the following parameter: \code{maxL} an integer giving the maximum number of mixtures of beta distributions. The function itself provides all default priors. }
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{selection}{flag to indicate whether variable selection is performed, where \code{FALSE} indicates that no variable selection will be performed.}
   
   \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters/regions. If CAR frailty model is specified in the formula, \code{Proximity} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID.}

   \item{truncation_time}{a vector of left-trucation times with length n.}
   
   \item{Knots}{an m by d matrix, where m is the number of selected knots for process convolution, and d is the dimension of each location. If \code{Knots} is not specified, the space-filling algorithm will be used to find the knots.}
   
   \item{InitParamMCMC}{flag to indicate wheter an initial MCMC will be run based on the centering parametric model, where \code{TRUE} indicates yes. }
}


\value{
  The results include the MCMC chains for the parameters; use \code{names} to find out what they are.

}

\seealso{
\code{\link{frailtyprior}}
}

\references{
Zhou, H. and Hanson, T. (2016). Bayesian semiparametric models for spatially correlated arbitrarily censored data. In preparation.

Higdon, D. (2002). Space and space-time modeling using process convolutions. In \emph{Quantitative methods for current environmental issues} (pp. 37-56). Springer London.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
\dontrun{
rm(list=ls())
library(coda)
library(survival)
library(spBayesSurv)
library(fields)

## True coeffs
betaT = c(1,1); 
## Baseline Survival
f0oft = function(t) 0.5*dlnorm(t, -1, 0.5)+0.5*dlnorm(t,1,0.5);
S0oft = function(t) (0.5*plnorm(t, -1, 0.5, lower.tail=FALSE)+
                       0.5*plnorm(t, 1, 0.5, lower.tail=FALSE))
## The Survival function:
Sioft = function(t,x,v=0)  exp( log(S0oft(t))*exp(sum(x*betaT)+v) ) ;
Fioft = function(t,x,v=0) 1-Sioft(t,x,v);
## The inverse for Fioft
Finv = function(u, x,v=0) uniroot(function (t) Fioft(t,x,v)-u, lower=1e-100, 
                                  upper=1e100, extendInt ="yes", tol=1e-6)$root
## kernel function
kern = function(dis, phi) exp(-0.5*(dis/phi)^2);

###############################################################################
########################### Start to simulation ###############################
###############################################################################
n = 500; ## sample size
nknots = 50; phiT=5; ## phiT is the kernel range parameter phi. 
tau2T = 1; ## true latent process variance; 
s1 = runif(n, 0, 40);
s2 = runif(n, 0, 100);
ss = cbind(s1, s2); ### the locations. 
Dnn = .Call("DistMat", t(ss), t(ss), PACKAGE = "spBayesSurv");
s0 = as.matrix(fields::cover.design(ss, nd=nknots)$design); ## knots selection
Dnm = .Call("DistMat", t(ss), t(s0), PACKAGE = "spBayesSurv");
ZZ = kern(Dnm, phi=phiT);
v = rnorm(nknots, 0, sqrt(tau2T)); ## generate laten process at each knot. 
vn = as.vector(crossprod(t(ZZ), v)); ## generate frailties at each location
## generate x 
x1 = rbinom(n, 1, 0.5); x2 = rnorm(n, 0, 1); X = cbind(x1, x2);
## generate survival times
u = runif(n);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = Finv(u[i], X[i,], vn[i]);
}

### ----------- right censored -------------###
t1=tT;t2=tT;
## right censored
Centime = runif(n, 2,6);
delta = (tT<=Centime) +0 ; length(which(delta==0))/n;
rcen = which(delta==0);
t1[rcen] = Centime[rcen];
t2[rcen] = NA;
## make a data frame
## Method 1: in the interval-censoring notation: 
## t1 is the left endpoint and t2 is the right endpoint.
## This way we could use Surv(t1, t2, type="interval2")
## Method 2: Because we have right-censored data, 
## we could use t1 as the observed survival times and delta as the indicator. 
## This way we could use Surv(t1, delta). This is the same as above. 
## (s1, s2) are the locations. 
d = data.frame(t1=t1, t2=t2, x1=x1, x2=x2, delta=delta, s1=s1, s2=s2); 
table(d$delta)/n;

##-------------spBayesSurv-------------------##
# MCMC parameters
nburn=3000; nsave=3000; nskip=0; niter = nburn+nsave
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(maxL=15, a0=1, b0=1); 
state <- list(cpar=1);
ptm<-proc.time()
res1 = survregbayes2(formula = Surv(t1, delta)~x1+x2+
                       frailtyprior("kriging", s1, s2), data=d, 
                     survmodel="PH", prior=prior, mcmc=mcmc, 
                     state=state, dist="loglogistic", Knots=s0);
## Or equivalently
#res1 = survregbayes2(formula = Surv(t1, t2, type="interval2")~x1+x2+
#                       frailtyprior("kriging", s1, s2), data=d, 
#                     survmodel="PH", prior=prior, mcmc=mcmc, 
#                     state=state, dist="loglogistic", Knots=s0);
sfit=summary(res1); sfit
systime1=proc.time()-ptm; systime1;

############################################
## Results
############################################
## acceptance rate of frailties
res1$ratev[1]
## traceplots;
par(mfrow=c(2,3));
traceplot(mcmc(res1$beta[1,]), main="beta1");
traceplot(mcmc(res1$beta[2,]), main="beta2");
traceplot(mcmc(res1$v[1,]), main="frailty");
traceplot(mcmc(res1$v[2,]), main="frailty");
traceplot(mcmc(res1$v[3,]), main="frailty");
#traceplot(mcmc(res1$v[4,]), main="frailty");
traceplot(mcmc(res1$phi), main="phi");

############################################
## Curves
############################################
wide=0.02;
tgrid = seq(1e-10,4,wide);
ngrid = length(tgrid);
p = length(betaT); # number of covariates
xpred = rbind(c(0,0), c(0,1)); 
estimates=plot(res1, xpred=xpred, tgrid=tgrid);
Shat = estimates$Shat;

## plot
plot(tgrid, Sioft(tgrid, xpred[2,]), "l", lwd=3);
lines(tgrid, Sioft(tgrid, xpred[1,]), "l", lwd=3);
lines(estimates$tgrid, estimates$Shat[,1], lty=2, lwd=3)
lines(estimates$tgrid, estimates$Shat[,2], lty=2, lwd=3)
}
}

