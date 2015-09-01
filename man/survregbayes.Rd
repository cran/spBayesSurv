\name{survregbayes}
\alias{survregbayes}
\alias{print.survregbayes}
\alias{plot.survregbayes}
\alias{summary.survregbayes}
\alias{print.summary.survregbayes}



\title{Bayesian Semiparametric Survival Models}
\description{
This function fits mixtures of Polya trees proportional hazards, proportional odds, accelerated failture time, and accelerated hazards models. It also allows to include exchangeable and CAR frailties for fitting clusterd survival data. The function can fit both Case I and Case II interval censored data, as well as standard right-censored, uncensored, and mixtures of these.
}
   
\usage{
survregbayes(formula, data, na.action, survmodel="PH", 
             mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
             prior=NULL, state=NULL, frailty=NULL, ID=NULL, Proximity=NULL)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{survmodel}{a character string for the assumed survival model. The options include \code{"PH"} for proportional hazards, \code{"PO"} for proportional odds, \code{"AFT"} for accelerated failture time, and \code{"AH"} for accelerated hazards.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The list includes the following parameter: \code{maxL} an integer giving the maximum level of LDTFP. The function itself provides all default priors. }
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{frailty}{a charater string with \code{NULL} fitting a non-frailty model, \code{CAR} fitting a CAR frailty model and \code{IID} fitting a independent Gaussian frailty model}
   
   \item{ID}{an n dimensional vector of cluster ID numbers. If \code{frailty=NULL}, \code{ID} is ignored; otherwise, \code{ID} is required. }
   
   \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters. If \code{frailty=CAR}, \code{Proximity} is required; otherwise it is ignored.}
}


\value{
  The results include the MCMC chains for the parameters; use \code{names} to find out what they are.

}

\seealso{
\code{\link{frailtyGAFT}}
}

\references{
Zhou, H. and Hanson, T. (2015). Bayesian spatial survival models. In \emph{Nonparametric Bayesian Inference in Biostatistics} (pp. 215-246). Springer International Publishing.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
\dontrun{
library(coda)
library(survival)
library(spBayesSurv)
## True coeffs
betaT = c(-1,1); 
## True cumulative Baseline hazard 
Lambda0 = function(t) {t+log(1+t)};
## The Survival function:
Sioft = function(t,x)  exp( -Lambda0(t)*exp(sum(x*betaT)) ) ;

tt=seq(0, 6, 0.01);
plot(tt, Sioft(tt, c(0,1)), "l");
lines(tt, Sioft(tt, c(0,0)), "l");

## The inverse for Sioft
Sinv = function(u, x) uniroot(function (t) Sioft(t,x)-u, lower=0, 
                              upper=1e100, tol=.Machine$double.eps^0.5)$root

##-------------Generate data-------------------##
## Generate x 
n = 200; 
x1 = rbinom(n, 1, 0.5);
x2 = rnorm(n, 0, 0.5)
X = cbind(x1, x2);
p = ncol(X); # number of covariates
## Generate survival times
u = runif(n, 0, 1);
tT = rep(0, n);
for (i in 1:n){
  tT[i] = Sinv(u[i], X[i,]) + 1e-100;
}

### interval-censored data:
t1 = tT; t2 = tT;
obs_num = 1+rpois(n, 2);
for(i in 1:n){
  obs_time = c( NA, cumsum( rexp(obs_num[i]) ), NA);
  jj=2;
  while( (jj<(obs_num[i]+2)) & (obs_time[jj]<tT[i]) ){
    jj=jj+1;
  }
  t1[i] = obs_time[jj-1];
  t2[i] = obs_time[jj];
}
sum(is.na(t2))/n

## make a data frame
d = data.frame(t1=t1, t2=t2, x1=x1, x2=x2, tT=tT);

# MCMC parameters
nburn <- 5000
nsave <- 5000
nskip <- 0
ndisplay <- 1000
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)

# Prior information
prior = list(a0 = 5, b0 = 1);
# current state values
state <- list(cpar=5);

##-------------PH-------------------##
fit0=survreg(formula = Surv(t1, t2, type="interval2")~x1+x2, data=d, 
             dist="weibull");
summary(fit0)
ptm<-proc.time()
res = survregbayes(formula = Surv(t1, t2, type="interval2")~x1+x2, data=d, 
                   survmodel="PH", prior=prior, mcmc=mcmc, state=state);
systime=proc.time()-ptm
summary(res);

### trace plots
par(mfrow = c(3,1));
traceplot(mcmc(res$beta[1,]), main="beta1");
traceplot(mcmc(res$beta[2,]), main="beta2");
traceplot(mcmc(res$cpar));

####################################################################
## Get curves
####################################################################
par(mfrow = c(1,2));
ygrid = seq(-10,1,length=300); tgrid = exp(ygrid);
ngrid = length(ygrid);
xnew = c(0,1)
xpred = cbind(c(0,0), xnew); 
nxpred = nrow(xpred);
estimates=plot(res, xpred, tgrid);
## plots
## survival function when x=(0,0)
i=1
par(cex=1.5,mar=c(4.1,4.1,1,1),cex.lab=1.4,cex.axis=1.1)
plot(tgrid, Sioft(tgrid, c(0,xnew[i])), "l", lwd=3, 
     xlim=c(0,3), xlab="time", ylab="survival");
polygon(x=c(rev(tgrid),tgrid),
        y=c(rev(estimates$Shatlow[,i]),estimates$Shatup[,i]),
        border=NA,col="lightgray");
lines(tgrid, Sioft(tgrid, c(0,xnew[i])), "l", lwd=3);
lines(tgrid, estimates$Shat[,i], lty=3, lwd=3, col=1);
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
}
}

