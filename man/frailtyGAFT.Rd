\name{frailtyGAFT}
\alias{frailtyGAFT}
\alias{frailtyGAFT.default}
\alias{print.frailtyGAFT}
\alias{summary.frailtyGAFT}
\alias{print.summaryfrailtyGAFT}


\title{Fit Generalized Accelerated Failure Time Frailty Model}
\description{
This function fits a a generalized accelerated failure time frailty model for clustered and/or areal-level time-to-event data.
}
   
\usage{
frailtyGAFT(tobs, x=NULL, xtf=NULL, prior, mcmc, state, frailty=NULL, 
            Proximity=NULL, data=sys.frame(sys.parent()), 
            na.action=na.fail, work.dir=NULL)
}

\arguments{
   \item{tobs}{   an n by 2 matrix giving the survival times as follows. If survival time is uncensored, two values
                  in the correspong row are the same; if survival time is right censored, the second value is needed to be \code{Inf};
                  if surival time is left censored, the first value is needed to be \code{0}. }
   
   \item{x}{      an n by p-1 matrix of covariates without intercept for median regression . The default is NULL, 
					 indicating no covariates included.} 
           
  \item{xtf}{     an n by q-1 matrix of covariates without intercept for LDTFP regression. The default is NULL, 
  				 indicating no covariates included.} 
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{L} an integer giving the maximum level of LDTFP. See Zhou, Hanson and Zhang (2014+) 
                     for more detailed prior specifications.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following elements: \code{nburn} an integer giving the number of burn-in 
                     scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer 
                     giving the total number of scans to be saved, \code{ndisplay} an integer giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.} 
    \item{frailty}{  an n by 1 vector of cluster ID numbers. If it is specified, a frailty term is included; otherwise, 
                    nonfrailty GAFT model is fitted. The default is NULL.}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{frailtyGAFT} to print an error message and terminate if there are any
                     incomplete observations.}       
    \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters. If \code{Proximity=NULL}, then 
                       independent Gaussian frailty is applied.}

    \item{work.dir}{working directory.} 

}

\details{
This function fits a a generalized accelerated failure time frailty model (Zhou, Hanson and Zhang, 2014+) for clustered and/or areal-level time-to-event data.
}

\value{
  The results include the MCMC chains for the parameters discussed in Zhou, Hanson and Zhang (2014+). 
  Use \code{names} to find out what they are.

}

\seealso{
\code{\link{GetCurves}}
}

\references{
Zhou, H., Hanson, T. and Zhang, J. (2014+). Generalized accelerated failure time spatial frailty model for arbitrarily censored data. Submitted
}


\author{
Haiming Zhou \email{<zhouh@email.sc.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
\dontrun{

###############################################################
# A simulated data: GAFT spatial frailty model
###############################################################
library(survival)
library(spBayesSurv)
library(coda)
library(mvtnorm)

## True densities
Finvsingle = function(u, F) {
  res = uniroot(function (x) F(x)-u, lower=-500, upper=500);
  res$root
}
Finv = function(u, F) {sapply(u, Finvsingle, F)};
f0 = function(x) dnorm(x, 0, 0.8);
F0 = function(x) pnorm(x, 0, 0.8);
shift=1
f1 = function(x) 0.5*dnorm(x, -shift, 0.5) + 0.5*dnorm(x, shift, 0.5)
F1 = function(x) 0.5*pnorm(x, -shift, 0.5) + 0.5*pnorm(x, shift, 0.5);
ff = function(x, xtf=0) {
  if(xtf==0) {res=f0(x)} else res=f1(x)
  res
}
FF = function(x, xtf=0){
  if(xtf==0) {res=F0(x)} else res=F1(x)
  res
}

# Simulation settings;
betaT = c(-1, 1, -0.5);
tau2T = 0.1;
m = 50; # blocks
mi = 5;
mis = rep(mi, m);
id = rep(1:m,mis);
n = length(id);        # Total number of subjects
# Generate symmetric adjaceny matrix, W 
wi = rep(0, m)
while(any(wi==0)){
  W = matrix(0,m,m)
  W[upper.tri(W,diag=FALSE)]<-rbinom(m*(m-1)/2,1,.1)
  W = W+t(W) 
  wi = apply(W,1,sum)  # No. neighbors
}

# Spatial effects, v
Wstar = matrix(0, m-1, m-1);
Dstar = diag(wi[-m]);
for(i in 1:(m-1)){
  for(j in 1:(m-1)){
    Wstar[i,j] = W[j,i]-W[j,m]-W[m,i]-wi[m]
  }
}
Qstar = Dstar-Wstar;
covT = tau2T*solve(Qstar);
v0 = c(rmvnorm(1,sigma=covT));
v = c(v0,-sum(v0));
vn = rep(v, mis);

# responses
x1 = rnorm(n, 0, 1);
x2 = rbinom(n, 1, 0.5);
xtf = x2; ptf = 2;
X = cbind(1,x1,x2); pce = ncol(X);
u = runif(n, 0, 1)
y = rep(0, n);
for(i in 1:n) {
  if(x2[i]==1) {
    y[i] = sum(betaT*X[i,]) + vn[i] + Finv(u[i], F1)
  } else {
    y[i] = sum(betaT*X[i,]) + vn[i] + Finv(u[i], F0)
  } 
}

# generate responses
Cen = runif(n, 0.5, 1)
delta = (exp(y)<=Cen)+0;
sum(delta)/n
tTrue = exp(y);
tobs = cbind(tTrue, tTrue);
tobs[which(delta==0),] = cbind(Cen[which(delta==0)], Inf);
dtotal = data.frame(tleft=tobs[,1], tright=tobs[,2], 
                    x1=x1, x2=x2, xtf=x2, ID=id, tTrue=tTrue);
## sort the data by ID
d = dtotal[order(dtotal$ID),];

# Prior information
prior = list(maxL = 4,
             a0 = 5,
             b0 = 1,
             siga0 = 2.001,
             sigb0 = 2.001,
             taua0 = 0.1,
             taub0 = 0.1);

# current state values
state <- list(alpha=5);

# MCMC parameters
nburn <- 1000
nsave <- 500
nskip <- 0
ndisplay <- 100
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)

# Fit the model
res = frailtyGAFT(tobs=d[,1:2], x=d[,c("x1","x2")], xtf=d[,c("x1","x2")], 
                  prior=prior, mcmc=mcmc, state=state, frailty=NULL, Proximity=W);
sfit = summary(res);


####################################################################
## Get curves
####################################################################
par(mfrow = c(2,1));
ygrid = seq(-3,3,0.03); tgrid = exp(ygrid);
ngrid = length(ygrid);
xnew = c(0, 1)
xpred = cbind(c(1,1.5), xnew); 
xtfpred = xpred;
nxpred = nrow(xpred);
IDpred = rep(as.character(which.max(wi)),nxpred);
estimates = GetCurves(res, xpred=xpred, xtfpred=xtfpred, 
                      ygrid=ygrid, CI=c(0.05, 0.95), IDpred=IDpred);
fhat = estimates$fhat; 
Shat = estimates$Shat;
## density in y
xlim=c(-3,3)
shiftpred = rep(0, nxpred);
for(i in 1:nxpred){
  shiftpred[i] = (sum(c(1,xpred[i,])*betaT) +v[as.numeric(IDpred)[i]]);
}
## plots
## density in y
plot(ygrid, ff(ygrid-shiftpred[1],xnew[1]), "l", lwd=2,  
     xlim=xlim, main="density in y")
for(i in 1:nxpred){
  lines(ygrid, ff(ygrid-shiftpred[i],xnew[i]), "l", lwd=2,  
        xlim=xlim, main="density in y")
  lines(ygrid, fhat[,i], lty=2, lwd=2, col=4);
}
## survival in y
plot(ygrid, 1-FF(ygrid-shiftpred[1],xnew[1]), "l", lwd=2,  
     xlim=xlim, main="survival in y")
for(i in 1:nxpred){
  lines(ygrid, 1-FF(ygrid-shiftpred[i],xnew[i]), "l", lwd=2,  
        xlim=xlim, main="survival in y")
  lines(ygrid, Shat[,i], lty=2, lwd=2, col=4);
}

### trace plots
par(mfrow = c(2,2))
rowMeans(res$beta);
traceplot(mcmc(res$beta[1,]), main="beta1");
traceplot(mcmc(res$beta[2,]), main="beta2");
traceplot(mcmc(res$beta[3,]), main="beta3");
mean(res$sigma2);
traceplot(mcmc(res$sigma2), main="sigma2");
}
}

