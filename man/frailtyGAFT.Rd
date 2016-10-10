\name{frailtyGAFT}
\alias{frailtyGAFT}
\alias{print.frailtyGAFT}
\alias{plot.frailtyGAFT}
\alias{summary.frailtyGAFT}
\alias{print.summary.frailtyGAFT}

\title{Generalized Accelerated Failure Time Frailty Model}
\description{
This function fits a generalized accelerated failure time frailty model for clustered and/or areal-level time-to-event data.
}
   
\usage{
frailtyGAFT(formula, data, na.action, 
            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
            prior=NULL, state=NULL, Proximity=NULL, Coordinates=NULL,
            DIST=NULL)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them. To include CAR frailties, add \code{frailtyprior("car",ID)} to the formula, where \code{ID} is an n dimensional vector of cluster ID numbers. Furthermore, use \code{frailtyprior("iid",ID)} for Gaussian exchangeable frailties, use \code{frailtyprior("grf",ID)} for Gaussian random fields (GRF) frailties, and exclude the term \code{frailtyprior()} for non-frailty models. Note: the data need to be sorted by ID. }
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The list includes the following parameter: \code{maxL} an integer giving the maximum level of LDTFP. See Zhou, et al. (2016) for more detailed prior specifications.}
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters/regions. If CAR frailty model is specified in the formula, \code{Proximity} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID. }
   
   \item{Coordinates}{an m by d coordinates matrix, where m is the number of clusters/regions, d is the dimension of coordiantes. If GRF frailty model is specified in the formula, \code{Coordinates} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID.}
   \item{DIST}{This is a function argument, used to calculate the distance. The default is Euclidean distance (\code{fields::rdist}). This function should have two arguments (X1,X2), where X1 and X2 are matrices with coordinates as the rows. The returned value of this function should be the pairwise distance matrix. If nrow(X1)=m and nrow(X2)=n then the function should return an m by n matrix of all distances between these two sets of points.}
}

\details{
This function fits a a generalized accelerated failure time frailty model (Zhou, et al., 2016) for clustered and/or areal-level time-to-event data. Note that the function arguments are slightly different with those presented in the original paper of Zhou, et al. (2016); see below examples. 
}

\value{
  The results include the MCMC chains for the parameters discussed in Zhou, et al. (2016). 
  Use \code{names} to find out what they are.

}

\seealso{
\code{\link{baseline}}, \code{\link{frailtyprior}}, \code{\link{survregbayes}}, \code{\link{rdist}}
}

\references{
Zhou, H., Hanson, T., and Zhang, J. (2016). Generalized accelerated failure time spatial frailty model for arbitrarily censored data. \emph{Lifetime Data Analysis}, in press.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
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
set.seed(1234567890)
Finvsingle = function(u, F) {
  res = uniroot(function (x) F(x)-u, lower=-1000, upper=1000, 
                tol=.Machine$double.eps^0.5);
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
n = length(id); # Total number of subjects
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
  }else{
    y[i] = sum(betaT*X[i,]) + vn[i] + Finv(u[i], F0)
  } 
}

# generate responses
Cen = runif(n, 0.5, 1)
delta = (exp(y)<=Cen)+0;
sum(delta)/n
tTrue = exp(y);
tobs = cbind(tTrue, tTrue);
tobs[which(delta==0),] = cbind(Cen[which(delta==0)], NA);
dtotal = data.frame(tleft=tobs[,1], tright=tobs[,2], x1=x1, 
                    x2=x2, xtf=x2, ID=id, tTrue=tTrue, censor=delta);
## sort the data by ID
d = dtotal[order(dtotal$ID),];

# Prior information and MCMC
fit0 <- survival::survreg(Surv(tleft, censor)~x1+x2, dist="lognormal", data=d);
prior = list(maxL = 4, a0 = 5, b0 = 1);
mcmc=list(nburn=100, nsave=500, nskip=1, ndisplay=100)

# Fit the model
res = frailtyGAFT(Surv(tleft, tright, type="interval2")~x1+x2+baseline(x1, x2)+
                    frailtyprior(prior="car", ID),  data=d, mcmc=mcmc, prior=prior, 
                  Proximity=W);
summary(res);


####################################################################
## Get curves
####################################################################
par(mfrow = c(2,1))
ygrid = seq(-3,3,0.03); tgrid = exp(ygrid);
ngrid = length(ygrid);
xnew = c(0, 1)
xpred = cbind(c(1,1.5), xnew); 
xtfpred = xpred;
nxpred = nrow(xpred);
estimates = plot(res, xpred, ygrid, xtfpred=xtfpred, CI=0.9, PLOT=F);
fhat = estimates$fhat; 
Shat = estimates$Shat;
## density in y
xlim=c(-3,3)
shiftpred = rep(0, nxpred);
for(i in 1:nxpred) shiftpred[i] = (sum(c(1,xpred[i,])*betaT));

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
par(mfrow = c(3,1))
traceplot(mcmc(res$beta[1,]), main="beta1");
traceplot(mcmc(res$beta[2,]), main="beta2");
traceplot(mcmc(res$beta[3,]), main="beta3");
}
}

