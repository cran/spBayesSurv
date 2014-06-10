\name{anovaDDP}
\alias{anovaDDP}
\alias{anovaDDP.default}

\title{Bayesian Nonparametric Nonproportional Hazards Survival Modeling}
\description{
    This function fits an ANOVA DDP model for survival analysis of right censored time-to-event data.
}
   
\usage{
anovaDDP(y, delta, x=NULL, prediction, prior, mcmc, state, 
         data=sys.frame(sys.parent()), na.action=na.fail, work.dir=NULL)
}

\arguments{
   \item{y}{         an n by 1 vector giving the log survival times.}
   
   \item{delta}{	   an n by 1 vector indicating whether it is right censored (=0) or not (=1).}
   
   \item{x}{         an n by p matrix of covariates without intercept. The default is NULL, 
                     indicating no covariates included.} 

   \item{prediction}{a list giving the information used to obtain conditional 
                     inferences. The list includes the following elements: 
                     \code{xpred} giving the n by p covariates matrix, used for prediction.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{N} an integer giving the truncation of the Dirichlet 
                     process, \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter \code{alpha}. See Zhou, Hanson and Knapp (2014+) 
                     for more detailed hyperprior specifications.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following elements: \code{nburn} an integer giving the number of burn-in 
                     scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer 
                     giving the total number of scans to be saved, \code{ndisplay} an integer giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}     

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{anovaDDP} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This function fits an ANOVA DDP model (De Iorio et al., 2009) for survival analysis of right 
  censored time-to-event data.
}

\value{
  The results include the MCMC chains for the parameters discussed in De Iorio et al. (2009) and 
  Zhou, Hanson and Knapp (2014+). Use \code{names} to find out what they are.

}

\seealso{
\code{\link{spCopulaDDP}}
}

\references{
De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. 
Biometrics, 65(3), 762-771.

Zhou, H., Hanson, T. and Knapp, R. (2014+). A Spatial Copula Approach to Fully Bayesian Nonparametric Survival Analysis. 
In preparation for publication.
}

\author{
Haiming Zhou \email{<zhouh@email.sc.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{ANOVA DDP}
\keyword{Bayesian nonparametric}

\examples{
\dontrun{

###############################################################
# A simulated data: mixture of two normals
###############################################################
rm(list=ls())
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(survival)
library(spBayesSurv)

## True parameters 
betaT = cbind(c(3.5, 0.5), c(2.5, -1)); 
wT = c(0.4, 0.6); 
sig2T = c(1^2, 0.5^2);
theta1 = 0.98; theta2 = 100000;

## generate coordinates: 
## npred is the # of locations for prediction
n = 300; npred = 30; ntot = n + npred;
ldist = 100; wdist = 40;
s1 = runif(ntot, 0, wdist); s2 = runif(ntot, 0, ldist);
s = rbind(s1,s2);
#plot(s[1,], s[2,]);

## Covariance matrix
corT = matrix(1, ntot, ntot);
for (i in 1:(ntot-1)){
  for (j in (i+1):ntot){
    dij = sqrt(sum( (s[,i]-s[,j])^2 ));
    corT[i,j] = theta1*exp(-theta2*dij);
    corT[j,i] = theta1*exp(-theta2*dij);
  }
}

## Generate x 
x = runif(ntot,-1.5,1.5);
X = cbind(rep(1,ntot), x);
p = ncol(X); # number of covariates + 1
## Generate transformed log of survival times
z = mvrnorm(1, rep(0, ntot), corT);
## The pdf of Ti:
fi = function(y, xi, w=wT){
  nw = length(w);
  ny = length(y);
  res = matrix(0, ny, nw);
  Xi = c(1,xi);
  for (k in 1:nw){
    res[,k] = w[k]*dnorm(y, sum(Xi*betaT[,k]), sqrt(sig2T[k]) )
  }
  apply(res, 1, sum)
}
## true plot
xx = seq(-2, 7, 0.01)
plot(xx, fi(xx, -1), "l", lwd=2, col=2)
lines(xx, fi(xx, 1), "l", lwd=2, col=3)
## The CDF of Ti:
Fi = function(y, xi, w=wT){
  nw = length(w);
  ny = length(y);
  res = matrix(0, ny, nw);
  Xi = c(1,xi);
  for (k in 1:nw){
    res[,k] = w[k]*pnorm(y, sum(Xi*betaT[,k]), sqrt(sig2T[k]) )
  }
  apply(res, 1, sum)
}
## The inverse for CDF of Ti
Finvsingle = function(u, xi) {
  res = uniroot(function (x) Fi(x, xi)-u, lower=-500, upper=500);
  res$root
}
Finv = function(u, xi) {sapply(u, Finvsingle, xi)};
## Generate log of survival times y
u = pnorm(z);
y = rep(0, ntot);
for (i in 1:ntot){
  y[i] = Finv(u[i], x[i]);
}
#plot(x,y);
yTrue = y;

## Censoring scheme
Centime = runif(ntot, 3.5,5);
Centime = 10000;
delta = (y<=Centime) +0 ;
sum(delta)/ntot;
cen = which(delta==0);
y[cen] = Centime[cen];

## make a data frame
dtotal = data.frame(s1=s1, s2=s2, y=y, x=x, delta=delta, yTrue=yTrue);
## Hold out npred=30 for prediction purpose
predindex = sample(1:ntot, npred);
dpred = dtotal[predindex,];
dtrain = dtotal[-predindex,];

# rename the variables 
d = dtrain;n=nrow(d); n;
s = cbind(d$s1, d$s2);
y = d$y;
x = d$x;
delta =d$delta;
 
# Prediction settings 
xpred = dpred$x;
s0 = cbind( dpred$s1, dpred$s2 );
prediction = list(spred=s0, xpred=xpred);

###############################################################
# ANOVA DDP 
###############################################################

# Prior information
prior = list(N = 10,
             m0 = rep(0,p), 
             S0 = diag(rep(1e5,p), nrow=p, ncol=p), 
             Sig0 = diag(rep(1e5,p), nrow=p, ncol=p), 
             k0 = 7,
             nua = 2, nub = 1, 
             a0 = 1, b0 = 1);

# current state values
state <- NULL

# MCMC parameters
nburn <- 1000
nsave <- 500
nskip <- 4
ndisplay <- 100
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)

# Fit the model
res1 = anovaDDP( y = y,
              delta =delta, 
              x = x, 
              prediction=prediction, 
              prior=prior, 
              mcmc=mcmc,
              state=state);

par(mfrow = c(2,2))
w.save = res1$w;
Kindex = which.max(rowMeans(w.save));
traceplot(mcmc(w.save[Kindex,]), main="w")
sig2.save = res1$sigma2;
traceplot(mcmc(sig2.save[Kindex,]), main="sig2")
beta.save = res1$beta;
traceplot(mcmc(beta.save[2,Kindex,]), main="beta")
alpha.save = res1$alpha;
traceplot(mcmc(alpha.save), main="alpha")

## LPML
#cpo=CPOanovaDDP(y, delta, X, beta.save, sig2.save, w.save )$cpo;
LPML1 = sum(log(res1$cpo)); LPML1;

## MSPE
mean((dpred$yTrue-apply(res1$Ypred, 1, median))^2); 

## plots
par(mfrow = c(2,2))
xnew = c(-1, 1)
xpred = cbind(xnew); 
nxpred = nrow(xpred);
ygrid = seq(0,6.0,0.05); tgrid = exp(ygrid);
ngrid = length(ygrid);
estimates = GetCurves(res1, xpred, ygrid, CI=c(0.05, 0.95));
fhat = estimates$fhat; 
Shat = estimates$Shat;
## density in y
plot(ygrid, fi(ygrid, xnew[1]), "l", lwd=2, ylim=c(0, 0.8), 
     xlim=c(0,6), main="density in y")
for(i in 1:nxpred){
  lines(ygrid, fi(ygrid, xnew[i]), lwd=2)
  lines(ygrid, fhat[,i], lty=2, lwd=2, col=4);
}
## survival in y
plot(ygrid, 1-Fi(ygrid, xnew[1]), "l", lwd=2, ylim=c(0, 1), 
     xlim=c(0,6), main="survival in y")
for(i in 1:nxpred){
  lines(ygrid, 1-Fi(ygrid, xnew[i]), lwd=2)
  lines(ygrid, Shat[,i], lty=2, lwd=2, col=4);
}
## density in t
plot(tgrid, fi(ygrid, xnew[1])/tgrid, "l", lwd=2, ylim=c(0, 0.15), 
     xlim=c(0,100), main="density in t")
for(i in 1:nxpred){
  lines(tgrid, fi(ygrid, xnew[i])/tgrid, lwd=2)
  lines(tgrid, fhat[,i]/tgrid, lty=2, lwd=2, col=4);
}
## survival in t
plot(tgrid, 1-Fi(ygrid, xnew[1]), "l", lwd=2, ylim=c(0, 1), 
     xlim=c(0,100), main="survival in t")
for(i in 1:nxpred){
  lines(tgrid, 1-Fi(ygrid, xnew[i]), lwd=2)
  lines(tgrid, Shat[,i], lty=2, lwd=2, col=4);
}

}
}
