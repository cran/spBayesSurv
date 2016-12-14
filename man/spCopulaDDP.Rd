\name{spCopulaDDP}
\alias{spCopulaDDP}
\alias{spCopulaDDP.default}

\title{Fit Marginal Bayesian Nonparametric Survival Model}
\description{
This function fits a marginal Bayesian Nonparametric model (Zhou, Hanson and Knapp, 2015) for point-referenced right censored time-to-event data.
}
   
\usage{
spCopulaDDP(y, delta, x=NULL, s, prediction, prior, mcmc, state, status=TRUE,
            FSA = TRUE, knots, data=sys.frame(sys.parent()), 
            na.action=na.fail, work.dir=NULL)
}

\arguments{
   \item{y}{         an n by 1 vector giving the log survival times.}
   
   \item{delta}{	   an n by 1 vector indicating whether it is right censored (=0) or not (=1).}
   
   \item{x}{         an n by p matrix of covariates without intercept. The default is NULL, 
					           indicating no covariates included.} 
					 
   \item{s}{         an n by d matrix of UMT coordinates, where d is the dimension of space.}

   \item{prediction}{a list giving the information used to obtain conditional
                     inferences. The list includes the following
                     elements: \code{spred} and \code{xpred} giving the n by 2 new locations 
					           and corresponding n by p covariates matrix, respectively, used for prediction. 
                     In addition, \code{predid} needs to be specified if \code{FSA=TRUE} .}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{N} an integer giving the truncation of the Dirichlet 
                     process, \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter \code{alpha}. See Zhou, Hanson and Knapp (2015) 
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
                     
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the continuation of 
                     a previous analysis (\code{FALSE}). In the latter case the current value of the parameters 
                     must be specified in the object state.}
    
    \item{FSA}{      indicate if the full scale approximation is need. The default is FALSE.} 
    
    \item{knots}{    a list giving the knots and block ids when \code{FSA=TRUE}. This list includes the
                     following parameter: \code{ss} an m by d matrix of UMT coordinates, where d is the
                     dimension of space, \code{blockid} an n by 1 id vector indicating which block that 
                     each observation is in.} 

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{spCopulaDDP} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
This function fits a marginal Bayesian Nonparametric model (Zhou, Hanson and Knapp, 2015) for point-referenced right censored time-to-event data. 
}

\value{
  The results include the MCMC chains for the parameters discussed in Zhou, Hanson and Knapp (2015). 
  Use \code{names} to find out what they are.

}

\seealso{
\code{\link{anovaDDP}}
}

\references{
Zhou, H., Hanson, T., and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, 71(4): 1101-1110. 
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{Spatial copula}
\keyword{Bayesian nonparametric}


\examples{
###############################################################
# A simulated data: mixture of two normals
###############################################################
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
library(MASS)

## True parameters 
set.seed(100)
betaT = cbind(c(3.5, 0.5), c(2.5, -1)); 
wT = c(0.4, 0.6); 
sig2T = c(1^2, 0.5^2);
theta1 = 0.98; theta2 = 0.1;

## generate coordinates: 
## npred is the # of locations for prediction
n = 50; npred = 3; ntot = n + npred;
ldist = 100; wdist = 40;
s1 = runif(ntot, 0, wdist); s2 = runif(ntot, 0, ldist);
s = rbind(s1,s2);
#plot(s[1,], s[2,]);
## divide them into blocks
nldist=5; nwdist=2;
nb=nldist*nwdist; nb; # number of blocks;
coor = matrix(0, nb, 4); ## four edges for each block;
tempindex=1; lstep=ldist/nldist; wstep=wdist/nwdist;
for(i in 1:nwdist){
  for(j in 1:nldist){
    coor[tempindex,] = c((i-1)*wstep, i*wstep, (j-1)*lstep, j*lstep );  
    tempindex = tempindex + 1;
  }
}
## Assign block id for each location
blockid = rep(NA,ntot);
for(i in 1:nb){
  blockid[((s1>coor[i,1])*(s1<=coor[i,2])*(s2>coor[i,3])*(s2<=coor[i,4]))==1]=i;
}
## Choose knots S*
nldist=10; nwdist=4;
m=nldist*nwdist; m; # number of knots;
ss = matrix(0, m, 2); 
tempindex=1; lstep=ldist/nldist; wstep=wdist/nwdist;
for(i in 1:nwdist){
  for(j in 1:nldist){
    ss[tempindex,] = c( (i-1)*wstep+wstep/2, (j-1)*lstep+lstep/2);  
    tempindex = tempindex + 1;
  }
}
## Covariance matrix
dnn = .Call("DistMat", s, s, PACKAGE = "spBayesSurv");
corT = theta1*exp(-theta2*dnn)+(1-theta1)*diag(ntot);

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
#plot(xx, fi(xx, -1), "l", lwd=2, col=2)
#lines(xx, fi(xx, 1), "l", lwd=2, col=3)
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
dtotal = data.frame(s1=s1, s2=s2, y=y, x=x, delta=delta, yTrue=yTrue, id=blockid);
## Hold out npred=3 for prediction purpose
predindex = sample(1:ntot, npred);
dpred = dtotal[predindex,];
dtrain = dtotal[-predindex,];

# rename the variables 
d = dtrain; n=nrow(d); n;
d = d[order(d$id), ];
s = cbind(d$s1, d$s2);
y = d$y;
x = d$x;
delta =d$delta;

# FSA settings
knots = list(ss=ss, blockid=d$id);

# Prediction settings 
xpred = dpred$x;
s0 = cbind( dpred$s1, dpred$s2 );
prediction = list(spred=s0, xpred=xpred, predid=dpred$id);

###############################################################
# spatial copula DDP 
###############################################################
# MCMC parameters
nburn <- 100
nsave <- 200
nskip <- 0
ndisplay <- 500
mcmc <- list(nburn=nburn,
             nsave=nsave,
             nskip=nskip,
             ndisplay=ndisplay)
# Note larger nburn, nsave and nskip should be used in practice.
# Prior information
prior = list(N = 10, 
             a0 = 2, b0 = 2);

# current state values
state <- NULL;

# Fit the model
ptm<-proc.time()
res = spCopulaDDP( y = y,
              delta =delta, 
              x = x, 
              s = s, 
              prediction=prediction, 
              prior=prior, 
              mcmc=mcmc,
              state=state,
              FSA=FALSE,status=TRUE,
              knots=knots);
systime1=proc.time()-ptm; systime1;
# trace plots
par(mfrow = c(3,2))
w.save2 = res$w;
Kindex = which.max(rowMeans(w.save2));
traceplot(mcmc(w.save2[Kindex,]), main="w")
sig2.save2 = res$sigma2;
traceplot(mcmc(sig2.save2[Kindex,]), main="sig2")
beta.save2 = res$beta;
alpha.save2 = res$alpha;
traceplot(mcmc(beta.save2[2,Kindex,]), main="beta")
traceplot(mcmc(alpha.save2), main="alpha")
theta1.save2 = res$theta1;
theta2.save2 = res$theta2
traceplot(mcmc(theta1.save2), main="theta1")
traceplot(mcmc(theta2.save2), main="theta2")

## LPML
LPML2 = sum(log(res$cpo)); LPML2;
## MSPE
mean((dpred$yTrue-apply(res$Ypred, 1, median))^2); 

## number of non-negligible components
quantile(colSums(res$w>0.05))

## plots
par(mfrow = c(2,2));
xnew = c(-1, 1);
xpred = cbind(xnew); 
nxpred = nrow(xpred);
ygrid = seq(0,6.0,0.3); tgrid = exp(ygrid);
ngrid = length(ygrid);
estimates = GetCurves(res, xpred, ygrid, CI=c(0.05, 0.95));
fhat = estimates$fhat; 
Shat = estimates$Shat;
## density in y
plot(ygrid, fi(ygrid, xnew[1]), "l", lwd=2, ylim=c(0, 0.8), 
      xlim=c(0,6), main="density in y")
for(i in 1:nxpred){
  lines(ygrid, fi(ygrid, xnew[i]), lwd=2)
  lines(ygrid, fhat[,i], lty=2, lwd=2, col=4);
  lines(ygrid, estimates$fhatup[,i], lty=2, lwd=1, col=4);
  lines(ygrid, estimates$fhatlow[,i], lty=2, lwd=1, col=4);
}
## survival in y
plot(ygrid, 1-Fi(ygrid, xnew[1]), "l", lwd=2, ylim=c(0, 1), 
      xlim=c(0,6), main="survival in y")
for(i in 1:nxpred){
  lines(ygrid, 1-Fi(ygrid, xnew[i]), lwd=2)
  lines(ygrid, Shat[,i], lty=2, lwd=2, col=4);
  lines(ygrid, estimates$Shatup[,i], lty=2, lwd=1, col=4);
  lines(ygrid, estimates$Shatlow[,i], lty=2, lwd=1, col=4);
}
## density in t
plot(tgrid, fi(ygrid, xnew[1])/tgrid, "l", lwd=2, ylim=c(0, 0.15), 
      xlim=c(0,100), main="density in t")
for(i in 1:nxpred){
  lines(tgrid, fi(ygrid, xnew[i])/tgrid, lwd=2)
  lines(tgrid, fhat[,i]/tgrid, lty=2, lwd=2, col=4);
  lines(tgrid, estimates$fhatup[,i]/tgrid, lty=2, lwd=1, col=4);
  lines(tgrid, estimates$fhatlow[,i]/tgrid, lty=2, lwd=1, col=4);
}
## survival in t
plot(tgrid, 1-Fi(ygrid, xnew[1]), "l", lwd=2, ylim=c(0, 1), 
      xlim=c(0,100), main="survival in t")
for(i in 1:nxpred){
  lines(tgrid, 1-Fi(ygrid, xnew[i]), lwd=2)
  lines(tgrid, Shat[,i], lty=2, lwd=2, col=4);
  lines(tgrid, estimates$Shatup[,i], lty=2, lwd=1, col=4);
  lines(tgrid, estimates$Shatlow[,i], lty=2, lwd=1, col=4);
}
}
