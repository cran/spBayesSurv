\name{frailtyGAFT}
\alias{frailtyGAFT}
\alias{print.frailtyGAFT}
\alias{summary.frailtyGAFT}
\alias{print.summary.frailtyGAFT}

\title{Generalized Accelerated Failure Time Frailty Model}
\description{
This function fits a generalized accelerated failure time frailty model (Zhou, et al., 2017) for clustered and/or areal-level time-to-event data. Note that the function arguments are slightly different with those presented in the original paper; see Zhou, Hanson and Zhang (2018) for new examples. 
}
   
\usage{
frailtyGAFT(formula, data, na.action, 
            mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
            prior=NULL, state=NULL, Proximity=NULL, Coordinates=NULL,
            DIST=NULL, scale.designX=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them. To include CAR frailties, add \code{frailtyprior("car",ID)} to the formula, where \code{ID} is an n dimensional vector of cluster ID numbers. Furthermore, use \code{frailtyprior("iid",ID)} for Gaussian exchangeable frailties, use \code{frailtyprior("grf",ID)} for Gaussian random fields (GRF) frailties, and exclude the term \code{frailtyprior()} for non-frailty models. Note: the data need to be sorted by ID. }
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. See Zhou, Hanson and Zhang (2018) for more detailed hyperprior specifications.}
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters/regions. If CAR frailty model is specified in the formula, \code{Proximity} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID. }
   
   \item{Coordinates}{an m by d coordinates matrix, where m is the number of clusters/regions, d is the dimension of coordiantes. If GRF frailty model is specified in the formula, \code{Coordinates} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID.}
   \item{DIST}{This is a function argument, used to calculate the distance. The default is Euclidean distance (\code{fields::rdist}). This function should have two arguments (X1,X2), where X1 and X2 are matrices with coordinates as the rows. The returned value of this function should be the pairwise distance matrix. If nrow(X1)=m and nrow(X2)=n then the function should return an m by n matrix of all distances between these two sets of points.}
   \item{scale.designX}{flag to indicate whether the design matrix X and Xtf will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
}

\details{
This function fits a a generalized accelerated failure time frailty model (Zhou, et al., 2017) for clustered and/or areal-level time-to-event data. Note that the function arguments are slightly different with those presented in the original paper of Zhou, et al. (2017); see Zhou, Hanson and Zhang (2018) for new examples. 
}

\value{
The \code{frailtyGAFT} object is a list containing at least the following components:
\item{modelname}{
the name of the fitted model
}
\item{terms}{
the \code{\link{terms}} object used
}
\item{coefficients}{
a named vector of coefficients. The last two elements are the estimates of scale parameter sigma and precision parameter alpha involved in the LDTFP prior. 
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
\item{pce}{
the number of columns in the model matrix including the intercept
}
\item{ptf}{
the number of columns in the model matrix used in the LDTFP baseline including the intercept
}
\item{Surv}{
the \code{\link{Surv}} object used 
}
\item{X.scaled}{
the n by pce-1 scaled design matrix
}
\item{X}{
the n by pce-1 orginal design matrix
}
\item{Xtf.scaled}{
the n by ptf-1 scaled design matrix used in the LDTFP baseline
}
\item{Xtf}{
the n by ptf-1 orginal design matrix used in the LDTFP baseline
}
\item{sigma2}{
the vector of posterior samples for the variance parameter used in the LDTFP prior. 
}
\item{beta}{
the pce by nsave matrix of posterior samples for the coefficients in the \code{linear.predictors} which includes the intercept
}
\item{beta.scaled}{
the pce by nsave matrix of posterior samples for the coefficients in the \code{linear.predictors}. Note that these posterior samples are based scaled design matrix.
}
\item{alpha}{
the vector of posterior samples for the precision parameter alpha in the LDTFP prior. 
}
\item{maxL}{
the truncation level used in the LDTFP prior. 
}
\item{logt}{
the n by nsave matrix of posterior samples for log survival times. 
}
\item{cpo}{
the length n vector of the stabilized estiamte of CPO; used for calculating LPML
}
\item{accept_beta}{
the acceptance rate in the posterior sampling of beta coefficient vector
}
\item{frail.prior}{
the frailty prior used in \code{\link{frailtyprior}} 
}
\item{BF}{
the Bayes factors for testing necessariness of each stratification covariate.  
}

The object will also have the following components when frailty models are fit:
\item{v}{
the nID by nsave matrix of posterior samples for frailties, where nID is the number of clusters considered.
}
\item{tau2}{
the vector of posterior samples for tau2 involved in the IID, GRF or CAR frailty prior. 
}
\item{ID}{
the cluster ID used in \code{\link{frailtyprior}} 
}
If GRF frailties are used, the object will also have:
\item{Coordinates}{
the \code{Coordinates} matrix used in \code{\link{survregbayes}} 
}
\item{ratephi}{
the acceptance rates in the posterior sampling of phi involved in the GRF prior 
}
\item{phi}{
the vector of posterior samples for phi involved in the GRF prior 
}
}

\seealso{
\code{\link{baseline}}, \code{\link{frailtyprior}}, \code{\link{survregbayes}}, \code{\link{rdist}}
}

\references{
Zhou, H., Hanson, T., and Zhang, J. (2020). spBayesSurv: Fitting Bayesian Spatial Survival Models Using R. \emph{Journal of Statistical Software}, 92(9): 1-33.

Zhou, H., Hanson, T., and Zhang, J. (2017). Generalized accelerated failure time spatial frailty model for arbitrarily censored data. \emph{Lifetime Data Analysis}, 23(3): 495-515.
}


\author{
Haiming Zhou and Timothy Hanson
}

\examples{
###############################################################
# A simulated data: GAFT spatial frailty model
###############################################################
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
library(MASS)

## True densities
set.seed(1)
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
mi = 2;
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
v0 = mvrnorm(1, mu=rep(0,m-1), Sigma=covT);
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
mcmc=list(nburn=200, nsave=200, nskip=0, ndisplay=100)
# Note larger nburn, nsave and nskip should be used in practice.

# Fit the model
ptm<-proc.time()
res = frailtyGAFT(Surv(tleft, tright, type="interval2")~x1+x2+baseline(x1, x2)+
                    frailtyprior(prior="car", ID),  data=d, mcmc=mcmc, prior=prior, 
                  Proximity=W);
summary(res);
systime1=proc.time()-ptm; systime1;

### trace plots
par(mfrow = c(3,1))
traceplot(mcmc(res$beta[1,]), main="beta1");
traceplot(mcmc(res$beta[2,]), main="beta2");
traceplot(mcmc(res$beta[3,]), main="beta3");

####################################################################
## Get curves
####################################################################
par(mfrow = c(1,1))
xpred = data.frame(x1=c(1,1.5), x2=c(0,1))
xtfpred = xpred;
plot(res, xnewdata=xpred, xtfnewdata=xtfpred, CI=0.9);
}

