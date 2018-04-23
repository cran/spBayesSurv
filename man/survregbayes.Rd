\name{survregbayes}
\alias{survregbayes}
\alias{print.survregbayes}
\alias{summary.survregbayes}
\alias{print.summary.survregbayes}

\title{Bayesian Semiparametric Survival Models}
\description{
This function fits semiparametric proportional hazards (PH), proportional odds (PO), accelerated failture time (AFT) and accelerated hazards (AH) models. Both georeferenced (location observed exactly) and areally observed (location known up to a geographic unit such as a county) spatial locations can be handled. Georeferenced data are modeled with Gaussian random field (GRF) frailties whereas areal data are modeled with a conditional autoregressive (CAR) prior on frailties. For non-spatial clustered data, an IID Gaussian frailties are assumed. Variable selection is also incorporated. The function can fit both Case I and Case II interval censored data, as well as standard right-censored, uncensored, and mixtures of these. The transformed Bernstein Polynomial (TBP) prior is used for fitting the baseline survival function. The full scale approximation (FSA) (Sang and Huang, 2012) could be used to inverse the spatial correlation matrix for georeferenced data. The function also fits all these models without frailties. The logarithm of the pseudo marginal likelihood (LPML), the deviance information criterion (DIC), and the Watanabe-Akaike information criterion (WAIC) are provided for model comparison.
}
   
\usage{
survregbayes(formula, data, na.action, survmodel="PH", dist="loglogistic", 
             mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
             prior=NULL, state=NULL, selection=FALSE, Proximity=NULL, 
             truncation_time=NULL, subject.num=NULL, Knots=NULL,
             Coordinates=NULL, DIST=NULL, InitParamMCMC=TRUE, 
             scale.designX=TRUE)
}

\arguments{
   \item{formula}{a formula expression with the response returned by the \code{Surv} function in the \code{survival} package. It supports right-censoring, left-censoring, interval-censoring, and mixtures of them. To include CAR frailties, add \code{frailtyprior("car",ID)} to the formula, where \code{ID} is an n dimensional vector of cluster ID numbers. Furthermore, use \code{frailtyprior("iid",ID)} for Gaussian exchangeable frailties, use \code{frailtyprior("grf",ID)} for Gaussian random fields (GRF) frailties, and exclude the term \code{frailtyprior()} for non-frailty models. Note: the data need to be sorted by ID.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{survmodel}{a character string for the assumed survival model. The options include \code{"PH"} for proportional hazards, \code{"PO"} for proportional odds, and \code{"AFT"} for accelerated failture time.}
   
   \item{dist}{centering distribution for TBP. Choices include \code{"loglogistic"}, \code{"lognormal"}, and \code{"weibull"}.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The function itself provides all default priors. Note if FSA is used, the number of knots \code{knots} and the number of blocks \code{nblock} are specified here; see examples below. See Zhou, Hanson and Zhang (2018) for more detailed hyperprior specifications.}
   
   \item{state}{a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis.} 
   
   \item{selection}{flag to indicate whether variable selection is performed, where \code{FALSE} indicates that no variable selection will be performed.}
   
   \item{Proximity}{an m by m symetric adjacency matrix, where m is the number of clusters/regions. If CAR frailty model is specified in the formula, \code{Proximity} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID.}

   \item{truncation_time}{a vector of left-trucation times with length n.}
   
   \item{subject.num}{a vector of subject id numbers when time dependent covariates are considered. For example, for subject 1 time dependent covariates are recorded over [0,1), [1,3), and for subject 2 covariates are recorded over [0,2), [2,3), [3,4). Suppose we only have two subjects, i.e. n=2. In this case, we save the data in the long format, set truncation_time=c(0,1,0,2,3) and subject.num=c(1,1,2,2,2).}
   
   \item{Knots}{an nknots by d matrix, where nknots is the number of selected knots for FSA, and d is the dimension of each location. If \code{Knots} is not specified, the space-filling algorithm will be used to find the knots.}
   
   \item{Coordinates}{an m by d coordinates matrix, where m is the number of clusters/regions, d is the dimension of coordiantes. If GRF frailty model is specified in the formula, \code{Coordinates} is required; otherwise it is ignored. Note: this matrix should be specified according to the data that have been sorted by ID.}
   
   \item{DIST}{This is a function argument, used to calculate the distance. The default is Euclidean distance (\code{fields::rdist}). This function should have two arguments (X1,X2), where X1 and X2 are matrices with coordinates as the rows. The returned value of this function should be the pairwise distance matrix. If nrow(X1)=m and nrow(X2)=n then the function should return an m by n matrix of all distances between these two sets of points.}
   
   \item{InitParamMCMC}{flag to indicate whether an initial MCMC will be run based on the centering parametric model, where \code{TRUE} indicates yes. }
   
   \item{scale.designX}{flag to indicate whether the design matrix X will be centered by column means and scaled by column standard deviations, where \code{TRUE} indicates yes. The default is \code{TRUE} for improving numerical stability. Even when it is scaled, the reported regression coefficients are in original scales. Note if we want to specify informative priors for regression coefficients, these priors should correspond to scaled predictors when \code{scale.designX=TRUE}.}
}


\value{
This class of objects is returned by the \code{survregbayes} function to represent a fitted Bayesian semiparametric survival model. Objects of this class have methods for the functions \code{print}, \code{summary} and \code{plot}.

The \code{survregbayes} object is a list containing at least the following components:
\item{modelname}{
the name of the fitted model
}
\item{terms}{
the \code{\link{terms}} object used
}
\item{dist}{
the centering distribution used in the TBP prior on baseline survival function
}
\item{survmodel}{
the model fitted 
}
\item{coefficients}{
a named vector of coefficients. The last two elements are the estimates of theta1 and theta2 involved in the centering baseline survival function. 
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
\item{nsubject}{
the number of subjects/individuals, which is equal to n in the absence of time-dependent covariates
}
\item{subject.num}{
the vector of subject id numbers when time dependent covariates are considered
}
\item{truncation_time}{
the vector of left-trucation times
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
\item{theta.scaled}{
the 2 by nsave matrix of posterior samples for theta1 and theta2 involved in the centering baseline survival function. Note that these posterior samples are based scaled design matrix.
}
\item{beta.scaled}{
the p by nsave matrix of posterior samples for the coefficients in the \code{linear.predictors}. Note that these posterior samples are based scaled design matrix.
}
\item{alpha}{
the vector of posterior samples for the precision parameter alpha in the TBP prior. 
}
\item{maxL}{
the truncation level used in the TBP prior. 
}
\item{weight}{
the maxL by nsave matrix of posterior samples for the weights in the TBP prior. 
}
\item{cpo}{
the length n vector of the stabilized estiamte of CPO; used for calculating LPML
}
\item{pD}{
the effective number of parameters involved in DIC
}
\item{DIC}{
the deviance information criterion (DIC)
}
\item{pW}{
the effective number of parameters involved in WAIC
}
\item{WAIC}{
the Watanabe-Akaike information criterion (WAIC)
}
\item{Surv.cox.snell}{
the \code{\link{Surv}} object used for Cox-Snell residual plot. This is not recommended for frailty models, for which please use the function \code{\link{cox.snell.survregbayes}}. 
}
\item{ratetheta}{
the acceptance rate in the posterior sampling of theta vector involved in the centering baseline survival function
}
\item{ratebeta}{
the acceptance rate in the posterior sampling of beta coefficient vector
}
\item{rateYs}{
the acceptance rate in the posterior sampling of weight vector involved in the TBP prior
}
\item{ratec}{
the acceptance rate in the posterior sampling of precision parameter alpha involved in the TBP prior
}
\item{frail.prior}{
the frailty prior used in \code{\link{frailtyprior}} 
}
\item{selection}{
whether the variable selection was performed 
}
\item{initial.values}{
the list of initial values used for the parameters 
}
\item{BF.baseline}{
the Bayes factor for comparing the parametric baseline vs. the TBP baseline 
}
\item{BF.bs}{
Bayes factors for testing linearity when \code{\link{bspline}} is added to the \code{linear.predictors}
}

The object will also have the following components when frailty models are fit:
\item{v}{
the nID by nsave matrix of posterior samples for frailties, where nID is the number of clusters considered.
}
\item{ratev}{
the vector of acceptance rates in the posterior sampling of frailties 
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
\item{Knots}{
the \code{Knots} matrix used in \code{\link{survregbayes}} 
}
If the variable selection is performed, the object will also include:
\item{gamma}{
the p by nsave matrix of posterior samples for gamma involved in the variable selection
}

}

\seealso{
\code{\link{frailtyprior}}, \code{\link{cox.snell.survregbayes}}, \code{\link{rdist}}, \code{\link{rdist.earth}}
}

\references{
Zhou, H., Hanson, T., and Zhang, J. (2018). spBayesSurv: Fitting Bayesian Spatial Survival Models Using R. \emph{Journal of Statistical Software}, accept minor.

Zhou, H. and Hanson, T. (2018). A unified framework for fitting Bayesian semiparametric models to arbitrarily censored survival data, including spatially-referenced data. \emph{Journal of the American Statistical Association}, accepted.

Sang, H. and Huang, J. Z. (2012). A full scale approximation of covariance functions for large spatial data sets. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 74(1), 111-132.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
rm(list=ls())
library(survival)
library(spBayesSurv)
library(coda)
library(MASS)
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
## correlation function
rho_Exp = function(dis, phi) exp(-(dis*phi));

###############################################################################
########################### Start to simulation ###############################
###############################################################################
phiT=1; sill=0.9999; ## phiT is the range parameter phi. 
tau2T = 1; ## true frailty variance; 
m = 50; mi=2
id=rep(1:m, each=mi)
mseq = rep(mi, m);
n = sum(mseq);
s1 = runif(m, 0, 10); s2 = runif(m, 0, 10);
locs = cbind(s1, s2); 
ss = cbind(rep(locs[,1],each=mi), rep(locs[,2],each=mi)); ### the locations. 
Dmm = .Call("DistMat", t(locs), t(locs), PACKAGE = "spBayesSurv");
Rmm = sill*rho_Exp(Dmm, phiT)+diag(1-sill, m, m);
v = mvrnorm(1, mu=rep(0,m), Sigma=tau2T*Rmm); 
vn = rep(v, each=mi)
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
d = data.frame(t1=t1, t2=t2, x1=x1, x2=x2, delta=delta, 
               s1=ss[,1], s2=ss[,2], id=id); 
table(d$delta)/n;

##-------------spBayesSurv-------------------##
# MCMC parameters
nburn=200; nsave=200; nskip=0;
# Note larger nburn, nsave and nskip should be used in practice.
mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
prior = list(maxL=15, a0=1, b0=1, nknots=m, nblock=m, nu=1);
# here if nknots<m, FSA will be used with nblock=m.
cor.dist = function(x1, x2) rdist(x1,x2)
ptm<-proc.time()
res1 = survregbayes(formula = Surv(t1, delta)~x1+x2+
                       frailtyprior("grf", id), data=d, InitParamMCMC=FALSE,
                     survmodel="PH", prior=prior, mcmc=mcmc, DIST=cor.dist,
                     dist="loglogistic", Coordinates = locs);
## Or equivalently formula = Surv(t1, t2, type="interval2") can also be used.
## Note InitParamMCMC=FALSE is used for speeding,
## InitParamMCMC=TRUE is recommended in general. 
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
par(mfrow=c(1,1));
wide=0.2;
tgrid = seq(1e-10,4,wide);
ngrid = length(tgrid);
p = length(betaT); # number of covariates
newdata = data.frame(x1=c(0,0), x2=c(0,1))
plot(res1, xnewdata=newdata, tgrid=tgrid, PLOT=TRUE);

## Cox-Snell plot
set.seed(1)
cox.snell.survregbayes(res1, ncurves=2, PLOT=TRUE);
}

