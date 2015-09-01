\name{GetCurves}
\alias{GetCurves}

\title{Density, Survival, and Hazard Estimates}
\description{
This function estimates mariginal density, survival, and hazard functions given covariates.
}
   
\usage{
GetCurves(fit, xpred, ygrid, CI=c(0.05, 0.95))
}

\arguments{
   \item{fit}{ an object obtained from functions \code{anovaDDP}, \code{spCopulaDDP},  \code{indeptCoxph} and \code{spCopulaCoxph}.}
   
   \item{xpred}{ an npred by p vector of covariates matrix used for curve estimates.}
   
   \item{ygrid}{ a vector of grid points (in log survival time scale) indicating where the curves will be estimated.} 

   \item{CI}{ a vector to indicate the level of credible interval, where first element indicate the lower probability.}
}

\details{
  This function estimates mariginal density, survival, and hazard functions given covariates.
}

\value{
  Use \code{names} to find out what they are, where \code{fhat} represents density, \code{Shat} represents survival and 
  \code{Hhat} represents hazard.
}

\seealso{
\code{\link{anovaDDP}, \link{spCopulaDDP}, \link{indeptCoxph}, \link{spCopulaCoxph}}
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\references{
De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. \emph{Biometrics}, 65(3): 762-771. 

Zhou, H., Hanson, T. and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, accepted. 
}