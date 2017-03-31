\name{GetCurves}
\alias{GetCurves}
\alias{plot.indeptCoxph}
\alias{plot.spCopulaCoxph}
\alias{plot.anovaDDP}
\alias{plot.spCopulaDDP}

\title{Density, Survival, and Hazard Estimates}
\description{
This function estimates density, survival, and hazard functions given covariates.
}

\usage{
GetCurves(x, xpred=NULL, tgrid=NULL, CI=0.95, PLOT=FALSE, ...)
## An alternative useage is 
## plot(x, xpred=NULL, tgrid=NULL, CI=0.95, PLOT=FALSE, ...)
}

\arguments{
   \item{x}{ an object obtained from the functions \code{anovaDDP}, \code{spCopulaDDP},  \code{indeptCoxph} and \code{spCopulaCoxph}.}
   
   \item{xpred}{an npred by p vector of covariates matrix used for curve estimates.}
   
   \item{tgrid}{a vector of grid points indicating where the curves will be estimated.} 

   \item{CI}{a numeric value indicating the level of credible interval.}
   
   \item{PLOT}{a logical value indicating whether the estimated survival curves will be plotted.}
   
   \item{...}{further arguments to be passed to or from other methods.}
}

\details{
  This function estimates density, survival, and hazard functions given covariates.
}

\value{
  Use \code{names} to find out what they are, where \code{fhat} represents density, \code{Shat} represents survival, \code{hhat} represents hazard. The credible bands are also provided, e.g., \code{Shatlow} represents the lower band and \code{Shatup} represents the upper band. 
}

\seealso{
\code{\link{anovaDDP}, \link{spCopulaDDP}, \link{indeptCoxph}, \link{spCopulaCoxph}}
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\references{
Zhou, H., Hanson, T., and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, 71(4): 1101-1110. 

De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. \emph{Biometrics}, 65(3): 762-771. 

}