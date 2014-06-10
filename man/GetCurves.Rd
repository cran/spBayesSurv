\name{GetCurves}
\alias{GetCurves}

\title{Obtain estimated mariginal density, survival, and hazard functions given covariates}
\description{
    This function estimates mariginal density, survival, and hazard functions given covariates.
}
   
\usage{
GetCurves(fit, xpred, ygrid, CI=c(0.05, 0.95))
}

\arguments{
   \item{fit}{       an object obtained from functions including \code{anovaDDP}, \code{spCopulaDDP}, 
                     \code{indeptCoxph} and \code{spCopulaCoxph}.}
   
   \item{xpred}{	   an n by p vector of covariates matrix used for curve estimates.}
   
   \item{ygrid}{     a vector of grid points (in log survival time scale) indicating where the curves will be estimated.} 

   \item{CI}{        a vector to indicate the level of credible interval, where first element indicate the lower probability. }
}

\details{
  This function fits an ANOVA DDP model (De Iorio et al., 2009) for survival analysis of right 
  censored time-to-event data.
}

\value{
  Use \code{names} to find out what they are, where \code{fhat} represents density, \code{Shat} represents survival and 
  \code{Hhat} represents hazard.
}

\seealso{
\code{\link{anovaDDP}, \link{spCopulaDDP}, \link{indeptCoxph}, \link{spCopulaCoxph}}
}

\author{
Haiming Zhou \email{<zhouh@email.sc.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}
