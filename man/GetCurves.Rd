\name{GetCurves}
\alias{GetCurves}

\title{Density, Survival, and Hazard Estimates}
\description{
This function estimates mariginal density, survival, and hazard functions given covariates.
}
   
\usage{
GetCurves(fit, xpred, ygrid, CI=c(0.05, 0.95), xtfpred=NULL, 
          IDpred=NULL, AFTbaseline=FALSE)
}

\arguments{
   \item{fit}{       an object obtained from functions including \code{anovaDDP}, \code{spCopulaDDP}, 
                     \code{indeptCoxph}, \code{spCopulaCoxph}, and \code{frailtyGAFT}.}
   
   \item{xpred}{	   an npred by p vector of covariates matrix used for curve estimates.}
   
   \item{ygrid}{     a vector of grid points (in log survival time scale) indicating where the curves will be estimated.} 

   \item{CI}{        a vector to indicate the level of credible interval, where first element indicate the lower probability. }
   \item{xtfpred}{   an npred by q vector of covariates matrix if \code{frailtyGAFT} is fitted. }
   \item{IDpred}{    an npred by 1 vector indicates the cluster ID of each predicted subject if frailty model is fitted.}
   \item{AFTbaseline}{ a logical value indicate if only baseline curves are needed. The default is FALSE.}
}

\details{
  This function estimates mariginal density, survival, and hazard functions given covariates.
}

\value{
  Use \code{names} to find out what they are, where \code{fhat} represents density, \code{Shat} represents survival and 
  \code{Hhat} represents hazard.
}

\seealso{
\code{\link{anovaDDP}, \link{spCopulaDDP}, \link{indeptCoxph}, \link{spCopulaCoxph}, \link{frailtyGAFT}}
}

\author{
Haiming Zhou \email{<zhouh@email.sc.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\references{
De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. Biometrics, 65(3): 762-771. 

Zhou, H., Hanson, T. and Knapp, R. (2014+). Marginal Bayesian nonparametric model for the time-to-extinction of the mountain yellow-legged frog. Biometrics. In revision. 

Zhou, H., Hanson, T. and Zhang, J. (2014+). Generalized accelerated failure time spatial frailty model for arbitrarily censored data. Submitted
}