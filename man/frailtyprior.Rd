\name{frailtyprior}
\alias{frailtyprior}

\title{Frailty prior specification}
\description{
This function allows one to add a frailty term to the linear predictor of semiparametric PH, PO and AFT models.
}
   
\usage{
frailtyprior(prior="car", ...)
}

\arguments{
   \item{prior}{a name string to be entered, e.g, \code{"car"}, \code{"iid"}, \code{"kriging"} or \code{"grf"}; see the example in \code{\link{frailtyGAFT}} and \code{\link{survregbayes2}}.}
   
   \item{...}{Cluster ID to be entered for clustered data or locations for point-referenced data; see the example in \code{\link{frailtyGAFT}} and \code{\link{survregbayes2}}.}
}

\references{
Zhou, H. and Hanson, T. (2016). Bayesian semiparametric models for spatially correlated arbitrarily censored data. In preparation.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}


