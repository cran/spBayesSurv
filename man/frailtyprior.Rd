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
   \item{prior}{a name string to be entered, e.g, \code{"car"}, \code{"iid"}, or \code{"grf"}; see the example in \code{\link{frailtyGAFT}} and \code{\link{survregbayes2}}.}
   
   \item{...}{Cluster ID to be entered for clustered data or locations for point-referenced data; see the example in \code{\link{frailtyGAFT}} and \code{\link{survregbayes2}}.}
}

\author{
Haiming Zhou and Timothy Hanson
}


