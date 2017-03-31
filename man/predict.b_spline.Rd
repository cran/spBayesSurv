\name{predict.b_spline}
\alias{predict.b_spline}

\title{Evaluate a Cubic Spline Basis}
\description{
Evaluate a predefined spline basis at given values.
}
   
\usage{
\method{predict}{b_spline}(object, newx, ...)
}

\arguments{
   \item{object}{the result of a call to \code{\link{b_spline}} having attributes describing knots, degree, etc..}
   
   \item{newx}{the \code{x} values at which evaluations are required.}
   
   \item{...}{Optional additional arguments. At present no additional arguments are used.}
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
require(stats)
basis <- b_spline(women$height, df = 5)
newX <- seq(58, 72, length.out = 51)
# evaluate the basis at the new data
predict(basis, newX)
}
