\name{predict.bspline}
\alias{predict.bspline}

\title{Evaluate a Cubic Spline Basis}
\description{
Evaluate a predefined spline basis at given values.
}
   
\usage{
\method{predict}{bspline}(object, newx, ...)
}

\arguments{
   \item{object}{the result of a call to \code{\link{bspline}} having attributes describing knots, df, etc..}
   
   \item{newx}{the \code{x} values at which evaluations are required.}
   
   \item{...}{Optional additional arguments. At present no additional arguments are used.}
}

\seealso{
\code{\link{bspline}}
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}

\examples{
require(stats)
basis <- bspline(women$height, df = 5)
newX <- seq(58, 72, length.out = 51)
# evaluate the basis at the new data
predict(basis, newX)
}
