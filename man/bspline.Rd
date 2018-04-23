\name{bspline}
\alias{bspline}
\alias{makepredictcall.bspline}

\title{Generate a Cubic B-Spline Basis Matrix}
\description{
Generate the B-spline basis matrix for a cubic spline with the first and last columns dropped. 
}
   
\usage{
bspline(x, df=NULL, knots=NULL, Boundary.knots = range(x))
}

\arguments{
   \item{x}{the predictor variable. Missing values are allowed.}
   
   \item{df}{degrees of freedom; one can specify \code{df} rather than \code{knots}; the function chooses \code{df-2} inner knots at suitable quantile of \code{x} (which will ignore missing values). The default, \code{NULL}, corresponds to one inner knots, i.e. \code{df=3}.}
   
   \item{knots}{the internal breakpoints that define the spline. The default is \code{NULL}, which corresponds the median for one knot, quantiles for more knots.}
   
   \item{Boundary.knots}{boundary points at which to anchor the B-spline basis (default the range of the non-NA data).}
}

\seealso{
\code{\link{predict.bspline}}
}

\references{
Hastie, T. J. (1992) Generalized additive models. Chapter 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
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
