\name{spBayesSurv-package}
\alias{spBayesSurv-package}
\alias{spBayesSurv}
\docType{package}
\title{
Bayesian Modeling and Analysis of Spatially Correlated Survival Data
}
\description{
This package provides several Bayesian survival models for spatial/non-spatial survival data: marginal Bayesian Nonparametric models, marginal Bayesian proportional hazards models, and generalized accelerated failure time frailty models.
}
\details{
\tabular{ll}{
Package: \tab spBayesSurv\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2014-10-24\cr
License: \tab GPL (>= 2)\cr
}
This package provides several Bayesian survival models for spatial/non-spatial survival data, including marginal Bayesian Nonparametric models, where \code{spCopulaDDP} is for point-referenced data and \code{anovaDDP} is for non-spatial data; marginal Bayesian proportional hazards models, where \code{spCopulaCox} is for point-referenced data and \code{indeptCoxph} is for non-spatial data; and generalized accelerated failure time frailty models via \code{frailtyGAFT} for both spatial and non-spatial survival data.
}
\author{
Haiming Zhou <zhouh@email.sc.edu> and Tim Hanson <hansont@stat.sc.edu>
}
\references{
De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. Biometrics, 65(3): 762-771. 

Zhou, H., Hanson, T. and Knapp, R. (2014+). Marginal Bayesian nonparametric model for the time-to-extinction of the mountain yellow-legged frog. Biometrics. In revision. 

Zhou, H., Hanson, T. and Zhang, J. (2014+). Generalized accelerated failure time spatial frailty model for arbitrarily censored data. Submitted
}
\keyword{ package }
