\name{spBayesSurv-package}
\alias{spBayesSurv-package}
\alias{spBayesSurv}
\docType{package}
\title{
Bayesian Modeling and Analysis of Spatially Correlated Survival Data
}
\description{
This package provides several Bayesian survival models for spatial/non-spatial survival data: marginal Bayesian Nonparametric models, marginal Bayesian proportional hazards models, generalized accelerated failure time frailty models, and standard semiparametric frailty models within the context of proportional hazards, proportional odds and accelerated failure time.
}
\details{
\tabular{ll}{
Package: \tab spBayesSurv\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2015-08-31\cr
License: \tab GPL (>= 2)\cr
}
This package provides several Bayesian survival models for spatial/non-spatial survival data, including marginal Bayesian Nonparametric models, where \code{spCopulaDDP} is for point-referenced data and \code{anovaDDP} is for non-spatial data; marginal Bayesian proportional hazards models, where \code{spCopulaCox} is for point-referenced data and \code{indeptCoxph} is for non-spatial data; generalized accelerated failure time frailty models via \code{frailtyGAFT} for both spatial and non-spatial survival data; and standard semiparametric frailty models via \code{survregbayes} within the context of proportional hazards, proportional odds and accelerated failure time.
}
\author{
Haiming Zhou <zhouh@niu.edu> and Tim Hanson <hansont@stat.sc.edu>
}
\references{
De Iorio, M., Johnson, W. O., Mueller, P., and Rosner, G. L. (2009). Bayesian nonparametric nonproportional hazards survival modeling. \emph{Biometrics}, 65(3): 762-771. 

Zhou, H., Hanson, T. and Knapp, R. (2015). Marginal Bayesian nonparametric model for time to disease arrival of threatened amphibian populations. \emph{Biometrics}, accepted.  

Zhou, H., Hanson, T. and Zhang, J. (2015). Generalized accelerated failure time spatial frailty model for arbitrarily censored data. \emph{Lifetime Data Analysis}. Revision submitted.

Zhou, H. and Hanson, T. (2015). Bayesian spatial survival models. In \emph{Nonparametric Bayesian Inference in Biostatistics} (pp. 215-246). Springer International Publishing.

}
\keyword{ package }
