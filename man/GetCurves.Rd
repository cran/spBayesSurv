\name{GetCurves}
\alias{GetCurves}
\alias{plot.survregbayes}
\alias{plot.frailtyGAFT}
\alias{plot.SuperSurvRegBayes}
\alias{plot.indeptCoxph}
\alias{plot.anovaDDP}
\alias{plot.spCopulaCoxph}
\alias{plot.spCopulaDDP}
\alias{plot.SpatDensReg}

\title{Density, Survival, and Hazard Estimates}
\description{
This function estimates density, survival, and hazard functions given covariates.
}

\usage{
GetCurves(x, xnewdata, xtfnewdata, tgrid = NULL, ygrid = NULL, 
          frail = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{survregbayes}(x, xnewdata, tgrid = NULL, 
     frail = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{frailtyGAFT}(x, xnewdata, xtfnewdata, tgrid = NULL, 
     frail = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{SuperSurvRegBayes}(x, xnewdata, tgrid = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{indeptCoxph}(x, xnewdata, tgrid = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{anovaDDP}(x, xnewdata, tgrid = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{spCopulaCoxph}(x, xnewdata, tgrid = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{spCopulaDDP}(x, xnewdata, tgrid = NULL, CI = 0.95, PLOT = TRUE, ...)
\method{plot}{SpatDensReg}(x, xnewdata, ygrid = NULL, CI = 0.95, PLOT = TRUE, ...)
}

\arguments{
   \item{x}{ an object obtained from the functions \code{\link{survregbayes}}, \code{\link{frailtyGAFT}}, \code{\link{SuperSurvRegBayes}}, \code{\link{indeptCoxph}}, \code{\link{anovaDDP}}, \code{\link{spCopulaCoxph}}, \code{\link{spCopulaDDP}} and \code{\link{SpatDensReg}}.}
   
   \item{xnewdata}{A data frame in which to look for variables with which to obtain estimated curves.}
   
   \item{xtfnewdata}{A data frame in which to look for variables with which to obtain estimated curves, used only for \code{\link{frailtyGAFT}}.}
   
   \item{tgrid}{a vector of grid points indicating where the curves will be estimated.} 
   
   \item{ygrid}{a vector of grid points indicating where the curves will be estimated, used only for \code{\link{SpatDensReg}}.} 
   
   \item{frail}{an optional matrix of posterior frailty values for \code{\link{survregbayes}} and \code{\link{frailtyGAFT}}, where the rows refer to clusters/regions and the number of columns is the length of thined MCMC chain. The default is to set frailties to be zero.} 

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
\code{\link{survregbayes}}, \code{\link{frailtyGAFT}}, \code{\link{SuperSurvRegBayes}}, \code{\link{indeptCoxph}}, \code{\link{anovaDDP}}, \code{\link{spCopulaCoxph}}, \code{\link{spCopulaDDP}} and \code{\link{SpatDensReg}}
}

\author{
Haiming Zhou and Timothy Hanson
}
