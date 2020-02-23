\name{cox.snell.survregbayes}
\alias{cox.snell.survregbayes}

\title{Cox-Snell Diagnostic Plot}
\description{
This function provides the Cox-Snell diagnostic plot (Zhou and Hanson, 2018) for fitting  for Bayesian semiparametric survival models.
}
   
\usage{
cox.snell.survregbayes(x, ncurves = 10, PLOT = TRUE)
}

\arguments{
   \item{x}{an object obtained from the function \code{\link{survregbayes}}.}
   
   \item{ncurves}{the number of posterior draws.}
   
   \item{PLOT}{a logical value indicating whether the estimated survival curves will be plotted.}
}


\value{
The function returns the plot (if \code{PLOT = TRUE}) and a list with the following components:
\item{resid}{
the \code{\link{Surv}} object using the everaged residuals
}
\item{resid*}{
the \code{*} goes from 1 to \code{ncurves}, and each is the \code{\link{Surv}} object using the residuals obtained from each posterior draw. 
}
\item{St1}{
the \code{n} by \code{ncurves} matrix of survival functions evaluated at left endpoints of the time interval.
}
\item{St2}{
the \code{n} by \code{ncurves} matrix of survival functions evaluated at right endpoints of the time interval.
}
\item{Delta}{
The status indicator: 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
}
}

\seealso{
\code{\link{survregbayes}} 
}

\references{
Zhou, H. and Hanson, T. (2018). A unified framework for fitting Bayesian semiparametric models to arbitrarily censored survival data, including spatially-referenced data. \emph{Journal of the American Statistical Association}, 113(522): 571-581.
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Tim Hanson \email{<hansont@stat.sc.edu>}
}


