\name{LeukSurv}
\alias{LeukSurv}

\title{The Leukemia Survival Data}
\description{
A dataset on the survival of acute myeloid leukemia in 1,043 pateietns, first analyzed by Henderson et al. (2002). It is of interest to investigate possible spatial variation in survival after accounting for known subject-specific prognostic factors, which include age, sex, white blood cell count (wbc) at diagnosis, and the Townsend score (tpi) for which higher values indicates less affluent areas. Both exact residential locations of all patients and their administrative districts (24 districts that make up the whole region) are available. 
}
   
\usage{data(LeukSurv)}

\format{
  \tabular{ll}{
    time:\tab survival time in days\cr
    cens:\tab right censoring status 0=censored, 1=dead\cr
    xcoord: \tab coordinates in x-axis of residence \cr
    ycoord: \tab coordinates in y-axis of residence \cr
    age:\tab age in years\cr
    sex:\tab  male=1 female=0\cr
    wbc:\tab white blood cell count at diagnosis, truncated at 500 \cr
    tpi:\tab the Townsend score for which higher values indicates less affluent areas\cr
    district:\tab administrative district of residence \cr
  }
}

\source{
Henderson, R., Shimakura, S., and Gorst, D. (2002), Modeling spatial variation in leukemia survival data, \emph{Journal of the American Statistical Association}, 97(460), 965-972.
}

\examples{
data(LeukSurv)
head(LeukSurv)
}

