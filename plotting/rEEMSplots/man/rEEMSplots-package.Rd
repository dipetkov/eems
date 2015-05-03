\name{rEEMSplots-package}
\alias{rEEMSplots-package}
\alias{rEEMSplots}
\docType{package}
\title{
EEMS graphics output
}
\description{
  Provides one function, eems.plots, to visualize the results of EEMS.
  It produces contour plots of the effective migration and diversity
  rates, as well as scatterplots of observed vs fitted genetic
  dissimilarities.
}
\details{
\tabular{ll}{
Package: \tab rEEMSplots\cr
Version: \tab 0.0.0.9000\cr
Date: \tab 2015-03-04\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Desislava Petkova <desislavka@gmail.com>
}
\references{
D. Petkova, J. Novembre, M. Stephens. Visualizing spatial population structure with estimated effective migration surfaces. biorXiv DOI: 10.1101/011809 (2014).
}
\keyword{ package }
\examples{
eems.plots(mcmcpath,plotpath,longlat,
           add.map=FALSE,
           add.grid=TRUE,
           add.demes=TRUE)
}
