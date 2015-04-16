\name{rEEMSplots-package}
\alias{rEEMSplots-package}
\alias{rEEMSplots}
\docType{package}
\title{
Generate figures from EEMS output
}
\description{
Provides one function, eems.plots, to visualize EEMS results.
}
\details{
\tabular{ll}{
Package: \tab rEEMSplots\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-04\cr
License: \tab GPL-3\cr
}
}
\author{
Desislava Petkova <desislavka@gmail.com>
}
\references{
D. Petkova, J. Novembre, M. Stephens. Visualizing spatial population structure with estimated effective migration surfaces. biorXiv DOI: 10.1101/011809 (2014).
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\examples{
eems.plots(mcmcpath,plotpath,longlat,
           add.map=FALSE,
           add.grid=TRUE,
           add.samples=TRUE)
}
