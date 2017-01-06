\name{rEEMSplots-package}
\alias{rEEMSplots-package}
\alias{rEEMSplots}
\docType{package}
\title{
EEMS graphics output
}
\description{
    Provides several functions to visualize the results of EEMS.
    
    \code{eems.plots} produces contour plots of the effective migration and diversity rates,
    as well as scatterplots of observed vs fitted genetic dissimilarities.
    
    \code{eems.population.grid} plots the population grid, with all edges in the same color,
    to visualize the grid before estimating migration and diversity rates.
    
    \code{eems.voronoi.samples} and \code{eems.posterior.draws} both take random samples
    from the posterior distribution of the migration and diversity rates and plot
    them to visualize the posterior variance (in a slightly different way).
    
    \code{eems.resid.heatmap} computes the matrix of residuals (differences between the
    observed and fitted dissimilarities) and plots it as a heat map.
}
\details{
\tabular{ll}{
Package: \tab rEEMSplots\cr
Version: \tab 0.0.1.0000\cr
Date: \tab 2016-04-24\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Desislava Petkova <desislavka@gmail.com>
}
\references{
D. Petkova, J. Novembre, M. Stephens. Visualizing spatial population
structure with estimated effective migration surfaces. Nature Genetics
48, 94 -- 100 (2016).
}
\keyword{ package }
