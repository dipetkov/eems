
# EEMS Shiny app

This Shiny app loads an .RData dataset and displays an interactive scatter plot of the observed and fitted dissimilarities between demes. The plot is interactive because you can click and select ("brush") points in it. Then the selected points are highlighted in two visualizations of the population grid.

## Launch the app

Get the latest version of `shiny` on GitHub. Make sure that `ggplot2` and `dplyr` are installed as well.

```r
devtools::install_github("rstudio/shiny")
install.packages(c("ggplot2", "dplyr"))
```

Then download the app directory `dev/rEEMSshiny` in your working directory and launch the app in R or RStudio.

```r
library("shiny")
runApp("rEEMSshiny")
```

Or use the version uploaded to the shinyapps.io server hosted by RStudio.

https://dipetkov.shinyapps.io/rEEMSshiny/

## Brush dissimilarities to highlight demes

* Each point in the dissimilarity scatter plot corresponds to a pair of demes, let's call them alpha and beta. The point indicates the observed and the fitted (average) dissimilarity between two individuals, one from alpha and another from beta.
* Each sampled deme is paired with every other sampled deme. So there are more points than demes in the *between* scatter plot: there are `d choose 2 = d(d-1)/2` points where `d` is the number of demes with at least two sampled individuals.

Since each point in the *between* scatter plot is associated with two demes, alpha and beta, the app displays the population grid twice and it highlights alpha in one plot and beta in the other. (For now) The edges connecting adjacent demes are not shown but the demes are plotted at their (longitude, latitude) locations. So the app should make it easier to identify which demes are associated with interesting patterns in the *between* scatter plot.
