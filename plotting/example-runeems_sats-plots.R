

## Check that the 'rEEMSplots' package, which provides the
## eems.plots function, is installed
available = suppressMessages(suppressWarnings(require('rEEMSplots',
    quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE)))
if (!available) {
    stop(paste("First install the 'rEEMSplots' package\n",
               "install.packages('rEEMSplots', repos = NULL, type='source')"))
}


library(rEEMSplots)


## mcmcpath is a list of three output directories; the results will be averaged
mcmcpath <- paste('../runeems_sats/data/sat-barrier-nIndiv150-nSites16-EEMS-nDemes200-simno',1:3,sep='')
plotpath <- '../runeems_sats/plot/sat-barrier-nIndiv150-nSites16-EEMS-nDemes200-simno1_3'
longlat <- TRUE


## eems.plots takes three required arguments:
##   mcmcpath: one or several output directories (for the same dataset)
##   plotpath: filename of figures to generate
##     The following figures are created by default:
##     * plotpath-mrates01.png: effective migration rates
##     * plotpath-qrates01.png: effective diversity rates
##     * plotpath-rist01/02.png: fitted vs observed distances
##     * plotpath-pilogl01.png: trace of posterior probability
##   longlat (TRUE or FALSE): are the coordinates ordered longitude/latitude or not?
##                          longlat = FALSE basically "transposes" the x and y axes
##
## eems.plots takes a variety of optional arguments:
##   plot.width and plot.height: width and height of the graphics region (in inches)
##   plot.voronoi (TRUE or FALSE): plot a few posterior Voronoi diagrams?
##   add.map: add 'worldHires' map (using the mapdata package)?
##   add.grid: add triangular population grid?
##   add.samples: add samples to their assigned location in the grid?
##   add.outline: add the habitat ring (as declared in the .outer file)?
##   col.map/col.grid/col.samples/col.outline: specify the colors
##   lwd.max/lwd.grid/lwd.outline: specify the line width
##   pch.samples: specify the character
##   cex.samples: specify the character size
##   max.cex.samples: some demes might be assigned more samples than others.
##     If max.cex.samples>0, then demes with more samples will also have bigger size.
##     If the sampling is uneven, then max.cex.samples>0 will underline this fact.


## The options projection.in/projection.out allow to specify the input/output projection.
## For example, use projection.in = "+proj=merc" for mercator
##               or projection.in = "+proj=longlat +datum=WGS84" if the coordinates are longitude and latitude
## Or use another projection by specifying the corresponding PROJ.4 string.


eems.plots(mcmcpath,plotpath,longlat,
           add.map=FALSE,
           add.grid=TRUE,
           add.samples=TRUE)

