
## rEEMSplots

This R package provides the function `eems.plots` to visualize the estimated effective migration and diversity surfaces. It is not on CRAN, so install it from source instead. In the R console:

```
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
if (file.exists("./rEEMSplots/")) {
  install.packages("rEEMSplots", repos=NULL, type="source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}
```

`rEEMSplots` requires a few other packages in order to work. If these are not already installed, an error message would appear.

```
ERROR: dependencies ‘Rcpp’, ‘RcppEigen’, ‘raster’, ‘rgeos’, ‘sp’ are not available for package ‘rEEMSplots’
* removing ‘/Library/Frameworks/R.framework/Versions/3.2/Resources/library/rEEMSplots’
Warning message:
In install.packages("rEEMSplots", repos = NULL, type = "source") :
  installation of package ‘rEEMSplots’ had non-zero exit status
```

`Rcpp` (as well as `rEEMSplots`) needs compilation and thus a compiler, as explained here:  [Advanced R by Hadley Wickham](http://adv-r.had.co.nz/Rcpp.html) and [Rcpp-FAQ](http://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf). `rgeos` requires the [Geometry Engine, Open Source](http://trac.osgeo.org/geos/) (GEOS).

Once the compiler and GEOS are both available, install the packages required by `rEEMSplots`.

```
install.packages(c("Rcpp","RcppEigen","raster","rgeos","sp"))
```

Finally, try to install `rEEMSplots` again:

```
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
if (file.exists("./rEEMSplots/")) {
  install.packages("rEEMSplots", repos=NULL, type="source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}
```

There are some optional packages, most importantly:

* `rgdal` to specify the projection (for example, Mercator) using the `projection.in` and `projection.out` options. Similarly to `rgeos`, `rgdal` requires to install the [Geospatial Data Abstraction Library](http://www.gdal.org) (GDAL).
* `rworldmap` and `rworldxtra` to add a geographic map using the `add.map` option.
