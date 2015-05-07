
## rEEMSplots

This R package provides the function `eems.plots`. It is not on CRAN, so install it from source instead. In the R console:

```
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
if (file.exists("./rEEMSplots/")) {
  install.packages("rEEMSplots", repos=NULL, type="source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}
```
