
#' @useDynLib rEEMSplots
#' @import Rcpp RcppEigen
#' @import raster rgeos sp
#' @import graphics grDevices
#' @importFrom stats lm na.omit
#' @importFrom utils read.table write.table
NULL

load.required.package <- function(package, required.by) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(paste0(
      "'", required.by, "' requires the '", package, "' package. ",
      "Please install '", package, "' first."
    ))
  } else {
    message(paste0("Loading ", package, " (required by ", required.by, ")"))
  }
}

# By default, all figures are saved as bitmap PNG images. However,
# it is straightforward to use another format (Here the alternative is PDF)
save.graphics <- function(plotpath, plot.params, ...) {
  plotpath <- path.expand(plotpath)
  if (plot.params$out.png) {
    bitmap(paste0(plotpath, "%02d.png"),
      type = "png16m",
      res = plot.params$res, units = "in",
      height = plot.params$height,
      width = plot.params$width, ...
    )
  } else {
    pdf(paste0(plotpath, "%02d.pdf"),
      height = plot.params$height,
      width = plot.params$width,
      onefile = FALSE, ...
    )
  }
}
