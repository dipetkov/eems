
myheatmap <- function(x, col = NULL, colscale = NULL) {
  if (!is.matrix(x)) x <- as.matrix(x)
  r <- nrow(x)
  c <- ncol(x)
  if (c == 1) {
    x <- matrix(rev(x[, 1]), nrow = 1)
  } else if (r == 1) {
    x <- matrix(rev(x[1, ]), ncol = 1)
  } else {
    x <- t(x)[rev(seq(r)), ]
  }
  if (is.null(colscale)) {
    colscale <- range(x, na.rm = TRUE)
  }
  palette <- mypalette(colors = col, colscale = colscale)
  image(x, col = palette$colors, breaks = palette$levels, axes = FALSE)
  palette
}

heatmap.resid <- function(datapath, mcmcpath) {
  mcmcpath <- mcmcpath
  nchains <- length(mcmcpath)
  message("Heatmap of n-by-n matrix of residuals (observed - fitted):")
  Diffs <- as.matrix(read.table(paste0(datapath, ".diffs"), header = FALSE))
  nIndiv <- nrow(Diffs)
  nChains <- 0
  Delta <- matrix(0, nIndiv, nIndiv)
  for (path in mcmcpath) {
    if (file.exists(file.path(path, "rdistJtDobsJ.txt")) &&
      file.exists(file.path(path, "rdistJtDhatJ.txt")) &&
      file.exists(file.path(path, "rdistoDemes.txt"))) {
      message(path)
      nChains <- nChains + 1
      ipmap <- scan(file.path(path, "ipmap.txt"), what = numeric(), quiet = TRUE)
      Sizes <- as.numeric(table(ipmap))
      n <- length(ipmap) # number of samples
      o <- length(Sizes) # number of observed demes
      J <- Matrix::spMatrix(n, o, i = seq(n), j = ipmap, x = rep(1, n)) # indicator matrix
      J <- as.matrix(J)
      JtDhatJ <- as.matrix(read.table(file.path(path, "rdistJtDhatJ.txt"),
        header = FALSE
      ))
      Delta <- Delta + J %*% JtDhatJ %*% t(J)
    }
  }
  if (nchains == 0) {
    return(NULL)
  }
  Delta <- Delta / nchains
  Delta <- Delta - diag(diag(Delta))
  resid <- Diffs - Delta
  diag(resid) <- NA # The diagonal entries would always be zero
  resid
}

#' A function to plot a heatmap of the residual pairwise dissimilarities, abs(observed - fitted)
#'
#' Given a set of EEMS output directories, this function generates a heatmap of the n-by-n matrix
#' of residual dissimilarities between pairs of individuals. The residuals are the differences
#' between the observed and the fitted dissimilarities.
#' The function also saves the residual matrix to a file called \code{plotpath-eems-resid.RData}.
#' In both the residual matrix and the corresponding heat map, individuals are are in the same
#' order as in the input files \code{datapath.coord} and \code{datapath.diffs}. Applicable only
#' in the case of SNP data when the observed dissimilarity matrix is computed explicitly.
#' @param datapath The full path and the file name of the input Diffs matrix, which is not copied
#' by \code{runeems} to the output directory.
#' @param heatmap.cols The heatmap color palette as a vector of colors, ordered from low to high.
#' Defaults to the "Reds" divergent palette in the RColorBrewer package.
#' @param heatmap.colscale A fixed range for the heatmap colors. The default is NULL, so the color
#' space is the observed range of the residuals.
#' @inheritParams eems.plots
#' @examples
#' # Use the provided example or supply the path to your own EEMS run
#' extdata_path <- system.file("extdata", package = "rEEMSplots")
#' eems_dataset <- file.path(extdata_path, "EEMS-barrier")
#' eems_results <- file.path(extdata_path, "EEMS-barrier")
#' name_figures <- file.path(path.expand("~"), "EEMS-barrier")
#'
#' eems.resid.heatmap(
#'   datapath = eems_dataset,
#'   mcmcpath = eems_results,
#'   plotpath = name_figures,
#'   heatmap.cols = c("gray99", "red")
#' )
#' @seealso \code{\link{eems.plots}, \link{eems.voronoi.samples}, \link{eems.posterior.draws},
#' \link{eems.population.grid}}
#' @export

eems.resid.heatmap <- function(datapath,
                               mcmcpath,
                               plotpath,
                               plot.width = 7,
                               plot.height = 7,
                               out.png = TRUE,
                               res = 600,
                               heatmap.cols = NULL,
                               heatmap.colscale = NULL) {
  load.required.package(
    package = "Matrix",
    required.by = "eems.resid.heatmap"
  )

  save.params <- list(height = plot.height, width = plot.width, res = res, out.png = out.png)

  # A vector of EEMS output directories, for the same dataset.
  # Assume that if eemsrun.txt exists, then all EEMS output files exist.
  mcmcpath <- mcmcpath[file.exists(file.path(mcmcpath, "eemsrun.txt"))]
  if (!length(mcmcpath)) {
    stop("Please provide one existing EEMS output directory, mcmcpath")
  }

  message("Processing the following EEMS output directories :")
  message(mcmcpath)

  eems.resid <- heatmap.resid(datapath, mcmcpath)
  save(eems.resid, file = paste0(plotpath, "-eems-resid.RData"))
  save.graphics(paste0(plotpath, "-eems-resid"), save.params)
  par(las = 1, font.main = 1, mar = c(0, 0, 0, 0) + 0.1)
  key <- myheatmap(abs(eems.resid), col = heatmap.cols, colscale = heatmap.colscale)
  dev.off()

  if (out.png) {
    save.params$height <- 6
    save.params$width <- 1.5
  } else {
    save.params$height <- 12
    save.params$width <- 3
  }

  save.graphics(paste0(plotpath, "-eems-resid-key"), save.params)
  par(las = 1, font.main = 1, mar = c(1, 1, 5, 8))
  myfilled.legend(
    levels = key$levels, col = key$colors,
    key.axes = axis(4, tick = FALSE, hadj = 1, line = 4, cex.axis = 2),
    key.title = mtext("abs(r)", side = 3, cex = 2.5, line = 1.5, font = 1)
  )
  dev.off()
}
