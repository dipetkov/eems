
#' A function to plot the constructed population grid and optionally the initial sampling locations
#'
#' Given an EEMS output directory, this function generates one figure to visualize the EEMS
#' population grid. All edges are shown in the same color to visualize the grid before estimating
#' migration and diversity rates. This can be helpful if EEMS exits with the error message
#' "The population grid is not connected".
#' @param add.coord A logical value indicating whether to add the original sampling locations to
#' the plot or not.
#' @param col.coord The color of the sampling locations. Defaults to \code{red}.
#' @param pch.coord The symbol, specified as an integer, or the character to be used for plotting
#' the sampling locations. Defaults to 3.
#' @param datapath The full path and the file name of the input dataset (the three files
#' datapath.coord, datapath.diffs, datapath.outer). Must be specified if \code{add_coord = TRUE}.
#' @inheritParams eems.plots
#' @examples
#' # Use the provided example or supply the path to your own EEMS run.
#' extdata_path <- system.file("extdata", package = "rEEMSplots")
#' eems_results <- file.path(extdata_path, "EEMS-example")
#' name_figures <- file.path(path.expand("~"), "EEMS-grid_connected")
#'
#' eems.population.grid(eems_results,
#'   name_figures,
#'   longlat = TRUE,
#'   add.outline = TRUE, col.outline = "purple", lwd.outline = 3,
#'   add.grid = TRUE, col.grid = "green", lwd.grid = 2
#' )
#'
#' # It is more interesting to see an example where the grid is unconnected
#' # due to the unusual shape of the habitat.
#' eems_results <- file.path(extdata_path, "EEMS-popgrid")
#' name_figures <- path.expand(file.path("~", "EEMS-grid_not_connected"))
#'
#' eems.population.grid(eems_results,
#'   name_figures,
#'   longlat = FALSE,
#'   add.outline = TRUE, col.outline = "purple", lwd.outline = 3,
#'   add.grid = TRUE, col.grid = "green", lwd.grid = 2
#' )
#' @seealso \code{\link{eems.plots}, \link{eems.voronoi.samples}, \link{eems.posterior.draws},
#' \link{eems.resid.heatmap}}
#' @export

eems.population.grid <- function(mcmcpath,
                                 plotpath,
                                 longlat,
                                 plot.width = 7,
                                 plot.height = 7,
                                 out.png = TRUE,
                                 res = 600,
                                 add.grid = TRUE,
                                 col.grid = "gray80",
                                 lwd.grid = 1,
                                 add.demes = FALSE,
                                 col.demes = "black",
                                 pch.demes = 19,
                                 min.cex.demes = 1,
                                 max.cex.demes = 3,
                                 add.outline = TRUE,
                                 col.outline = "gray90",
                                 lwd.outline = 2,
                                 add.coord = FALSE,
                                 col.coord = "red",
                                 pch.coord = 3,
                                 datapath = NULL) {
  plot.params <- list(
    add.grid = add.grid, add.outline = add.outline, add.demes = add.demes,
    col.grid = col.grid, col.outline = col.outline, col.demes = col.demes,
    lwd.grid = lwd.grid, lwd.outline = lwd.outline, pch.demes = pch.demes,
    min.cex.demes = min.cex.demes, max.cex.demes = max.cex.demes
  )
  plot.params <- check.plot.params(plot.params)

  # A vector of EEMS output directories, for the same dataset.
  # Assume that if eemsrun.txt exists, then all EEMS output files exist.
  mcmcpath <- mcmcpath[file.exists(file.path(mcmcpath, "eemsrun.txt"))]
  if (!length(mcmcpath)) {
    stop("Please provide at least one existing EEMS output directory, mcmcpath")
  }

  save.params <- list(height = plot.height, width = plot.width, res = res, out.png = out.png)

  message("Processing the following EEMS output directories :")
  message(mcmcpath)

  save.graphics(paste0(plotpath, "-popgrid"), save.params)

  # Read the habitat outline from the mcmcpath (output) directory
  graph <- read.graph(mcmcpath, longlat)
  outer <- graph$outer

  # Read the sampling locations from the datapath (input) directory
  coordpath <- paste0(datapath, ".coord")
  coord <- NULL

  if (!is.null(datapath) && file.exists(coordpath)) {
    message(paste0("Read the original sampling locations from ", coordpath))
    coord <- scan(coordpath, what = numeric(), quiet = TRUE)
    coord <- matrix(coord, ncol = 2, byrow = TRUE)
    if (!longlat) {
      coord <- coord[, c(2, 1)]
    }
  }

  # Make the plot canvas large enough to fit both the sampling locations and the habitat outline
  plot(rbind(outer, coord), xlab = "", ylab = "", axes = FALSE, type = "n", asp = 1)

  filled.contour.outline(mcmcpath, longlat, plot.params)
  filled.contour.graph(mcmcpath, longlat, plot.params)

  # Plot the sampling locations
  if (add.coord) {
    if (is.null(coord)) {
      message(
        "\nTo add the sampling locations with option add.coord = TRUE,\n",
        "Specify the full path to the input dataset (coord/diffs/outer files)\n\n"
      )
    } else {
      points(coord, pch = pch.coord, col = col.coord)
    }
  }
  dev.off()
}
