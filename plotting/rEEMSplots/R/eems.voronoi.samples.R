
# This function is mainly for testing purposes and will create on Voronoi diagram
# for each saved MCMC iteration
voronoi.diagram <- function(mcmcpath, dimns, longlat, plot.params,
                            post.draws = 1, is.mrates = TRUE) {
  message("Plotting Voronoi tessellation of estimated effective rates")
  if (is.mrates) {
    log_scale <- plot.params$m.log_scale
  } else {
    log_scale <- plot.params$q.log_scale
  }
  voronoi <- read.voronoi(mcmcpath, longlat, is.mrates, log_scale)
  rates <- voronoi$rates
  tiles <- voronoi$tiles
  xseed <- voronoi$xseed
  yseed <- voronoi$yseed

  # Choose one saved posterior draw at random
  niters <- length(tiles)
  riter <- sample(seq(niters), 1)
  message(mcmcpath)
  message("Draw ", riter, " (out of ", niters, ")")

  eems.colors <- plot.params$eems.colors
  num.levels <- length(eems.colors)
  if (is.mrates) {
    eems.levels <- eems.colscale(rates, num.levels, plot.params$m.colscale)
  } else {
    eems.levels <- eems.colscale(rates, num.levels, plot.params$q.colscale)
  }
  n.levels <- length(eems.levels)
  max.levels <- max(eems.levels)
  min.levels <- min(eems.levels)
  # Jump over stored parameters for draws 1 to (riter - 1)
  skip <- sum(tiles[riter:1][-1])
  now.tiles <- tiles[riter]
  now.rates <- rates[(skip + 1):(skip + now.tiles)]
  now.xseed <- xseed[(skip + 1):(skip + now.tiles)]
  now.yseed <- yseed[(skip + 1):(skip + now.tiles)]
  # Standardize the log-transformed rates, without taking into account
  # the relative size of the tiles (this is hard to do without a grid)
  now.rates <- now.rates - mean(now.rates)
  now.rates <- ifelse(now.rates > max.levels, max.levels, now.rates)
  now.rates <- ifelse(now.rates < min.levels, min.levels, now.rates)
  par(mar = c(0, 0, 0, 0) + 0.1)
  plot(0, 0,
    type = "n", xlim = dimns$xlim, ylim = dimns$ylim, asp = 1,
    axes = FALSE, xlab = "", ylab = "", main = ""
  )
  if (now.tiles == 1) {
    # There is only one tile
    tile.color <- eems.colors[round(n.levels / 2)]
    polygon(dimns$xlim, dimns$ylim, col = tile.color, border = "white")
  } else {
    # Plot each tile in turn (as a polygon)
    Voronoi <- deldir::deldir(now.xseed, now.yseed, rw = c(dimns$xlim, dimns$ylim))
    tilelist <- deldir::tile.list(Voronoi)
    for (c in 1:now.tiles) {
      tile.color <- eems.colors[findInterval(now.rates[c], eems.levels, all.inside = TRUE)]
      polygon(tilelist[[c]]$x, tilelist[[c]]$y, col = tile.color, border = "white")
    }
    filled.contour.graph(mcmcpath, longlat, plot.params)
  }
  if (plot.params$add.seeds) {
    points(now.xseed, now.yseed,
      pch = plot.params$pch.seeds, cex = plot.params$cex.seeds,
      col = plot.params$col.seeds, lwd = 3
    )
  }
  list(colors = eems.colors, levels = eems.levels)
}

#' A function to plot Voronoi diagrams of effective migration and diversity rates
#'
#' Given a set of EEMS output directories, this function takes random draws from the posterior
#' distribution of the migration and diversity rate parameters. Each draw is visualized as two
#' Voronoi diagrams; the migration diagram is saved to a file ending in \code{mvoronoiXX},
#' the diversity diagram is saved to a file ending in \code{qvoronoiXX} where \code{XX} is
#' a numeric id. Specify the number of times to draw from the posterior with the argument
#' \code{post.draws}. If \code{post.draws = 10}, then \code{eems.voronoi.samples} will generate
#' plots with id \code{XX = 1} to \code{XX = 10}.
#'
#' Note about the implementation: \code{eems.voronoi.samples} samples randomly from the posterior
#' draws saved during the execution of EEMS, after burn-in and thinning.
#' @param post.draws Number of times to sample from the posterior. The default is 1.
#' @param add.seeds A logical value indicating whether to add the Voronoi seeds or not.
#' @param col.seeds The color of the Voronoi seeds. Defaults to \code{green}.
#' @param pch.seeds The symbol, specified as an integer, or the character to be used for
#' plotting the Voronoi seeds. Defaults to 4.
#' @param cex.seeds The size of the symbol/character used for plotting the Voronoi seeds.
#' Defaults to 1.
#' @param cex.demes The size of the symbol/character used for plotting observed demes.
#' Defaults to 1.
#' @inheritParams eems.plots
#' @references Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for
#' Improved Data Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
#' @examples
#' # Use the provided example or supply the path to your own EEMS run.
#' extdata_path <- system.file("extdata", package = "rEEMSplots")
#' eems_results <- file.path(extdata_path, "EEMS-example")
#' name_figures <- file.path(path.expand("~"), "EEMS-example")
#'
#' library("deldir")
#'
#' # Plot a series of Voronoi diagrams for the EEMS model parameters:
#' # the effective migration rates (m) and the effective diversity rates (q).
#' eems.voronoi.samples(
#'   mcmcpath = eems_results,
#'   plotpath = paste0(name_figures, "-voronoi-diagrams"),
#'   longlat = TRUE, post.draws = 10
#' )
#' @seealso \code{\link{eems.plots}, \link{eems.posterior.draws}, \link{eems.resid.heatmap},
#' \link{eems.population.grid}}
#' @export

eems.voronoi.samples <- function(mcmcpath,
                                 plotpath,
                                 longlat,
                                 post.draws = 1,
                                 plot.width = 7,
                                 plot.height = 7,
                                 out.png = TRUE,
                                 res = 600,
                                 add.grid = FALSE,
                                 col.grid = "gray80",
                                 lwd.grid = 1,
                                 add.outline = TRUE,
                                 col.outline = "gray80",
                                 lwd.outline = 2,
                                 add.demes = FALSE,
                                 col.demes = "gray80",
                                 pch.demes = 19,
                                 cex.demes = 1,
                                 add.seeds = TRUE,
                                 col.seeds = "#8AE234",
                                 pch.seeds = 4,
                                 cex.seeds = 1,
                                 eems.colors = NULL,
                                 m.colscale = NULL,
                                 q.colscale = NULL,
                                 add.title = FALSE) {
  load.required.package(package = "deldir", required.by = "eems.voronoi.samples")

  plot.params <- list(
    eems.colors = eems.colors, m.colscale = m.colscale, q.colscale = q.colscale,
    add.grid = add.grid, add.outline = add.outline, add.demes = add.demes, add.seeds = add.seeds,
    col.grid = col.grid, col.outline = col.outline, col.demes = col.demes, col.seeds = col.seeds,
    lwd.grid = lwd.grid, lwd.outline = lwd.outline, pch.demes = pch.demes, pch.seeds = pch.seeds,
    cex.seeds = cex.seeds, min.cex.demes = cex.demes, max.cex.demes = cex.demes,
    add.title = add.title
  )
  plot.params <- check.plot.params(plot.params)

  # A vector of EEMS output directories, for the same dataset.
  # Assume that if eemsrun.txt exists, then all EEMS output files exist.
  mcmcpath <- mcmcpath[file.exists(file.path(mcmcpath, "eemsrun.txt"))]
  if (!length(mcmcpath)) {
    stop("Please provide at least one existing EEMS output directory, mcmcpath")
  }

  dimns <- read.dimns(mcmcpath, longlat)
  save.params <- list(height = plot.height, width = plot.width, res = res, out.png = out.png)

  message("Processing the following EEMS output directory :")
  message(mcmcpath)

  plot.params$add.grid <- add.grid
  plot.params$all.demes <- FALSE
  plot.params$add.demes <- FALSE

  save.graphics(paste0(plotpath, "-mvoronoi"), save.params)
  par(las = 1, font.main = 1)
  for (draw in seq(post.draws)) {
    # Choose one output directory at random
    voronoi.diagram(sample(mcmcpath, 1),
      dimns, longlat, plot.params,
      post.draws = post.draws, is.mrates = TRUE
    )
  }
  dev.off()

  plot.params$add.grid <- FALSE
  plot.params$all.demes <- add.demes

  save.graphics(paste0(plotpath, "-qvoronoi"), save.params)
  par(las = 1, font.main = 1)
  for (draw in seq(post.draws)) {
    # Choose one output directory at random
    voronoi.diagram(sample(mcmcpath, 1),
      dimns, longlat, plot.params,
      post.draws = post.draws, is.mrates = FALSE
    )
  }
  dev.off()
}
