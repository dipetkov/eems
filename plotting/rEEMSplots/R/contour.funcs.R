
# Redefine the 'filledContour' function from the 'raster' package
# since I don't like how the legend looks.
# For now, I have only changed how the legend/bar is plotted,
# the change itself is in the 'myfilled.contour' function
myfilledContour <- function(x, y = 1, maxpixels = 1e+05, add.key = TRUE, ...) {
  if (nlayers(x) > 1) {
    y <- min(max(1, y), nlayers(x))
    x <- raster(x, y)
  }
  x <- sampleRegular(x, maxpixels, asRaster = TRUE, useGDAL = TRUE)
  seq_cols <- seq_len(ncol(x))
  seq_rows <- rev(seq_len(nrow(x)))
  X <- xFromCol(x, seq_cols)
  Y <- yFromRow(x, seq_rows)
  Z <- t(matrix(getValues(x), ncol = x@ncols, byrow = TRUE)[seq_rows, ])
  myfilled.contour(
    x = X, y = Y, z = Z, add.key = add.key,
    axes = FALSE, frame.plot = FALSE, ...
  )
}

# I have only changed the legend/bar to remove the black border.
myfilled.contour <- function(x = seq(0, 1, length.out = nrow(z)),
                             y = seq(0, 1, length.out = ncol(z)), z,
                             xlim = range(x, finite = TRUE),
                             ylim = range(y, finite = TRUE),
                             zlim = range(z, finite = TRUE),
                             levels = pretty(zlim, nlevels), nlevels = 20,
                             color.palette = cm.colors, col = color.palette(length(levels) - 1),
                             plot.title, plot.axes, key.title, key.axes, asp = NA,
                             xaxs = "i", yaxs = "i", las = 1,
                             axes = TRUE, frame.plot = axes, add.title = FALSE, add.key = TRUE,
                             ...) {
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    } else {
      stop("no 'z' matrix specified")
    }
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
    stop("increasing 'x' and 'y' values expected")
  }
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  if (add.key) {
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  }
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  if (add.key) {
    plot.new()
    par(plt = c(0, 0.3, 0.1, 0.9))
    # Change: add 'border = NA' and remove 'box()'
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
    if (missing(key.axes)) {
      if (axes) {
        axis(4)
      }
    } else {
      key.axes
    }
    if (!missing(key.title)) {
      key.title
    }
  }
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  # Change: make the plot region almost as big as the figure region
  #         but leave space for the main title
  if (add.title) {
    par(plt = c(0.02, 0.98, 0.02, 0.95))
  } else {
    par(plt = c(0.02, 0.98, 0.02, 0.98))
  }
  plot.window(xlim, ylim, "", asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  } else {
    plot.axes
  }
  if (frame.plot) {
    box()
  }
  if (add.title) {
    if (missing(plot.title)) {
      title(...)
    } else {
      plot.title
    }
  }
  invisible()
}

myfilled.legend <-
  function(levels, color.palette = cm.colors,
           col = color.palette(length(levels) - 1), plot.title, plot.axes,
           key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
           axes = TRUE, frame.plot = axes, ...) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
    if (!missing(key.axes)) {
      key.axes
    } else {
      axis(4, tick = FALSE)
    }
    if (!missing(key.title)) {
      key.title
    }
  }

filled.contour.points <- function(mcmcpath, longlat, plot.params, highlight) {
  if (is.null(highlight)) {
    return(NULL)
  }
  if (is.null(highlight$index)) {
    message("Specify the indices of points to highlight")
    return(NULL)
  }
  if (is.null(highlight$col)) highlight$col <- "red"
  if (is.null(highlight$cex)) highlight$cex <- 1
  if (is.null(highlight$pch)) highlight$pch <- 4
  if (is.null(highlight$lwd)) highlight$lwd <- 2

  graph <- read.graph(mcmcpath, longlat)
  coord <- graph$demes[graph$ipmap, ]

  index2coord <- match(highlight$index, seq_len(nrow(coord)))
  index2highlight <- which(!is.na(index2coord))
  index2coord <- index2coord[index2highlight]

  # Check that there is at least one point to highlight
  if (!length(index2highlight)) {
    return(NULL)
  }

  coord <- coord[index2coord, ]
  coord <- matrix(coord, ncol = 2)
  # Transform the coordinates to a different projection if necessary
  if (!is.null(plot.params$proj.out)) {
    coord <- sp::SpatialPoints(coord, proj4string = CRS(plot.params$proj.in))
    coord <- sp::spTransform(coord, CRSobj = CRS(plot.params$proj.out))
  }
  points(coord,
    cex = highlight$cex[index2highlight],
    col = highlight$col[index2highlight],
    pch = highlight$pch[index2highlight],
    lwd = highlight$lwd[index2highlight]
  )
}

transform.rates <- function(dimns, tiles, rates, xseed, yseed, zero_mean) {
  # Bind together the points which make up the raster image
  # and the extra points provided by the user, if any
  marks <- rbind(dimns$marks, dimns$coord)
  coord_raster <- nrow(dimns$marks)
  coord_extra <- nrow(dimns$coord)

  # Compute the rates at all the points
  if (zero_mean) {
    rslts <- tiles2contours_standardize(
      tiles, rates, cbind(xseed, yseed),
      marks, dimns$dist.metric
    )
  } else {
    rslts <- tiles2contours(
      tiles, rates, cbind(xseed, yseed),
      marks, dimns$dist.metric
    )
  }

  Zvals <- matrix(rslts$zvals[seq_len(coord_raster)],
    nrow = dimns$nxmrks, ncol = dimns$nymrks
  )
  if (coord_extra > 0) {
    Zvals_extra <- rslts$zvals[seq_len(coord_extra) + coord_raster]
  } else {
    Zvals_extra <- NULL
  }
  if (zero_mean) {
    # Without normalizing the migration rates m to have mean 0,
    # it doesn't make sense to keep track of the number of times (m > 0) and (m < 0).
    PrGT0 <- matrix(rslts$prgt0[seq_len(coord_raster)],
      nrow = dimns$nxmrks, ncol = dimns$nymrks
    )
    PrLT0 <- matrix(rslts$prlt0[seq_len(coord_raster)],
      nrow = dimns$nxmrks, ncol = dimns$nymrks
    )
  } else {
    PrGT0 <- NULL
    PrLT0 <- NULL
  }
  list(
    Zvals = Zvals, PrGT0 = PrGT0, PrLT0 = PrLT0,
    Zvals_extra = Zvals_extra, niters = length(tiles)
  )
}

check.habitat.outer.is.valid <- function(mcmcpath, longlat) {
  graph <- read.graph(mcmcpath, longlat)
  outer <- graph$outer
  # Check that the geometry of the habitat outline is both Simple (gIsSimple)
  # and Closed (end points overlap).
  x <- outer[, 1]
  y <- outer[, 2]
  l <- readWKT(paste("LINESTRING(", paste(paste(x, y), collapse = ", "), ")"))
  if (!rgeos::gIsRing(l, byid = FALSE)) {
    message(
      "Error (rgeos):\n",
      "The habitat geometry is not a valid ring (a ring is both simple and closed)\n",
      "Let the habitat be the rectangle defined by (xmin, ymin) and (xmax, yxmax)\n",
      "where xmin, xmax = range(longitude) and ymin, ymax = range(latitude)."
    )
    xmin <- min(outer[, 1])
    xmax <- max(outer[, 1])
    ymin <- min(outer[, 2])
    ymax <- max(outer[, 2])
    outer <- cbind(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin))
  }
  outer
}

filled.contour.map <- function(mcmcpath, longlat, plot.params) {
  if (!is.null(plot.params$proj.in) && plot.params$add.map) {
    map <- rworldmap::getMap(resolution = "high")
    map <- sp::spTransform(map, CRSobj = CRS(plot.params$proj.out))
    plot(map, col = NA, border = plot.params$col.map, lwd = plot.params$lwd.map, add = TRUE)
  }
}

filled.contour.outline <- function(mcmcpath, longlat, plot.params) {
  outer <- check.habitat.outer.is.valid(mcmcpath, longlat)
  boundary <- sp::SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)), "1")))
  if (!is.null(plot.params$proj.in)) {
    proj4string(boundary) <- CRS(plot.params$proj.in)
    boundary <- sp::spTransform(boundary, CRSobj = CRS(plot.params$proj.out))
  }
  # The filledContour fills a rectangular plot; now color the habitat exterior white.
  exterior <- rgeos::gDifference(rgeos::gEnvelope(boundary), boundary)
  if (!is.null(exterior)) {
    plot(exterior, col = "white", border = "white", add = TRUE)
  }
  if (plot.params$add.outline) {
    plot(boundary,
      col = NA, border = plot.params$col.outline,
      lwd = plot.params$lwd.outline, add = TRUE
    )
  }
}

filled.contour.graph <- function(mcmcpath, longlat, plot.params) {
  graph <- read.graph(mcmcpath, longlat)
  if (plot.params$add.grid) {
    segments <- list()
    for (e in seq_len(nrow(graph$edges))) {
      segments[[e]] <- sp::Line(graph$demes[graph$edges[e, ], ])
    }
    segments <- sp::SpatialLines(list(Lines(segments, ID = "a")))
    if (!is.null(plot.params$proj.in)) {
      proj4string(segments) <- CRS(plot.params$proj.in)
      segments <- sp::spTransform(segments, CRSobj = CRS(plot.params$proj.out))
    }
    lines(segments, col = plot.params$col.grid, lwd = plot.params$lwd.grid)
  }
  if (plot.params$all.demes) {
    all.demes <- sp::SpatialPoints(matrix(graph$demes[graph$alpha, ],
      ncol = 2
    ))
    if (!is.null(plot.params$proj.in)) {
      proj4string(all.demes) <- CRS(plot.params$proj.in)
      all.demes <- sp::spTransform(all.demes, CRSobj = CRS(plot.params$proj.out))
    }
    points(all.demes,
      col = plot.params$col.demes,
      pch = plot.params$pch.demes, cex = plot.params$min.cex.demes
    )
  } else if (plot.params$add.demes) {
    # Make sure that observed.demes is a matrix even if there is one alpha
    observed.demes <- sp::SpatialPoints(matrix(graph$demes[graph$alpha, ],
      ncol = 2
    ))
    if (!is.null(plot.params$proj.in)) {
      proj4string(observed.demes) <- CRS(plot.params$proj.in)
      observed.demes <- sp::spTransform(observed.demes, CRSobj = CRS(plot.params$proj.out))
    }
    cex.points <- plot.params$min.cex.demes
    if (min(graph$sizes) < max(graph$sizes)) {
      cex.points <- cex.points +
        (plot.params$max.cex.demes - plot.params$min.cex.demes) *
          (graph$sizes - min(graph$sizes)) / (max(graph$sizes) - min(graph$sizes))
    }
    points(observed.demes,
      col = plot.params$col.demes,
      pch = plot.params$pch.demes, cex = cex.points
    )
  }
}

project_raster <- function(surface, proj.in, proj.out) {
  raster::projection(surface) <- CRS(proj.in)
  raster::projectRaster(surface, crs = CRS(proj.out))
}

points_to_raster <- function(points, dimns, params) {
  surface <- raster::raster(t(points),
    xmn = dimns$xlim[1], xmx = dimns$xlim[2],
    ymn = dimns$ylim[1], ymx = dimns$ylim[2]
  )
  surface <- raster::flip(surface, direction = "y")
  if (!is.null(params$proj.in)) {
    surface <- project_raster(surface, params$proj.in, params$proj.out)
  }
  surface
}

null.eems.contour <- function(mcmcpath, longlat, plot.params,
                              plot.xy = NULL, highlight.samples = NULL) {
  dimns <- read.dimns(mcmcpath, longlat)
  eems.colors <- plot.params$eems.colors
  num.levels <- length(eems.colors)
  eems.levels <- seq(from = -2.5, to = +2.5, length = num.levels + 1)
  main.title <- ""
  key.title <- ""
  Zeros <- matrix(0, nrow = dimns$nxmrks, ncol = dimns$nymrks)
  myfilledContour(
    points_to_raster(Zeros, dimns, plot.params),
    col = eems.colors, levels = eems.levels, asp = 1,
    add.key = plot.params$add.colbar,
    key.axes = axis(4, tick = FALSE, hadj = 1, line = 3, cex.axis = 1.5),
    key.title = mtext(key.title, side = 3, cex = 1.5, line = 1.5, font = 1),
    add.title = plot.params$add.title,
    plot.title = mtext(text = main.title, side = 3, line = 0, cex = 1.5),
    plot.axes = {
      filled.contour.outline(mcmcpath, longlat, plot.params)
      filled.contour.map(mcmcpath, longlat, plot.params)
      plot.xy
      filled.contour.graph(mcmcpath, longlat, plot.params)
      filled.contour.points(mcmcpath, longlat, plot.params, highlight.samples)
    }
  )
  list(eems.colors = eems.colors, eems.levels = eems.levels)
}

one.eems.contour <- function(mcmcpath, dimns, Zmean, longlat, plot.params, is.mrates,
                             plot.xy = NULL, highlight.samples = NULL) {
  eems.colors <- plot.params$eems.colors
  num.levels <- length(eems.colors)
  if (is.mrates) {
    eems.levels <- eems.colscale(Zmean, num.levels, plot.params$m.colscale)
    main.title <- "Posterior mean migration rates m (on the log10 scale)"
    key.title <- "log(m)"
  } else {
    eems.levels <- eems.colscale(Zmean, num.levels, plot.params$q.colscale)
    main.title <- "Posterior mean diversity rates q (on the log10 scale)"
    key.title <- "log(q)"
  }
  myfilledContour(
    points_to_raster(Zmean, dimns, plot.params),
    col = eems.colors, levels = eems.levels, asp = 1,
    add.key = plot.params$add.colbar,
    key.axes = axis(4, tick = FALSE, hadj = 1, line = 3, cex.axis = 1.5),
    key.title = mtext(key.title, side = 3, cex = 1.5, line = 1.5, font = 1),
    add.title = plot.params$add.title,
    plot.title = mtext(text = main.title, side = 3, line = 0, cex = 1.5),
    plot.axes = {
      filled.contour.outline(mcmcpath, longlat, plot.params)
      filled.contour.map(mcmcpath, longlat, plot.params)
      plot.xy
      filled.contour.graph(mcmcpath, longlat, plot.params)
      filled.contour.points(mcmcpath, longlat, plot.params, highlight.samples)
    }
  )
  list(eems.colors = eems.colors, eems.levels = eems.levels)
}

one.prob.contour <- function(mcmcpath, dimns, Props, longlat, plot.params, is.mrates,
                             plot.xy = NULL, highlight.samples = NULL) {
  Props <- (Props + 1) / 2
  Props[Props < 0] <- 0
  Props[Props > 1] <- 1

  if (is.mrates) v <- "m" else v <- "q"
  main.title <- paste0("Posterior probabilities P{log(", v, ") <> 0 | diffs}")
  key.greaterthan0 <- paste0("P{log(", v, ") > 0}")
  key.lessthan0 <- paste0("P{log(", v, ") < 0}")

  # Probabilities are lie in the range [0, 1] but I can experiment with different ways
  # to split the range [0, 1] into bins.
  alpha <- plot.params$prob.levels
  alpha <- alpha[alpha > 0.5 & alpha < 1]
  alpha <- sort(unique(alpha))
  prob.labels <- c("", as.character(rev(alpha)), as.character(alpha), "")
  prob.levels <- c(0, 1 - rev(alpha), alpha, 1)
  prob.colors <- default.eems.colors()
  prob.colors <- colorRampPalette(prob.colors)(length(prob.levels) - 1)

  myfilledContour(
    points_to_raster(Props, dimns, plot.params),
    col = prob.colors, levels = prob.levels, asp = 1,
    add.key = plot.params$add.colbar,
    key.axes = axis(4,
      at = prob.levels, labels = prob.labels,
      tick = FALSE, hadj = 0, line = 1, cex.axis = 1.5
    ),
    key.title = {
      mtext(key.greaterthan0, side = 3, cex = 1.5, line = 1.5, font = 1)
      mtext(key.lessthan0, side = 1, cex = 1.5, line = 1.5, font = 1)
    },
    add.title = plot.params$add.title,
    plot.title = mtext(text = main.title, side = 3, line = 0, cex = 1.5),
    plot.axes = {
      filled.contour.outline(mcmcpath, longlat, plot.params)
      filled.contour.map(mcmcpath, longlat, plot.params)
      # Including `plot.xy` has no effect because it has already been
      # evaluated once to add extra points to the rates plot
      plot.xy
      filled.contour.graph(mcmcpath, longlat, plot.params)
      filled.contour.points(mcmcpath, longlat, plot.params, highlight.samples)
    }
  )
  list(prob.colors = prob.colors, prob.levels = prob.levels)
}

random.eems.contour <- function(mcmcpath, dimns, longlat, plot.params, is.mrates,
                                plot.xy = NULL, highlight.samples = NULL) {
  if (is.mrates) {
    message("Plotting effective migration surface (one posterior draw of m rates)")
    zero_mean <- plot.params$m.zero_mean
    log_scale <- plot.params$m.log_scale
  } else {
    message("Plotting effective diversity surface (one posterior draw of q rates)")
    zero_mean <- plot.params$q.zero_mean
    log_scale <- plot.params$q.log_scale
  }

  voronoi <- read.voronoi(mcmcpath, longlat, is.mrates, log_scale)
  tiles <- voronoi$tiles
  rates <- voronoi$rates
  xseed <- voronoi$xseed
  yseed <- voronoi$yseed

  # Choose one saved posterior draw at random
  niters <- length(tiles)
  riter <- sample(seq_len(niters), 1)
  message(mcmcpath)
  message("Draw #", riter, " (out of ", niters, ")")

  # Jump over stored parameters for draws 1 to (riter - 1)
  skip <- sum(tiles[riter:1][-1])
  now.tiles <- tiles[riter]
  now.rates <- rates[(skip + 1):(skip + now.tiles)]
  now.xseed <- xseed[(skip + 1):(skip + now.tiles)]
  now.yseed <- yseed[(skip + 1):(skip + now.tiles)]

  rslt <- transform.rates(dimns, now.tiles, now.rates, now.xseed, now.yseed, zero_mean)
  Zvals <- rslt$Zvals
  # Actually plot the colored contour plot
  one.eems.contour(mcmcpath[1], dimns, Zvals, longlat, plot.params, is.mrates,
    plot.xy = plot.xy, highlight.samples = highlight.samples
  )
}

average.eems.contours <- function(mcmcpath, dimns, longlat, plot.params, is.mrates,
                                  plot.xy = NULL, highlight.samples = NULL) {
  if (is.mrates) {
    message("Plotting effective migration surface (posterior mean of m rates)")
    zero_mean <- plot.params$m.zero_mean
    log_scale <- plot.params$m.log_scale
  } else {
    message("Plotting effective diversity surface (posterior mean of q rates)")
    zero_mean <- plot.params$q.zero_mean
    log_scale <- plot.params$q.log_scale
  }
  Zmean <- matrix(0, dimns$nxmrks, dimns$nymrks)
  PrGT0 <- matrix(0, dimns$nxmrks, dimns$nymrks)
  PrLT0 <- matrix(0, dimns$nxmrks, dimns$nymrks)
  coord_extra <- nrow(dimns$coord)
  Zmean_extra <- numeric(coord_extra)
  niters <- 0
  # Loop over each output directory in mcmcpath to average the colored contour plots
  for (path in mcmcpath) {
    message(path)
    stopifnot(all(file.exists(file.path(path, c(
      "mcmcmtiles.txt",
      "mcmcmrates.txt",
      "mcmcxcoord.txt",
      "mcmcycoord.txt",
      "mcmcqtiles.txt",
      "mcmcqrates.txt",
      "mcmcwcoord.txt",
      "mcmczcoord.txt"
    )))))
    voronoi <- read.voronoi(path, longlat, is.mrates, log_scale)
    tiles <- voronoi$tiles
    rates <- voronoi$rates
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    rslt <- transform.rates(dimns, tiles, rates, xseed, yseed, zero_mean)
    Zmean <- Zmean + rslt$Zvals
    Zmean_extra <- Zmean_extra + rslt$Zvals_extra
    niters <- niters + rslt$niters
    if (zero_mean) {
      PrGT0 <- PrGT0 + rslt$PrGT0
      PrLT0 <- PrLT0 + rslt$PrLT0
    }
  }
  Zmean <- Zmean / niters
  PrGT0 <- PrGT0 / niters
  PrLT0 <- PrLT0 / niters
  Zmean_extra <- Zmean_extra / niters

  # Actually plot the colored contour plot
  # Pass one mcmcpath (shouldn't matter which one) in case adding extra information
  # (demes, edges, etc.) on top of the contour plot.
  rates.raster <- one.eems.contour(
    mcmcpath[1], dimns, Zmean, longlat, plot.params, is.mrates,
    plot.xy = plot.xy, highlight.samples = highlight.samples
  )
  xyz.values <- cbind(dimns$coord, Zmean_extra)
  rates.raster$xyz.values <- xyz.values

  if (zero_mean) {
    raster.prob <- one.prob.contour(
      mcmcpath[1], dimns, PrGT0 - PrLT0, longlat, plot.params, is.mrates,
      plot.xy = plot.xy, highlight.samples = highlight.samples
    )
  }
  rates.raster
}
