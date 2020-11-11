
check.plot.params <- function(params) {

  ## Is there a way to check if a string is a valid PROJ.4 string?
  if (is.character(params$proj.in)) {
    params$proj.in <- params$proj.in[1]
  } else {
    params$proj.in <- NULL
  }
  if (is.character(params$proj.out)) {
    params$proj.out <- params$proj.out[1]
  } else {
    params$proj.out <- params$proj.in
  }
  message(
    "Input projection: ", params$proj.in, "\n",
    "Output projection: ", params$proj.out, "\n\n\n"
  )
  if (is.null(params$proj.in)) {
    if (!is.null(params$proj.out)) {
      stop("Specify the input projection, projection.in, as a PROJ.4 string")
    }
  } else {
    load.required.package(package = "rgdal", required.by = "projection.in")
  }

  if (is.logical(params$add.map)) {
    params$add.map <- params$add.map[1]
  } else {
    params$add.map <- FALSE
  }
  if (is.color(params$col.map)) {
    params$col.map <- params$col.map[1]
  } else {
    params$col.map <- "#AAAAAA"
  }
  if (is.numeric(params$lwd.map)) {
    params$lwd.map <- params$lwd.map[1]
  } else {
    params$lwd.map <- 2
  }

  if (is.logical(params$add.grid)) {
    params$add.grid <- params$add.grid[1]
  } else {
    params$add.grid <- FALSE
  }
  if (is.color(params$col.grid)) {
    params$col.grid <- params$col.grid[1]
  } else {
    params$col.grid <- "#BBBBBB"
  }
  if (is.numeric(params$lwd.grid)) {
    params$lwd.grid <- params$lwd.grid[1]
  } else {
    params$lwd.grid <- 1
  }

  if (is.logical(params$add.outline)) {
    params$add.outline <- params$add.outline[1]
  } else {
    params$add.outline <- FALSE
  }
  if (is.color(params$col.outline)) {
    params$col.outline <- params$col.outline[1]
  } else {
    params$col.outline <- "#EEEEEE"
  }
  if (is.numeric(params$lwd.outline)) {
    params$lwd.outline <- params$lwd.outline[1]
  } else {
    params$lwd.outline <- 2
  }

  if (is.logical(params$add.demes)) {
    params$add.demes <- params$add.demes[1]
  } else {
    params$add.demes <- FALSE
  }
  if (is.logical(params$all.demes)) {
    params$all.demes <- params$all.demes[1]
  } else {
    params$all.demes <- FALSE
  }
  if (is.color(params$col.demes)) {
    params$col.demes <- params$col.demes[1]
  } else {
    params$col.demes <- "#000000"
  }
  if (is.numeric(params$pch.demes)) {
    params$pch.demes <- params$pch.demes[1]
  } else {
    params$pch.demes <- 19
  }
  if (is.numeric(params$min.cex.demes)) {
    params$min.cex.demes <- params$min.cex.demes[1]
  } else {
    params$min.cex.demes <- 1
  }
  if (is.numeric(params$max.cex.demes)) {
    params$max.cex.demes <- params$max.cex.demes[1]
  } else {
    params$max.cex.demes <- 3
  }
  if (params$max.cex.demes < params$min.cex.demes) {
    params$max.cex.demes <- params$min.cex.demes
  }

  if (is.logical(params$add.seeds)) {
    params$add.seeds <- params$add.seeds[1]
  } else {
    params$add.seeds <- FALSE
  }
  if (is.color(params$col.seeds)) {
    params$col.seeds <- params$col.seeds[1]
  } else {
    params$col.seeds <- "#8AE234"
  }
  if (is.numeric(params$pch.seeds)) {
    params$pch.seeds <- params$pch.seeds[1]
  } else {
    params$pch.seeds <- 19
  }
  if (is.numeric(params$cex.seeds)) {
    params$cex.seeds <- params$cex.seeds[1]
  } else {
    params$cex.seeds <- 1
  }

  if (params$add.map) {
    load.required.package(package = "rworldmap", required.by = "add.map")
    load.required.package(package = "rworldxtra", required.by = "add.map")
    if (is.null(params$proj.in)) {
      stop(paste0(
        "To add a geographical map, specify the input and output projections.\n",
        "For example, if the coordinates are longitude and latitude, \n",
        "you can use the PROJ.4 string '+proj=longlat +datum=WGS84'\n\n\n"
      ))
    }
  }

  if (is.numeric(params$m.colscale)) {
    params$m.colscale <- set.colscale(params$m.colscale)
  } else {
    params$m.colscale <- c(-2.5, 2.5)
  }
  if (is.numeric(params$q.colscale)) {
    params$q.colscale <- set.colscale(params$q.colscale)
  } else {
    params$q.colscale <- c(-0.1, 0.1)
  }
  if (length(params$eems.colors) < 2 || !all(is.color(params$eems.colors))) {
    params$eems.colors <- default.eems.colors()
  }

  if (is.logical(params$add.title)) {
    params$add.title <- params$add.title[1]
  } else {
    params$add.title <- TRUE
  }
  if (is.logical(params$add.colbar)) {
    params$add.colbar <- params$add.colbar[1]
  } else {
    params$add.colbar <- TRUE
  }

  ## Additional options: By default, both the effectove migration rates and
  ## the effective diversity rates are log10-transformed and mean-zero standardized
  if (is.logical(params$m.zero_mean)) {
    params$m.zero_mean <- params$m.zero_mean[1]
  } else {
    params$m.zero_mean <- TRUE
  }
  if (is.logical(params$q.zero_mean)) {
    params$q.zero_mean <- params$q.zero_mean[1]
  } else {
    params$q.zero_mean <- TRUE
  }
  if (is.logical(params$m.log_scale)) {
    params$m.log_scale <- params$m.log_scale[1]
  } else {
    params$m.log_scale <- TRUE
  }
  if (is.logical(params$q.log_scale)) {
    params$q.log_scale <- params$q.log_scale[1]
  } else {
    params$q.log_scale <- TRUE
  }

  params$prob.levels <- params$prob.levels[params$prob.levels >= 0]
  params$prob.levels <- params$prob.levels[params$prob.levels <= 1]
  params$prob.levels <- sort(unique(c(0, params$prob.levels, 1)))

  params
}

sub.axes.labels <- function() {
  JtDJ <- list(
    xlab = expression(paste("Fitted dissimilarity between individuals  ", Delta[i * j])),
    ylab = expression(paste("Observed dissimilarity between individuals  ", D[i * j])),
    mainTRUE = "There should be at least two observed demes to plot pairwise dissimilarities"
  )
  Between <- list(
    xlab = expression(paste(
      "Fitted dissimilarity between demes:   ",
      Delta[alpha * beta], " - (", Delta[alpha * alpha], "+", Delta[beta * beta], ") / 2"
    )),
    ylab = expression(paste(
      "Observed dissimilarity between demes:   ",
      D[alpha * beta], " - (", D[alpha * alpha], "+", D[beta * beta], ") / 2"
    )),
    mainTRUE = expression(paste(
      "Dissimilarities between pairs of sampled demes (",
      alpha, ", ", beta, ")"
    )),
    subTRUE = "Singleton demes, if any, are excluded from this plot (but not from EEMS)",
    mainFALSE = expression(paste(
      "Dissimilarities between pairs of sampled demes (",
      alpha, ", ", beta, ")"
    )),
    subFALSE = expression(paste(
      "Gray means that a single individual is sampled from either ",
      alpha, " or ", beta
    ))
  )
  Within <- list(
    xlab = expression(paste("Fitted dissimilarity within demes:   ", Delta[alpha * alpha])),
    ylab = expression(paste("Observed dissimilarity within demes:   ", D[alpha * alpha])),
    mainTRUE = expression(paste("Dissimilarities within sampled demes ", alpha)),
    subTRUE = "Singleton demes, if any, are excluded from this plot (but not from EEMS)",
    mainFALSE = expression(paste("Dissimilarities within sampled demes ", alpha))
  )
  GeoDist <- Between
  GeoDist$xlab <- "Great circle distance between demes (km)"
  list(JtDJ = JtDJ, Between = Between, Within = Within, GeoDist = GeoDist)
}

# dist.type takes the values: "JtDJ", "Between", "Within", "GeoDist"
dist.axes.labels <- function(dist.type, remove.singletons = TRUE, subtitle = NULL) {
  labels <- sub.axes.labels()[[dist.type]]
  title(xlab = labels$xlab)
  title(ylab = labels$ylab)
  if (remove.singletons) {
    mtext(side = 3, line = 2, cex = 1.3, text = labels$mainTRUE)
    if (is.null(subtitle)) {
      mtext(side = 3, line = 0.5, cex = 1, text = labels$subTRUE)
    } else {
      mtext(side = 3, line = 0.5, cex = 1, text = subtitle)
    }
  } else {
    mtext(side = 3, line = 2, cex = 1.3, text = labels$mainFALSE)
    if (is.null(subtitle)) {
      mtext(side = 3, line = 0.5, cex = 1, text = labels$subFALSE)
    } else {
      mtext(side = 3, line = 0.5, cex = 1, text = subtitle)
    }
  }
}

sub.scattercols <- function(sizes) {
  c("black", "gray60")[1 + (sizes < 2)]
}

sub.scatterplot <- function(dist.type, dist.data, remove.singletons, add.abline, add.r.squared,
                            subtitle = NULL, add = FALSE) {
  if (remove.singletons) {
    dist.data <- dist.data[dist.data$size > 1, ]
  }
  if (is.null(dist.data$size)) dist.data$size <- 2
  if (is.null(dist.data$pch)) dist.data$pch <- 1
  if (is.null(dist.data$cex)) dist.data$cex <- 1
  if (is.null(dist.data$col)) dist.data$col <- 1

  group <- dist.data$col
  # It turns out sort.list sorts alphabetically by group name
  # The following works to sort by number of occurrences:
  ord <- sort.list(table(group)[group], decreasing = TRUE)
  dist.data <- dist.data[ord, ]

  if (!add) {
    plot(dist.data$fitted, dist.data$obsrvd, type = "n", xlab = "", ylab = "")
    dist.axes.labels(dist.type, remove.singletons, subtitle)
    if (add.abline) abline(a = 0, b = 1, col = "red", lwd = 2)
    if (add.r.squared) {
      # Fit a linear model for the observed dissimilarities as a function of the fitted ones
      r.squared <- summary(lm(dist.data$obsrvd ~ dist.data$fitted))$r.squared
      r.squared <- round(r.squared, digits = 3)
      legend("topleft", legend = substitute(
        paste(R^2, " = ", val),
        list(val = r.squared)
      ), bty = "n")
    }
  }
  points(dist.data$fitted,
    dist.data$obsrvd,
    col = dist.data$col,
    pch = dist.data$pch,
    cex = dist.data$cex
  )
}

# The distance metric used by EEMS
# Euclidean distance by default; great circle (haversine) distance as an alternative
which.dist.metric <- function(mcmcpath) {
  dist.metric <- "euclidean"
  Lines <- readLines(file.path(mcmcpath, "eemsrun.txt"))
  for (i in seq(Lines)) {
    s <- gsub("\\s", "", Lines[i]) # Remove any empty space
    x <- strsplit(s, "distance = ")[[1]] # Is there a line 'distance = xxx'?
    if (length(x) == 2) dist.metric <- tolower(x[2]) # What is the distance metric?
  }
  if ((dist.metric == "euclidean") || (dist.metric == "greatcirc")) {
    message(
      "Using '", dist.metric,
      "' distance to assign interpolation points to Voronoi tiles.\n\n\n"
    )
  } else {
    stop("Specify either 'euclidean' or 'greatcirc' distance metric in eemsrun.txt.")
  }
  dist.metric
}

read.dimns <- function(path, longlat, nxmrks = NULL, nymrks = NULL, coord = NULL) {
  eems.output <- NULL
  datapath.outer <- paste0(path, ".outer")
  mcmcpath.outer <- file.path(path, "outer.txt")
  if (file.exists(mcmcpath.outer)) {
    eems.output <- TRUE
  } else if (file.exists(datapath.outer)) {
    eems.output <- FALSE
  }
  if (is.null(eems.output)) {
    stop(paste0(path, " is neither a datapath nor a mcmcpath."))
  }
  if (eems.output) {
    outer <- scan(mcmcpath.outer, what = numeric(), quiet = TRUE)
  } else {
    outer <- scan(datapath.outer, what = numeric(), quiet = TRUE)
  }
  outer <- matrix(outer, ncol = 2, byrow = TRUE)
  if (is.null(coord)) {
    coord <- matrix(0, ncol = 2, nrow = 0)
  } else {
    coord <- as.matrix(coord, ncol = 2, byrow = TRUE)
  }
  # "Close" the outline if the first row is not the same as the last row
  if (sum(head(outer, 1) != tail(outer, 1))) {
    outer <- rbind(outer, head(outer, 1))
  }
  if (!longlat) {
    outer <- outer[, c(2, 1)]
    coord <- coord[, c(2, 1)]
  }
  xlim <- range(outer[, 1])
  ylim <- range(outer[, 2])
  aspect <- abs((diff(ylim) / diff(xlim)) / cos(mean(ylim) * pi / 180))
  # Choose the number of interpolation in each direction
  if (is.null(nxmrks) && is.null(nymrks)) {
    if (aspect > 1) {
      nxmrks <- 100
      nymrks <- round(nxmrks * aspect)
    } else {
      nymrks <- 100
      nxmrks <- round(nymrks / aspect)
    }
  }
  # The interpolation points are equally spaced
  xmrks <- seq(from = xlim[1], to = xlim[2], length = nxmrks)
  ymrks <- seq(from = ylim[1], to = ylim[2], length = nymrks)
  marks <- cbind(rep(xmrks, times = nymrks), rep(ymrks, each = nxmrks))
  if (eems.output) {
    dist.metric <- which.dist.metric(path)
  } else {
    dist.metric <- "euclidean"
  }
  # If no additional coordinates at which to estimate m and q rates
  # are specified, report the m and q estimates at the raster marks
  if (!nrow(coord)) {
    coord <- marks
  }
  list(
    nxmrks = nxmrks, xmrks = xmrks, xlim = xlim, xspan = diff(xlim),
    nymrks = nymrks, ymrks = ymrks, ylim = ylim, yspan = diff(ylim),
    marks = marks, nmrks = c(nxmrks, nymrks), aspect = aspect,
    outer = outer, coord = coord, dist.metric = dist.metric
  )
}

read.edges <- function(mcmcpath) {
  edges <- read.table(file.path(mcmcpath, "edges.txt"), colClasses = numeric())
  edges <- as.matrix(edges)
  # Previously EEMS output the edges with one vertex per line and
  # the six neighbors of each vertex listed in order
  # Currently EEMS outputs the edges with one edge per line (as a
  # pair of vertices) which is more general
  if (ncol(edges) == 6) {
    # Convert the old format to the new format
    edges0 <- edges
    edges <- matrix(0, nrow = sum(edges0 > 0), ncol = 2)
    nv <- nrow(edges0)
    nn <- ncol(edges0)
    nodes <- 1:nv
    e <- 0
    for (a in 1:nv) {
      for (i in 1:nn) {
        b <- edges0[a, i]
        if (b %in% nodes) {
          e <- e + 1
          edges[e, 1] <- a
          edges[e, 2] <- b
        }
      }
    }
  }
  edges
}

read.graph <- function(path, longlat) {
  eems.output <- NULL
  if (file.exists(file.path(path, "demes.txt")) &&
    file.exists(file.path(path, "ipmap.txt")) &&
    file.exists(file.path(path, "outer.txt"))) {
    eems.output <- TRUE
  } else if (file.exists(paste0(path, ".coord")) &&
    file.exists(paste0(path, ".diffs")) &&
    file.exists(paste0(path, ".outer"))) {
    eems.output <- FALSE
  }
  if (is.null(eems.output)) {
    stop(paste0(path, " is neither a datapath nor a mcmcpath."))
  }
  if (eems.output) {
    # Read the assigned sample coordinates
    ipmap <- scan(file.path(path, "ipmap.txt"), what = numeric(), quiet = TRUE)
    demes <- scan(file.path(path, "demes.txt"), what = numeric(), quiet = TRUE)
    outer <- scan(file.path(path, "outer.txt"), what = numeric(), quiet = TRUE)
    demes <- matrix(demes, ncol = 2, byrow = TRUE)
    outer <- matrix(outer, ncol = 2, byrow = TRUE)
    edges <- read.edges(path)
  } else {
    # Read the original sample coordinates
    coord <- scan(paste0(path, ".coord"), what = numeric(), quiet = TRUE)
    outer <- scan(paste0(path, ".outer"), what = numeric(), quiet = TRUE)
    coord <- matrix(coord, ncol = 2, byrow = TRUE)
    outer <- matrix(outer, ncol = 2, byrow = TRUE)
    edges <- NULL
    # The two sampling coordinates are combined to define a "deme",
    # with the maximum possible precision
    ipmap <- factor(paste(coord[, 1], coord[, 2], sep = "x"))
    demes <- levels(ipmap)
    index <- match(demes, ipmap)
    o <- length(index)
    demes <- coord[index, ]
    ipmap <- (1:o)[ipmap]
    # "Close" the outline if the first row is not the same as the last row
    if (sum(head(outer, 1) != tail(outer, 1))) {
      outer <- rbind(outer, head(outer, 1))
    }
  }
  if (!longlat) {
    demes <- demes[, c(2, 1)]
    outer <- outer[, c(2, 1)]
  }
  sizes <- table(ipmap)
  alpha <- as.numeric(names(sizes))
  sizes <- as.numeric(sizes)
  list(ipmap = ipmap, demes = demes, edges = edges, alpha = alpha, sizes = sizes, outer = outer)
}

read.voronoi <- function(mcmcpath, longlat, is.mrates, log_scale) {
  if (is.mrates) {
    rates <- scan(file.path(mcmcpath, "mcmcmrates.txt"),
      what = numeric(), quiet = TRUE
    )
    tiles <- scan(file.path(mcmcpath, "mcmcmtiles.txt"),
      what = numeric(), quiet = TRUE
    )
    xseed <- scan(file.path(mcmcpath, "mcmcxcoord.txt"),
      what = numeric(), quiet = TRUE
    )
    yseed <- scan(file.path(mcmcpath, "mcmcycoord.txt"),
      what = numeric(), quiet = TRUE
    )
  } else {
    rates <- scan(file.path(mcmcpath, "mcmcqrates.txt"),
      what = numeric(), quiet = TRUE
    )
    tiles <- scan(file.path(mcmcpath, "mcmcqtiles.txt"),
      what = numeric(), quiet = TRUE
    )
    xseed <- scan(file.path(mcmcpath, "mcmcwcoord.txt"),
      what = numeric(), quiet = TRUE
    )
    yseed <- scan(file.path(mcmcpath, "mcmczcoord.txt"),
      what = numeric(), quiet = TRUE
    )
  }
  if (!longlat) {
    tempi <- xseed
    xseed <- yseed
    yseed <- tempi
  }
  if (log_scale) {
    rates <- log10(rates)
  }
  list(rates = rates, tiles = tiles, xseed = xseed, yseed = yseed)
}
