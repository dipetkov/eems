
#' @useDynLib rEEMSplots
#' @import raster rgeos sp geosphere RcppEigen
#' @importFrom Rcpp evalCpp

default.eems.colors <- function( ) {
    writeLines(paste0("Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.\n", 
                      "It combines two color schemes from the 'dichromat' package, which itself is based on\n", 
                      "a collection of color schemes for scientific data graphics:\n", 
                      "\tLight A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data\n", 
                      "\tGraphics. EOS Transactions of the American Geophysical Union, 85(40), 385.\n", 
                      "See also http://geog.uoregon.edu/datagraphics/color_scales.htm\n\n\n"))
    ## To reproduce the default eems colors:
    ## Oranges <- dichromat::colorschemes$BluetoDarkOrange.12[12:7]
    ## Blues <- dichromat::colorschemes$BrowntoBlue.12[7:12]
    ## eems.colors <- c(Oranges, "#FBFBFB", Blues)
    eems.colors <- c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", ## orange sequence
                     "#FBFBFB", ## very slightly off-white
                     "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99") ## blue sequence
    return (eems.colors)
}
sub.axes.labels <- function() {
    JtDJ <- list(
        xlab = expression(paste("Fitted dissimilarity between individuals  ", Delta[i * j])), 
        ylab = expression(paste("Observed dissimilarity between individuals  ", D[i * j])), 
        mainTRUE = "There should be at least two observed demes to plot pairwise dissimilarities")
    Between <- list(
        xlab = expression(paste("Fitted dissimilarity between demes:   ", 
            Delta[alpha * beta], " - (", Delta[alpha*alpha], "+", Delta[beta * beta], ") / 2")), 
        ylab = expression(paste("Observed dissimilarity between demes:   ", 
            D[alpha * beta], " - (", D[alpha * alpha], "+", D[beta * beta], ") / 2")), 
        mainTRUE = expression(paste("Dissimilarities between pairs of sampled demes (", alpha , ", ", beta, ")")), 
        subTRUE = "Singleton demes, if any, are excluded from this plot (but not from EEMS)", 
        mainFALSE = expression(paste("Dissimilarities between pairs of sampled demes (", alpha , ", ", beta, ")")), 
        subFALSE = expression(paste("Gray means that a single individual is sampled from either ", alpha, " or ", beta)))
    Within <- list(
        xlab = expression(paste("Fitted dissimilarity within demes:   ", Delta[alpha * alpha])), 
        ylab = expression(paste("Observed dissimilarity within demes:   ", D[alpha * alpha])), 
        mainTRUE = expression(paste("Dissimilarities within sampled demes ", alpha)), 
        subTRUE = "Singleton demes, if any, are excluded from this plot (but not from EEMS)", 
        mainFALSE = expression(paste("Dissimilarities within sampled demes ", alpha)))
    GeoDist <- Between
    GeoDist$xlab <- "Great circle distance between demes (km)"    
    return (list(JtDJ = JtDJ, Between = Between, Within = Within, GeoDist = GeoDist))
}
## dist.type takes the values: "JtDJ", "Between", "Within", "GeoDist"
dist.axes.labels <- function(dist.type, remove.singletons = TRUE, subtitle = NULL) {
    labels = sub.axes.labels()[[dist.type]]
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
    return (c("black", "gray60")[1 + (sizes < 2)])
}
sub.scatterplot <- function(dist.type, dist.data, remove.singletons, add.abline, subtitle = NULL, add = FALSE) {
    if (remove.singletons) {
        dist.data <- dist.data[dist.data$size > 1, ]
    }
    if (is.null(dist.data$size)) dist.data$size <- 2
    if (is.null(dist.data$pch))  dist.data$pch <- 1
    if (is.null(dist.data$cex))  dist.data$cex <- 1
    if (is.null(dist.data$col))  dist.data$col <- 1

    group <- dist.data$col
    ## It turns out sort.list sorts alphabetically by group name:
    ## ord = sort.list(group, decreasing = TRUE)
    ## The following works to sort by number of occurrences:
    ord <- sort.list(table(group)[group], decreasing = TRUE)
    dist.data <- dist.data[ord, ]
    
    if (!add) {
        plot(dist.data$fitted, dist.data$obsrvd, type = "n", xlab = "", ylab = "")
        dist.axes.labels(dist.type, remove.singletons, subtitle)
        if (add.abline) abline(a = 0, b = 1, col = "red", lwd = 2)
    }
    points(dist.data$fitted, 
           dist.data$obsrvd, 
           col = dist.data$col, 
           pch = dist.data$pch, 
           cex = dist.data$cex)
}
JtDJ2BandW <- function(JtDJ, sizes = NULL) {
    ## There should be NAs on the main diagonal of JtDobsJ -- these elements correspond to demes with
    ## a single observed sample. And since there is a single sample taken, the average dissimilarity
    ## between two distinct individuals from such demes cannot be observed.
    ## Instead, I will use the average within dissimilarity, W, across the demes with multiple samples.
    ## This is make a difference only if remove.singletons = FALSE, which is not the default option.
    JtDJ <- as.matrix(JtDJ)
    if (!is.null(sizes)) { diag(JtDJ)[sizes < 2] = mean(diag(JtDJ)[sizes >= 2]) }
    n <- nrow(JtDJ)
    W <- diag(JtDJ)
    S <- matrix(W, n, n)
    B <- JtDJ - (S + t(S)) / 2
    B <- B[upper.tri(B, diag = FALSE)]
    return (list(W = W, B = B))
}
## Check that a string, or a vector of strings, represents a hex color
is.color <- function(x) {
    if (is.null(x)) { return(FALSE) }
    ## grepl("^#[0-9A-F]{6}$", x)
    sapply(x, function(x) { tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })
}
set.colscale <- function(x, is.mrate) {
    colscale <- c(-1, 1)
    if (is.numeric(x) && min(x) < max(x)) {
        colscale <- c(min(x), max(x))
    }
    return(colscale)
}
check.plot.params <- function(params) {

    ## Is there a way to check if a string is a valid PROJ.4 string?
    if (is.character(params$proj.in)) params$proj.in <- params$proj.in[1]
    else params$proj.in <- NULL
    if (is.character(params$proj.out)) params$proj.out <- params$proj.out[1]
    else params$proj.out <- params$proj.in
    writeLines(paste0("Input projection: ", params$proj.in, "\n", 
                      "Output projection: ", params$proj.out, "\n\n\n"))
    if (is.null(params$proj.in)) {
        if (!is.null(params$proj.out))
            stop("Specify the input projection, projection.in, as a PROJ.4 string")
    } else {
        load.required.package(package = "rgdal", required.by = "projection.in")
    }
    
    if (is.logical(params$add.map)) params$add.map <- params$add.map[1]
    else params$add.map <- FALSE
    if (is.color(params$col.map))   params$col.map <- params$col.map[1]
    else params$col.map <- "#AAAAAA"
    if (is.numeric(params$lwd.map)) params$lwd.map <- params$lwd.map[1]
    else params$lwd.map <- 2

    if (is.logical(params$add.grid)) params$add.grid <- params$add.grid[1]
    else params$add.grid <- FALSE
    if (is.color(params$col.grid))   params$col.grid <- params$col.grid[1]
    else params$col.grid <- "#BBBBBB"
    if (is.numeric(params$lwd.grid)) params$lwd.grid <- params$lwd.grid[1]
    else params$lwd.grid <- 1

    if (is.logical(params$add.outline)) params$add.outline <- params$add.outline[1]
    else params$add.outline <- FALSE
    if (is.color(params$col.outline))   params$col.outline <- params$col.outline[1]
    else params$col.outline <- "#EEEEEE"
    if (is.numeric(params$lwd.outline)) params$lwd.outline <- params$lwd.outline[1]
    else params$lwd.outline <- 2

    if (is.logical(params$add.demes)) params$add.demes <- params$add.demes[1]
    else params$add.demes <- FALSE
    if (is.logical(params$all.demes)) params$all.demes <- params$all.demes[1]
    else params$all.demes <- FALSE
    if (is.color(params$col.demes))   params$col.demes <- params$col.demes[1]
    else params$col.demes <- "#000000"
    if (is.numeric(params$pch.demes)) params$pch.demes <- params$pch.demes[1]
    else params$pch.demes <- 19
    if (is.numeric(params$min.cex.demes)) params$min.cex.demes <- params$min.cex.demes[1]
    else params$min.cex.demes <- 1
    if (is.numeric(params$max.cex.demes)) params$max.cex.demes <- params$max.cex.demes[1]
    else params$max.cex.demes <- 3
    if (params$max.cex.demes < params$min.cex.demes)
        params$max.cex.demes <- params$min.cex.demes

    if (is.logical(params$add.seeds)) params$add.seeds <- params$add.seeds[1]
    else params$add.seeds <- FALSE
    if (is.color(params$col.seeds))   params$col.seeds <- params$col.seeds[1]
    else params$col.seeds <- "#8AE234"
    if (is.numeric(params$pch.seeds)) params$pch.seeds <- params$pch.seeds[1]
    else params$pch.seeds <- 19
    if (is.numeric(params$cex.seeds)) params$cex.seeds <- params$cex.seeds[1]
    else params$cex.seeds <- 1

    if (params$add.map) {
        load.required.package(package = 'rworldmap', required.by = 'add.map')
        load.required.package(package = 'rworldxtra', required.by = 'add.map')
        if (is.null(params$proj.in)) {
            stop(paste0("To add a geographical map, specify the input and output projections.\n", 
                        "For example, if the coordinates are longitude and latitude, \n", 
                        "you can use the PROJ.4 string '+proj=longlat +datum=WGS84'\n\n\n"))
        }
    }

    if (is.numeric(params$m.colscale)) params$m.colscale <- set.colscale(params$m.colscale)
    else params$m.colscale <- c(-2.5, 2.5)
    if (is.numeric(params$q.colscale)) params$q.colscale <- set.colscale(params$q.colscale)
    else params$q.colscale <- c(-0.1, 0.1)
    if (length(params$eems.colors) < 2 || sum(!is.color(params$eems.colors)))
        params$eems.colors <- default.eems.colors( )
    
    if (is.logical(params$add.title))  params$add.title <- params$add.title[1]
    else params$add.title <- TRUE
    if (is.logical(params$add.colbar)) params$add.colbar <- params$add.colbar[1]
    else params$add.colbar <- TRUE

    if (is.logical(params$m.zero_mean)) params$m.zero_mean <- params$m.zero_mean[1]
    else params$m.zero_mean <- TRUE
    if (is.logical(params$q.zero_mean)) params$q.zero_mean <- params$q.zero_mean[1]
    else params$q.zero_mean <- TRUE
    if (is.logical(params$m.log_scale)) params$m.log_scale <- params$m.log_scale[1]
    else params$m.log_scale <- TRUE
    if (is.logical(params$q.log_scale)) params$q.log_scale <- params$q.log_scale[1]
    else params$q.log_scale <- TRUE

    return(params)
}
eems.colscale <- function(Zvals, num.levels, colscale) {
    maxZ <- max(Zvals)
    minZ <- min(Zvals)
    ## Check whether the range of values should be fixed so that
    ## migration rates can be plotted on the same scale
    if (is.numeric(colscale) && minZ >= min(colscale) && maxZ <= max(colscale)) {
        minZ <- colscale[1]
        maxZ <- colscale[2]
        step <- (maxZ - minZ) / num.levels
        eems.levels <- seq(from = minZ, to = maxZ, by = step)
    } else {
        eems.levels <- seq(from = minZ, to = maxZ, length = num.levels + 1)
    }
    return(eems.levels)
}
## The distance metric used by EEMS
## Euclidean distance by default; great circle (haversine) distance as an alternative
which.dist.metric <- function(mcmcpath) {
  dist.metric <- "euclidean"
  Lines <- readLines(paste0(mcmcpath, "/eemsrun.txt"))
  for (i in seq(Lines)) {
    s <- gsub("\\s", "", Lines[i])       ## Remove any empty space
    x <- strsplit(s, "distance = ")[[1]]  ## Is there a line 'distance = xxx'?
    if (length(x) == 2) dist.metric <- tolower(x[2]) ## What is the distance metric, in all lower case?
  }
  if ((dist.metric == "euclidean") || (dist.metric == "greatcirc")) {
      writeLines(paste0("Using '", dist.metric, "' distance to assign interpolation points to Voronoi tiles.\n\n\n"))
  } else {
      stop("Specify either 'euclidean' or 'greatcirc' distance metric in eemsrun.txt.")
  }
  return (dist.metric)
}
read.dimns <- function(path, longlat, nxmrks = NULL, nymrks = NULL) {
    eems.output <- NULL
    datapath.outer <- paste0(path, '.outer')
    mcmcpath.outer <- paste0(path, '/outer.txt')
    if (file.exists(mcmcpath.outer)) {
        eems.output <- TRUE
    } else if (file.exists(datapath.outer)) {
        eems.output <- FALSE
    }
    if (is.null(eems.output)) {
        stop(paste0(path, ' is neither a datapath nor a mcmcpath.'))
    }
    if (eems.output) {
        outer <- scan(mcmcpath.outer, what = numeric(), quiet = TRUE)
    } else {
        outer <- scan(datapath.outer, what = numeric(), quiet = TRUE)
    }
    outer <- matrix(outer, ncol = 2, byrow = TRUE)
    ## "Close" the outline if the first row is not the same as the last row
    if (sum(head(outer, 1) != tail(outer, 1))) {
        outer <- rbind(outer, head(outer, 1))
    }
    if (!longlat) { outer <- outer[, c(2, 1)] }
    xlim <- range(outer[, 1])
    ylim <- range(outer[, 2])
    aspect <- abs((diff(ylim) / diff(xlim)) / cos(mean(ylim) * pi / 180))
    ## Choose the number of interpolation in each direction
    if (is.null(nxmrks) && is.null(nymrks)) {
        if (aspect > 1) {
            nxmrks <- 100
            nymrks <- round(nxmrks * aspect)
        } else {
            nymrks <- 100
            nxmrks <- round(nymrks / aspect)
        }
    }
    ## The interpolation points are equally spaced
    xmrks <- seq(from = xlim[1], to = xlim[2], length = nxmrks)
    ymrks <- seq(from = ylim[1], to = ylim[2], length = nymrks)
    marks <- cbind(rep(xmrks, times = nymrks), rep(ymrks, each = nxmrks))
    if (eems.output) {
        dist.metric <- which.dist.metric(path)
    } else {
        dist.metric <- 'euclidean'
    }
    return(list(nxmrks = nxmrks, xmrks = xmrks, xlim = xlim, xspan = diff(xlim), 
                nymrks = nymrks, ymrks = ymrks, ylim = ylim, yspan = diff(ylim), 
                marks = marks, nmrks = c(nxmrks, nymrks), aspect = aspect, 
                outer = outer, dist.metric = dist.metric))
}
read.edges <- function(mcmcpath) {
    edges <- read.table(paste0(mcmcpath, '/edges.txt'), colClasses = numeric())
    edges <- as.matrix(edges)
    ## Previously EEMS output the edges with one vertex per line and
    ## the six neighbors of each vertex listed in order
    ## Currently EEMS outputs the edges with one edge per line (as a
    ## pair of vertices) which is more general
    if (ncol(edges) == 6) {
        ## Convert the old format to the new format
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
        } }
    }
    return(edges)
}
read.graph <- function(path, longlat) {
    eems.output <- NULL
    if (file.exists(paste0(path, '/demes.txt')) &&
        file.exists(paste0(path, '/ipmap.txt')) &&
        file.exists(paste0(path, '/outer.txt'))) {
        eems.output <- TRUE
    } else if (file.exists(paste0(path, '.coord')) &&
               file.exists(paste0(path, '.diffs')) &&
               file.exists(paste0(path, '.outer'))) {
        eems.output <- FALSE
    }
    if (is.null(eems.output)) {
        stop(paste0(path, ' is neither a datapath nor a mcmcpath.'))
    }
    if (eems.output) {
        ## Read the assigned sample coordinates
        ipmap <- scan(paste0(path, '/ipmap.txt'), what = numeric(), quiet = TRUE)
        demes <- scan(paste0(path, '/demes.txt'), what = numeric(), quiet = TRUE)
        outer <- scan(paste0(path, '/outer.txt'), what = numeric(), quiet = TRUE)
        demes <- matrix(demes, ncol = 2, byrow = TRUE)
        outer <- matrix(outer, ncol = 2, byrow = TRUE)
        edges <- read.edges(path)
    } else {
        ## Read the original sample coordinates
        coord <- scan(paste0(path, '.coord'), what = numeric(), quiet = TRUE)
        outer <- scan(paste0(path, '.outer'), what = numeric(), quiet = TRUE)
        coord <- matrix(coord, ncol = 2, byrow = TRUE)
        outer <- matrix(outer, ncol = 2, byrow = TRUE)
        edges <- NULL
        ## Case 1: Each sample is its own "deme"
        ## Samples with exactly the same location are overplotted
        #ipmap <- seq(nrow(coord))
        #demes <- coord
        ## Case 2: The two sampling coordinates are combined
        ## to define a "deme", with the maximum possible precision
        ipmap <- factor(paste(coord[, 1], coord[, 2], sep = "x"))
        demes <- levels(ipmap)
        index <- match(demes, ipmap)
        o <- length(index)
        demes <- coord[index, ]
        ipmap <- (1:o)[ipmap]
        ## "Close" the outline if the first row is not the same as the last row
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
    return(list(ipmap = ipmap, demes = demes, edges = edges, alpha = alpha, sizes = sizes, outer = outer))
}
read.voronoi <- function(mcmcpath, longlat, is.mrates, log_scale) {
    if (is.mrates) {
        rates <- scan(paste0(mcmcpath, '/mcmcmrates.txt'), what = numeric(), quiet = TRUE)
        tiles <- scan(paste0(mcmcpath, '/mcmcmtiles.txt'), what = numeric(), quiet = TRUE)
        xseed <- scan(paste0(mcmcpath, '/mcmcxcoord.txt'), what = numeric(), quiet = TRUE)
        yseed <- scan(paste0(mcmcpath, '/mcmcycoord.txt'), what = numeric(), quiet = TRUE)
    } else {
        rates <- scan(paste0(mcmcpath, '/mcmcqrates.txt'), what = numeric(), quiet = TRUE)
        tiles <- scan(paste0(mcmcpath, '/mcmcqtiles.txt'), what = numeric(), quiet = TRUE)
        xseed <- scan(paste0(mcmcpath, '/mcmcwcoord.txt'), what = numeric(), quiet = TRUE)
        yseed <- scan(paste0(mcmcpath, '/mcmczcoord.txt'), what = numeric(), quiet = TRUE)
    }
    if (!longlat) {
        tempi <- xseed
        xseed <- yseed
        yseed <- tempi
    }
    if (log_scale) {
        rates <- log10(rates)
    }
    return(list(rates = rates, tiles = tiles, xseed = xseed, yseed = yseed))
}
transform.rates <- function(mcmcpath, dimns, longlat, is.mrates, log_scale, zero_mean) {
    voronoi <- read.voronoi(mcmcpath, longlat, is.mrates, log_scale)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    if (zero_mean) {
        Zvals <- .Call("rEEMSplots__rcppstandardize_rates", PACKAGE = "rEEMSplots", 
                       tiles, rates, xseed, yseed, dimns$marks, dimns$nmrks, dimns$dist.metric)
    } else {
        Zvals <- .Call("rEEMSplots__rcppdont_standardize_rates", PACKAGE = "rEEMSplots", 
                       tiles, rates, xseed, yseed, dimns$marks, dimns$nmrks, dimns$dist.metric)
    }
    return(list(Zvals = Zvals, niters = length(tiles)))
}
filled.contour.points <- function(mcmcpath, longlat, plot.params, highlight) {
    if (is.null(highlight)) {
        return (NULL)
    }
    if (is.null(highlight$index)) {
        writeLines('Specify the indices of points to highlight')
        return (NULL)
    }            
    if (is.null(highlight$col)) highlight$col <- "red"
    if (is.null(highlight$cex)) highlight$cex <- 1
    if (is.null(highlight$pch)) highlight$pch <- 4
    if (is.null(highlight$lwd)) highlight$lwd <- 2

    graph <- read.graph(mcmcpath, longlat)
    coord <- graph$demes[graph$ipmap, ]
    
    index2coord <- match(highlight$index, seq(nrow(coord)))
    index2highlight <- which(!is.na(index2coord))
    index2coord <- index2coord[index2highlight]

    ## Check that there is at least one point to highlight
    if (!length(index2highlight)) { return(NULL) }
    
    coord <- coord[index2coord, ]
    coord <- matrix(coord, ncol = 2)
    ## Tranform the coordinates to a different projection if necessary
    if (!is.null(plot.params$proj.out)) {
        coord <- sp::SpatialPoints(coord, proj4string = CRS(plot.params$proj.in))
        coord <- sp::spTransform(coord, CRSobj = CRS(plot.params$proj.out))
    }                
    points(coord, 
           cex = highlight$cex[index2highlight], 
           col = highlight$col[index2highlight], 
           pch = highlight$pch[index2highlight], 
           lwd = highlight$lwd[index2highlight])
}
filled.contour.map <- function(mcmcpath, longlat, plot.params) {
    if (!is.null(plot.params$proj.in) && plot.params$add.map) {
        map <- rworldmap::getMap(resolution = "high")
        map <- sp::spTransform(map, CRSobj = CRS(plot.params$proj.out))
        plot(map, col = NA, border = plot.params$col.map, lwd = plot.params$lwd.map, add = TRUE)
    }
}
filled.contour.exterior <- function(mcmcpath, longlat, plot.params) {
    if (is.null(plot.params$proj.in)) {
        filled.contour.exterior.proj.unknown(mcmcpath, longlat, plot.params)
    } else {
        filled.contour.exterior.proj.known(mcmcpath, longlat, plot.params)
    }
}
filled.contour.axes <- function(mcmcpath, longlat, plot.params) {
    if (is.null(plot.params$proj.in)) {
        filled.contour.axes.proj.unknown(mcmcpath, longlat, plot.params)
    } else {
        filled.contour.axes.proj.known(mcmcpath, longlat, plot.params)
    }
}
check.habitat.outer.is.vald <- function(mcmcpath, longlat) {
    graph <- read.graph(mcmcpath, longlat)
    outer <- graph$outer
    ## Check that the geometry of the habitat outline is both Simple (gIsSimple)
    ## and Closed (end points intersect).
    x <- outer[, 1]
    y <- outer[, 2]
    l <- readWKT( paste("LINESTRING(", paste(paste(x, y), collapse = ", "), ")") )
    if (!rgeos::gIsRing(l, byid = FALSE)) {
        writeLines(paste0("Error (rgeos):\n", 
                          "The habitat geometry is not a valid ring (a ring is both simple and closed)\n", 
                          "Let the habitat be the rectangle defined by (xmin, ymin) and (xmax, yxmax)\n", 
                          "where xmin, xmax = range(longitude) and ymin, ymax = range(latitude)."))
        xmin <- min(outer[, 1])
        xmax <- max(outer[, 1])
        ymin <- min(outer[, 2])
        ymax <- max(outer[, 2])
        outer <- cbind(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin))
    }
    return(outer)
}
filled.contour.exterior.proj.unknown <- function(mcmcpath, longlat, plot.params) {
    outer <- check.habitat.outer.is.vald(mcmcpath, longlat)
    ## The filledContour fills a rectangular plot; now color the habitat exterior white.
    boundary <- sp::SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)), "1")))
    exterior <- rgeos::gDifference(rgeos::gEnvelope(boundary), boundary)    
    if (!is.null(exterior)) {
        plot(exterior, col = "white", border = "white", add = TRUE)
    }
    if (plot.params$add.outline) {
        plot(boundary, col = NA, border = plot.params$col.outline, lwd = plot.params$lwd.outline, add = TRUE)
    }
}
filled.contour.exterior.proj.known <- function(mcmcpath, longlat, plot.params) {
    outer <- check.habitat.outer.is.vald(mcmcpath, longlat)
    ## The filledContour fills a rectangular plot; now color the habitat exterior white.
    boundary <- sp::SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)), "1")), 
                                    proj4string = CRS(plot.params$proj.in))
    boundary <- sp::spTransform(boundary, CRSobj = CRS(plot.params$proj.out))
    exterior <- rgeos::gDifference(rgeos::gEnvelope(boundary), boundary)
    if (!is.null(exterior)) {
        plot(exterior, col = "white", border = "white", add = TRUE)
    }
    if (plot.params$add.outline) {
        plot(boundary, col = NA, border = plot.params$col.outline, lwd = plot.params$lwd.outline, add = TRUE)
    }
}
filled.contour.axes.proj.unknown <- function(mcmcpath, longlat, plot.params) {
    graph <- read.graph(mcmcpath, longlat)
    if (plot.params$add.grid) {
        segments <- list()
        for (e in 1:nrow(graph$edges)) {
            segments[[e]] <- sp::Line(graph$demes[graph$edges[e, ], ])
        }
        segments <- sp::SpatialLines(list(Lines(segments, ID = "a")))
        lines(segments, col = plot.params$col.grid, lwd = plot.params$lwd.grid)
    }
    if (plot.params$all.demes) {
        all.demes <- sp::SpatialPoints(graph$demes)
        cex.points <- plot.params$min.cex.demes
        points(all.demes, col = plot.params$col.demes, pch = plot.params$pch.demes, cex = cex.points)
    } else if (plot.params$add.demes) {
        observed.demes <- sp::SpatialPoints(graph$demes[graph$alpha, ])
        cex.points <- plot.params$min.cex.demes
        if (min(graph$sizes) < max(graph$sizes)) {
            cex.points <- cex.points +
                (plot.params$max.cex.demes - plot.params$min.cex.demes) *
                    (graph$sizes - min(graph$sizes)) / (max(graph$sizes) - min(graph$sizes))
        }
        points(observed.demes, col = plot.params$col.demes, pch = plot.params$pch.demes, cex = cex.points)
    }
}
filled.contour.axes.proj.known <- function(mcmcpath, longlat, plot.params) {
    graph <- read.graph(mcmcpath, longlat)
    if (plot.params$add.grid) {
        segments <- list()
        for (e in 1:nrow(graph$edges)) {
            segments[[e]] <- sp::Line(graph$demes[graph$edges[e, ], ])
        }
        segments <- sp::SpatialLines(list(Lines(segments, ID = "a")), proj4string = CRS(plot.params$proj.in))
        segments <- sp::spTransform(segments, CRSobj = CRS(plot.params$proj.out))
        lines(segments, col = plot.params$col.grid, lwd = plot.params$lwd.grid)
    }
    if (plot.params$all.demes) {
        all.demes <- sp::SpatialPoints(graph$demes[graph$alpha, ], proj4string = CRS(plot.params$proj.in))
        all.demes <- sp::spTransform(all.demes, CRSobj = CRS(plot.params$proj.out))
        cex.points <- plot.params$min.cex.demes
        points(all.demes, col = plot.params$col.demes, pch = plot.params$pch.demes, cex = cex.points)
    } else if (plot.params$add.demes) {
        observed.demes <- sp::SpatialPoints(graph$demes[graph$alpha, ], proj4string = CRS(plot.params$proj.in))
        observed.demes <- sp::spTransform(observed.demes, CRSobj = CRS(plot.params$proj.out))
        cex.points <- plot.params$min.cex.demes
        if (min(graph$sizes) < max(graph$sizes)) {
            cex.points <- cex.points +
                (plot.params$max.cex.demes - plot.params$min.cex.demes) *
                    (graph$sizes - min(graph$sizes)) / (max(graph$sizes) - min(graph$sizes))
        }
        points(observed.demes, col = plot.params$col.demes, pch = plot.params$pch.demes, cex = cex.points)
    }
}
null.eems.contour <- function(mcmcpath, longlat, plot.params, 
                              plot.xy = NULL, highlight.samples = NULL) {
    dimns <- read.dimns(mcmcpath, longlat)
    eems.colors <- plot.params$eems.colors
    num.levels <- length(eems.colors)
    eems.levels <- seq(from = -2.5, to = +2.5, length = num.levels+1)
    main.title <- ""
    key.title <- ""
    Z <- matrix(0, nrow = dimns$nxmrks, ncol = dimns$nymrks)
    rr <- flip(raster::raster(t(Z), 
                              xmn = dimns$xlim[1], xmx = dimns$xlim[2], 
                              ymn = dimns$ylim[1], ymx = dimns$ylim[2]), direction = 'y')
    if (!is.null(plot.params$proj.in)) {
        raster::projection(rr) <- CRS(plot.params$proj.in)
        rr <- raster::projectRaster(rr, crs = CRS(plot.params$proj.out))
    }
    myfilledContour(rr, col = eems.colors, levels = eems.levels, asp = 1, 
                    add.key = plot.params$add.colbar, 
                    key.axes = axis(4, tick = FALSE, hadj = 1, line = 3, cex.axis = 1.5), 
                    key.title = mtext(key.title, side = 3, cex = 1.5, line = 1.5, font = 1), 
                    add.title = plot.params$add.title, 
                    plot.title = mtext(text = main.title, side = 3, line = 0, cex = 1.5), 
                    plot.axes = {
                        filled.contour.exterior(mcmcpath, longlat, plot.params);
                        filled.contour.map(mcmcpath, longlat, plot.params);
                        plot.xy;
                        filled.contour.axes(mcmcpath, longlat, plot.params);
                        filled.contour.points(mcmcpath, longlat, plot.params, highlight.samples);
                    })
    return (list(eems.colors = eems.colors, eems.levels = eems.levels))
}
one.eems.contour <- function(mcmcpath, dimns, Zmean, longlat, plot.params, is.mrates, 
                             plot.xy = NULL, highlight.samples = NULL) {
    eems.colors <- plot.params$eems.colors
    num.levels <- length(eems.colors)
    if (is.mrates) {
        eems.levels <- eems.colscale(Zmean, num.levels, plot.params$m.colscale)
        main.title <- "Effective migration rates m : posterior mean"
        key.title <- expression(paste(log, "(", italic(m), ")", sep = ""))
    } else {
        eems.levels <- eems.colscale(Zmean, num.levels, plot.params$q.colscale)
        main.title <- "Effective diversity rates q : posterior mean"
        key.title <- expression(paste(log, "(", italic(q), ")", sep = ""))
    }
    rr <- flip(raster::raster(t(Zmean), 
                              xmn = dimns$xlim[1], xmx = dimns$xlim[2], 
                              ymn = dimns$ylim[1], ymx = dimns$ylim[2]), direction = 'y')
    if (!is.null(plot.params$proj.in)) {
        raster::projection(rr) <- CRS(plot.params$proj.in)
        rr <- raster::projectRaster(rr, crs = CRS(plot.params$proj.out))
    }
    myfilledContour(rr, col = eems.colors, levels = eems.levels, asp = 1, 
                    add.key = plot.params$add.colbar, 
                    key.axes = axis(4, tick = FALSE, hadj = 1, line = 3, cex.axis = 1.5), 
                    key.title = mtext(key.title, side = 3, cex = 1.5, line = 1.5, font = 1), 
                    add.title = plot.params$add.title, 
                    plot.title = mtext(text = main.title, side = 3, line = 0, cex = 1.5), 
                    plot.axes = {
                        filled.contour.exterior(mcmcpath, longlat, plot.params);
                        filled.contour.map(mcmcpath, longlat, plot.params);
                        plot.xy;
                        filled.contour.axes(mcmcpath, longlat, plot.params);
                        filled.contour.points(mcmcpath, longlat, plot.params, highlight.samples);
                    })
    return (list(eems.colors = eems.colors, eems.levels = eems.levels))
}
average.eems.contours <- function(mcmcpath, dimns, longlat, plot.params, is.mrates, 
                                  plot.xy = NULL, highlight.samples = NULL) {
    mcmcpath1 <- character( )
    if (is.mrates) {
        writeLines('Plotting effective migration surface (posterior mean of m rates)')
        for (path in mcmcpath) {
            Files <- paste0(path, c('/mcmcmtiles.txt', '/mcmcmrates.txt', 
                                    '/mcmcxcoord.txt', '/mcmcycoord.txt'))
            if (min(file.exists(Files))) { mcmcpath1 <- c(mcmcpath1, path) }
        }
        zero_mean <- plot.params$m.zero_mean
        log_scale <- plot.params$m.log_scale
    } else {
        writeLines('Plotting effective diversity surface (posterior mean of q rates)')
        for (path in mcmcpath) {
            Files <- paste0(path, c('/mcmcqtiles.txt', '/mcmcqrates.txt', 
                                    '/mcmcwcoord.txt', '/mcmczcoord.txt'))
            if (min(file.exists(Files))) { mcmcpath1 <- c(mcmcpath1, path) }
        }
        zero_mean <- plot.params$q.zero_mean
        log_scale <- plot.params$q.log_scale
    }
    mcmcpath <- mcmcpath1
    if (length(mcmcpath) == 0) { return(NULL) }
    Zmean <- matrix(0, dimns$nxmrks, dimns$nymrks)
    niters <- 0
    for (path in mcmcpath) {
        writeLines(path)
        rslt <- transform.rates(path, dimns, longlat, is.mrates, log_scale, zero_mean)
        Zmean <- Zmean + rslt$Zvals
        niters <- niters + rslt$niters
    }
    Zmean <- Zmean / niters
    rates.raster <- one.eems.contour(mcmcpath[1], dimns, Zmean, longlat, plot.params, is.mrates, 
                                     plot.xy = plot.xy, highlight.samples = highlight.samples)
    return(rates.raster)
}
## This function is mainly for testing purposes and will create on Voronoi diagram for each saved MCMC iteration
voronoi.diagram <- function(mcmcpath, dimns, longlat, plot.params, mcmc.iters = NULL, is.mrates = TRUE) {
    mcmcpath <- mcmcpath[1]
    writeLines('Plotting Voronoi tessellation of estimated effective rates')
    writeLines(mcmcpath)
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
    niters <- length(tiles)
    eems.colors <- plot.params$eems.colors 
    num.levels <- length(eems.colors)
    if (is.mrates) {
        main.title <- 'Effective migration rates m'
        eems.levels <- eems.colscale(rates, num.levels, plot.params$m.colscale)
    } else {
        main.title <- 'Effective diversity rates q'
        eems.levels <- eems.colscale(rates, num.levels, plot.params$q.colscale)
    }
    if (is.null(mcmc.iters)) {
        mcmc.iters <- seq(niters)
    }
    L <- length(eems.levels)
    max.levels <- max(eems.levels)
    min.levels <- min(eems.levels)
    count <- 0
    for (i in 1:niters) {
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        count <- count + now.tiles
        if (!(i %in% mcmc.iters)) { next }
        ## Standardize the log-transformed rates, without taking into account
        ## the relative size of the tiles (this is hard to do without a grid)
        now.rates <- now.rates - mean(now.rates)
        now.rates <- ifelse(now.rates > max.levels, max.levels, now.rates)
        now.rates <- ifelse(now.rates < min.levels, min.levels, now.rates)
        now.main.title <- paste0(main.title, " : iteration ", i, " (after burn-in and thinning)")
        if (!plot.params$add.title) { now.main.title <- "" }
        plot(0, 0, type = "n", xlim = dimns$xlim, ylim = dimns$ylim, asp = 1, 
             axes = FALSE, xlab = "", ylab = "", main = now.main.title)
        if (now.tiles == 1) {
            ## There is only one tile
            tile.color <- eems.colors[round(L / 2)]
            polygon(dimns$xlim, dimns$ylim, col = tile.color, border = "white")
        } else {
            ## Plot each tile in turn (as a polygon)
            Voronoi <- deldir::deldir(now.xseed, now.yseed, rw = c(dimns$xlim, dimns$ylim))
            tilelist <- deldir::tile.list(Voronoi)
            for (c in 1:now.tiles) {
                tile.color <- eems.colors[ findInterval(now.rates[c], eems.levels, all.inside = TRUE) ]
                polygon(tilelist[[c]]$x, tilelist[[c]]$y, col = tile.color, border = "white")
            }
            filled.contour.axes(mcmcpath, longlat, plot.params)
        }
        if (plot.params$add.seeds) {
            points(now.xseed, now.yseed,
                   pch = plot.params$pch.seeds, cex = plot.params$cex.seeds,
                   col = plot.params$col.seeds, lwd = 3)
        }
    }
    return(list(colors = eems.colors, levels = eems.levels))    
}
plot.logposterior <- function(mcmcpath) {
    writeLines('Plotting posterior probability trace')
    mcmcpath1 <- character()
    for (path in mcmcpath) {
        if (file.exists(paste0(path, '/mcmcpilogl.txt'))) {
            mcmcpath1 <- c(mcmcpath1, path)
        }
    }
    mcmcpath <- mcmcpath1
    nchains <- length(mcmcpath)
    if (nchains == 0) { return(NULL) }
    ## Here colors = RColorBrewer::brewer.pal(12, "Paired") + "black"
    colors <- rep_len(c("#000000", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                        "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"), length.out = nchains)
    ltypes <- rep_len(1:4, length.out = nchains)
    posteriors <- list()
    yrange <- NULL
    niters <- NULL
    for (i in 1:nchains) {
        path <- mcmcpath[i]; writeLines(path)
        pilogl <- scan(paste0(path, '/mcmcpilogl.txt'), quiet = TRUE)
        pilogl <- matrix(pilogl, ncol = 2, byrow = TRUE)
        posterior <- pilogl[, 1] + pilogl[, 2]
        posteriors[[i]] <- posterior
        yrange <- range(c(yrange, posterior))
        niters <- max(niters, length(posterior))
    }
    plot(c(1, niters), yrange, type = "n", xlab = "MCMC iteration  (after burn-in and thinning)", ylab = "log posterior")
    if (nchains == 1) {
        mtext(side = 3, line = 2, cex = 1.3, text = "Has the MCMC chain converged?")
        mtext(side = 3, line = 0.5, cex = 1, text = "If not, restart EEMS and/or increase numMCMCIter, numBurnIter, numThinIter")
    } else {
        mtext(side = 3, line = 2, cex = 1.3, text = "Have the MCMC chains converged?")
        mtext(side = 3, line = 0.5, cex = 1, text = "If not, restart EEMS and/or increase numMCMCIter, numBurnIter, numThinIter")
    }
    for (i in 1:nchains) {
        posterior <- posteriors[[i]]
        niters <- length(posterior)
        lines(1:niters, posterior, col = colors[i], lty = ltypes[i], lwd = 2)
    }
    legend("topright", legend = 1:nchains, col = colors[1:nchains], 
           lty = ltypes[1:nchains], lwd = 2, bty = "n", inset = c(-0.12, 0), cex = 0.5)
}
geo.distm <- function(coord, longlat, plot.params) {
    if (!longlat) {
        coord <- coord[, c(2, 1)]
    }
    if (!is.null(plot.params$proj.in)) {
        coord <- sp::SpatialPoints(coord, proj4string = CRS(plot.params$proj.in))
        coord <- sp::spTransform(coord, CRSobj = CRS("+proj=longlat +datum=WGS84"))
    }
    Dist <- geosphere::distm(coord, fun = distHaversine) / 1000
    Dist <- Dist[upper.tri(Dist, diag = FALSE)]
    return(Dist)
}
dist.scatterplot <- function(mcmcpath, longlat, plot.params, 
                             remove.singletons = TRUE, add.abline = FALSE, 
                             highlight = NULL) {
    writeLines('Plotting average dissimilarities within and between demes')
    mcmcpath1 <- character()
    for (path in mcmcpath) {
        if (file.exists(paste0(path, '/rdistJtDobsJ.txt'))&&
            file.exists(paste0(path, '/rdistJtDhatJ.txt'))&&
            file.exists(paste0(path, '/rdistoDemes.txt'))) {
            mcmcpath1 <- c(mcmcpath1, path)
        }
    }
    mcmcpath <- mcmcpath1
    nchains <- length(mcmcpath)
    if (nchains == 0) { return(NULL) }
    ## List of observed demes, with number of samples taken collected
    ## Each row specifies: x coordinate, y coordinate, n samples
    oDemes <- scan(paste0(mcmcpath[1], '/rdistoDemes.txt'), quiet = TRUE)
    oDemes <- matrix(oDemes, ncol = 3, byrow = TRUE)
    Sizes <- oDemes[, 3]
    nPops <- nrow(oDemes)
    Demes <- seq(nPops)
    JtDobsJ <- matrix(0, nPops, nPops)
    JtDhatJ <- matrix(0, nPops, nPops)
    for (path in mcmcpath) {
        writeLines(path)
        oDemes1 <- scan(paste0(path, '/rdistoDemes.txt'), quiet = TRUE)
        oDemes1 <- matrix(oDemes1, ncol = 3, byrow = TRUE)
        if ( sum(dim(oDemes) != dim(oDemes1)) || sum(oDemes != oDemes1)) {
            plot.params$add.demes <- TRUE
            plot.params$add.grid <- TRUE
            plot.params$add.title <- TRUE
            plot.params$add.colbar <- FALSE
            null.eems.contour(mcmcpath[1], longlat, plot.params)
            mtext(side = 3, line = 2, cex = 1.3, text = 'EEMS results for at least two different population grids')
            null.eems.contour(    path   , longlat, plot.params)
            mtext(side = 3, line = 2, cex = 1.3, text = 'EEMS results for at least two different population grids')
            writeLines('EEMS results for at least two different population grids')
            return(NULL)
        }
        JtDobsJ <- JtDobsJ + as.matrix(read.table(paste0(path, '/rdistJtDobsJ.txt'), header = FALSE))
        JtDhatJ <- JtDhatJ + as.matrix(read.table(paste0(path, '/rdistJtDhatJ.txt'), header = FALSE))
    }
    JtDobsJ <- JtDobsJ/nchains
    JtDhatJ <- JtDhatJ/nchains
    colnames(JtDobsJ) <- Demes
    rownames(JtDobsJ) <- Demes
    colnames(JtDhatJ) <- Demes
    rownames(JtDhatJ) <- Demes
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
    Deme1 <- matrix(Demes, nrow = nPops, ncol = nPops)
    Deme2 <- t(Deme1)
    tempi <- matrix(Sizes, nPops, nPops)
    Small <- pmin(tempi, t(tempi))
    Small <- Small[upper.tri(Small, diag = FALSE)]    
    alpha <- Deme1[upper.tri(Deme1, diag = FALSE)]
    beta  <- Deme2[upper.tri(Deme2, diag = FALSE)]
    ## Under pure isolation by distance, we expect the genetic dissimilarities
    ## between demes increase with the geographic distance separating them
    Dist = geo.distm(oDemes[, 1:2], longlat, plot.params)
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##     
    if (sum(Sizes > 1) < 2) {
        ## Sizes(alpha) > 1 means that there are at least two individuals assigned to deme alpha
        JtDJ.component <- data.frame(alpha.x = oDemes[, 1][alpha], 
                                     alpha.y = oDemes[, 2][alpha], 
                                     beta.x = oDemes[, 1][beta], 
                                     beta.y = oDemes[, 2][beta], 
                                     size = Small, 
                                     col = sub.scattercols(Small), 
                                     fitted = JtDhatJ[upper.tri(JtDhatJ, diag = FALSE)], 
                                     obsrvd = JtDobsJ[upper.tri(JtDobsJ, diag = FALSE)], 
                                     stringsAsFactors = FALSE)
        sub.scatterplot("JtDJ", JtDJ.component, remove.singletons, add.abline)
        writeLines('There should be at least two observed demes to plot pairwise dissimilarities')
        return (NULL)
    }
    out1 <- JtDJ2BandW(JtDobsJ, Sizes)
    Wobs <- out1$W
    Bobs <- out1$B
    out2 <- JtDJ2BandW(JtDhatJ)
    What <- out2$W
    Bhat <- out2$B
    B.component <- data.frame(alpha.x = oDemes[, 1][alpha], 
                              alpha.y = oDemes[, 2][alpha], 
                              beta.x = oDemes[, 1][beta], 
                              beta.y = oDemes[, 2][beta], 
                              fitted = Bhat, 
                              obsrvd = Bobs, 
                              size = Small, 
                              col = sub.scattercols(Small), 
                              stringsAsFactors = FALSE)
    W.component <- data.frame(alpha.x = oDemes[, 1][Demes], 
                              alpha.y = oDemes[, 2][Demes], 
                              fitted = What, 
                              obsrvd = Wobs, 
                              size = Sizes, 
                              col = sub.scattercols(Sizes), 
                              stringsAsFactors = FALSE)
    G.component <- data.frame(alpha.x = oDemes[, 1][alpha], 
                              alpha.y = oDemes[, 2][alpha], 
                              beta.x = oDemes[, 1][beta], 
                              beta.y = oDemes[, 2][beta], 
                              fitted = Dist, 
                              obsrvd = Bobs, 
                              size = Small, 
                              col = sub.scattercols(Small), 
                              stringsAsFactors = FALSE)
    highlight.demes <- FALSE
    if (!is.null(highlight) &&
        !is.null(highlight$x) &&
        !is.null(highlight$y)) {
        ## Highligh some demes
        if (is.null(highlight$col)) highlight$col <- "red"
        if (is.null(highlight$cex)) highlight$cex <- 1 
        if (is.null(highlight$pch)) highlight$pch <- 4 
        if (is.null(highlight$lwd)) highlight$lwd <- 2 
        ## In the between-deme component, each point is associated with two demes:
        ## alpha and beta
        highlight <- dplyr::select(highlight, x, y, col, cex, pch, lwd)
        highlight <- dplyr::rename(highlight, alpha.col = col)
        highlight <- dplyr::rename(highlight, alpha.cex = cex)
        highlight <- dplyr::rename(highlight, alpha.pch = pch)
        highlight <- dplyr::rename(highlight, alpha.lwd = lwd)
        B.component <- dplyr::left_join(B.component, highlight, by = c("alpha.x" = "x", "alpha.y" = "y"))
        W.component <- dplyr::left_join(W.component, highlight, by = c("alpha.x" = "x", "alpha.y" = "y"))
        G.component <- dplyr::left_join(G.component, highlight, by = c("alpha.x" = "x", "alpha.y" = "y"))
        highlight <- dplyr::rename(highlight, beta.col = alpha.col)
        highlight <- dplyr::rename(highlight, beta.cex = alpha.cex)
        highlight <- dplyr::rename(highlight, beta.pch = alpha.pch)
        highlight <- dplyr::rename(highlight, beta.lwd = alpha.lwd)
        B.component <- dplyr::left_join(B.component, highlight, by = c("beta.x" = "x", "beta.y" = "y"))
        G.component <- dplyr::left_join(G.component, highlight, by = c("beta.x" = "x", "beta.y" = "y"))
        if (sum(!is.na(B.component$alpha.col)) ||
            sum(!is.na(B.component$beta.col))) { highlight.demes <- TRUE }
    }    
    if (highlight.demes) {
        B.component.alpha <- dplyr::select(B.component, fitted, obsrvd, size, starts_with("alpha"))
        B.component.beta <- dplyr::select(B.component, fitted, obsrvd, size, starts_with("beta"))
        B.component.alpha <- dplyr::rename(B.component.alpha, col = alpha.col, cex = alpha.cex, pch = alpha.pch, lwd = alpha.lwd)
        B.component.beta <- dplyr::rename(B.component.beta, col = beta.col, cex = beta.cex, pch = beta.pch, lwd = beta.lwd)
        sub.scatterplot("Between", B.component, remove.singletons, add.abline, 
                        subtitle = expression(paste("Highlight ", alpha, " in (", alpha, ", ", beta, ")")))
        sub.scatterplot("Between", B.component.alpha, remove.singletons, add.abline, add = TRUE)
        sub.scatterplot("Between", B.component, remove.singletons, add.abline, 
                        subtitle = expression(paste("Highlight ", beta, " in (", alpha, ", ", beta, ")")))
        sub.scatterplot("Between", B.component.beta, remove.singletons, add.abline, add = TRUE)
    } else {
        sub.scatterplot("Between", B.component, remove.singletons, add.abline)
    }
    if (highlight.demes) {
        W.component.alpha <- dplyr::select(W.component, fitted, obsrvd, size, starts_with("alpha"))
        W.component.alpha <- dplyr::rename(W.component.alpha, col = alpha.col, cex = alpha.cex, pch = alpha.pch, lwd = alpha.lwd)
        sub.scatterplot("Within", W.component, remove.singletons, add.abline)
        sub.scatterplot("Within", W.component.alpha, remove.singletons, add.abline, add = TRUE)
    } else {    
        sub.scatterplot("Within", W.component, remove.singletons, add.abline)
    }
    if (highlight.demes) {
        G.component.alpha <- dplyr::select(G.component, fitted, obsrvd, size, starts_with("alpha"))
        G.component.beta <- dplyr::select(G.component, fitted, obsrvd, size, starts_with("beta"))
        G.component.alpha <- dplyr::rename(G.component.alpha, col = alpha.col, cex = alpha.cex, pch = alpha.pch, lwd = alpha.lwd)
        G.component.beta <- dplyr::rename(G.component.beta, col = beta.col, cex = beta.cex, pch = beta.pch, lwd = beta.lwd)
        sub.scatterplot("GeoDist", G.component, remove.singletons, add.abline = FALSE, 
                        subtitle = expression(paste("Highlight ", alpha, " in (", alpha, ", ", beta, ")")))
        sub.scatterplot("GeoDist", G.component.alpha, remove.singletons, add = TRUE)
        sub.scatterplot("GeoDist", G.component, remove.singletons, add.abline = FALSE, 
                        subtitle = expression(paste("Highlight ", beta, " in (", alpha, ", ", beta, ")")))
        sub.scatterplot("GeoDist", G.component.beta, remove.singletons, add = TRUE)
    } else {
        sub.scatterplot("GeoDist", G.component, remove.singletons, add.abline = FALSE)
    }
    return (list(B.component = B.component, W.component = W.component))
}
heatmap.resid <- function(datapath, mcmcpath) {
    mcmcpath <- mcmcpath
    nchains <- length(mcmcpath)
    writeLines('Heatmap of n-by-n matrix of residuals (observed - fitted):')
    Diffs <- as.matrix(read.table(paste0(datapath, ".diffs"), header = FALSE))
    nIndiv <- nrow(Diffs)
    nChains <- 0
    Delta <- matrix(0, nIndiv, nIndiv)
    for (path in mcmcpath) {
        if (file.exists(paste0(path, '/rdistJtDobsJ.txt')) &&
            file.exists(paste0(path, '/rdistJtDhatJ.txt')) &&
            file.exists(paste0(path, '/rdistoDemes.txt'))) {
            writeLines(path)
            nChains <- nChains + 1
            ipmap <- scan(paste0(path, "/ipmap.txt"), what = numeric(), quiet = TRUE)
            Sizes <- as.numeric(table(ipmap))
            n <- length(ipmap) ## number of samples
            o <- length(Sizes) ## number of observed demes
            J <- Matrix::spMatrix(n, o, i = seq(n), j = ipmap, x = rep(1, n)) ## indicator matrix
            J <- as.matrix(J)
            JtDhatJ <- as.matrix(read.table(paste0(path, '/rdistJtDhatJ.txt'), header = FALSE))
            Delta <- Delta + J %*% JtDhatJ %*% t(J)
        }
    }
    if (nChains == 0) { return(NULL) }
    Delta <- Delta / nChains
    Delta <- Delta - diag(diag(Delta))
    resid <- Diffs - Delta
    diag(resid) <- NA ## The diagonal entries would always be zero
    return (resid)
}
mypalette <- function(colors = NULL, colscale = NULL, n = 299) {
    if (!is.null(colors)) {
        colors <- colorRampPalette(colors)(n = n)
    } else {
        ## Here colors = RColorBrewer::brewer.pal(9, "Reds")
        colors <- colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272",
                                     "#FB6A4A",
                                     "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(n = n)
    }
    if (!is.null(colscale)) {
        if (length(colscale) != (n + 1)) {
            levels <- seq(from = min(colscale), to = max(colscale), length.out = n + 1)
        } else {
            levels <- sort(colscale)
        }
    }
    return(list(colors = colors, levels = levels))
}
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
    return (palette)
}
## By default, all figures are saved as bitmap PNG images. However, 
## it is straightforward to use another format (Here the alternative is PDF)
save.graphics <- function(plotpath, plot.params) {
    if (plot.params$out.png) {
        bitmap(paste0(plotpath, '%02d.png'), type = 'png16m', res = plot.params$res, 
               height = plot.params$height, width = plot.params$width, units = 'in')
    } else {
        pdf(paste0(plotpath, '%02d.pdf'), 
              height = plot.params$height, width = plot.params$width, onefile = FALSE)
    }
}

## Redefine the 'filledContour' function from the 'raster' package
## since I don't like how the legend looks.
## For now, I have only changed how the legend/bar is plotted, 
## the change itself is in the 'myfilled.contour' function
myfilledContour <- function (x, y = 1, maxpixels = 1e+05, add.key = TRUE, ...) {
    if (nlayers(x) > 1) {
        y <- min(max(1, y), nlayers(x))
        x <- raster(x, y)
    }
    x <- sampleRegular(x, maxpixels, asRaster = TRUE, useGDAL = TRUE)
    X <- xFromCol(x, 1:ncol(x))
    Y <- yFromRow(x, nrow(x):1)
    Z <- t(matrix(getValues(x), ncol = x@ncols, byrow = TRUE)[nrow(x):1, ])
    myfilled.contour(x = X, y = Y, z = Z, add.key = add.key, 
                     axes = FALSE, frame.plot = FALSE, ...)
}
## I have only changed the legend/bar to remove the black border.
myfilled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, 
                              xlim = range(x, finite = TRUE), 
                              ylim = range(y, finite = TRUE), 
                              zlim = range(z, finite = TRUE), 
                              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                              axes = TRUE, frame.plot = axes, add.title = FALSE, add.key = TRUE, 
                              ...) {

    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")
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
        ## Change: add 'border = NA' and remove 'box()'
        plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
        rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
        if (missing(key.axes)) {
            if (axes)
                axis(4)
        }
        else key.axes
        ## box()
        if (!missing(key.title))
            key.title
    }
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    ## Change: make the plot region almost as big as the figure region
    ##         but leave space for the main title
    if (add.title) {
        par(plt = c(0.02, 0.98, 0.02, 0.95))
    } else {
        par(plt = c(0.02, 0.98, 0.02, 0.98))
    }
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (add.title) {
        if (missing(plot.title))
            title(...)
        else plot.title
    }
    invisible()
}
myfilled.legend <-
    function (levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes, ...) 
{

    plot.new( )
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
    if (!missing(key.axes)) {
        key.axes
    } else {
        axis(4, col.axis = fg.col, tick = FALSE)
    }
    if (!missing(key.title)) {
        key.title
    }
}
load.required.package <- function(package, required.by) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(paste0("'", required.by, "' requires the '", package, "' package. Please install it first."))
    } else {
        writeLines(paste0("Loading ", package, " (required by ", required.by, ")"))
    }
}

#' A function to plot effective migration and diversity surfaces from EEMS output
#'
#' Given a vector of EEMS output directories, this function generates five figures to visualize EEMS results:
#' \itemize{
#'  \item \code{plotpath}-mrates01 (effective migration surface)
#'  \item \code{plotpath}-qrates01 (effective diversity surface)
#'  \item \code{plotpath}-rdist01 (between-demes component of genetic dissimilarity)
#'  \item \code{plotpath}-rdist02 (within-demes component of genetic dissimilarity)
#'  \item \code{plotpath}-rdist03 (between-demes dissimilarities vs geographic distances)
#'  \item \code{plotpath}-pilogl01 (posterior probability trace)
#' }
#' The latter figures are helpful in checking that the MCMC sampler has converged (the trace plot \code{pilogl01}) and that the EEMS model fits the data well (the scatter plots of genetic dissimilarities \code{rdist0*}). \code{eems.plots} will work with a single EEMS output directory but it is better to run EEMS several times, randomly initializing the MCMC chain each time. In other words, it is a good idea to simulate several realizations of the Markov chain, each realization starting with a different value of the EEMS parameters.
#' @param mcmcpath A vector of EEMS output directories, for the same dataset. Warning: There is minimal checking that the given  directories are for the same dataset.
#' @param plotpath The full path and the file name for the graphics to be generated. 
#' @param longlat A logical value indicating whether the coordinates are given as pairs (longitude, latitude) or (latitude, longitude).
#' @param plot.width, plot.height The width and height of the graphics region for the two rate contour plots, in inches. The default values are both 7.
#' @param out.png A logical value indicating whether to generate output graphics as PNGs (the default) or PDFs.
#' @param res Resolution, in dots per inch; used only if \code{out.png} is set to TRUE. The default is 600.
#' @param xpd A logical value indicating whether to clip plotting to the figure region (\code{xpd} = TRUE, which is the default) or clip plotting to the plot region (\code{xpd} = FALSE).
#' @param add.grid A logical value indicating whether to add the population grid or not.
#' @param col.grid The color of the population grid. Defaults to \code{gray80}.
#' @param lwd.grid The line width of the population grid. Defaults to 1.
#' @param add.outline A logical value indicating whether to add the habitat outline or not.
#' @param col.outline The color of the habitat outline. Defaults to \code{white}.
#' @param lwd.outline The line width of the habitat outline. Defaults to 2.
#' @param add.demes A logical value indicating whether to add the observed demes or not.
#' @param col.demes The color of the demes. Defaults to \code{black}.
#' @param pch.demes The symbol, specified as an integer, or the character to be used for plotting the demes. Defaults to 19.
#' @param min.cex.demes, max.cex.demes The minimum and the maximum size of the deme symbol/character. Defaults to 1 and 3, respectively. If \code{max.cex.demes} > \code{min.cex.demes}, then demes with more samples also have bigger size: the deme with the fewest samples has size \code{min.cex.demes} and the deme with the most samples has size \code{max.cex.demes}.
#' @param projection.in, projection.out The input and the output cartographic projections, each specified as a PROJ.4 string. Requires the \code{rgdal} package.
#' @param add.map A logical value indicating whether to add a high-resolution geographic map. Requires the \code{rworldmap} and \code{rworldxtra} packages. It also requires that \code{projection.in} is specified.
#' @param col.map The color of the geographic map. Default is \code{gray60}.
#' @param lwd.map The line width of the geographic map. Defaults to 2.
#' @param eems.colors The EEMS color scheme as a vector of colors, ordered from low to high. Defaults to a DarkOrange to Blue divergent palette with six orange shades, white in the middle, six blue shades. Acknowledgement: The default color scheme is adapted from the \code{dichromat} package.
#' @param m.colscale, q.colscale A fixed range for log10-transformed migration and diversity rates, respectively. If the estimated rates fall outside the specified range, then the color scale is ignored. By default, no range is specified for either type of rates.
#' @param add.colbar A logical value indicating whether to add the color bar (the key that shows how colors map to rates) to the right of the plot. Defaults to TRUE.
#' @param remove.singletons Remove demes with a single observation from the diagnostic scatterplots. Defaults to TRUE.
#' @param add.abline Add the line \code{y = x} to the diagnostic scatterplots of observed vs fitted genetic dissimilarities.
#' @param add.title A logical value indicating whether to add the main title in the contour plots. Defaults to TRUE.
#' @references Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
#' @export
#' @examples
#'
#' ## Use the provided example or supply the path to your own EEMS run.
#' eems.results.to.plot = paste0(path.package("rEEMSplots"), "/extdata/EEMS-example")
#' name.figures.to.save = "EEMS-example-rEEMSplots"
#' 
#' ## Produce the five EEMS figures, with default values for all optional parameters.
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-default"), 
#'            longlat = TRUE)
#'
#' ## Flip the x and y axis, i.e., assume that the x coordinate is the latitude
#' ## and the y coordinate is the longitude.
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-axes-flipped"), 
#'            longlat = FALSE)
#'
#' ## Generate PNG figures with height 9 inches, width 8 inches
#' ## and resolution 600 dots per inch.
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-demes-and-edges"), 
#'            longlat = TRUE, 
#'            plot.height = 8, 
#'            plot.width = 7, 
#'            res = 600, 
#'            out.png = TRUE)
#'
#' ## Generate PDF figures with height 9 inches and width 8 inches.
#' ## The resolution option, res, is ignored.
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-output-PDFs"), 
#'            longlat = TRUE, 
#'            plot.height = 8, 
#'            plot.width = 7, 
#'            res = 600, 
#'            out.png = FALSE)
#' 
#' ## Choose somewhat impractical colors and shapes for the outline, the grid and the demes.
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-demes-and-edges"), 
#'            longlat = TRUE, 
#'            add.grid = TRUE, 
#'            col.grid = "gray90", 
#'            lwd.grid = 2, 
#'            add.outline = TRUE, 
#'            col.outline = "blue", 
#'            lwd.outline = 5, 
#'            add.demes = TRUE, 
#'            col.demes = "red", 
#'            pch.demes = 5, 
#'            min.cex.demes = 0.5, 
#'            max.cex.demes = 1.5)
#'
#' library(rgdal)
#' 
#' ## Produce contour plots in the Mercator projection (used by Google Maps)
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-cartographic-projections"), 
#'            longlat = TRUE, 
#'            projection.in = "+proj=longlat +datum=WGS84", 
#'            projection.out = "+proj=merc +datum=WGS84")
#'
#' library(rworldmap)
#' library(rworldxtra)
#' 
#' ## Add a high-resolution geographic map
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-geographic-map"), 
#'            longlat = TRUE, 
#'            projection.in = "+proj=longlat +datum=WGS84", 
#'            projection.out = "+proj=merc +datum=WGS84", 
#'            add.map = TRUE, 
#'            col.map = "black", 
#'            lwd.map = 5)
#'
#' library(RColorBrewer)
#'
#' ## Use a divergent Red to Blue color scheme from the RColorBrewer package
#' ## instead of the default DarkOrange to Blue color scheme.
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-new-eems-colors"), 
#'            longlat = TRUE, 
#'            projection.in = "+proj=longlat +datum=WGS84", 
#'            projection.out = "+proj=merc +datum=WGS84", 
#'            eems.colors = brewer.pal(11, "RdBu"))
#'
#' ## Specify the color scales for the migration and diversity rates 
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-fix-colscales"), 
#'            longlat = TRUE, 
#'            add.outline = TRUE, 
#'            col.outline = "gray", 
#'            m.colscale = c(-3, 3), 
#'            q.colscale = c(-0.3, +0.3))
#'
#' ## In the diagnostic scatter plots, highlight demes in the South in red
#' ## All demes in the population grid, as (longitude, latitude) pairs
#' demes = read.table(paste0(eems.results.to.plot, "/demes.txt"))
#' demes.long = demes[, 1]
#' demes.lat = demes[, 2]
#' ## The demes in the South, defined as the region below the 10th parallel south
#' south = demes[demes.lat < -10, ]
#' highlight.demes = data.frame(
#'     x = south[, 1], 
#'     y = south[, 2], 
#'     col = "red", 
#'     stringsAsFactors = FALSE)
#' eems.plots(mcmcpath = eems.results.to.plot, 
#'            plotpath = paste0(name.figures.to.save, "-highlight-demes"), 
#'            longlat = TRUE, 
#'            highlight.demes = highlight.demes)

eems.plots <- function(mcmcpath, 
                       plotpath, 
                       longlat, 
                       
                       ## Properties of the figures to generate
                       plot.width = 7, 
                       plot.height = 7, 
                       out.png = TRUE, 
                       res = 600, 
                       xpd = TRUE, 

                       ## Properties of the population grid
                       add.grid = FALSE, 
                       col.grid = "gray80", 
                       lwd.grid = 1, 

                       ## Properties of the observed demes
                       add.demes = FALSE, 
                       col.demes = "black", 
                       pch.demes = 19, 
                       min.cex.demes = 1, 
                       max.cex.demes = 3, 

                       ## Properties of the habitat outline
                       add.outline = FALSE, 
                       col.outline = "gray90", 
                       lwd.outline = 2, 

                       ## Properties of the cartographic projection
                       projection.in = NULL, 
                       projection.out = NULL, 

                       ## Properties of the geographic map
                       add.map = FALSE, 
                       col.map = "gray60", 
                       lwd.map = 2, 

                       ## Color palette
                       eems.colors = NULL, 

                       ## Properties of the color key
                       add.colbar = TRUE, 
                       m.colscale = NULL, 
                       q.colscale = NULL, 

                       ## Remove demes with a single observation
                       remove.singletons = TRUE, 

                       ## Add the line y = x to scatter plots
                       add.abline = FALSE, 

                       ## Highlight some points
                       highlight.demes = NULL, 
                       
                       ## Extra options
                       add.title = TRUE, 
                       m.plot.xy = NULL, 
                       q.plot.xy = NULL,

                       m.zero_mean = TRUE,
                       q.zero_mean = TRUE,
                       m.log_scale = TRUE,
                       q.log_scale = TRUE) {
    
    plot.params <- list(eems.colors = eems.colors, m.colscale = m.colscale, q.colscale = q.colscale, 
                        add.map = add.map, add.grid = add.grid, add.outline = add.outline, add.demes = add.demes, 
                        col.map = col.map, col.grid = col.grid, col.outline = col.outline, col.demes = col.demes, 
                        lwd.map = lwd.map, lwd.grid = lwd.grid, lwd.outline = lwd.outline, pch.demes = pch.demes, 
                        min.cex.demes = min.cex.demes, proj.in = projection.in, add.colbar = add.colbar,
                        max.cex.demes = max.cex.demes, proj.out = projection.out, add.title = add.title,
                        m.zero_mean = m.zero_mean, m.log_scale = m.log_scale, 
                        q.zero_mean = q.zero_mean, q.log_scale = q.log_scale)
    plot.params <- check.plot.params(plot.params)

    ## A vector of EEMS output directories, for the same dataset.
    ## Assume that if eemsrun.txt exists, then all EEMS output files exist.
    mcmcpath <- mcmcpath[file.exists(paste0(mcmcpath, '/eemsrun.txt'))]
    if (!length(mcmcpath))
        stop('Please provide at least one existing EEMS output directory, mcmcpath')

    dimns <- read.dimns(mcmcpath[1], longlat)
    
    save.params <- list(height = plot.height, width = plot.width, res = res, out.png = out.png)

    writeLines('Processing the following EEMS output directories :')
    writeLines(mcmcpath)
    
    ## Plot filled contour of estimated effective migration rates
    save.graphics(paste0(plotpath, '-mrates'), save.params)
    par(las = 1, font.main = 1, xpd = xpd)
    mrates.raster <- average.eems.contours(mcmcpath, dimns, longlat, plot.params,
                                           is.mrates = TRUE, plot.xy = m.plot.xy)
    dev.off( )

    ## Plot filled contour of estimated effective diversity rates
    save.graphics(paste0(plotpath, '-qrates'), save.params)
    par(las = 1, font.main = 1, xpd = xpd)
    qrates.raster <- average.eems.contours(mcmcpath, dimns, longlat, plot.params,
                                           is.mrates = FALSE, plot.xy = q.plot.xy)
    dev.off( )

    if (!add.colbar) {

        if (out.png) {
            save.params$height <- 6
            save.params$width <- 1.5
        } else {
            save.params$height <- 12
            save.params$width <- 3
        }
        
        save.graphics(paste0(plotpath, '-mkey'), save.params)
        par(las = 1, font.main = 1, xpd = xpd, mar = c(0, 1, 5, 8))
        myfilled.legend(col = mrates.raster$eems.colors, 
                        levels = mrates.raster$eems.levels, 
                        key.axes = axis(4, tick = FALSE, hadj = 1, line = 4, cex.axis = 2), 
                        key.title = mtext(expression(paste(log, "(", italic(m), ")", sep = "")), 
                            side = 3, cex = 2.5, line = 1.5, font = 1))
        dev.off( )
        save.graphics(paste0(plotpath, '-qkey'), save.params)
        par(las = 1, font.main = 1, xpd = xpd, mar = c(0, 1, 5, 8))
        myfilled.legend(col = qrates.raster$eems.colors, 
                        levels = qrates.raster$eems.levels, 
                        key.axes = axis(4, tick = FALSE, hadj = 1, line = 6, cex.axis = 2), 
                        key.title = mtext(expression(paste(log, "(", italic(q), ")", sep = "")), 
                            side = 3, cex = 2.5, line = 1.5, font = 1))
        dev.off( )

    }
    
    save.params$height <- 6
    save.params$width <- 6.5

    if (!is.null(highlight.demes)) {
        load.required.package(package = 'dplyr', required.by = 'dist.scatterplot')
    }

    ## Plot scatter plots of observed vs fitted genetic differences
    save.graphics(paste0(plotpath, '-rdist'), save.params)
    par(las = 1, font.main = 1, mar = c(5, 5, 4, 2) + 0.1)
    B.and.W <- dist.scatterplot(mcmcpath, longlat, plot.params, 
                                remove.singletons = remove.singletons,
                                add.abline = add.abline, highlight = highlight.demes)
    dev.off( )

    B.component <- B.and.W$B.component
    W.component <- B.and.W$W.component
    save(B.component, W.component, file = paste0(plotpath, '-rdist.RData'))
    
    ## Plot trace plot of posterior probability to check convergence
    save.graphics(paste0(plotpath, '-pilogl'), save.params)
    par(las = 0, font.main = 1, mar = c(5, 5, 4, 5) + 0.1, xpd = TRUE)
    plot.logposterior(mcmcpath)
    dev.off( )
    
}

#' A function to plot Voronoi diagrams of effective migration and diversity rates
#'
#' Given one EEMS output directory, this function produces a series of Voronoi diagrams drawn from the posterior distribution of the rate parameters. There is one series of tessellations for the effective migration rates and another series for the effective diversity rates. Both series contain one figure for each saved MCMC iteration, after burn-in and thinning. Warning: \code{eems.voronoi} potentially generates a large number of figures.
#' @param mcmcpath An EEMS output directory. If mcmcpath is a vector of directories, then only the first (existing) directory is used.
#' @param plotpath The full path and the file name for the graphics to be generated.
#' @param longlat A logical value indicating whether the coordinates are given as pairs (longitude, latitude) or (latitude, longitude).
#' @param plot.width, plot.height The width and height of the graphics region for the two rate contour plots, in inches. The default values are both 7.
#' @param out.png A logical value indicating whether to generate output graphics as PNGs (the default) or PDFs.
#' @param res Resolution, in dots per inch; used only if \code{out.png} is set to TRUE. The default is 600.
#' @param eems.colors The EEMS color scheme as a vector of colors, ordered from low to high. Defaults to a DarkOrange to Blue divergent palette with six orange shades, white in the middle, six blue shades. Acknowledgement: The default color scheme is adapted from the \code{dichromat} package.
#' @param mcmc.iters The indices of the Voronoi diagrams to be plotted. By default, \code{eems.voronoi} generates two diagrams -- one for the migration rates and the other for the diversity rates -- for each save MCMC iteration, after burn-in and thinning.
#' @param add.title A logical value indicating whether to add the main title in the contour plots. Defaults to TRUE.
#' @param add.grid A logical value indicating whether to add the population grid or not.
#' @param col.grid The color of the population grid. Defaults to \code{gray80}.
#' @param lwd.grid The line width of the population grid. Defaults to 1.
#' @param add.outline A logical value indicating whether to add the habitat outline or not.
#' @param col.outline The color of the habitat outline. Defaults to \code{gray80}.
#' @param lwd.outline The line width of the habitat outline. Defaults to 2.
#' @param add.seeds A logical value indicating whether to add the Voronoi seeds or not.
#' @param col.seeds The color of the Voronoi seeds. Defaults to \code{green}.
#' @param pch.seeds The symbol, specified as an integer, or the character to be used for plotting the Voronoi seeds. Defaults to 4.
#' @param cex.seeds The size of the symbol/character used for plotting the Voronoi seeds. Defaults to 1.
#' @param add.demes A logical value indicating whether to add the graph demes or not.
#' @param col.demes The color of the demes. Defaults to \code{gray80}.
#' @param pch.demes The symbol, specified as an integer, or the character to be used for plotting the demes. Defaults to 1.
#' @param cex.demes The size of the deme symbol/character. Defaults to 1.
#' @param m.colscale, q.colscale A fixed range for log10-transformed migration and diversity rates, respectively. If the estimated rates fall outside the specified range, then the color scale is ignored. By default, no range is specified for either type of rates.
#' @param add.colbar A logical value indicating whether to add the color bar (the key that shows how colors map to rates) to the right of the plot. Defaults to TRUE.
#' @param add.title A logical value indicating whether to add the main title in the contour plots. Defaults to TRUE.
#' @references Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
#' @export
#' @examples
#' 
#' ## Use the provided example or supply the path to your own EEMS run.
#' eems.results.to.plot = paste0(path.package("rEEMSplots"), "/extdata/EEMS-example", sep = ""w)
#' name.figures.to.save = "EEMS-example-rEEMSplots"
#'
#' library(deldir)
#'
#' ## Plot a series of Voronoi diagrams for the EEMS model parameters:
#' ## the effective migration rates (m) and the effective diversity rates (q).
#' eems.voronoi(mcmcpath = eems.results.to.plot, 
#'             plotpath = paste0(name.figures.to.save, "-voronoi-diagrams"), 
#'             longlat = TRUE, mcmc.iters = seq(9))

eems.voronoi <- function(mcmcpath, 
                         plotpath, 
                         longlat, 

                         ## Properties of the figures to generate
                         plot.width = 7, 
                         plot.height = 7, 
                         out.png = TRUE, 
                         res = 600, 

                         ## Properties of the population grid
                         add.grid = FALSE, 
                         col.grid = "gray80", 
                         lwd.grid = 1, 
                         
                         ## Properties of the habitat outline
                         add.outline = TRUE, 
                         col.outline = "gray80",
                         lwd.outline = 2, 
                         
                         ## Properties of the graph vertices
                         add.demes = FALSE, 
                         col.demes = "gray80", 
                         pch.demes = 19, 
                         cex.demes = 1, 

                         ## Properties of the Voronoi seeds
                         add.seeds = TRUE, 
                         col.seeds = "#8AE234", 
                         pch.seeds = 4, 
                         cex.seeds = 1, 
                         
                         ## Color palette
                         eems.colors = NULL, 
                         
                         ## MCMC iterations to plot
                         mcmc.iters = NULL, 

                         ## Properties of the color key
                         add.colbar = TRUE, 
                         m.colscale = NULL, 
                         q.colscale = NULL, 
                         
                         ## Extra options
                         add.title = TRUE) {
    
    plot.params <- list(eems.colors = eems.colors, m.colscale = m.colscale, q.colscale = q.colscale, 
                        add.grid = add.grid, add.outline = add.outline, add.demes = add.demes, add.seeds = add.seeds, 
                        col.grid = col.grid, col.outline = col.outline, col.demes = col.demes, col.seeds = col.seeds, 
                        lwd.grid = lwd.grid, lwd.outline = lwd.outline, pch.demes = pch.demes, pch.seeds = pch.seeds, 
                        cex.seeds = cex.seeds, min.cex.demes = cex.demes, max.cex.demes = cex.demes, 
                        add.colbar = add.colbar, add.title = add.title)    
    plot.params <- check.plot.params(plot.params)

    load.required.package(package = 'deldir', required.by = 'eems.voronoi')

    ## A vector of EEMS output directories, for the same dataset.
    ## Assume that if eemsrun.txt exists, then all EEMS output files exist.
    mcmcpath <- mcmcpath[file.exists(paste0(mcmcpath, '/eemsrun.txt'))]
    if (!length(mcmcpath))
        stop('Please provide one existing EEMS output directory, mcmcpath')

    ## If mcmcpath is a vector of directories, use only the first existing directory.
    mcmcpath <- mcmcpath[1]

    dimns <- read.dimns(mcmcpath, longlat)
    save.params <- list(height = plot.height, width = plot.width, res = res, out.png = out.png)

    writeLines('Processing the following EEMS output directory :')
    writeLines(mcmcpath)

    plot.params$add.grid <- add.grid
    plot.params$all.demes <- FALSE
    plot.params$add.demes <- FALSE    

    save.graphics(paste0(plotpath, '-mvoronoi'), save.params)
    par(las = 1, font.main = 1)
    voronoi.diagram(mcmcpath, dimns, longlat, plot.params, mcmc.iters = mcmc.iters, is.mrates = TRUE)
    dev.off( )    

    plot.params$add.grid = FALSE
    plot.params$all.demes = add.demes

    save.graphics(paste0(plotpath, '-qvoronoi'), save.params)
    par(las = 1, font.main = 1)
    voronoi.diagram(mcmcpath, dimns, longlat, plot.params, mcmc.iters = mcmc.iters, is.mrates = FALSE)
    dev.off( )
    
}

#' A function to plot a heatmap of the residual pairwise dissimilarities, abs(observed - fitted)
#'
#' Given a vector of EEMS output directories, this function generates a heatmap of the n-by-n matrix of residual dissimilarities between pairs of individuals. It also saves the matrix of residuals to a file called \code{plotpath}-eems-resid.RData. In both the residual matrix and the corresponding heatmap, individuals are ordered according to the input \code{datapath}.coord and \code{datapath}.diffs. Applies only to the case of SNP data (when a Diffs matrix is computed).
#' @param datapath The full path and the file name of the input Diffs matrix, which is not copied by \code{runeems} to the output directory.
#' @param mcmcpath An EEMS output directory. If mcmcpath is a vector of directories, then only the first (existing) directory is used.
#' @param plotpath The full path and the file name for the graphics to be generated.
#' @param plot.width, plot.height The width and height of the graphics region for the two rate contour plots, in inches. The default values are both 7.
#' @param out.png A logical value indicating whether to generate output graphics as PNGs (the default) or PDFs.
#' @param res Resolution, in dots per inch; used only if \code{out.png} is set to TRUE. The default is 600.
#' @param heatmap.cols The heatmap color palette as a vector of colors, ordered from low to high. Defaults to the "Reds" divergent palette in the RColorBrewer package.
#' @export
#' @examples
#'
#' ## Use the provided example or supply the path to your own EEMS run
#' eems.dataset.to.plot = paste0(path.package("rEEMSplots"), "/extdata/EEMS-barrier")
#' eems.results.to.plot = paste0(path.package("rEEMSplots"), "/extdata/EEMS-barrier")
#' name.figures.to.save = "EEMS-barrier"
#'
#' eems.resid.heatmap(datapath = eems.dataset.to.plot, 
#'                    mcmcpath = eems.results.to.plot, 
#'                    plotpath = name.figures.to.save, 
#'                    heatmap.cols = c("gray99", "red"))

eems.resid.heatmap <- function(datapath, 
                               mcmcpath, 
                               plotpath, 

                               ## Properties of the figure to generate
                               plot.width = 7, 
                               plot.height = 7, 
                               out.png = TRUE, 
                               res = 600, 

                               ## (Divergent) color palette
                               heatmap.cols = NULL, 
                               heatmap.colscale = NULL) {

    save.params = list(height = plot.height, width = plot.width, res = res, out.png = out.png)
    
    ## A vector of EEMS output directories, for the same dataset.
    ## Assume that if eemsrun.txt exists, then all EEMS output files exist.
    mcmcpath <- mcmcpath[file.exists(paste0(mcmcpath, '/eemsrun.txt'))]
    if (!length(mcmcpath))
        stop('Please provide one existing EEMS output directory, mcmcpath')

    writeLines('Processing the following EEMS output directories :')
    writeLines(mcmcpath)

    eems.resid <- heatmap.resid(datapath, mcmcpath)
    save(eems.resid, file = paste0(plotpath, '-eems-resid.RData'))
    save.graphics(paste0(plotpath, '-eems-resid'), save.params)
    par(las = 1, font.main = 1, mar = c(0, 0, 0, 0)+0.1)
    key <- myheatmap(abs(eems.resid), col = heatmap.cols, colscale = heatmap.colscale)
    dev.off( )

    if (out.png) {
        save.params$height <- 6
        save.params$width <- 1.5
    } else {
        save.params$height <- 12
        save.params$width <- 3
    }
    
    save.graphics(paste0(plotpath, '-eems-resid-key'), save.params)
    par(las = 1, font.main = 1, mar = c(1, 1, 5, 8))
    myfilled.legend(col = key$colors, levels = key$levels, 
                    key.axes = axis(4, tick = FALSE, hadj = 1, line = 4, cex.axis = 2), 
                    key.title = mtext("abs(r)", side = 3, cex = 2.5, line = 1.5, font = 1))
    dev.off( )
}
