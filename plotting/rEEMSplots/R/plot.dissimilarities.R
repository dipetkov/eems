
JtDJ2BandW <- function(JtDJ, sizes = NULL) {
  # There should be NAs on the main diagonal of JtDobsJ -- these elements correspond to demes with
  # a single observed sample. And since there is a single sample taken, the average dissimilarity
  # between two distinct individuals from such demes cannot be observed.
  # Instead, I will use the average within dissimilarity, W, across the demes with multiple samples.
  # This is make a difference only if remove.singletons = FALSE, which is not the default option.
  JtDJ <- as.matrix(JtDJ)
  if (!is.null(sizes)) {
    diag(JtDJ)[sizes < 2] <- mean(diag(JtDJ)[sizes >= 2])
  }
  n <- nrow(JtDJ)
  W <- diag(JtDJ)
  S <- matrix(W, n, n)
  B <- JtDJ - (S + t(S)) / 2
  B <- B[upper.tri(B, diag = FALSE)]
  list(W = W, B = B)
}

geo.distm <- function(coord, longlat, plot.params) {
  if (!longlat) {
    coord <- coord[, c(2, 1)]
  }
  if (!is.null(plot.params$proj.in)) {
    # If the locations are projected, convert them to longitude/latitude pairs
    coord <- sp::SpatialPoints(coord, proj4string = CRS(plot.params$proj.in))
    coord <- sp::spTransform(coord, CRSobj = CRS("+proj=longlat +datum=WGS84"))
    coord <- sp::coordinates(coord)
  }
  Dist <- sp::spDists(coord, longlat = TRUE)
  Dist <- Dist[upper.tri(Dist, diag = FALSE)]
  Dist
}

dist.scatterplot <- function(mcmcpath, longlat, plot.params,
                             remove.singletons = TRUE, add.abline = FALSE,
                             add.r.squared = FALSE, highlight = NULL) {
  message("Plotting average dissimilarities within and between demes")
  for (path in mcmcpath) {
    stopifnot(all(file.exists(file.path(path, c(
      "rdistJtDhatJ.txt",
      "rdistJtDobsJ.txt",
      "rdistoDemes.txt"
    )))))
  }
  nchains <- length(mcmcpath)
  # List of observed demes, with number of samples taken collected
  # Each row specifies: x coordinate, y coordinate, n samples
  oDemes <- scan(file.path(mcmcpath[1], "rdistoDemes.txt"), quiet = TRUE)
  oDemes <- matrix(oDemes, ncol = 3, byrow = TRUE)
  Sizes <- oDemes[, 3]
  nPops <- nrow(oDemes)
  Demes <- seq(nPops)
  if (nPops < 2) {
    message(
      "All individuals sampled from the same deme. ",
      "Check that individual and habitat coordinates are given ",
      "in consistent order (either latitude/longitude or ",
      "longitude/latitude) in the *.outer and *.coord files."
    )
    return(NULL)
  }
  JtDobsJ <- matrix(0, nPops, nPops)
  JtDhatJ <- matrix(0, nPops, nPops)
  for (path in mcmcpath) {
    message(path)
    oDemes1 <- scan(file.path(path, "rdistoDemes.txt"), quiet = TRUE)
    oDemes1 <- matrix(oDemes1, ncol = 3, byrow = TRUE)
    if (sum(dim(oDemes) != dim(oDemes1)) || sum(oDemes != oDemes1)) {
      plot.params$add.demes <- TRUE
      plot.params$add.grid <- TRUE
      plot.params$add.title <- TRUE
      plot.params$add.colbar <- FALSE
      null.eems.contour(mcmcpath[1], longlat, plot.params)
      mtext(
        side = 3, line = 2, cex = 1.3, text =
          "EEMS results for at least two different population grids"
      )
      null.eems.contour(path, longlat, plot.params)
      mtext(
        side = 3, line = 2, cex = 1.3, text =
          "EEMS results for at least two different population grids"
      )
      message("EEMS results for at least two different population grids")
      return(NULL)
    }
    JtDobsJ <- JtDobsJ +
      as.matrix(read.table(file.path(path, "rdistJtDobsJ.txt"), header = FALSE))
    JtDhatJ <- JtDhatJ +
      as.matrix(read.table(file.path(path, "rdistJtDhatJ.txt"), header = FALSE))
  }
  JtDobsJ <- JtDobsJ / nchains
  JtDhatJ <- JtDhatJ / nchains
  colnames(JtDobsJ) <- Demes
  rownames(JtDobsJ) <- Demes
  colnames(JtDhatJ) <- Demes
  rownames(JtDhatJ) <- Demes
  Deme1 <- matrix(Demes, nrow = nPops, ncol = nPops)
  Deme2 <- t(Deme1)
  tempi <- matrix(Sizes, nPops, nPops)
  Small <- pmin(tempi, t(tempi))
  Small <- Small[upper.tri(Small, diag = FALSE)]
  alpha <- Deme1[upper.tri(Deme1, diag = FALSE)]
  beta <- Deme2[upper.tri(Deme2, diag = FALSE)]
  # Under pure isolation by distance, we expect the genetic dissimilarities
  # between demes increase with the geographic distance separating them
  Dist <- geo.distm(oDemes[, 1:2], longlat, plot.params)
  if (sum(Sizes > 1) < 2) {
    # Sizes(alpha) > 1 means that there are at least two individuals assigned to deme alpha
    JtDJ.component <- data.frame(
      alpha.x = oDemes[, 1][alpha],
      alpha.y = oDemes[, 2][alpha],
      beta.x = oDemes[, 1][beta],
      beta.y = oDemes[, 2][beta],
      size = Small,
      col = sub.scattercols(Small),
      fitted = JtDhatJ[upper.tri(JtDhatJ, diag = FALSE)],
      obsrvd = JtDobsJ[upper.tri(JtDobsJ, diag = FALSE)],
      stringsAsFactors = FALSE
    )
    sub.scatterplot("JtDJ", JtDJ.component,
      remove.singletons = FALSE,
      add.abline, add.r.squared
    )
    message(
      "There is one or zero demes with multiple observed individuals; ",
      "plotting dissimilarities between singletons instead."
    )
    return(NULL)
  }
  out1 <- JtDJ2BandW(JtDobsJ, Sizes)
  Wobs <- out1$W
  Bobs <- out1$B
  out2 <- JtDJ2BandW(JtDhatJ)
  What <- out2$W
  Bhat <- out2$B
  B.component <- data.frame(
    alpha.x = oDemes[, 1][alpha],
    alpha.y = oDemes[, 2][alpha],
    beta.x = oDemes[, 1][beta],
    beta.y = oDemes[, 2][beta],
    fitted = Bhat,
    obsrvd = Bobs,
    size = Small,
    col = sub.scattercols(Small),
    stringsAsFactors = FALSE
  )
  W.component <- data.frame(
    alpha.x = oDemes[, 1][Demes],
    alpha.y = oDemes[, 2][Demes],
    fitted = What,
    obsrvd = Wobs,
    size = Sizes,
    col = sub.scattercols(Sizes),
    stringsAsFactors = FALSE
  )
  G.component <- data.frame(
    alpha.x = oDemes[, 1][alpha],
    alpha.y = oDemes[, 2][alpha],
    beta.x = oDemes[, 1][beta],
    beta.y = oDemes[, 2][beta],
    fitted = Dist,
    obsrvd = Bobs,
    size = Small,
    col = sub.scattercols(Small),
    stringsAsFactors = FALSE
  )
  # Code to highlight some demes in the scatterplots of observed vs fitted distances
  # This can help to identify demes that are of special interest among all the other points
  # Since `highlights` is not fully working yet and to avoid notes about "no visible
  # binding for global variable XXX", I have commented the code out.
  sub.scatterplot("Between", B.component, remove.singletons, add.abline, add.r.squared)
  sub.scatterplot("Within", W.component, remove.singletons, add.abline, add.r.squared)
  sub.scatterplot("GeoDist", G.component, remove.singletons, add.abline = FALSE, add.r.squared)
  list(B.component = B.component, W.component = W.component, G.component = G.component)
}
