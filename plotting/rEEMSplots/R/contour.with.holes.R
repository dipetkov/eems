
#' Create a habitat with holes
#' @param outer A two-column matrix of coordinates (longitude, latitude).
#' @param holes List of holes. Each hole is a two-column matrix of coordinates
#' (longitude, latitude).
#' @export
create_habitat_with_holes <- function(outer, holes) {
  polys <- list(Polygon(outer, FALSE))

  for (hole in holes) {
    polys <- c(polys, Polygon(hole, TRUE))
  }

  SpatialPolygons(list(Polygons(polys, "1")))
}

point_is_in_habitat <- function(point, habitat) {
  point <- matrix(point, nrow = 1)
  point <- SpatialPoints(point)
  point <- point[habitat]
  length(point) == 1
}

reindex_demes_sequentially <- function(in_habitat) {
  indices_in_habitat <- cumsum(in_habitat)
  indices_in_habitat[!in_habitat] <- NA_integer_
  indices_in_habitat
}

reindex_demes <- function(demes, new_indices) {
  demes[!is.na(new_indices), ]
}

reindex_edges <- function(edges, new_indices) {
  edges[, 1] <- new_indices[edges[, 1]]
  edges[, 2] <- new_indices[edges[, 2]]
  na.omit(edges)
}

plot_basic_popgrid <- function(gridpath, col.grid = "green", col.demes = "red", ...) {
  demes <- read.table(paste0(gridpath, ".demes"))
  edges <- read.table(paste0(gridpath, ".edges"))

  bitmap(
    paste0(gridpath, ".png"),
    type = "png16m", res = 600, units = "in", height = 6, width = 8, ...
  )

  plot(demes, col = col.demes, pch = 19, xlab = "", ylab = "")

  for (i in seq_len(nrow(edges))) {
    a <- edges[i, 1]
    b <- edges[i, 2]
    lines(demes[c(a, b), ], col = col.grid)
  }
  dev.off()
}

#' Remove demes that lie outside of the habitat
#' @param habitat A SpatialPolygons habitat (possibly with holes).
#' @param mcmcpath An EEMS output directory which contains demes.txt and edges.txt.
#' @param gridpath Output filename.
#' @export
#' @examples
#' extdata <- system.file("extdata", package = "rEEMSplots")
#'
#' # The input is an EEMS run with the default regular triangular grid without holes
#' eems_run_with_default_habitat<- file.path(extdata, "EEMS-barrier")
#' # The output filepath; three files with the same name but different extension to be created
#' custom_grid_path_with_holes <- file.path(path.expand("~"), "gridpath-with-holes")
#'
#' # Load the population habitat
#' outer <- read.table(file.path(eems_run_with_default_habitat, "outer.txt"))
#' # Each hole is a ring (simple closed polygon) and the holes don't overlap
#' hole1 <- data.frame(V1 = c(2., 5., 5., 2., 2.), V2 = c(2., 2., 5., 5., 2.))
#' hole2 <- data.frame(V1 = c(6.5, 10., 8., 6.5), V2 = c(2.5, 5., 5., 2.5))
#' 
#' # Create the new habitat with holes
#' new_habitat_with_holes <- create_habitat_with_holes(outer, list(hole1, hole2))
#'
#' # plot(new_habitat_with_holes, col = "gray")
#'
#' # This function creates three files:
#' # * The pair gridpath.demes and gridpath.edges which specify a custom population grid with holes
#' # * The plot gridpath.png which visualizes the custom grid with demes in red and edges in green
#' remove_demes_outside_habitat(
#'   habitat=new_habitat_with_holes,
#'   mcmcpath=eems_run_with_default_habitat,
#'   gridpath=custom_grid_path_with_holes)

remove_demes_outside_habitat <- function(habitat, mcmcpath, gridpath) {
  demes <- read.table(file.path(mcmcpath, "demes.txt"))
  edges <- read.table(file.path(mcmcpath, "edges.txt"))

  in_habitat <- apply(demes, 1, point_is_in_habitat, habitat)
  new_indices <- reindex_demes_sequentially(in_habitat)

  demes <- reindex_demes(demes, new_indices)
  edges <- reindex_edges(edges, new_indices)

  write.table(demes, paste0(gridpath, ".demes"),
    row.names = FALSE, col.names = FALSE
  )
  write.table(edges, paste0(gridpath, ".edges"),
    row.names = FALSE, col.names = FALSE
  )

  plot_basic_popgrid(gridpath)
}
