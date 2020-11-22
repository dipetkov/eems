
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

reindex_ipmap <- function(ipmap, new_indices) {
  new_indices[ipmap[, 1]]
}

#' Remove demes that lie outside of the habitat
#' @param habitat A SpatialPolygons habitat (possibly with holes).
#' @param mcmcpath An EEMS output directory which contrains at a minimum the files
#' demes.txt, edges.txt, ipmap.txt.
#' @param gridpath Output filename.
#' @export
remove_demes_outside_habitat <- function(habitat, mcmcpath, gridpath) {
  demes <- read.table(file.path(mcmcpath, "demes.txt"))
  edges <- read.table(file.path(mcmcpath, "edges.txt"))
  ipmap <- read.table(file.path(mcmcpath, "ipmap.txt"))

  in_habitat <- apply(demes, 1, point_is_in_habitat, habitat)
  new_indices <- reindex_demes_sequentially(in_habitat)

  demes <- reindex_demes(demes, new_indices)
  edges <- reindex_edges(edges, new_indices)
  ipmap <- reindex_ipmap(ipmap, new_indices)

  write.table(demes, paste0(gridpath, ".demes"),
    row.names = FALSE, col.names = FALSE
  )
  write.table(edges, paste0(gridpath, ".edges"),
    row.names = FALSE, col.names = FALSE
  )
  write.table(ipmap, paste0(gridpath, ".ipmap"),
    row.names = FALSE, col.names = FALSE
  )
}
