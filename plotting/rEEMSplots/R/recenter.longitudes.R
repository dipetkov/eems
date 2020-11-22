
#' Shift spatial object longitudinally
#'
#' The function `sp::recenter` shifts coordinates for a Pacific view by shifting longitudes
#' from [-180, 180] to [0, 360]. I modified `sp::recenter` to shift longitudes by any degree.
#'
#' @param obj A SpatialPolygons or a SpatialPolygonsDataFrame object.
#' @param degree Number of degrees (positive or negative) to shift.
#' @export
setGeneric("recenter_by", function(obj, degree) {
  standardGeneric("recenter_by")
})

recenter_by.SpatialPolygons <- function(obj, degree) {
  proj <- is.projected(obj)
  if (is.na(proj)) stop("unknown coordinate reference system")
  if (proj) stop("cannot shift projected coordinate reference system")
  projargs <- methods::slot(obj, "proj4string")
  pls <- methods::slot(obj, "polygons")
  Srl <- lapply(pls, recenter_by.Polygons, degree)
  res <- SpatialPolygons(Srl, proj4string = projargs)
  res
}

setMethod("recenter_by", "SpatialPolygons", recenter_by.SpatialPolygons)
setMethod("recenter_by", "SpatialPolygonsDataFrame", recenter_by.SpatialPolygons)

recenter_by.Polygons <- function(obj, degree) {
  ID <- methods::slot(obj, "ID")
  rings <- methods::slot(obj, "Polygons")
  srl <- lapply(rings, recenter_by.Polygon, degree)
  res <- Polygons(srl, ID = ID)
  res
}

recenter_by.Polygon <- function(obj, degree) {
  crds <- methods::slot(obj, "coords")
  hole <- methods::slot(obj, "hole")
  crds[, 1] <- crds[, 1] + degree
  res <- Polygon(crds, hole)
  res
}
