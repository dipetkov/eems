
default.eems.colors <- function() {
  c(
    # orange sequence, from dark to light
    "#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC",
    # very slightly off-white
    "#FBFBFB",
    # blue sequence, from light to dark
    "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99"
  )
}

is.color <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  # Check that a string, or a vector of strings, represents a valid color
  sapply(x, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
  })
}

set.colscale <- function(x) {
  if (min(x) < max(x)) range(x) else c(-1, 1)
}

mypalette <- function(colors = NULL, colscale = NULL, n = 299) {
  if (is.null(colors)) {
    # Use RColorBrewer's 9-color "Reds" palette.
    colors <- c(
      "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
      "#EF3B2C", "#CB181D", "#A50F15", "#67000D"
    )
  }
  colors <- colorRampPalette(colors)(n = n)

  if (!is.null(colscale)) {
    if (length(colscale) != (n + 1)) {
      levels <- seq(from = min(colscale), to = max(colscale), length.out = n + 1)
    } else {
      levels <- sort(colscale)
    }
  }
  list(colors = colors, levels = levels)
}

eems.colscale <- function(values, nlevels, colscale) {
  max_val <- max(values + 0.001, colscale)
  min_val <- min(values - 0.001, colscale)
  seq(from = min_val, to = max_val, length = nlevels + 1)
}
