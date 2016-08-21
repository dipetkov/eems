
## Check whether the required packages are already installed
## and installed them if necessary.

req_packages <- c("shiny", "ggplot2", "dplyr")
new_packages <- setdiff(req_packages, installed.packages()[, "Package"])
if(length(new_packages)) { install.packages(new_packages) }
suppressWarnings(suppressMessages(lapply(req_packages, library, 
                                         character.only = TRUE)))

emptyDF <- function() {
  data_frame(
    alpha.x = numeric(),
    alpha.y = numeric(),
    beta.x = numeric(),
    beta.y = numeric(),
    fitted = numeric(),
    obsrvd = numeric(),
    selected_ = logical())
}

c2Palette <- c("#000000", "#66D65C") ## simple two-color scale (black and green)

scatterPlot <- function(df, cols, xlab, ylab) {
  if (nrow(df) > 0) {
    ggplot(df, aes_string(x = cols[1], y = cols[2])) +
      ## Plot the *selected* points on top of the *unselected* ones.
      geom_point(data = df %>% filter(selected_ == FALSE),
                 color = c2Palette[1], alpha = 0.3) +
      geom_point(data = df %>% filter(selected_ == TRUE),
                 color = c2Palette[2], alpha = 1) +
      coord_fixed(ratio = 1) + 
      labs(x = xlab, y = ylab) + 
      theme_minimal()
  }
}
scatterPlotDist <- function(df, cols) {
  scatterPlot(df %>% mutate(selected_ = FALSE), cols,
              xlab = "Fitted dissimilarity",
              ylab = "Observed dissimilarity")
}
scatterPlotDeme <- function(df, cols) {
  scatterPlot(df, cols,
              xlab = "longitude",
              ylab = "latitude")
}
