
scatterPlot <- function(df, cols, xlab, ylab) {
  if (nrow(df) > 0) {
    ggplot(df, aes_string(x = cols[1], y = cols[2])) +
      ## Plot the *selected* points on top of the *unselected* ones.
      geom_point(data = df %>% filter(selected_ == FALSE),
                 color = "#000000", alpha = 0.3) +
      geom_point(data = df %>% filter(selected_ == TRUE),
                 color = "#66D65C", alpha = 1) +
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
