
plot.logposterior <- function(mcmcpath) {
  message("Plotting posterior probability trace")
  nchains <- length(mcmcpath)
  # Here colors = RColorBrewer::brewer.pal(12, "Paired") + "black"
  colors <- rep_len(c(
    "#000000", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
    "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
  ), length.out = nchains)
  ltypes <- rep_len(1:4, length.out = nchains)
  posteriors <- list()
  yrange <- NULL
  niters <- NULL
  for (i in 1:nchains) {
    path <- mcmcpath[i]
    message(path)
    stopifnot(all(file.exists(file.path(path, "mcmcpilogl.txt"))))
    pilogl <- scan(file.path(path, "mcmcpilogl.txt"), quiet = TRUE)
    pilogl <- matrix(pilogl, ncol = 2, byrow = TRUE)
    posterior <- pilogl[, 1] + pilogl[, 2]
    posteriors[[i]] <- posterior
    yrange <- range(c(yrange, posterior))
    niters <- max(niters, length(posterior))
  }
  plot(c(1, niters), yrange,
    type = "n",
    xlab = "MCMC iteration  (after burn-in and thinning)",
    ylab = "log posterior"
  )
  if (nchains == 1) {
    mtext(
      side = 3, line = 2, cex = 1.3,
      text = "Has the MCMC chain converged?"
    )
    mtext(
      side = 3, line = 0.5, cex = 1,
      text = "If not, restart EEMS and/or increase numMCMCIter, numBurnIter, numThinIter"
    )
  } else {
    mtext(
      side = 3, line = 2, cex = 1.3,
      text = "Have the MCMC chains converged?"
    )
    mtext(
      side = 3, line = 0.5, cex = 1,
      text = "If not, restart EEMS and/or increase numMCMCIter, numBurnIter, numThinIter"
    )
  }
  for (i in 1:nchains) {
    posterior <- posteriors[[i]]
    niters <- length(posterior)
    lines(1:niters, posterior, col = colors[i], lty = ltypes[i], lwd = 2)
  }
  legend("topright",
    legend = 1:nchains, col = colors[1:nchains],
    lty = ltypes[1:nchains], lwd = 2, bty = "n", inset = c(-0.12, 0), cex = 0.5
  )
}
