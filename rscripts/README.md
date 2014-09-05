

## Visualize EEMS results

These scripts assume that you have run EEMS to generate output files that all have name `mcmcpath` (with different extensions).

### Required R packages 

The plotting scripts use two R packages, `fields` and `deldir`, that should be installed first. Additional color schemes can be generated with the `dichromat` package.

### Plot mean posterior surfaces

There are two mean posterior surfaces: the estimated effective migration surface and the estimated effective diversity surface. The former characterizes the differences between two difference demes and the latter characterizes the differences between two different samples from the same deme. (In population genetics, a deme is a group of randomly mating individuals, in the same location.)

```
# Read in information about the habitat
dimns <- read.dimns(datapath)

# Plot the two estimated effective surfaces; return information about the scaling
mlegend <- mcmc.mrates(mcmcpath,dimns)
qlegend <- mcmc.qrates(mcmcpath,dimns)

# Plot the legends to know what values the colors stand for
mcmc.mrates.legend(datapath,mcmcpath,dimns,mmrks,mlegend)
mcmc.qrates.legend(datapath,mcmcpath,dimns,qmrks,qlegend)
```

### Plot scatterplots of the observed and fitted distances

EEMS estimates the model parameters theta, which include the migration rates and the diversity rates, so that the fitted distances Delta(theta) are approximatedly equal to the observed distances D. Plot Delta against D to check whether the relationship seems to hold. Otherwise the EEMS might not have converged or the underlying model might not fit the data very well. 

```
dist.scatterplot(mcmcpath)
```

### Plot a trace of the posterior

Plot a trace of the log posterior, after burn-in and thinning. This might also help to check whether EEMS appears to have converged.

```
pilogl <- read.table(paste(mcmcpath,'.mcmcpilogl',sep=''))
# log(posterior) = log(prior) + log(likelihood)
post <- pilogl[,1] + pilogl[,2]
plot(post)
```

### Plot draws from the posterior Voronoi tessellation 

Plot draws from the posterior distribution on the two colored Voronoi tessellations, to illustrate the uncertainty in the estimated effective rate parameters. (This will generate a large number of plots.)

```
mlegend <- mcmc.mrates.voronoi(mcmcpath,dimns)
qlegend <- mcmc.qrates.voronoi(mcmcpath,dimns)
```

### Acknowledgements

* The default Dark Orange to Blue color scheme, combines two color schemes from the `dichromat` package, which itself is based on a collection of [color schemes for scientific data graphics] (http://geog.uoregon.edu/datagraphics/color_scales.htm). See also Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
* We use a modified version of the [`filled.contour`] (http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots) function, which has been modified by Ian Taylor to remove the key and facilitate overplotting with contour( ), and further modified by Carey McGilliard and Bridget Ferris to allow multiple plots on one page.
