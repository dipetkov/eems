

## Visualize EEMS results

This script assumes that you have run EEMS three times, with three different random seeds, and the results are stored in the following directories:

```
barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes-simno1/
barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes-simno2/
barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes-simno3/
```

### Required R packages 

The plotting scripts use two R packages, `fields` and `deldir`, that should be installed first. Additional color schemes can be generated with the `dichromat` package. The geographical map can be drawn with the `maps` package.

One function generates several figures that represent the EEMS results:

```
source('default.eems.plots.R')

mcmcpath <-
paste('barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes-simno',1:3,sep='')
plotpath <- 'barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes-simno1_3'

## mcmcpath is a list of three output directories; the results are averaged
eemsplots(mcmcpath,plotpath,add.map=FALSE)
```

There are two mean posterior surfaces: the estimated effective migration surface and the estimated effective diversity surface. The former characterizes the differences between two difference demes and the latter characterizes the differences between two different samples from the same deme. (In population genetics, a deme is a group of randomly mating individuals, in the same location.)

EEMS estimates the model parameters theta, which include the migration rates and the diversity rates, so that the fitted distances Delta(theta) are approximatedly equal to the observed distances D. Plot Delta against D to check whether the relationship seems to hold. Otherwise the EEMS might not have converged or the underlying model might not fit the data very well. 

Plot a trace of the log posterior, after burn-in and thinning. This might also help to check whether EEMS appears to have converged.

### Acknowledgements

* The default Dark Orange to Blue color scheme, combines two color schemes from the `dichromat` package, which itself is based on a collection of [color schemes for scientific data graphics] (http://geog.uoregon.edu/datagraphics/color_scales.htm). See also Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
* We use a modified version of the [`filled.contour`] (http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots) function, which has been modified by Ian Taylor to remove the key and facilitate overplotting with contour( ), and further modified by Carey McGilliard and Bridget Ferris to allow multiple plots on one page.
