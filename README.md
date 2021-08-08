## EEMS

This repository contains an implementation of the EEMS method for analyzing and visualizing spatial population structure from geo-referenced genetic samples. EEMS uses the concept of *effective migration* to model the relationship between genetics and geography, and it outputs an **e**stimated **e**ffective **m**igration **s**urface (hence, EEMS) - a visual representation of population structure that can highlight potential regions of higher-than-average and lower-than-average historic gene flow.

Please consider reading the paper `EEMS-article.pdf` (or the published version [here](http://www.nature.com/ng/journal/v48/n1/full/ng.3464.html)) and the documentation `EEMS-doc.pdf` first.

* The directory `bed2diffs` contains a small program that reads genotypes in Plink binary format and computes the matrix of average pairwise differences, to pass to EEMS.
* The directory `pipeline` contains a pipeline with the goal to combine all scripts.
* The directory `plotting` contains an R package, rEEMSplots, to visualize EEMS results.
* The directory `runeems_sats` contains a C++ implementation of EEMS to use with microsatellite data.
* The directory `runeems_snps` contains a C++ implementation of EEMS to use with SNP data.

### C++ implementation

The C++ implementation uses the [Eigen](http://eigen.tuxfamily.org) template library for linear algebra computations and the [Boost](http://www.boost.org) libraries for random number generation and the habitat geometry. EEMS has been tested with Eigen ~~3.2.2~~ 3.3.9 and Boost ~~1.57~~ 1.76.0 and might not be compatible with newer versions of Boost/Eigen. After downloading Eigen (which does not need installation) and installing Boost, update the variables `EIGEN_INC`, `BOOST_LIB`, `BOOST_INC` in the Makefile. The dynamic Boost libraries are linked to slightly differently on a Mac and a Linux machine, so run `make darwin` on a Mac or `make linux` on a Linux machine. This C++ implementation has not been tested on Windows.

### EEMS for SNPs and microsatellites

There are two versions of EEMS: `runeems_snps` for SNP data and `runeems_sats` for microsatellite data. The data input format and the EEMS model are somewhat different for SNPs and microsatellites, hence the two versions. The source code can be found in `runeems_snps/src` and `runeems_sats/src`, respectively. The directories `runeems_snps/data` and `runeems_sats/data` contain data, simulated with `ms`, to illustrate the input file format and how EEMS is run.

For SNP data, the input files are:

* datapath.**diffs**: the matrix of average pairwise genetic differences (which can be computed with the program `bed2diffs`)

* datapath.**coord**: the sample coordinates (two coordinates per sample, one sample per line)

* datapath.**outer**: the habitat coordinates (as a sequence of vertices that outline a closed polygon)

Here `datapath` is the full path + the file name (but without the extension).

And for microsatellite data, the input files are:

* datapath.**sites**: genotype data (one sample per line, a negative number indicates the allele is missing)

* datapath.**coord**: as described above for SNP data.

* datapath.**outer**: as described above for SNP data.

EEMS also requires a configuration file where various program options can be specified. For example, consider `runeems_snps/src/params-simno1.ini` which contains the following information:

```
datapath = ./data/barrier-schemeZ-nIndiv300-nSites3000
mcmcpath = ./data/barrier-schemeZ-nIndiv300-nSites3000-EEMS-nDemes200-simno1
nIndiv = 300
nSites = 3000
nDemes = 200
diploid = false
numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 9999
```

This file specifies the following *required* input arguments: the path to the input data (**datapath**), the path to the output data (**mcmcpath**), the number of samples (**nIndiv**), the number of markers (**nSites**), the density of the population grid (**nDemes**), is the species diploid or haploid (**diploid**), the number of MCMC and burn-in iterations (**numMCMCIter**, **numBurnIter**), and the thinning interval (**numThinIter**).

```
./runeems_snps --params params-simno1.ini --seed 123
```

(Specifying the random seed is optional.)

Finally, the EEMS results can be visualized with the function `eems.plots` defined in the R package `rEEMSplot`. The package is not on CRAN, so install it from source instead. (The code is in the directory `plotting`.)

```
## Part 1: Install rEEMSplots
## Check that the current directory contains the rEEMSplots source directory
if (file.exists("./rEEMSplots")) {
  install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}


## Possibly change the working directory with setwd()


## Part 2: Generate graphics
library(rEEMSplots)

mcmcpath = "./data/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes200-simno1"
plotpath = "./plot/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes200-simno1-rEEMSplots"

eems.plots(mcmcpath, plotpath, longlat = TRUE)
```

The function `eems.plots` generates several figures automatically (to encourage looking at all the figures). There are examples in `EEMS-doc.pdf`, with captions that explain each figure. `eems.plots` also saves several objects to an RData file, which can be read back from the file with `load`.

```
load(paste0(plotpath, "-rdist.RData"))

ls()
#> [1] "B.component" "G.component" "W.component" "xym.values"  "xyq.values"

library("ggplot2")
library("dplyr")

## Reproduce plotpath-rdist01.png,
## which plots observed vs fitted dissimilarities between demes
ggplot(B.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point()
```

The data frame `B.component` contains information about the dissimilarities between observed pairs of demes (alpha, beta). "Observed" means that each deme has at least one sampled individual assigned to it. Furthermore, each deme has two coordinates (x and y, longitude and latitude) and these are labeled `alpha.x, alpha.y` for deme alpha and `beta.x, beta.y` for deme beta.

```
## Reproduce plotpath-rdist02.png,
## which plots observed vs fitted dissimilarities within demes
ggplot(W.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point()

## Reproduce plotpath-rdist03.png,
## which plots observed dissimilarities against great circle distances between demes
ggplot(W.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point()

## Empty matrices unless you have specified additional coordinates
## at which to estimate the migration and diversity rates
## with the `xy.coords` argument to `eems.plots`
xym.values
xyq.values
```

### Shiny app

It is straightforward to make the scatter plots which plot fitted vs observed dissimilarities between demes. But it is even more useful to link the dissimilarities to the corresponding locations in the habitat. The between-demes dissimilarities saved in `plotpath-rdist.RData` can be visualized with a little interactive Shiny app available [here](https://dipetkov.shinyapps.io/rEEMSshiny/). It can help to identify outliers in the pairwise scatter plots, if any, with specific geographic locations on the map.

### The MATLAB/Octave implementation

This version is less efficient and it is provided here only for completeness, as it is the original implementation.
