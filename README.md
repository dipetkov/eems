
This repository contains an implementation of the EEMS method for analyzing and visualizing spatial population structure from geo-referenced genetic samples. EEMS uses the concept of *effective migration* to model the relationship between genetics and geography, and it outputs an **e**stimated **e**ffective **m**igration **s**urface (hence, EEMS) - a visual representation of population structure that can highlight potential regions of higher-than-average and lower-than-average historic gene flow.

## The C++ implementation

* `bed2diffs` contains a small program that reads genotypes in Plink binary format and computes the matrix of average pairwise differences, to pass to EEMS
* `examples` shows how to run EEMS using data simulated with `ms`
* `mscripts` contains the source code for MATLAB/OCTAVE
* `rscripts` contains scripts to visualize EEMS results
* `pipeline` contains a pipeline with the goal to combine all scripts

Please consider reading the paper, `EEMS-article.pdf`, and the documentation, `eems-doc.pdf`, first.

There are two programs: `runeems_snps` for SNP data and `runeems_sats` for microsatellite data. (The data input format and the EEMS model are somewhat different for SNPs and microsatellites, hence the two versions.) The source code can be found in `runeems_snps/src` and `runeems_sats/src`, respectively.

The C++ implementation uses the [Eigen](http://eigen.tuxfamily.org) template library for linear algebra computations and the [Boost](http://www.boost.org) libraries for random number generation and the habitat geometry. EEMS has been tested with Eigen 3.2.2 and Boost 1_57.

After obtaining Eigen and installing Boost, update the variables `EIGEN_INC`, `BOOST_LIB`, `BOOST_INC` in the Makefile and run `make` on the command line.

The folders `runeems_snps/data` and `runeems_sats/data` contain data, simulated with `ms`, to illustrate how input files are formatted and how EEMS is run.

For SNP data, the input files are:

* datapath.diffs: the matrix of average pairwise genetic differences (which can be computed with the program `bed2diffs`)

* datapath.coord: the sample coordinates (two coordinates per sample, one sample per line)

* datapath.outer: the habitat coordinates (as a sequence of vertices that outline a closed polygon)

Here `datapath` indicates the filename, e.g., `runeems_snps/data/barrier-schemeX-nIndiv300-nSites3000`.

And for microsatellite data, the input files are:

* datapath.sites: genotype data (one sample per line, a negative number indicates the allele is missing)

* datapath.coord: as described above for SNP data.

* datapath.outer: as described above for SNP data.

EEMS also requires an input parameter file. For example, consider `runeems_snps/src/params-simno1.ini` which contains the following information:

```
datapath = ./data/barrier-schemeZ-nIndiv300-nSites3000
mcmcpath = ./data/barrier-schemeZ-nIndiv300-nSites3000-EEMS-nDemes200-simno1
nIndiv = 300
nSites = 3000
nDemes = 200
diploid = false
numMCMCIter = 12000
numBurnIter = 6000
numThinIter = 99
```

This file specified the following *required* input arguments: the path to the input data, the path to the output data, the number of samples, the number of markers, the density of the grid, is the species diploid or haploid, the number of MCMC and burn-in iterations, and the thinning interval.

```
./runeems_snps --params params-simno1.ini --seed 123
```

(Specifying the random seed is optional.)

Finally, the EEMS results can be plotted with the R script `default.eems.plots.R` available in `runeems_snps/plot/`.

```
source('default.eems.plots.R')mcmcpath <- './data/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes153-simno1'plotpath <- './plot/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes153-simno1'eemsplots(mcmcpath,plotpath,add.map=FALSE)
```

## The MATLAB/Octave implementation

This version is less efficient and it is provided here only for completeness, as it is the original implementation.
