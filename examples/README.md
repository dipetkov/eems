

## How to run EEMS on a simulated dataset

This script uses a dataset simulated with `ms` to show how to run EEMS, in either MATLAB or Octave.

### EEMS input files

There are three input files:
* datapath.coord: sample coordinates (two coordinates per sample, one sample per line)
* datapath.diffs: matrix of average pairwise genetic differences (can be computed with the program bed2diffs if the dataset is already in plink binary format)
* datapath.dimns: 3x2 matrix which specifies the habitat range in the x direction (longitude), the habitat range in the y direction (latitude), and on the third line, the number of samples and the number of markers

Comment about the last file, `datapath.dimns`: It might be a good idea to add a "border" around the sampling locations, which is the reason for introducing `datapath.dimns`. Suppose that `x` is a vector of longitudes and `y` is a vector of latitudes. (Note that `datapath.coord` is simply the matrix `[x,y]`, without a header or sample IDs.) Most straightforwardly, the first two lines in `datapath.dimns` can be specified as

```
min(x) max(x)
min(y) max(y)
```

However, in practice it seems to be a good idea to extend the habitat range by a little bit in all four directions. (EEMS can be sensitive to border effects, so it is better to avoid placing samples at the border.) For example,

```
min(x)-1 max(x)+1
min(y)-1 max(y)+1
```

The values -1 and +1 might not be appropriate for a particular dataset; they should be chosen in view of the actual range of $x$ and $y$.

### Simulated datasets

The `data` folder contains six datasets simulated with `ms` and in EEMS format. (These datasets were used for Figure 2 in the EEMS paper.)

* barrier to migration, under three different sampling schemes
* uniform migration, under the same three sampling schemes

### Required input arguments

EEMS can be executed as follows:
```
MCMC_haploid(sourcepath,datapath,mcmcpath,xDemes,yDemes,...)
```

* `sourcepath` is the directory with the EEMS source scripts
* `datapath` is the path to the three EEMS input files (without file extensions)
* `mcmcpath` is the path to the EEMS output files
* `xDemes` and `yDemes` specify that the population graph is assumed to be a regular triangular grid of size xDemes-by-yDemes

### Optional input arguments

There are optional input arguments (that have default values but can be tuned to improve convergence). Here is a complete list:

Three parameters that specify the Markov Chain Monte Carlo sampling
* `numMCMCIter`: number of MCMC iterations
* `numBurnIter`: number of burn-in iterations to discard at the start
* `numThinIter`: number of iterations to thin between two writing steps

Variances for the proposal distributions:
* Proposal distributions for the parameters of the Voronoi tessellation which describes the effective migration rates (between demes) and thus the differences between two different demes
  * `mEffctProposalS2`:  variance for the cell effects
  * `mSeedsProposalS2`:  variance for the cell locations
  * `mrateMuProposalS2`: variance for the overall migration rate on the log10 scale
* Proposal distributions for the parameters of the Voronoi tessellation which describes the effective diversity rates (within demes) and thus the differences between two different samples from the same deme
  * `qEffctProposalS2`:  variance for the cell effects
  * `qSeedsProposalS2`:  variance for the cell locations
* `dfProposalS2`: proposal variance for the degrees of freedom parameter; the default is sqrt(nSNPs)

Hyperparameters:
* `mrateShape,mrateScale`: hyperparameters for the variance omega_m^2 of the migration cell effects
* `qrateShape,qrateScale`: hyperparameters for the variance omega_q^2 of the diversity cell effects
* `s2locShape,s2locShape`: hyperparameters for the scale parameter sigma^2
* `negBiSize`: number of failures for the Negative-Binomial prior on the number of Voronoi tiles
* `negBiProb`: success probability for the Negative-Binomial prior on the number of Voronoi tiles

### Acknowledgements

We have used `ms` to simulate genotypes under Kimura's stepping-stone model.

* R. R. Hudson. Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18(2):337–338, 2002.
