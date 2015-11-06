
## runeems_snps

The program `runeems_snps` implements the EEMS method for analyzing spatial population structure. This version uses the pairwise genetic dissimilarity matrix computed from SNP data.

* `src` contains the C++ source code, which needs to be compiled.
* `data` contains simulated data generated with [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) to illustrate the input file format.

For help on building `runeems_snps` see the instruction manual.

### Input data format (SNPs)

`runeems_snps` requires three data input files that have the same file name but different extension. The description below assumes that `datapath` is the full path + the file name (but without the extension).

#### datapath.diffs

`datapath.diffs` is the matrix of average pairwise genetic dissimilarities. This can be computed with [bed2diffs](bed2diffs/README.md) from genetic data in plink binary format.

The dissimilarity matrix is nonnegative, symmetric, with 0s on the main diagonal. These conditions are necessary but not sufficient for `diffs` to be a valid dissimilarity matrix. Mathematically, `diffs` should be conditionally negative definite.

#### datapath.coord

`datapath.coord` are the sample coordinates, two coordinates per sample, one sample per line. The sampling locations should be given in the same order as the rows and columns of the dissimilarity matrix.

#### datapath.outer

`datapath.outer` are the habitat coordinates, as a sequence of vertices that form a closed polygon. The habitat vertices should be listed *counterclockwise* and the first vertex should also be the last vertex, so that the outline is a *closed ring*. Otherwise, EEMS attempts to "correct" the polygon and prints a warning message.
