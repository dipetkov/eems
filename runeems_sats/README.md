
## runeems_sats

The program `runeems_sats` implements the EEMS method for analyzing spatial population structure. This version uses raw microsatellite data.

* `src` contains the C++ source code, which needs to be compiled.
* `data` contains simulated data generated with [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) to illustrate the input file format.

For help on building `runeems_sats` see the [documentation](README.md).

### Input data format (microsatellites)

`runeems_sats` requires three data input files that have the same file name but different extension. The description below assumes that `datapath` is the full path + the file name (but without the extension).

#### datapath.sites

`datapath.sites` is the matrix of allele copies; missing alleles are specified by any negative number. For n individuals at L loci, the `sites` matrix is n-by-L for L haploid markers and n-by-2L for L diploid markers.

#### datapath.coord

`datapath.coord` are the sample coordinates, two coordinates per sample, one sample per line. The sampling locations should be given in the same order as the rows and columns of the dissimilarity matrix.

#### datapath.outer

`datapath.outer` are the habitat coordinates, as a sequence of vertices that form a closed polygon. The habitat vertices should be listed *counterclockwise* and the first vertex should also be the last vertex, so that the outline is a *closed ring*. Otherwise, EEMS attempts to "correct" the polygon and prints a warning message.
