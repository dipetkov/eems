## bed2diffs

`bed2diffs` is a small program that reads the genotypes in plink binary format (i.e., from a set of bed/bim/fam files) and computes the average pairwise differences.

### Installation

`bed2diffs` uses the `libplinkio` C library to read plink binary data. The library can be downloaded here: https://github.com/fadern/libplinkio.

`bed2diffs` uses OpenMP to parallelize the computation of the pairwise differences `diffs`. Multithreading is useful if the data contains millions of SNPs.

After obtaining `libplinkio`, update the path to PLINKIO in the Makefile and run `make linux`. (I have tested `bed2diffs` only on a linux machine.)

### Usage

```
./bed2diffs_v1
# Usage: ./bed2diffs_v1 --bfile PlinkData 
# Options:
#   --nthreads K	Set number of OpenMP threads

./bed2diffs_v1 --bfile test/example
```

### Output

bed2diffs generates two files. The matrix of average pairwise differences is written to a text file without row names, column names or comments, with extension `diffs`. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the `diffs` matrix is the same as in the `fam` file, but in any case, the order is explicitly written to a text file with one sample per line, with extension `order`.


### Two versions

There are two versions, `bed2diffs_v1` and `bed2diffs_v2`. If there is missing data, the two versions will produce different results.

* In version 1, the average difference between samples i and j is computed across SNPs where both i and j are called. This version will also output the matrix of counts (number of SNPs) where both i and j are called.
* In version 2, the average difference between samples i and j is computed across all SNPs, and missing genotypes at a given SNP are "imputed" as the average genotype at that SNP.

We used `bed2diffs_v1` to compute genetic dissimilaritiy matrices for the EEMS paper. However, this version is not guaranteed to produce a matrix that is (mathematically) an Euclidean distance matrix, especially, if genotypes are not missing at random.

The alternative `bed2diffs_v2` always produces an Euclidean distance matrix, as a consequence of imputing missing genotypes with the observed average genotype. However, the imputation step might not be reasonable if genotypes are not missing at random.

Therefore, the data should be cleaned up beforehand, to remove SNPs with high missingness, as it is usually done before any analysis of population structure.
