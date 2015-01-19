## bed2diffs

bed2diffs is a small program that reads the genotypes in plink binary format (i.e., a set of bed/bim/fam files) and computes the average pairwise differences.

### Installation

bed2diffs uses the libplinkio C library to read a plink binary dataset. The library can be downloaded from the https://github.com/fadern/libplinkio.

bed2diffs uses OpenMP to parallelize the computation of `diffs`. Multithreading can be useful if the data contains millions of SNPs.

After obtaining libplinkio, update the path to PLINKIO in the Makefile and run `make`.

### Usage

```
./bed2diffs_v1
# Usage: ./bed2diffs_v1 --bfile PlinkData 
# Options:
#   --nthreads K	Set number of OpenMP threads

./bed2diffs_v1 --bfile test/example
```

### Output

bed2diffs generates two files, example.diffs and example.order. The matrix of average pairwise differences is written to a text file without a header or any other information, with extension `diffs`. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the `diffs` matrix is the same as in the `fam` file, but in any case, the order is explicitly written to a text file with one sample per line, with extension `order`.


### Two versions

There are two versions, `bed2diffs_v1` and `bed2diffs_v2`. If there is missing data, the two versions will produce different results because

* version 1: The average difference between samples i and j is computed across SNPs where both i and j are called. This version will also output the matrix of counts (number of SNPs) where both i and j are called.
* version 2: The average difference between samples i and j is computed across all SNPs, and missing genotypes at a given SNP are "imputed" as the average genotype at that SNP.

We used `bed2diffs_v1` to compute matrices of differences for the EEMS paper. However, this version is not guaranteed to produce a matrix that is (mathematically) an Euclidean distance matrix, especially, if genotypes are not missing at random.

The alternative `bed2diffs_v2` always produces an Euclidean distance matrix, as a consequence of imputing missing genotypes with the observed verage genotype. However, this imputation might not be reasonable if genotypes are not missing at random.

Therefore, the data should be cleaned up beforehand, to remove SNPs with high missingness, as it usually done before any analysis of population structure.
