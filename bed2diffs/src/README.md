## bed2diffs

bed2diffs is a small program that reads the genotypes in plink binary format (i.e., a set of bed/bim/fam files) and computes the average pairwise differences.

### Installation

bed2diffs uses the libplinkio C library to read a plink binary dataset. The library can be downloaded from the https://github.com/fadern/libplinkio.

bed2diffs uses OpenMP to parallelize the computation. Multithreading can be useful if the data contains millions of SNPs.

After obtaining libplinkio, update the path to PLINKIO in the Makefile and run `make`.

### Usage

```
./bed2diffs
# Usage: ./bed2diffs --bfile PlinkData 
# Options:
#   --nthreads K	Set number of OpenMP threads

./bed2diffs --bfile test/example
```

### Output

bed2diffs generates two files, example.diffs and example.order. The matrix of average pairwise differences is written to a text file without a header or any other information, with extension `diffs`. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the diffs matrix is the same as in the fam file, but in any case, the order is explicitly written to a text file with one sample per line, with extension `order`.
