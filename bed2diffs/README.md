## bed2diffs

bed2diffs is a small program that reads in plink files in bed/bim/fam format and computes the average pairwise differences.

### Installation

For reading the plink genotype files, bed2diffs requires the libplinkio C library, which has to be installed first and which can be downloaded from the https://github.com/fadern/libplinkio.

For the actual computation, bed2diffs requires the Eigen template library, which does not need to be installed and  which can be downloaded from http://eigen.tuxfamily.org.

After obtaining the required libraries, update the paths to PLINKIO and EIGEN in the Makefile and run `make`.

### Usage

```
./bed2diffs
# Usage: ./bed2diffs --bfile PlinkData 
./bed2diffs --bfile test/example
```

### Output

bed2diffs generates two files, example.diffs and example.order. The matrix of average pairwise differences is saved in a text file without header or any other information, with extension `diffs`. The size of the matrix is $n\times n$ where $n$ is the number of samples. The order of the samples in the diffs matrix is the same as in the fam file, but in any case, the order is explicitly saved in a text file with one sample per line, with extension `order`.
