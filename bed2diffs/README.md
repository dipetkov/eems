## bed2diffs

`bed2diffs` is a small program that reads the genotypes in plink binary format (i.e., from a set of bed/bim/fam files) and computes the average pairwise differences.

### Installation

`bed2diffs` uses the `libplinkio` C library to read plink binary data. The library can be downloaded here: https://github.com/fadern/libplinkio.

`bed2diffs` uses OpenMP to parallelize the computation of the pairwise differences `diffs`. Multithreading is useful if the data contains millions of SNPs.

After obtaining `libplinkio`, update the path to PLINKIO in the Makefile and run `make linux`. (I have tested `bed2diffs` only on a linux machine.)

### Usage

You must specify a plink dataset, and optionally, the number of threads.
```
./bed2diffs_v1
```
This produces a very, very short help message: 
```
## Usage: ./bed2diffs_v1 --bfile PlinkData 
## Options:
##   --nthreads K	Set number of OpenMP threads
```
There are two example datasets in plink binary format. These are actually the same genetic data, in two different modes.

The first dataset is in SNP-major mode. This means that the bed file stores the genotypes of all individuals for the first SNP, then the genotypes of all individuals for the second SNP, and so on. In other words, if the genetic data is represented as a matrix, then each row corresponds to a SNP.
```
./bed2diffs_v1 --bfile ../test/example-SNP-major-mode
## or
./bed2diffs_v1 --bfile ../test/example-SNP-major-mode --nthreads 2
```
`bed2diffs_v1` computes the genetic dissimilarity matrix for 6 samples at 6 SNPs:
```
## Detected plink dataset ../test/example-SNP-major-mode.[bed/bim/fam] with 6 samples and 6 SNPs
## Computed average pairwise differences across 6 SNPs

cat ../test/example-SNP-major-mode.diffs 
## 0 1.50 1.00 0.50 0.50 2.00
## 1.50 0 0.25 0.50 1.00 0.75
## 1.00 0.25 0 0.33 0.67 1.00
## 0.50 0.50 0.33 0 0.33 0.67
## 0.50 1.00 0.67 0.33 0 0.33
## 2.00 0.75 1.00 0.67 0.33 0
```
Actually, `bed2diffs_v1` reports the dissimilarities with higher precision: 12 digits after the decimal point. I have rounded to 2 digits after the decimal point for more readability. 

The second dataset is in sample-major mode. This means that the bed file stores the genotypes of the first individual at all SNPs, then the genotypes of the second individuals at all SNPs, an so on. In other words, if the genetic data is represented as a matrix, then each row corresponds to a sample.
```
./bed2diffs_v1 --bfile ../test/example-sample-major-mode
```
`bed2diffs_v1` expects that the data is in SNP-major mode and so the above produces an error message:
```
## [Data::getsize] This program requires plink files [bed/bim/fam] in SNP-major mode
```
`plink` can transpose the data into SNP-major mode, which is the default:
```
plink --bfile ../test/example-sample-major-mode --transpose --make-bed --out ../test/example-sample-major-mode-transposed
```
Finally, `bed2diffs_v1` produces the same genetic dissimilarity matrix from the transposed data.
```
./bed2diffs_v1 --bfile ../test/example-sample-major-mode-transposed

cat ../test/example-sample-major-mode-transposed.diffs 
## 0 1.50 1.00 0.50 0.50 2.00
## 1.50 0 0.25 0.50 1.00 0.75
## 1.00 0.25 0 0.33 0.67 1.00
## 0.50 0.50 0.33 0 0.33 0.67
## 0.50 1.00 0.67 0.33 0 0.33
## 2.00 0.75 1.00 0.67 0.33 0
```

### Output

bed2diffs generates two files. The matrix of average pairwise differences is written to a text file without row names, column names or comments, with extension `diffs`. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the `diffs` matrix is the same as in the `fam` file, but in any case, the order is explicitly written to a text file with one sample per line, with extension `order`.


### Two versions

There are two versions, `bed2diffs_v1` and `bed2diffs_v2`. If there is missing data, the two versions will produce different results.

* In version 1, the average difference between samples i and j is computed across SNPs where both i and j are called. This version will also output the matrix of counts (number of SNPs) where both i and j are called.
* In version 2, the average difference between samples i and j is computed across all SNPs, and missing genotypes at a given SNP are "imputed" as the average genotype at that SNP.

We used `bed2diffs_v1` to compute genetic dissimilarity matrices for the EEMS paper. However, this version is not guaranteed to produce a matrix that is (mathematically) an Euclidean distance matrix, especially, if genotypes are not missing at random.

Here is again the dissimilarity matrix for the example dataset, computed with `bed2diffs_v1`. The longish message reminds us the formula used by `bed2diffs_v1`.
```
./bed2diffs_v1 --bfile ../test/example-SNP-major-mode
## Set the number of OpenMP threads to 1
## 
## Compute the average genetic differences according to: 
##   Dij = (1/|Mij|) sum_{m in Mij} (z_{im} - z_{jm})^2
##   where Mij is the set of SNPs where both i and j are called
## 
## Detected plink dataset ../test/example-SNP-major-mode.[bed/bim/fam] with 6 samples and 3 SNPs
## Computed average pairwise differences across 3 SNPs

cat ../test/example-SNP-major-mode.diffs 
## 0 1.50 1.00 0.50 0.50 2.00
## 1.50 0 0.25 0.50 1.00 0.75
## 1.00 0.25 0 0.33 0.67 1.00
## 0.50 0.50 0.33 0 0.33 0.67
## 0.50 1.00 0.67 0.33 0 0.33
## 2.00 0.75 1.00 0.67 0.33 0
```
By design, there is a lot of missingness in this dataset. The `diffs` matrix (version 1) has the following eigenvalues:
```
-2.2136305, -1.4233232, -0.3067077, -0.2311668,  0.0991569,  4.0756712
```
So even though the matrix is symmetric, with 0s on the main diagonal and positive entries off the diagonal, it is not an Euclidean matrix because it has two positive eigenvalues.

The alternative `bed2diffs_v2` always produces an Euclidean distance matrix, because it "imputes" missing genotypes with the observed average genotype. However, the imputation step might not be reasonable if genotypes are not missing at random.

Here is the dissimilarity matrix for the example dataset, computed with `bed2diffs_v2`, which uses a different formula.
```
./bed2diffs_v2 --bfile ../test/example-SNP-major-mode
## Set the number of OpenMP threads to 1
## 
## Compute the average genetic differences according to: 
##   Dij = (1/|Mtot|) sum_{m in Mtot} (z*_{im} - z*_{jm})^2
##   where Mtot is the set of all SNPs and
##   z*_{im} = z_{im} if sample i is called at marker m
##           = zbar_m (the average genotype at m) otherwise
## 
## Detected plink dataset ../test/example-SNP-major-mode.[bed/bim/fam] with 6 samples and 6 SNPs
## Computed average pairwise differences across 6 SNPs

cat ../test/example-SNP-major-mode.diffs 
##  0 2.00 1.00 1.67 1.67 2.67
##  2.00 0 0.33 0.33 0.33 0.67
##  1.00 0.33 0 0.67 0.67 1.00
##  1.67 0.33 0.67 0 0.00 0.33
##  1.67 0.33 0.67 0.00 0 0.33
##  2.67 0.67 1.00 0.33 0.33 0
```
The `diffs` matrix (version 2) has the following eigenvalues:
```
-1.6992500, -1.0755689, -0.6971778, -0.1489801, -0.1047345,  3.7257113
```
This matrix is a valid Euclidean distance matrix as it has one positive eigenvalue and the rest of the eigenvalues are negative. (The eigenvalues sum up to 0, which is another property of a distance matrix.)

However, rather than using `bed2diffs_v2`, it might be better to clean the data beforehand and to remove SNPs with high missingness, as it is usually done before any analysis of population structure. And then use `bed2diffs_v1`.

See bed2diffs-doc.pdf for a slightly longer explanation about the difference between the two versions of `bed2diffs`.
