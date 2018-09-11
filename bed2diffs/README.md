bed2diffs
---------

`bed2diffs` is a small program that reads genetic data in plink binary format (i.e., from a set of three files with extensions bed/bim/fam) and computes the average genetic dissimilarity matrix.

### Compilation

`bed2diffs` uses the [libplinkio](https://github.com/fadern/libplinkio) library to read genotype data stored in plink binary format. To install libplinkio, first clone the GitHub repository and get the latest version (commit 781e9ee37076).

```
git clone https://github.com/mfranberg/libplinkio
cd libplinkio
git checkout 781e9ee37076
```

Then follow the instructions to install libplinkio to a custom location /path/to/plinkio. Finally, update `PLINKIO` in the Makefile in `src` and `src-without-openmp` directories.

```
mkdir build
cd build
../configure --prefix=/path/to/plinkio
make && make check && make install
```

Optionally, `bed2diffs` uses OpenMP to parallelize the computation of the pairwise differences. Multithreading is useful if the data contains millions of SNPs. Choose the `src` directory to compile `bed2diffs` with OpenMP support; otherwise, choose the `src-wout-openmp` directory. Finally, to compile, use `make linux` on a Linux machine or `make darwin` on a Mac.

### Usage

You must specify the name of a binary plink file, and optionally, the number of threads. On the command line, type:

    ./src/bed2diffs_v1

This produces a very, very short help message:

    ## Usage: bed2diffs_v1 --bfile PlinkData
    ## Options:
    ##   --nthreads K   Set the number of OpenMP threads

In the `test` directory there are two example datasets in plink binary format. These actually store the same genetic data, in two different modes.

The first dataset is in SNP-major mode. This means that the bed file stores the genotypes of all individuals for the first SNP, then the genotypes of all individuals for the second SNP, and so on. In other words, if the genetic data is represented as a matrix, then each row corresponds to a SNP. To compute the genetic dissimilarity matrix, type

    ./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode

or, to use two threads, type

    ./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode --nthreads 2

    ## Set the number of OpenMP threads to 1
    ##
    ## Compute the average genetic differences according to:
    ##   Dij = (1/|Mij|) sum_{m in Mij} (z_{im} - z_{jm})^2
    ##   where Mij is the set of SNPs where both i and j are called
    ##
    ## Detected plink dataset ./test/example-SNP-major-mode.[bed/bim/fam] with 6 samples and 6 SNPs
    ## Computed average pairwise differences across 6 SNPs

Let's load the dissimilarity matrix in `R`.

``` r
diffs <- read.table("./test/example-SNP-major-mode.diffs")
```

``` r
diffs
```

    ##      IID1 IID2 IID3 IID4 IID5 IID6
    ## IID1  0.0 1.50 1.00 0.50 0.50 2.00
    ## IID2  1.5 0.00 0.25 0.50 1.00 0.75
    ## IID3  1.0 0.25 0.00 0.33 0.67 1.00
    ## IID4  0.5 0.50 0.33 0.00 0.33 0.67
    ## IID5  0.5 1.00 0.67 0.33 0.00 0.33
    ## IID6  2.0 0.75 1.00 0.67 0.33 0.00

The matrix is nonnegative, symmetric with zeros on the diagonal because `diffs[i, j]` is the squared gentics difference between individuals i and j, average across the SNPs. (These conditions are are necessary but not sufficient for a matrix to be a valid distance matrix. See the note on the two versions of `bed2diffs` below.)

The second dataset is in sample-major mode. This means that the bed file stores the genotypes of the first individual at all SNPs, then the genotypes of the second individuals at all SNPs, and so on. In other words, if the genetic data is represented as a matrix, then each row corresponds to a sample. `bed2diffs` expects that the data is in SNP-major format.

    ./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode

    ## Set the number of OpenMP threads to 1
    ##
    ## Compute the average genetic differences according to:
    ##   Dij = (1/|Mij|) sum_{m in Mij} (z_{im} - z_{jm})^2
    ##   where Mij is the set of SNPs where both i and j are called
    ##
    ## [Data::getsize] bed2diffs requires plink files [bed/bim/fam] in SNP-major mode

We can use `plink` can transpose the data into SNP-major mode, which is the default:

    plink --bfile ./test/example-sample-major-mode --transpose --make-bed --out ./test/example-sample-major-mode-transposed

    ##
    ## @----------------------------------------------------------@
    ## |        PLINK!       |     v1.07      |   10/Aug/2009     |
    ## |----------------------------------------------------------|
    ## |  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
    ## |----------------------------------------------------------|
    ## |  For documentation, citation & bug-report instructions:  |
    ## |        http://pngu.mgh.harvard.edu/purcell/plink/        |
    ## @----------------------------------------------------------@
    ##
    ## Skipping web check... [ --noweb ]
    ## Writing this text to log file [ ./test/example-sample-major-mode-transposed.log ]
    ## Analysis started: Thu Aug  3 00:41:18 2017
    ##
    ## Options in effect:
    ##  --noweb
    ##  --bfile ./test/example-sample-major-mode
    ##  --transpose
    ##  --recode
    ##  --make-bed
    ##  --out ./test/example-sample-major-mode-transposed
    ##
    ## ** For gPLINK compatibility, do not use '.' in --out **
    ## Reading map (extended format) from [ ./test/example-sample-major-mode.bim ]
    ## 6 markers to be included from [ ./test/example-sample-major-mode.bim ]
    ## Reading pedigree information from [ ./test/example-sample-major-mode.fam ]
    ## 6 individuals read from [ ./test/example-sample-major-mode.fam ]
    ## 1 individuals with nonmissing phenotypes
    ## Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
    ## Missing phenotype value is also -9
    ## 1 cases, 0 controls and 5 missing
    ## 4 males, 0 females, and 2 of unspecified sex
    ## Warning, found 2 individuals with ambiguous sex codes
    ## Writing list of these individuals to [ ./test/example-sample-major-mode-transposed.nosex ]
    ## Reading genotype bitfile from [ ./test/example-sample-major-mode.bed ]
    ## Detected that binary PED file is v1.00 individual-major mode
    ## Before frequency and genotyping pruning, there are 6 SNPs
    ## Converting data to SNP-major format
    ## 4 founders and 2 non-founders found
    ## Total genotyping rate in remaining individuals is 0.722222
    ## 0 SNPs failed missingness test ( GENO > 1 )
    ## 0 SNPs failed frequency test ( MAF < 0 )
    ## Converting data to Individual-major format
    ## After frequency and genotyping pruning, there are 6 SNPs
    ## After filtering, 1 cases, 0 controls and 5 missing
    ## After filtering, 4 males, 0 females, and 2 of unspecified sex
    ## Writing pedigree information to [ ./test/example-sample-major-mode-transposed.fam ]
    ## Writing map (extended format) information to [ ./test/example-sample-major-mode-transposed.bim ]
    ## Writing genotype bitfile to [ ./test/example-sample-major-mode-transposed.bed ]
    ## Using (default) SNP-major mode
    ## Converting data to SNP-major format
    ##
    ## Analysis finished: Thu Aug  3 00:41:18 2017

Finally, `bed2diffs_v1` produces the same genetic dissimilarity matrix from the transposed data.

    ./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode-transposed

``` r
diffs <- read.table("./test/example-sample-major-mode-transposed.diffs")
```

``` r
diffs
```

    ##      IID1 IID2 IID3 IID4 IID5 IID6
    ## IID1  0.0 1.50 1.00 0.50 0.50 2.00
    ## IID2  1.5 0.00 0.25 0.50 1.00 0.75
    ## IID3  1.0 0.25 0.00 0.33 0.67 1.00
    ## IID4  0.5 0.50 0.33 0.00 0.33 0.67
    ## IID5  0.5 1.00 0.67 0.33 0.00 0.33
    ## IID6  2.0 0.75 1.00 0.67 0.33 0.00

### Output

`bed2diffs` generates two files. The matrix of average pairwise differences is written to a text file without row names, column names or comments, with extension `diffs`. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the `diffs` matrix is the same as in the `fam` file, but in any case, the order is explicitly written to a text file with one sample per line, with extension `order`.

### The two versions of bed2diffs

There are two versions, `bed2diffs_v1` and `bed2diffs_v2`. If there is missing data, the two versions will produce different results.

-   In version 1, the average difference between samples i and j is computed across SNPs where both i and j are called. This version will also output the matrix of counts (number of SNPs) where both i and j are called.
-   In version 2, the average difference between samples i and j is computed across all SNPs, and missing genotypes at a given SNP are "imputed" as the average genotype at that SNP.

We used `bed2diffs_v1` to compute genetic dissimilarity matrices for the EEMS paper. However, this version is not guaranteed to produce a matrix that is (mathematically) an Euclidean distance matrix, especially, if genotypes are not missing at random.

By design, there is a lot of missingness in the example dataset. For completeness, here is the genotype matrix in SNP-major mode, before dealing with the missing calls:

    ##      IID1 IID2 IID3 IID4 IID5 IID6
    ## snp1    2    0   NA   NA   NA   NA
    ## snp2    2   NA    1    2    2    2
    ## snp3    2    1    1   NA   NA    0
    ## snp4    1    1    2    2    2    1
    ## snp5   NA    2    2    2    1    1
    ## snp6    2    1    1   NA   NA    0

The `diffs` matrix (version 1) has the following eigenvalues:

``` r
eigenvals <- eigen(diffs)$values
```

    ## -2.21 -1.43 -0.31 -0.23 0.1 4.08

So even though the matrix is symmetric, with zeros on the main diagonal and positive entries off the diagonal, it is not an Euclidean matrix because it has two positive eigenvalues. An Euclidean distance matrix is conditionally negative definite and thus it has one positive eigenvalue and n-1 negative eigenvalues (1).

The alternative `bed2diffs_v2` produces an Euclidean distance matrix by construction, because it "imputes" missing genotypes with the observed average genotype (2). UPDATE: Fixed a bug in v2 on 18 May 2016.

Here is the completed genotype matrix after imputation.

    ##      IID1 IID2 IID3 IID4 IID5 IID6
    ## snp1  2.0  0.0    1    1    1    1
    ## snp2  2.0  1.8    1    2    2    2
    ## snp3  2.0  1.0    1    1    1    0
    ## snp4  1.0  1.0    2    2    2    1
    ## snp5  1.6  2.0    2    2    1    1
    ## snp6  2.0  1.0    1    1    1    0

Here is the dissimilarity matrix for the example dataset, computed with `bed2diffs_v2`.

    ./src/bed2diffs_v2 --bfile ./test/example-SNP-major-mode

    ## Set the number of OpenMP threads to 1
    ##
    ## Compute the average genetic differences according to:
    ##   Dij = (1/|Mtot|) sum_{m in Mtot} (z*_{im} - z*_{jm})^2
    ##   where Mtot is the set of all SNPs and
    ##   z*_{im} = z_{im} if sample i is called at marker m
    ##           = zbar_m (the average genotype at m) otherwise
    ##
    ## Detected plink dataset ./test/example-SNP-major-mode.[bed/bim/fam] with 6 samples and 6 SNPs
    ## Computed average pairwise differences across 6 SNPs

``` r
diffs <- read.table("./test/example-SNP-major-mode.diffs")
```

``` r
diffs
```

    ##      IID1 IID2 IID3 IID4 IID5 IID6
    ## IID1 0.00 1.03 0.86 0.69 0.73 1.56
    ## IID2 1.03 0.00 0.44 0.34 0.51 0.67
    ## IID3 0.86 0.44 0.00 0.17 0.33 0.83
    ## IID4 0.69 0.34 0.17 0.00 0.17 0.67
    ## IID5 0.73 0.51 0.33 0.17 0.00 0.50
    ## IID6 1.56 0.67 0.83 0.67 0.50 0.00

The `diffs` matrix (version 2) has the following eigenvalues:

``` r
eigenvals <- eigen(diffs)$values
```

    ## -1.63 -0.9 -0.55 -0.26 -0.1 3.43

This matrix is a valid Euclidean distance matrix as it has one positive eigenvalue and the rest of the eigenvalues are negative. (The eigenvalues sum up to 0, which is another property of a distance matrix.)

However, the imputation performed by `bed2diffs_v2` would not be appropriate if genotypes are not missing at random. Therefore, rather than using `bed2diffs_v2`, it would be better to clean the data beforehand and to remove SNPs with high missingness, as it is usually done before any analysis of population structure. And then use `bed2diffs_v1`.

See Documentation/bed2diffs-doc.pdf for a slightly longer explanation about the difference between the two versions of `bed2diffs`.

#### References

1.  Balaji and Bapat. On Euclidean distance matrices. *Linear Algebra Appl*, 424:108-117, 2007.
2.  By plugging the observed genotype frequency for any missing calls, we construct a representation *X*<sub>*n* × *p*</sub> of *n* individuals in *p*-dimensional space such that *D*<sub>*i*, *j*</sub> = ||*x*<sub>*i*</sub> − *x*<sub>*j*</sub>||<sub>2</sub><sup>2</sup> is the squared Euclidean distance between *i* and *j*.

#### bed2diffs in R

``` r
# Use the "pairwise.complete.obs" method to compute pairwise dissimilarities
# This straightforward implementation
# uses a double loop, so would be slow if the sample size is large.
bed2diffs_v1 <- function(genotypes) {

  nIndiv <- nrow(genotypes)
  nSites <- ncol(genotypes)
  diffs <- matrix(0, nIndiv, nIndiv)

  for (i in seq(nIndiv - 1)) {
    for (j in seq(i + 1, nIndiv)) {
      x <- genotypes[i, ]
      y <- genotypes[j, ]
      diffs[i, j] <- mean((x - y)^2, na.rm = TRUE)
      diffs[j, i] <- diffs[i, j]
    }
  }

  diffs
}

# Compute the diffs matrix using the "mean allele frequency"
# imputation method
bed2diffs_v2 <- function(genotypes) {

  nIndiv <- nrow(genotypes)
  nSites <- ncol(genotypes)
  missing <- is.na(genotypes)

  ## Impute NAs with the column means (= twice the allele frequencies)
  geno_means <- colMeans(genotypes, na.rm = TRUE)
  # nIndiv rows of genotype means
  geno_means <- matrix(geno_means, nrow = nIndiv, ncol = nSites, byrow = TRUE)

  ## Set the means which correspond to observed genotypes to 0
  geno_means[missing == FALSE] <- 0
  ## Set the missing genotypes to 0 (used to be NA)
  genotypes[missing == TRUE] <- 0
  genotypes <- genotypes + geno_means

  similarities <- genotypes %*% t(genotypes) / nSites
  self_similarities <- diag(similarities)
  vector1s <- rep(1, nIndiv)

  diffs <-
    self_similarities %*% t(vector1s) +
    vector1s %*% t(self_similarities) - 2 * similarities
  diffs
}
```
