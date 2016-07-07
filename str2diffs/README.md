Read genetic data in `structure` format, compute `diffs` matrix
===============================================================

This script shows how to use the `adegenet` package to load data in `structure` format into `R` and then compute the matrix `diffs` of genetic dissimilarities.

First I looked at `adegenet::dist.genpop()`. This function computes genetic distances between populations and there are various metrics, which is interesting, but these are between populations, not between individuals, so not what is required for EEMS. Instead, let's compute the dissimilarity matrix explicitly.

First I need to load the genotypes from a `structure` file. Fortunately `adegenet` also has a `read.structure()` function.

``` r
if (!require("adegenet")) install.packages("adegenet")
library("adegenet")
```

Here is the raw data. It is the same small example with 6 individuals and 6 SNPs used to illustrate `bed2diffs`.

``` r
raw <- read.table("example-structure.str")
colnames(raw) <- c("indiv", "label", "pop", paste0("snp", 1:6))
raw
```

    ##    indiv label pop snp1 snp2 snp3 snp4 snp5 snp6
    ## 1   IID1     1   1    1    1    1    0   -9    1
    ## 2   IID1     1   1    1    1    1    1   -9    1
    ## 3   IID2     1   1    0   -9    0    0    1    0
    ## 4   IID2     1   1    0   -9    1    1    1    1
    ## 5   IID3     1   1   -9    0    0    1    1    0
    ## 6   IID3     1   1   -9    1    1    1    1    1
    ## 7   IID4     2   1   -9    1   -9    1    1   -9
    ## 8   IID4     2   1   -9    1   -9    1    1   -9
    ## 9   IID5     2   1   -9    1   -9    1    0   -9
    ## 10  IID5     2   1   -9    1   -9    1    1   -9
    ## 11  IID6     2   1   -9    1    0    0    0    0
    ## 12  IID6     2   1   -9    1    0    1    1    0

Specify that there are 6 individuals (`n.ind = 6`), 6 SNPs (`n.loc = 6`), that the first column contains genotype labels (`col.lab = 1`) and that the second column contains population labels (`col.pop = 2`).

``` r
data <- read.structure("example-structure.str",
                       onerowperind = FALSE, 
                       n.ind = 6, n.loc = 6, col.lab = 1, col.pop = 2, ask = FALSE)
```

    ## 
    ##  Converting data from a STRUCTURE .stru file to a genind object...

``` r
data
```

    ## /// GENIND OBJECT /////////
    ## 
    ##  // 6 individuals; 6 loci; 12 alleles; size: 8.8 Kb
    ## 
    ##  // Basic content
    ##    @tab:  6 x 12 matrix of allele counts
    ##    @loc.n.all: number of alleles per locus (range: 2-2)
    ##    @loc.fac: locus factor for the 12 columns of @tab
    ##    @all.names: list of allele names for each locus
    ##    @ploidy: ploidy of each individual  (range: 2-2)
    ##    @type:  codom
    ##    @call: read.structure(file = "example-structure.str", n.ind = 6, n.loc = 6, 
    ##     onerowperind = FALSE, col.lab = 1, col.pop = 2, ask = FALSE)
    ## 
    ##  // Optional content
    ##    @pop: population of each individual (group size range: 3-3)

The 6 loci in this dataset are observed to be polymorphic and bi-allelic. In general, it is a good idea to remove SNPs with high missingness or with some other indication of low genotyping quality. But this example only shows how to compute the genetic dissimilarity matrix, not how to clean the data for population structure analysis.

Start with the genotype matrix.

``` r
Geno <- data@tab
```

EEMS works with either haploid/diploid data but not multi-allelic loci. There are no multi-allelic loci in this example but the following chunk would have removed them, if there were any.

``` r
## Keep SNPs that are observed to be bi-allelic.
multi.loci <- names(which(data@loc.n.all != 2))
## Explanation: 
## Suppose we want to remove locus, `L100` with alleles, `L100.00`, `L100.01`, `L100.02`, 
## then detect columns whose names matches the regular expression `^L100\\.\\d+$`
multi.cols <- which(str_detect(colnames(Geno), 
                               paste0("^", multi.loci, "\\.\\d+$", collapse = "|")))
if (length(multi.cols)) Geno <- Geno[, - multi.cols]
dim(Geno)
```

    ## [1]  6 12

There are six individuals at six polymorphic bi-allelic SNPs (2 \* 6 alleles). So `Geno` has 12 columns, one for each of 12 alleles.

Let's convert the matrix to 0-1 labeling. I arbitrarily choose one allele to be the "derived" allele and, for each individual, count how many copies of the derived allele it carries. This is very easy if the tab matrix is of type "codom".

``` r
stopifnot(identical(data@type, 'codom'))
```

Since the labeling does not matter for computing differences, I pick the second allele to label as the "derived" allele. That is, I pick all loci whose name ends with `.01`.

``` r
Geno <- Geno[, str_detect(colnames(Geno), "\\.01$")]
```

Next compute the dissimilarity matrix, using the "mean allele frequency" imputation method, which corresponds to `bed2diffs-v2`. See the `bed2diffs` README for more information about the two methods of imputing missing that are implemented in `bed2diffs-v1` and `bed2diffs-v2`.

``` r
bed2diffs_v2 <- function(Geno) {
  nIndiv <- nrow(Geno)
  nSites <- ncol(Geno)
  Miss <- is.na(Geno)
  ## Impute NAs with the column means (= twice the allele frequencies)
  Mean <- matrix(colMeans(Geno, na.rm = TRUE), ## a row of means
                 nrow = nIndiv, ncol = nSites, byrow = TRUE) ## a matrix with nIndiv identical rows of means
  Mean[Miss == 0] <- 0 ## Set the means that correspond to observed genotypes to 0
  Geno[Miss == 1] <- 0 ## Set the missing genotypes to 0 (used to be NA) 
  Geno <- Geno + Mean
  ## Compute similarities
  Sim <- Geno %*% t(Geno) / nSites
  SelfSim <- diag(Sim) ## self-similarities
  vector1s <- rep(1, nIndiv) ## vector of 1s
  ## This chunk generates a `diffs` matrix
  Diffs <- SelfSim %*% t(vector1s) + vector1s %*% t(SelfSim) - 2 * Sim
  Diffs
}
```

Here is a function that implements the "pairwise.complete.obs" method, which corresponds to `bed2diffs-v1`. The straightforward implementation uses a double loop, so would be slow if the sample size is large.

``` r
bed2diffs_v1 <- function(Geno) {
  nIndiv <- nrow(Geno)
  nSites <- ncol(Geno)
  Diffs <- matrix(0, nIndiv, nIndiv)
  
  for (i in seq(nIndiv - 1)) {
    for (j in seq(i + 1, nIndiv)) {
      x <- Geno[i, ]
      y <- Geno[j, ]
      Diffs[i, j] <- mean((x - y)^2, na.rm = TRUE)
      Diffs[j, i] <- Diffs[i, j]
    }
  }
  Diffs
}
```

Let's compute both versions of the dissimilarity matrix and inspect the eigenvalues.

``` r
Diffs_v1 <- bed2diffs_v1(Geno)
Diffs_v2 <- bed2diffs_v2(Geno)
Diffs_v1 <- round(Diffs_v1, digits = 6)
Diffs_v2 <- round(Diffs_v2, digits = 6)
```

Check that the dissimilarity matrix has one positive eigenvalue and `nIndiv-1` negative eigenvalues, as required by a full-rank Euclidean distance matrix.

``` r
sort(round(eigen(Diffs_v1)$values, digits = 2))
```

    ## [1] -2.21 -1.42 -0.31 -0.23  0.10  4.08

``` r
sort(round(eigen(Diffs_v2)$values, digits = 2))
```

    ## [1] -1.63 -0.90 -0.55 -0.26 -0.10  3.44

The condition does not hold for version 1, so let's save version 2 and use it as input to `runeems`.

``` r
write.table(Diffs_v2, "example-structure.diffs", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
```