
#+ global_options, include = FALSE
suppressWarnings(suppressMessages(library(knitr)))
opts_chunk$set(echo = TRUE, warning = TRUE, message = TRUE)

#' ## bed2diffs
#' `bed2diffs` is a small program that reads genetic data in plink binary format (i.e., from a set of 
#' three files with extensions bed/bim/fam) and computes the average genetic dissimilarity matrix.
#' 
#' ### Compilation
#' 
#' `bed2diffs` uses the [libplinkio](https://github.com/fadern/libplinkio) library to read genotype
#' data stored in plink binary format. Use the instructions for installing to a custom location 
#' /path/to/plinkio. Then update `PLINKIO` in the Makefile in `src` and `src-without-openmp` directories.
#' 
#' Optionally, `bed2diffs` uses OpenMP to parallelize the computation of the pairwise differences.
#' Multithreading is useful if the data contains millions of SNPs. Choose the `src` directory to compile
#' `bed2diffs` with OpenMP support; otherwise, choose the `src-wout-openmp` directory. Finally, to 
#' compile, use `make linux` on a Linux machine or `make darwin` on a Mac.
#' 
#' ### Usage
#' 
#' You must specify the name of a binary plink file, and optionally, the number of threads. On the command
#' line, type:
#' ```
#' ./src/bed2diffs_v1
#' ```
#' This produces a very, very short help message:
#+ , echo = FALSE, warning = FALSE
out <- system("./src/bed2diffs_v1 2>&1", intern = TRUE)
cat(out, sep = "\n")

#' In the `test` directory there are two example datasets in plink binary format. These actually 
#' store the same genetic data, in two different modes.
#' 
#' The first dataset is in SNP-major mode. This means that the bed file stores the genotypes of all 
#' individuals for the first SNP, then the genotypes of all individuals for the second SNP, and so on. 
#' In other words, if the genetic data is represented as a matrix, then each row corresponds to a SNP.
#' To compute the genetic dissimilarity matrix, type
#' 
#' ```
#' ./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode
#' ```
#' or, to use two threads, type
#' ```
#' ./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode --nthreads 2
#' ```
#+ , echo = FALSE
out <- system("./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode", intern = TRUE)
cat(out, sep = "\n")

#' Let's load the dissimilarity matrix in `R`.
#+ 
diffs <- read.table("./test/example-SNP-major-mode.diffs")
#+ echo = FALSE
colnames(diffs) <- rownames(diffs) <- paste0("IID", 1:6)
diffs <- round(diffs, digits = 2)
#+ 
diffs

#' The matrix is nonnegative, symmetric with zeros on the diagonal because `diffs[i, j]` is the squared
#' gentics difference between individuals i and j, average across the SNPs. (These conditions are
#' are necessary but not sufficient for a matrix to be a valid distance matrix. See the note on the 
#' two versions of `bed2diffs` below.)
#' 
#' The second dataset is in sample-major mode. This means that the bed file stores the genotypes of 
#' the first individual at all SNPs, then the genotypes of the second individuals at all SNPs, and 
#' so on. In other words, if the genetic data is represented as a matrix, then each row corresponds 
#' to a sample. `bed2diffs` expects that the data is in SNP-major format.
#' ```
#' ./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode
#' ```
#+ , echo = FALSE, warning = FALSE
out <- system("./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode 2>&1", intern = TRUE)
cat(out, sep = "\n")

#' We can use `plink` can transpose the data into SNP-major mode, which is the default:
#' ```
#' plink --bfile ./test/example-sample-major-mode --transpose --make-bed --out ./test/example-sample-major-mode-transposed
#' ```
#+ , echo = FALSE
out <- system("plink --bfile ./test/example-sample-major-mode --transpose --make-bed --out ./test/example-sample-major-mode-transposed", intern = TRUE)
cat(out, sep = "\n")

#' Finally, `bed2diffs_v1` produces the same genetic dissimilarity matrix from the transposed data.
#' ```
#' ./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode-transposed
#' ```
#+ , echo = FALSE
system("./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode-transposed")
#+ 
diffs <- read.table("./test/example-sample-major-mode-transposed.diffs")
#+ echo = FALSE
colnames(diffs) <- rownames(diffs) <- paste0("IID", 1:6)
diffs <- round(diffs, digits = 2)
#+ 
diffs

#' ### Output
#' 
#' `bed2diffs` generates two files. The matrix of average pairwise differences is written to a text file 
#' without row names, column names or comments, with extension `diffs`. The size of the matrix is NxN 
#' where N is the number of samples; there are only 0s on the main diagonal. The order of the samples 
#' in the `diffs` matrix is the same as in the `fam` file, but in any case, the order is explicitly 
#' written to a text file with one sample per line, with extension `order`.
#' 
#' ### The two versions of bed2diffs
#' 
#' There are two versions, `bed2diffs_v1` and `bed2diffs_v2`. If there is missing data, the two versions 
#' will produce different results.
#' 
#' * In version 1, the average difference between samples i and j is computed across SNPs where both i 
#'   and j are called. This version will also output the matrix of counts (number of SNPs) where both i 
#'   and j are called.
#' * In version 2, the average difference between samples i and j is computed across all SNPs, and 
#'   missing genotypes at a given SNP are "imputed" as the average genotype at that SNP.
#'   
#' We used `bed2diffs_v1` to compute genetic dissimilarity matrices for the EEMS paper. However, 
#' this version is not guaranteed to produce a matrix that is (mathematically) an Euclidean distance
#' matrix, especially, if genotypes are not missing at random.
#' 
#' By design, there is a lot of missingness in the example dataset. For completeness, here is 
#' the genotype matrix in SNP-major mode, before dealing with the missing calls:
#+ , echo = FALSE
suppressWarnings(suppressMessages(library(MultiPhen)))
Geno <- read.plink("./test/example-SNP-major-mode")
attributes(Geno)$closeConnection <- NULL
colnames(Geno) <- paste0("snp", 1:6)
t(Geno)

#' The `diffs` matrix (version 1) has the following eigenvalues:
#+ 
eigenvals <- eigen(diffs)$values
#+ , echo = FALSE
eigenvals <- round(sort(eigenvals), digits = 2)
cat(eigenvals)

#' So even though the matrix is symmetric, with zeros on the main diagonal and positive entries off 
#' the diagonal, it is not an Euclidean matrix because it has two positive eigenvalues. (An
#' Euclidean distance matrix is conditionally negative definite and thus it has one positive 
#' eigenvalue and n-1 negative eigenvalues [1].)
#'
#' The alternative `bed2diffs_v2` produces an Euclidean distance matrix by construction, because it
#' "imputes" missing genotypes with the observed average genotype [2]. UPDATE: Fixed a bug in v2 on 
#' 18 May 2016. 
#' 
#' Here is the completed genotype matrix after imputation.
#+ , echo = FALSE
nIndiv <- nrow(Geno)
nSites <- ncol(Geno)
Miss <- is.na(Geno)
## Impute NAs with the column means (= twice the allele frequencies)
Mean <- matrix(colMeans(Geno, na.rm = TRUE), ## a row of means
               nrow = nIndiv, ncol = nSites, byrow = TRUE) ## a matrix with nIndiv identical rows of means
Mean[Miss == 0] <- 0 ## Set the means that correspond to observed genotypes to 0
Geno[Miss == 1] <- 0 ## Set the missing genotypes to 0 (used to be NA) 
Geno <- Geno + Mean
t(Geno)

#' Here is the dissimilarity matrix for the example dataset, computed with `bed2diffs_v2`.
#' ```
#' ./src/bed2diffs_v2 --bfile ./test/example-SNP-major-mode
#' ```
#+ , echo = FALSE
out <- system("./src/bed2diffs_v2 --bfile ./test/example-SNP-major-mode", intern = TRUE)
cat(out, sep = "\n")
#+ 
diffs <- read.table("./test/example-SNP-major-mode.diffs")
#+ echo = FALSE
colnames(diffs) <- rownames(diffs) <- paste0("IID", 1:6)
diffs <- round(diffs, digits = 2)
#+ 
diffs

#' The `diffs` matrix (version 2) has the following eigenvalues:
#+ 
eigenvals <- eigen(diffs)$values
#+ , echo = FALSE
eigenvals <- round(sort(eigenvals), digits = 2)
cat(eigenvals)

#' This matrix is a valid Euclidean distance matrix as it has one positive eigenvalue and the rest 
#' of the eigenvalues are negative. (The eigenvalues sum up to 0, which is another property of a 
#' distance matrix.)
#' 
#' However, the imputation performed by `bed2diffs_v2` would not be appropriate if genotypes are 
#' not missing at random. Therefore, rather than using `bed2diffs_v2`, it would be better to clean 
#' the data beforehand and to remove SNPs with high missingness, as it is usually done before 
#' any analysis of population structure. And then use `bed2diffs_v1`.
#' 
#' See Documentation/bed2diffs-doc.pdf for a slightly longer explanation about the difference between 
#' the two versions of `bed2diffs`.
#' 
#' #### References
#' 1. Balaji and Bapat. On Euclidean distance matrices. *Linear Algebra Appl*, 424:108-117, 2007.
#' 2. By plugging the observed genotype frequency for any missing calls, we construct a representation
#'    $X_{n \times p}$ of $n$ individuals in $p$-dimensional space such that 
#'    $D_{i, j} = ||x_i - x_j||_2^2$ is the squared Euclidean distance between $i$ and $j$.
