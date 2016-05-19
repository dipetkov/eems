


## bed2diffs
`bed2diffs` is a small program that reads genetic data in plink binary format (i.e., from a set of 
three files with extensions bed/bim/fam) and computes the average genetic dissimilarity matrix.

### Compilation

`bed2diffs` uses the [libplinkio](https://github.com/fadern/libplinkio) library to read genotype
data stored in plink binary format. Use the instructions for installing to a custom location 
/path/to/plinkio. Then update `PLINKIO` in the Makefile in `src` and `src-without-openmp` directories.

Optionally, `bed2diffs` uses OpenMP to parallelize the computation of the pairwise differences.
Multithreading is useful if the data contains millions of SNPs. Choose the `src` directory to compile
`bed2diffs` with OpenMP support; otherwise, choose the `src-wout-openmp` directory. Finally, to 
compile, use `make linux` on a Linux machine or `make darwin` on a Mac.

### Usage

You must specify the name of a binary plink file, and optionally, the number of threads. On the command
line, type:
```
./src/bed2diffs_v1
```
This produces a very, very short help message:


```
## Usage: bed2diffs_v1 --bfile PlinkData
## Options:
##   --nthreads K	Set the number of OpenMP threads
```

In the `test` directory there are two example datasets in plink binary format. These actually 
store the same genetic data, in two different modes.

The first dataset is in SNP-major mode. This means that the bed file stores the genotypes of all 
individuals for the first SNP, then the genotypes of all individuals for the second SNP, and so on. 
In other words, if the genetic data is represented as a matrix, then each row corresponds to a SNP.
To compute the genetic dissimilarity matrix, type

```
./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode
```
or, to use two threads, type
```
./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode --nthreads 2
```


```
## Set the number of OpenMP threads to 1
## 
## Compute the average genetic differences according to: 
##   Dij = (1/|Mij|) sum_{m in Mij} (z_{im} - z_{jm})^2
##   where Mij is the set of SNPs where both i and j are called
## 
## Detected plink dataset ./test/example-SNP-major-mode.[bed/bim/fam] with 6 samples and 6 SNPs
## Computed average pairwise differences across 6 SNPs
```

Let's load the dissimilarity matrix in `R`.


```r
diffs <- read.table("./test/example-SNP-major-mode.diffs")
```


```r
diffs
```

```
##      IID1 IID2 IID3 IID4 IID5 IID6
## IID1  0.0 1.50 1.00 0.50 0.50 2.00
## IID2  1.5 0.00 0.25 0.50 1.00 0.75
## IID3  1.0 0.25 0.00 0.33 0.67 1.00
## IID4  0.5 0.50 0.33 0.00 0.33 0.67
## IID5  0.5 1.00 0.67 0.33 0.00 0.33
## IID6  2.0 0.75 1.00 0.67 0.33 0.00
```

The matrix is nonnegative, symmetric with zeros on the diagonal because `diffs[i, j]` is the squared
gentics difference between individuals i and j, average across the SNPs. (These conditions are
are necessary but not sufficient for a matrix to be a valid distance matrix. See the note on the 
two versions of `bed2diffs` below.)

The second dataset is in sample-major mode. This means that the bed file stores the genotypes of 
the first individual at all SNPs, then the genotypes of the second individuals at all SNPs, and 
so on. In other words, if the genetic data is represented as a matrix, then each row corresponds 
to a sample. `bed2diffs` expects that the data is in SNP-major format.
```
./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode
```


```
## Set the number of OpenMP threads to 1
## 
## Compute the average genetic differences according to: 
##   Dij = (1/|Mij|) sum_{m in Mij} (z_{im} - z_{jm})^2
##   where Mij is the set of SNPs where both i and j are called
## 
## [Data::getsize] bed2diffs requires plink files [bed/bim/fam] in SNP-major mode
```

We can use `plink` can transpose the data into SNP-major mode, which is the default:
```
plink --bfile ./test/example-sample-major-mode --transpose --make-bed --out ./test/example-sample-major-mode-transposed
```


```
## PLINK v1.90b3.37 64-bit (16 May 2016)      https://www.cog-genomics.org/plink2
## (C) 2005-2016 Shaun Purcell, Christopher Chang   GNU General Public License v3
## Logging to ./test/example-sample-major-mode-transposed.log.
## Options in effect:
##   --bfile ./test/example-sample-major-mode
##   --make-bed
##   --out ./test/example-sample-major-mode-transposed
##   --transpose
## 
## Note: --transpose flag deprecated.  Use '--recode transpose ...'.
## 16384 MB RAM detected; reserving 8192 MB for main workspace.
## 6 variants loaded from .bim file.
## 6 people (4 males, 0 females, 2 ambiguous) loaded from .fam.
## Ambiguous sex IDs written to ./test/example-sample-major-mode-transposed.nosex
## .
## 1 phenotype value loaded from .fam.
## Sample-major .bed file detected.  Transposing...
## Variant-major .bed written to
## ./test/example-sample-major-mode-transposed.bed.vmaj .
## Using 1 thread (no multithreaded calculations invoked).
## Before main variant filters, 4 founders and 2 nonfounders present.
## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
## Total genotyping rate is 0.722222.
## 6 variants and 6 people pass filters and QC.
## Among remaining phenotypes, 1 is a case and 0 are controls.  (5 phenotypes are
## missing.)
## --make-bed to ./test/example-sample-major-mode-transposed.bed +
## ./test/example-sample-major-mode-transposed.bim +
## ./test/example-sample-major-mode-transposed.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
```

Finally, `bed2diffs_v1` produces the same genetic dissimilarity matrix from the transposed data.
```
./src/bed2diffs_v1 --bfile ./test/example-sample-major-mode-transposed
```



```r
diffs <- read.table("./test/example-sample-major-mode-transposed.diffs")
```


```r
diffs
```

```
##      IID1 IID2 IID3 IID4 IID5 IID6
## IID1  0.0 1.50 1.00 0.50 0.50 2.00
## IID2  1.5 0.00 0.25 0.50 1.00 0.75
## IID3  1.0 0.25 0.00 0.33 0.67 1.00
## IID4  0.5 0.50 0.33 0.00 0.33 0.67
## IID5  0.5 1.00 0.67 0.33 0.00 0.33
## IID6  2.0 0.75 1.00 0.67 0.33 0.00
```

### Output

`bed2diffs` generates two files. The matrix of average pairwise differences is written to a text file 
without row names, column names or comments, with extension `diffs`. The size of the matrix is NxN 
where N is the number of samples; there are only 0s on the main diagonal. The order of the samples 
in the `diffs` matrix is the same as in the `fam` file, but in any case, the order is explicitly 
written to a text file with one sample per line, with extension `order`.

### The two versions of bed2diffs

There are two versions, `bed2diffs_v1` and `bed2diffs_v2`. If there is missing data, the two versions 
will produce different results.

* In version 1, the average difference between samples i and j is computed across SNPs where both i 
  and j are called. This version will also output the matrix of counts (number of SNPs) where both i 
  and j are called.
* In version 2, the average difference between samples i and j is computed across all SNPs, and 
  missing genotypes at a given SNP are "imputed" as the average genotype at that SNP.
  
We used `bed2diffs_v1` to compute genetic dissimilarity matrices for the EEMS paper. However, 
this version is not guaranteed to produce a matrix that is (mathematically) an Euclidean distance
matrix, especially, if genotypes are not missing at random.

By design, there is a lot of missingness in the example dataset. For completeness, here is 
the genotype matrix in SNP-major mode, before dealing with the missing calls:


```
##      IID1 IID2 IID3 IID4 IID5 IID6
## snp1    2    0   NA   NA   NA   NA
## snp2    2   NA    1    2    2    2
## snp3    2    1    1   NA   NA    0
## snp4    1    1    2    2    2    1
## snp5   NA    2    2    2    1    1
## snp6    2    1    1   NA   NA    0
```

The `diffs` matrix (version 1) has the following eigenvalues:


```r
eigenvals <- eigen(diffs)$values
```

```
## -2.21 -1.43 -0.31 -0.23 0.1 4.08
```

So even though the matrix is symmetric, with zeros on the main diagonal and positive entries off 
the diagonal, it is not an Euclidean matrix because it has two positive eigenvalues. (An
Euclidean distance matrix is conditionally negative definite and thus it has one positive 
eigenvalue and n-1 negative eigenvalues [1].)

The alternative `bed2diffs_v2` produces an Euclidean distance matrix by construction, because it
"imputes" missing genotypes with the observed average genotype [2]. UPDATE: Fixed a bug in v2 on 
18 May 2016. 

Here is the completed genotype matrix after imputation.


```
##      IID1 IID2 IID3 IID4 IID5 IID6
## snp1  2.0  0.0    1    1    1    1
## snp2  2.0  1.8    1    2    2    2
## snp3  2.0  1.0    1    1    1    0
## snp4  1.0  1.0    2    2    2    1
## snp5  1.6  2.0    2    2    1    1
## snp6  2.0  1.0    1    1    1    0
```

Here is the dissimilarity matrix for the example dataset, computed with `bed2diffs_v2`.
```
./src/bed2diffs_v2 --bfile ./test/example-SNP-major-mode
```


```
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
```

```r
diffs <- read.table("./test/example-SNP-major-mode.diffs")
```


```r
diffs
```

```
##      IID1 IID2 IID3 IID4 IID5 IID6
## IID1 0.00 1.03 0.86 0.69 0.73 1.56
## IID2 1.03 0.00 0.44 0.34 0.51 0.67
## IID3 0.86 0.44 0.00 0.17 0.33 0.83
## IID4 0.69 0.34 0.17 0.00 0.17 0.67
## IID5 0.73 0.51 0.33 0.17 0.00 0.50
## IID6 1.56 0.67 0.83 0.67 0.50 0.00
```

The `diffs` matrix (version 2) has the following eigenvalues:


```r
eigenvals <- eigen(diffs)$values
```

```
## -1.63 -0.9 -0.55 -0.26 -0.1 3.43
```

This matrix is a valid Euclidean distance matrix as it has one positive eigenvalue and the rest 
of the eigenvalues are negative. (The eigenvalues sum up to 0, which is another property of a 
distance matrix.)

However, the imputation performed by `bed2diffs_v2` would not be appropriate if genotypes are 
not missing at random. Therefore, rather than using `bed2diffs_v2`, it would be better to clean 
the data beforehand and to remove SNPs with high missingness, as it is usually done before 
any analysis of population structure. And then use `bed2diffs_v1`.

See Documentation/bed2diffs-doc.pdf for a slightly longer explanation about the difference between 
the two versions of `bed2diffs`.

#### References
1. Balaji and Bapat. On Euclidean distance matrices. *Linear Algebra Appl*, 424:108-117, 2007.
2. By plugging the observed genotype frequency for any missing calls, we construct a representation
   $X_{n \times p}$ of $n$ individuals in $p$-dimensional space such that 
   $D_{i, j} = ||x_i - x_j||_2^2$ is the squared Euclidean distance between $i$ and $j$.
