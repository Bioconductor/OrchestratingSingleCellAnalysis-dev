---
output: html_document
bibliography: ref.bib
---

# Dealing with big data

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("aaron-collapse")) {
        event.target.classList.toggle("active");
        var content = event.target.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
    }
})
</script>

<style>
.aaron-collapse {
  background-color: #eee;
  color: #444;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

.aaron-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>

## Motivation 

Advances in scRNA-seq technologies have increased the number of cells that can be assayed in routine experiments.
Public databases such as GEO are continually expanding with more scRNA-seq studies, while large-scale projects such as the Human Cell Atlas are expected to generate data for billions of cells.
For effective data analysis, the computational methods need to scale with the increasing size of scRNA-seq data sets.
This section discusses how we can use various aspects of the Bioconductor ecosystem to tune our analysis pipelines for greater speed and efficiency.

## Fast approximations

### Nearest neighbor searching

Identification of neighbouring cells in PC or expression space is a common procedure that is used in many functions, e.g., `buildSNNGraph()`, `doubletCells()`.
The default is to favour accuracy over speed by using an exact nearest neighbour (NN) search, implemented with the $k$-means for $k$-nearest neighbours algorithm [@wang2012fast].
However, for large data sets, it may be preferable to use a faster approximate approach.
The *[BiocNeighbors](https://bioconductor.org/packages/3.12/BiocNeighbors)* framework makes it easy to switch between search options by simply changing the `BNPARAM=` argument in compatible functions.
To demonstrate, we will use the 10X PBMC data:

<button class="aaron-collapse">View history</button>
<div class="aaron-content">
   
```r
#--- loading ---#
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

#--- gene-annotation ---#
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
    rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
    column="SEQNAME", keytype="GENEID")

#--- cell-detection ---#
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

#--- quality-control ---#
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

#--- variance-modelling ---#
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

#--- dimensionality-reduction ---#
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

#--- clustering ---#
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.pbmc) <- factor(clust)
```

</div>


```r
sce.pbmc
```

```
## class: SingleCellExperiment 
## dim: 33694 3985 
## metadata(1): Samples
## assays(2): counts logcounts
## rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
## rowData names(2): ID Symbol
## colnames(3985): AAACCTGAGAAGGCCT-1 AAACCTGAGACAGACC-1 ...
##   TTTGTCAGTTAAGACA-1 TTTGTCATCCCAAGAT-1
## colData names(4): Sample Barcode sizeFactor label
## reducedDimNames(3): PCA TSNE UMAP
## altExpNames(0):
```

We had previously clustered on a shared nearest neighbor graph generated with an exact neighbour search (Section \@ref(clustering-graph)).
We repeat this below using an approximate search, implemented using the [Annoy](https://github.com/spotify/Annoy) algorithm.
This involves constructing a `AnnoyParam` object to specify the search algorithm and then passing it to the `buildSNNGraph()` function.
The results from the exact and approximate searches are consistent with most clusters from the former re-appearing in the latter.
This suggests that the inaccuracy from the approximation can be largely ignored.


```r
library(scran)
library(BiocNeighbors)
snn.gr <- buildSNNGraph(sce.pbmc, BNPARAM=AnnoyParam(), use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
table(Exact=colLabels(sce.pbmc), Approx=clusters$membership)
```

```
##      Approx
## Exact   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
##    1  205   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
##    2    0 479   0   0   2   0   0   0   0  27   0   0   0   0   0   0
##    3    0   0 540   0   1   0   0   0   0   0   0   0   0   0   0   0
##    4    0   0   0  55   0   0   0   0   0   1   0   0   0   0   0   0
##    5    0  25   0   0 349   0   0   0   0   0   0   0   0   0   0   0
##    6    0   0   0   0   0 125   0   0   0   0   0   0   0   0   0   0
##    7    0   0   0   0   0   0  46   0   0   0   0   0   0   0   0   0
##    8    0   0   0   0   0   0   0 432   0   0   0   0   0   0   0   0
##    9    0   0   0   1   0   0   0  10 291   0   0   0   0   0   0   0
##    10   0  28   0   0   0   0   0   0   0 839   0   0   0   0   0   0
##    11   0   0   0   0   0   0   0   0   0   0  47   0   0   0   0   0
##    12   0   0   0   0   0   0   0   0   0   0   0 155   0   0   0   0
##    13   0   0   0   0   0   0   0   0   0   0   0   0 166   0   0   0
##    14   0   0   0   0   0   0   0   0   0   0   0   0   0  61   0   0
##    15   0   0   0   0   0   0   0   0   0   0   0   0   0   0  84   0
##    16   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  16
```



Note that Annoy writes the NN index to disk prior to performing the search.
Thus, it may not actually be faster than the default exact algorithm for small datasets, depending on whether the overhead of disk write is offset by the computational complexity of the search.
It is also not difficult to find situations where the approximation deteriorates, especially at high dimensions, though this may not have an appreciable impact on the biological conclusions.


```r
set.seed(1000)
y1 <- matrix(rnorm(50000), nrow=1000)
y2 <- matrix(rnorm(50000), nrow=1000)
Y <- rbind(y1, y2)
exact <- findKNN(Y, k=20)
approx <- findKNN(Y, k=20, BNPARAM=AnnoyParam())
mean(exact$index!=approx$index)
```

```
## [1] 0.5619
```

### Singular value decomposition {#big-data-svd}

The singular value decomposition (SVD) underlies the PCA used throughout our analyses, e.g., in `denoisePCA()`, `fastMNN()`, `doubletCells()`.
(Briefly, the right singular vectors are the eigenvectors of the gene-gene covariance matrix, where each eigenvector represents the axis of maximum remaining variation in the PCA.)
The default `base::svd()` function performs an exact SVD that is not performant for large datasets.
Instead, we use fast approximate methods from the *[irlba](https://CRAN.R-project.org/package=irlba)* and *[rsvd](https://CRAN.R-project.org/package=rsvd)* packages, conveniently wrapped into the *[BiocSingular](https://bioconductor.org/packages/3.12/BiocSingular)* package for ease of use and package development.
Specifically, we can change the SVD algorithm used in any of these functions by simply specifying an alternative value for the `BSPARAM=` argument.


```r
library(scater)
library(BiocSingular)

# As the name suggests, it is random, so we need to set the seed.
set.seed(101000)
r.out <- runPCA(sce.pbmc, ncomponents=20, BSPARAM=RandomParam())
str(reducedDim(r.out))
```

```
##  num [1:3985, 1:20] 15.05 13.43 -8.67 -7.74 6.45 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:3985] "AAACCTGAGAAGGCCT-1" "AAACCTGAGACAGACC-1" "AAACCTGAGGCATGGT-1" "AAACCTGCAAGGTTCT-1" ...
##   ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
##  - attr(*, "varExplained")= num [1:20] 85.36 40.43 23.22 8.99 6.66 ...
##  - attr(*, "percentVar")= num [1:20] 19.85 9.4 5.4 2.09 1.55 ...
##  - attr(*, "rotation")= num [1:500, 1:20] 0.203 0.1834 0.1779 0.1063 0.0647 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:500] "LYZ" "S100A9" "S100A8" "HLA-DRA" ...
##   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
```

```r
set.seed(101001)
i.out <- runPCA(sce.pbmc, ncomponents=20, BSPARAM=IrlbaParam())
str(reducedDim(i.out))
```

```
##  num [1:3985, 1:20] 15.05 13.43 -8.67 -7.74 6.45 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:3985] "AAACCTGAGAAGGCCT-1" "AAACCTGAGACAGACC-1" "AAACCTGAGGCATGGT-1" "AAACCTGCAAGGTTCT-1" ...
##   ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
##  - attr(*, "varExplained")= num [1:20] 85.36 40.43 23.22 8.99 6.66 ...
##  - attr(*, "percentVar")= num [1:20] 19.85 9.4 5.4 2.09 1.55 ...
##  - attr(*, "rotation")= num [1:500, 1:20] 0.203 0.1834 0.1779 0.1063 0.0647 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:500] "LYZ" "S100A9" "S100A8" "HLA-DRA" ...
##   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
```

Both IRLBA and randomized SVD (RSVD) are much faster than the exact SVD with negligible loss of accuracy.
This motivates their default use in many *[scran](https://bioconductor.org/packages/3.12/scran)* and *[scater](https://bioconductor.org/packages/3.12/scater)* functions, at the cost of requiring users to set the seed to guarantee reproducibility.
IRLBA can occasionally fail to converge and require more iterations (passed via `maxit=` in `IrlbaParam()`), while RSVD involves an explicit trade-off between accuracy and speed based on its oversampling parameter (`p=`) and number of power iterations (`q=`).
We tend to prefer IRLBA as its default behavior is more accurate, though RSVD is much faster for file-backed matrices (Section \@ref(data-integration)).

## Parallelization

Parallelization of calculations across genes or cells is an obvious strategy for speeding up scRNA-seq analysis workflows.
The *[BiocParallel](https://bioconductor.org/packages/3.12/BiocParallel)* package provides a common interface for parallel computing throughout the Bioconductor ecosystem, manifesting as a `BPPARAM=` argument in compatible functions.
We can pick from a diverse range of parallelization backends depending on the available hardware and operating system.
For example, we might use forking across 2 cores to parallelize the variance calculations on a Unix system:


```r
library(BiocParallel)
dec.pbmc.mc <- modelGeneVar(sce.pbmc, BPPARAM=MulticoreParam(2))
dec.pbmc.mc
```

```
## DataFrame with 33694 rows and 6 columns
##                     mean       total        tech          bio   p.value
##                <numeric>   <numeric>   <numeric>    <numeric> <numeric>
## RP11-34P13.3 0.000000000 0.000000000 0.000000000  0.00000e+00       NaN
## FAM138A      0.000000000 0.000000000 0.000000000  0.00000e+00       NaN
## OR4F5        0.000000000 0.000000000 0.000000000  0.00000e+00       NaN
## RP11-34P13.7 0.002166050 0.002227438 0.002227466 -2.81875e-08  0.500035
## RP11-34P13.8 0.000522431 0.000549601 0.000537244  1.23571e-05  0.436506
## ...                  ...         ...         ...          ...       ...
## AC233755.2     0.0000000   0.0000000   0.0000000   0.00000000       NaN
## AC233755.1     0.0000000   0.0000000   0.0000000   0.00000000       NaN
## AC240274.1     0.0102893   0.0121099   0.0105809   0.00152898  0.157651
## AC213203.1     0.0000000   0.0000000   0.0000000   0.00000000       NaN
## FAM231B        0.0000000   0.0000000   0.0000000   0.00000000       NaN
##                    FDR
##              <numeric>
## RP11-34P13.3       NaN
## FAM138A            NaN
## OR4F5              NaN
## RP11-34P13.7  0.756451
## RP11-34P13.8  0.756451
## ...                ...
## AC233755.2         NaN
## AC233755.1         NaN
## AC240274.1    0.756451
## AC213203.1         NaN
## FAM231B            NaN
```

Another approach would be to distribute jobs across a network of computers, which yields the same result:


```r
dec.pbmc.snow <- modelGeneVar(sce.pbmc, BPPARAM=SnowParam(5))
dec.pbmc.snow
```

```
## DataFrame with 33694 rows and 6 columns
##                     mean       total        tech          bio   p.value
##                <numeric>   <numeric>   <numeric>    <numeric> <numeric>
## RP11-34P13.3 0.000000000 0.000000000 0.000000000  0.00000e+00       NaN
## FAM138A      0.000000000 0.000000000 0.000000000  0.00000e+00       NaN
## OR4F5        0.000000000 0.000000000 0.000000000  0.00000e+00       NaN
## RP11-34P13.7 0.002166050 0.002227438 0.002227466 -2.81875e-08  0.500035
## RP11-34P13.8 0.000522431 0.000549601 0.000537244  1.23571e-05  0.436506
## ...                  ...         ...         ...          ...       ...
## AC233755.2     0.0000000   0.0000000   0.0000000   0.00000000       NaN
## AC233755.1     0.0000000   0.0000000   0.0000000   0.00000000       NaN
## AC240274.1     0.0102893   0.0121099   0.0105809   0.00152898  0.157651
## AC213203.1     0.0000000   0.0000000   0.0000000   0.00000000       NaN
## FAM231B        0.0000000   0.0000000   0.0000000   0.00000000       NaN
##                    FDR
##              <numeric>
## RP11-34P13.3       NaN
## FAM138A            NaN
## OR4F5              NaN
## RP11-34P13.7  0.756451
## RP11-34P13.8  0.756451
## ...                ...
## AC233755.2         NaN
## AC233755.1         NaN
## AC240274.1    0.756451
## AC213203.1         NaN
## FAM231B            NaN
```



For high-performance computing (HPC) systems with a cluster of compute nodes, we can distribute jobs via the job scheduler using the `BatchtoolsParam` class.
The example below assumes a SLURM cluster, though the settings can be easily configured for a particular system (see [here](https://bioconductor.org/packages/3.12/BiocParallel/vignettes/BiocParallel_BatchtoolsParam.pdf) for details).


```r
# 2 hours, 8 GB, 1 CPU per task, for 10 tasks.
bpp <- BatchtoolsParam(10, cluster="slurm",
	resources=list(walltime=7200, memory=8000, ncpus=1))
```

Parallelization is best suited for CPU-intensive calculations where the division of labor results in a concomitant reduction in compute time.
It is not suited for tasks that are bounded by other compute resources, e.g., memory or file I/O (though the latter is less of an issue on HPC systems with parallel read/write).
In particular, R itself is inherently single-core, so many of the parallelization backends involve (i) setting up one or more separate R sessions, (ii) loading the relevant packages and (iii) transmitting the data to that session.
Depending on the nature and size of the task, this overhead may outweigh any benefit from parallel computing. 

## Out of memory representations

The count matrix is the central structure around which our analyses are based.
In most of the previous chapters, this has been held fully in memory as a dense `matrix` or as a sparse `dgCMatrix`.
Howevever, in-memory representations may not be feasible for very large data sets, especially on machines with limited memory.
For example, the 1.3 million brain cell data set from 10X Genomics [@zheng2017massively] would require over 100 GB of RAM to hold as a `matrix` and around 30 GB as a `dgCMatrix`.
This makes it challenging to explore the data on anything less than a HPC system.

The obvious solution is to use a file-backed matrix representation where the data are held on disk and subsets are retrieved into memory as requested.
While a number of implementations of file-backed matrices are available (e.g., *[bigmemory](https://CRAN.R-project.org/package=bigmemory)*, *[matter](https://bioconductor.org/packages/3.12/matter)*), we will be using the implementation from the *[HDF5Array](https://bioconductor.org/packages/3.12/HDF5Array)* package.
This uses the popular HDF5 format as the underlying data store, which provides a measure of standardization and portability across systems.
We demonstrate with a subset of 20,000 cells from the 1.3 million brain cell data set, as provided by the *[TENxBrainData](https://bioconductor.org/packages/3.12/TENxBrainData)* package.


```r
library(TENxBrainData)
sce.brain <- TENxBrainData20k() 
sce.brain
```

```
## class: SingleCellExperiment 
## dim: 27998 20000 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(2): Ensembl Symbol
## colnames: NULL
## colData names(4): Barcode Sequence Library Mouse
## reducedDimNames(0):
## altExpNames(0):
```

Examination of the `SingleCellExperiment` object indicates that the count matrix is a `HDF5Matrix`.
From a comparison of the memory usage, it is clear that this matrix object is simply a stub that points to the much larger HDF5 file that actually contains the data.
This avoids the need for large RAM availability during analyses.


```r
counts(sce.brain)
```

```
## <27998 x 20000> matrix of class HDF5Matrix and type "integer":
##              [,1]     [,2]     [,3]     [,4] ... [,19997] [,19998] [,19999]
##     [1,]        0        0        0        0   .        0        0        0
##     [2,]        0        0        0        0   .        0        0        0
##     [3,]        0        0        0        0   .        0        0        0
##     [4,]        0        0        0        0   .        0        0        0
##     [5,]        0        0        0        0   .        0        0        0
##      ...        .        .        .        .   .        .        .        .
## [27994,]        0        0        0        0   .        0        0        0
## [27995,]        0        0        0        1   .        0        2        0
## [27996,]        0        0        0        0   .        0        1        0
## [27997,]        0        0        0        0   .        0        0        0
## [27998,]        0        0        0        0   .        0        0        0
##          [,20000]
##     [1,]        0
##     [2,]        0
##     [3,]        0
##     [4,]        0
##     [5,]        0
##      ...        .
## [27994,]        0
## [27995,]        0
## [27996,]        0
## [27997,]        0
## [27998,]        0
```

```r
object.size(counts(sce.brain))
```

```
## 2328 bytes
```

```r
file.info(path(counts(sce.brain)))$size
```

```
## [1] 76264332
```

Manipulation of the count matrix will generally result in the creation of a `DelayedArray` object from the *[DelayedArray](https://bioconductor.org/packages/3.12/DelayedArray)* package.
This remembers the operations to be applied to the counts and stores them in the object, to be executed when the modified matrix values are realized for use in calculations.
The use of delayed operations avoids the need to write the modified values to a new file at every operation, which would unnecessarily require time-consuming disk I/O.


```r
tmp <- counts(sce.brain)
tmp <- log2(tmp + 1)
tmp
```

```
## <27998 x 20000> matrix of class DelayedMatrix and type "double":
##              [,1]     [,2]     [,3] ... [,19999] [,20000]
##     [1,]        0        0        0   .        0        0
##     [2,]        0        0        0   .        0        0
##     [3,]        0        0        0   .        0        0
##     [4,]        0        0        0   .        0        0
##     [5,]        0        0        0   .        0        0
##      ...        .        .        .   .        .        .
## [27994,]        0        0        0   .        0        0
## [27995,]        0        0        0   .        0        0
## [27996,]        0        0        0   .        0        0
## [27997,]        0        0        0   .        0        0
## [27998,]        0        0        0   .        0        0
```

Many functions described in the previous workflows are capable of accepting `HDF5Matrix` objects.
This is powered by the availability of common methods for all matrix representations (e.g., subsetting, combining, methods from *[DelayedMatrixStats](https://bioconductor.org/packages/3.12/DelayedMatrixStats)*) as well as representation-agnostic C++ code using *[beachmat](https://bioconductor.org/packages/3.12/beachmat)* [@lun2018beachmat].
For example, we compute QC metrics below with the same `calculateQCMetrics()` function that we used in the other workflows.


```r
library(scater)
is.mito <- grepl("^mt-", rowData(sce.brain)$Symbol)
qcstats <- perCellQCMetrics(sce.brain, subsets=list(Mt=is.mito))
qcstats
```

```
## DataFrame with 20000 rows and 6 columns
##             sum  detected subsets_Mt_sum subsets_Mt_detected subsets_Mt_percent
##       <integer> <integer>      <integer>           <integer>          <numeric>
## 1          3060      1546            123                  10            4.01961
## 2          3500      1694            118                  11            3.37143
## 3          3092      1613             58                   9            1.87581
## 4          4420      2050            131                  10            2.96380
## 5          3771      1813            100                   8            2.65182
## ...         ...       ...            ...                 ...                ...
## 19996      4431      2050            127                   9           2.866170
## 19997      6988      2704             60                   9           0.858615
## 19998      8749      2988            305                  11           3.486113
## 19999      3842      1711            129                   8           3.357626
## 20000      1775       945             26                   6           1.464789
##           total
##       <integer>
## 1          3060
## 2          3500
## 3          3092
## 4          4420
## 5          3771
## ...         ...
## 19996      4431
## 19997      6988
## 19998      8749
## 19999      3842
## 20000      1775
```

Needless to say, data access from file-backed representations is slower than that from in-memory representations.
The time spent retrieving data from disk is an unavoidable cost of reducing memory usage.
Whether this is tolerable depends on the application.
One example usage pattern involves performing the heavy computing quickly with in-memory representations on HPC systems with plentiful memory, and then distributing file-backed counterparts to individual users for exploration and visualization on their personal machines.

## Session Info {-}

<button class="aaron-collapse">View session info</button>
<div class="aaron-content">
```
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /home/biocbuild/bbs-3.12-bioc/R/lib/libRblas.so
LAPACK: /home/biocbuild/bbs-3.12-bioc/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=C              
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] TENxBrainData_1.9.0         HDF5Array_1.17.3           
 [3] rhdf5_2.33.7                BiocParallel_1.23.2        
 [5] BiocSingular_1.5.0          scater_1.17.4              
 [7] ggplot2_3.3.2               bluster_0.99.1             
 [9] BiocNeighbors_1.7.0         scran_1.17.15              
[11] SingleCellExperiment_1.11.6 SummarizedExperiment_1.19.6
[13] DelayedArray_0.15.7         matrixStats_0.56.0         
[15] Matrix_1.2-18               Biobase_2.49.0             
[17] GenomicRanges_1.41.6        GenomeInfoDb_1.25.10       
[19] IRanges_2.23.10             S4Vectors_0.27.12          
[21] BiocGenerics_0.35.4         BiocStyle_2.17.0           
[23] simpleSingleCell_1.13.16   

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0              colorspace_1.4-1             
 [3] ellipsis_0.3.1                scuttle_0.99.12              
 [5] XVector_0.29.3                bit64_4.0.2                  
 [7] interactiveDisplayBase_1.27.5 AnnotationDbi_1.51.3         
 [9] codetools_0.2-16              knitr_1.29                   
[11] dbplyr_1.4.4                  graph_1.67.1                 
[13] shiny_1.5.0                   BiocManager_1.30.10          
[15] compiler_4.0.2                httr_1.4.2                   
[17] dqrng_0.2.1                   assertthat_0.2.1             
[19] fastmap_1.0.1                 limma_3.45.10                
[21] later_1.1.0.1                 htmltools_0.5.0              
[23] tools_4.0.2                   rsvd_1.0.3                   
[25] igraph_1.2.5                  gtable_0.3.0                 
[27] glue_1.4.1                    GenomeInfoDbData_1.2.3       
[29] dplyr_1.0.1                   rappdirs_0.3.1               
[31] Rcpp_1.0.5                    vctrs_0.3.2                  
[33] rhdf5filters_1.1.2            ExperimentHub_1.15.1         
[35] DelayedMatrixStats_1.11.1     xfun_0.16                    
[37] stringr_1.4.0                 ps_1.3.4                     
[39] beachmat_2.5.1                mime_0.9                     
[41] lifecycle_0.2.0               irlba_2.3.3                  
[43] statmod_1.4.34                XML_3.99-0.5                 
[45] AnnotationHub_2.21.2          edgeR_3.31.4                 
[47] zlibbioc_1.35.0               scales_1.1.1                 
[49] promises_1.1.1                yaml_2.2.1                   
[51] curl_4.3                      memoise_1.1.0                
[53] gridExtra_2.3                 stringi_1.4.6                
[55] RSQLite_2.2.0                 BiocVersion_3.12.0           
[57] rlang_0.4.7                   pkgconfig_2.0.3              
[59] bitops_1.0-6                  evaluate_0.14                
[61] lattice_0.20-41               purrr_0.3.4                  
[63] Rhdf5lib_1.11.3               CodeDepends_0.6.5            
[65] bit_4.0.4                     processx_3.4.3               
[67] tidyselect_1.1.0              magrittr_1.5                 
[69] bookdown_0.20                 R6_2.4.1                     
[71] snow_0.4-3                    generics_0.0.2               
[73] DBI_1.1.0                     pillar_1.4.6                 
[75] withr_2.2.0                   RCurl_1.98-1.2               
[77] tibble_3.0.3                  crayon_1.3.4                 
[79] BiocFileCache_1.13.1          rmarkdown_2.3                
[81] viridis_0.5.1                 locfit_1.5-9.4               
[83] grid_4.0.2                    blob_1.2.1                   
[85] callr_3.4.3                   digest_0.6.25                
[87] xtable_1.8-4                  httpuv_1.5.4                 
[89] munsell_0.5.0                 beeswarm_0.2.3               
[91] viridisLite_0.3.0             vipor_0.4.5                  
```
</div>
