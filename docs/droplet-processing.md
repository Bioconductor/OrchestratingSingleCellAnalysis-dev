---
output:
  html_document
bibliography: ref.bib
---

# Droplet processing

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("rebook-collapse")) {
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
.rebook-collapse {
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

.rebook-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>

## Motivation

Droplet-based single-cell protocols aim to isolate each cell inside its own droplet in a water-in-oil emulsion, such that each droplet serves as a miniature reaction chamber for highly multiplexed library preparation [@macosko2015highly;@klein2015droplet].
Upon sequencing, reads are assigned to individual cells based on the presence of droplet-specific barcodes.
This enables a massive increase in the number of cells that can be processed in typical scRNA-seq experiments, contributing to the dominance^[As of time of writing.] of technologies such as the 10X Genomics platform [@zheng2017massively].
However, as the allocation of cells to droplets is not known in advance, the data analysis requires some special steps to determine what each droplet actually contains.
This chapter explores some of the more common preprocessing procedures that might be applied to the count matrices generated from droplet protocols.

## Calling cells from empty droplets {#qc-droplets}

### Background

An unique aspect of droplet-based data is that we have no prior knowledge about whether a particular library (i.e., cell barcode) corresponds to cell-containing or empty droplets.
Thus, we need to call cells from empty droplets based on the observed expression profiles.
This is not entirely straightforward as empty droplets can contain ambient (i.e., extracellular) RNA that can be captured and sequenced, resulting in non-zero counts for libraries that do not contain any cell.
To demonstrate, we obtain the **unfiltered** count matrix for the PBMC dataset from 10X Genomics.

<button class="rebook-collapse">View history</button>
<div class="rebook-content">
   
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
```

</div>


```r
sce.pbmc
```

```
## class: SingleCellExperiment 
## dim: 33694 737280 
## metadata(1): Samples
## assays(1): counts
## rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
##   ENSG00000268674
## rowData names(2): ID Symbol
## colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ...
##   TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
## colData names(2): Sample Barcode
## reducedDimNames(0):
## altExpNames(0):
```

The distribution of total counts exhibits a sharp transition between barcodes with large and small total counts (Figure \@ref(fig:rankplot)), probably corresponding to cell-containing and empty droplets respectively.
A simple approach would be to apply a threshold on the total count to only retain those barcodes with large totals.
However, this unnecessarily discards libraries derived from cell types with low RNA content.


```r
library(DropletUtils)
bcrank <- barcodeRanks(counts(sce.pbmc))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
        col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
```

<div class="figure">
<img src="droplet-processing_files/figure-html/rankplot-1.png" alt="Total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts). The inferred locations of the inflection and knee points are also shown." width="672" />
<p class="caption">(\#fig:rankplot)Total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts). The inferred locations of the inflection and knee points are also shown.</p>
</div>

### Testing for empty droplets

We use the `emptyDrops()` function to test whether the expression profile for each cell barcode is significantly different from the ambient RNA pool [@lun2018distinguishing].
Any significant deviation indicates that the barcode corresponds to a cell-containing droplet.
This allows us to discriminate between well-sequenced empty droplets and droplets derived from cells with little RNA, both of which would have similar total counts in Figure \@ref(fig:rankplot).
We call cells at a false discovery rate (FDR) of 0.1%, meaning that no more than 0.1% of our called barcodes should be empty droplets on average.


```r
# emptyDrops performs Monte Carlo simulations to compute p-values,
# so we need to set the seed to obtain reproducible results.
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

# See ?emptyDrops for an explanation of why there are NA values.
summary(e.out$FDR <= 0.001)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical     989    4300  731991
```

`emptyDrops()` uses Monte Carlo simulations to compute $p$-values for the multinomial sampling transcripts from the ambient pool.
The number of Monte Carlo iterations determines the lower bound for the $p$-values [@phipson2010permutation].
The `Limited` field in the output indicates whether or not the computed $p$-value for a particular barcode is bounded by the number of iterations.
If any non-significant barcodes are `TRUE` for `Limited`, we may need to increase the number of iterations.
A larger number of iterations will result in a lower $p$-value for these barcodes, which may allow them to be detected after correcting for multiple testing.


```r
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)
```

```
##        Limited
## Sig     FALSE TRUE
##   FALSE   989    0
##   TRUE   1728 2572
```

As mentioned above, `emptyDrops()` assumes that barcodes with low total UMI counts are empty droplets.
Thus, the null hypothesis should be true for all of these barcodes. 
We can check whether the hypothesis testing procedure holds its size by examining the distribution of $p$-values for low-total barcodes with `test.ambient=TRUE`.
Ideally, the distribution should be close to uniform (Figure \@ref(fig:ambientpvalhist)).
Large peaks near zero indicate that barcodes with total counts below `lower` are not all ambient in origin.
This can be resolved by decreasing `lower` further to ensure that barcodes corresponding to droplets with very small cells are not used to estimate the ambient profile.


```r
set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(sce.pbmc), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
    xlab="P-value", main="", col="grey80") 
```

<div class="figure">
<img src="droplet-processing_files/figure-html/ambientpvalhist-1.png" alt="Distribution of $p$-values for the assumed empty droplets." width="672" />
<p class="caption">(\#fig:ambientpvalhist)Distribution of $p$-values for the assumed empty droplets.</p>
</div>

Once we are satisfied with the performance of `emptyDrops()`, we subset our `SingleCellExperiment` object to retain only the detected cells.
Discerning readers will notice the use of `which()`, which conveniently removes the `NA`s prior to the subsetting.
 

```r
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]
```

It is worth pointing out that, at this point, we do not attempt to remove the ambient contamination from each library.
Accurate quantification of the contamination rate in each cell is difficult as it generally requires some prior biological knowledge about genes that are expected to have mutually exclusive expression profiles _and_ are highly abundant in the ambient solution [@young2018soupx].
Fortunately, ambient contamination usually has little effect on the downstream conclusions for routine analyses; cell type identities are usually easy enough to determine from the affected genes, notwithstanding a (mostly harmless) low background level of expression for marker genes that should be unique to a cell type.
However, more susceptible analyses may require specific remedies like those discussed in Section \@ref(ambient-problems).

### Relationship with other QC metrics

While `emptyDrops()` will distinguish cells from empty droplets, it makes no statement about the quality of the cells.
It is entirely possible for droplets to contain damaged or dying cells, which need to be removed prior to downstream analysis.
This is achieved using the same outlier-based strategy described in Section \@ref(quality-control-outlier).
Filtering on the mitochondrial proportion provides the most additional benefit in this situation, provided that we check that we are not removing a subpopulation of metabolically active cells (Figure \@ref(fig:qc-mito-pbmc)). 


```r
library(scuttle)
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
pbmc.qc <- perCellQCMetrics(sce.pbmc, subsets=list(MT=is.mito))
discard.mito <- isOutlier(pbmc.qc$subsets_MT_percent, type="higher")
summary(discard.mito)
```

```
##    Mode   FALSE    TRUE 
## logical    3985     315
```

```r
plot(pbmc.qc$sum, pbmc.qc$subsets_MT_percent, log="x",
    xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(discard.mito, "thresholds")["higher"], col="red")
```

<div class="figure">
<img src="droplet-processing_files/figure-html/qc-mito-pbmc-1.png" alt="Percentage of reads assigned to mitochondrial transcripts, plotted against the library size. The red line represents the upper threshold used for QC filtering." width="672" />
<p class="caption">(\#fig:qc-mito-pbmc)Percentage of reads assigned to mitochondrial transcripts, plotted against the library size. The red line represents the upper threshold used for QC filtering.</p>
</div>

`emptyDrops()` already removes cells with very low library sizes or (by association) low numbers of expressed genes.
Thus, further filtering on these metrics is not strictly necessary.
It may still be desirable to filter on both of these metrics to remove non-empty droplets containing cell fragments or stripped nuclei that were not caught by the mitochondrial filter.
However, this should be weighed against the risk of losing genuine cell types as discussed in Section \@ref(outlier-assumptions).

Note that _CellRanger_ version 3 automatically performs cell calling using an algorithm similar to `emptyDrops()`.
If we had started our analysis with the **filtered** count matrix, we could go straight to computing other QC metrics.
We would not need to run `emptyDrops()` manually as shown here, and indeed, attempting to do so would lead to nonsensical results if not outright software errors.
Nonetheless, it may still be desirable to load the **unfiltered** matrix and apply `emptyDrops()` ourselves, on occasions where more detailed inspection or control of the cell-calling statistics is desired.

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.0.0 Patched (2020-05-01 r78341)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /home/luna/Software/R/R-4-0-branch-dev/lib/libRblas.so
LAPACK: /home/luna/Software/R/R-4-0-branch-dev/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] scuttle_0.99.13             DropletUtils_1.9.10        
 [3] SingleCellExperiment_1.11.6 SummarizedExperiment_1.19.6
 [5] DelayedArray_0.15.7         matrixStats_0.56.0         
 [7] Matrix_1.2-18               Biobase_2.49.0             
 [9] GenomicRanges_1.41.6        GenomeInfoDb_1.25.10       
[11] IRanges_2.23.10             S4Vectors_0.27.12          
[13] BiocGenerics_0.35.4         BiocStyle_2.17.0           
[15] rebook_0.99.4              

loaded via a namespace (and not attached):
 [1] locfit_1.5-9.4            xfun_0.16                
 [3] HDF5Array_1.17.3          lattice_0.20-41          
 [5] rhdf5_2.33.7              htmltools_0.5.0          
 [7] yaml_2.2.1                XML_3.99-0.5             
 [9] rlang_0.4.7               R.oo_1.23.0              
[11] R.utils_2.9.2             BiocParallel_1.23.2      
[13] CodeDepends_0.6.5         dqrng_0.2.1              
[15] GenomeInfoDbData_1.2.3    stringr_1.4.0            
[17] zlibbioc_1.35.0           R.methodsS3_1.8.0        
[19] codetools_0.2-16          evaluate_0.14            
[21] knitr_1.29                callr_3.4.3              
[23] ps_1.3.4                  highr_0.8                
[25] Rcpp_1.0.5                edgeR_3.31.4             
[27] BiocManager_1.30.10       limma_3.45.10            
[29] graph_1.67.1              XVector_0.29.3           
[31] digest_0.6.25             stringi_1.4.6            
[33] bookdown_0.20             processx_3.4.3           
[35] grid_4.0.0                tools_4.0.0              
[37] bitops_1.0-6              rhdf5filters_1.1.2       
[39] magrittr_1.5              RCurl_1.98-1.2           
[41] DelayedMatrixStats_1.11.1 rmarkdown_2.3            
[43] Rhdf5lib_1.11.3           compiler_4.0.0           
```
</div>
