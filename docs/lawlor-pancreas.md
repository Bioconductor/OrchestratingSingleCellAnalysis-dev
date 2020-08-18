# Lawlor human pancreas (SMARTer)

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

## Introduction

This performs an analysis of the @lawlor2017singlecell dataset,
consisting of human pancreas cells from various donors.

## Data loading


```r
library(scRNAseq)
sce.lawlor <- LawlorPancreasData()
```


```r
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(edb, keys=rownames(sce.lawlor), keytype="GENEID", 
    columns=c("SYMBOL", "SEQNAME"))
rowData(sce.lawlor) <- anno[match(rownames(sce.lawlor), anno[,1]),-1]
```

## Quality control


```r
unfiltered <- sce.lawlor
```


```r
library(scater)
stats <- perCellQCMetrics(sce.lawlor, 
    subsets=list(Mito=which(rowData(sce.lawlor)$SEQNAME=="MT")))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
    batch=sce.lawlor$`islet unos id`)
sce.lawlor <- sce.lawlor[,!qc$discard]
```


```r
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, x="islet unos id", y="sum", colour_by="discard") +
        scale_y_log10() + ggtitle("Total count") +
        theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x="islet unos id", y="detected", 
        colour_by="discard") + scale_y_log10() + ggtitle("Detected features") +
        theme(axis.text.x = element_text(angle = 90)), 
    plotColData(unfiltered, x="islet unos id", y="subsets_Mito_percent",
        colour_by="discard") + ggtitle("Mito percent") +
        theme(axis.text.x = element_text(angle = 90)),
    ncol=2
)
```

<div class="figure">
<img src="lawlor-pancreas_files/figure-html/unref-lawlor-qc-dist-1.png" alt="Distribution of each QC metric across cells from each donor of the Lawlor pancreas dataset. Each point represents a cell and is colored according to whether that cell was discarded." width="960" />
<p class="caption">(\#fig:unref-lawlor-qc-dist)Distribution of each QC metric across cells from each donor of the Lawlor pancreas dataset. Each point represents a cell and is colored according to whether that cell was discarded.</p>
</div>


```r
plotColData(unfiltered, x="sum", y="subsets_Mito_percent",
    colour_by="discard") + scale_x_log10()
```

<div class="figure">
<img src="lawlor-pancreas_files/figure-html/unref-lawlor-qc-comp-1.png" alt="Percentage of mitochondrial reads in each cell in the 416B dataset compared to the total count. Each point represents a cell and is colored according to whether that cell was discarded." width="672" />
<p class="caption">(\#fig:unref-lawlor-qc-comp)Percentage of mitochondrial reads in each cell in the 416B dataset compared to the total count. Each point represents a cell and is colored according to whether that cell was discarded.</p>
</div>


```r
colSums(as.matrix(qc))
```

```
##              low_lib_size            low_n_features high_subsets_Mito_percent 
##                         9                         5                        25 
##                   discard 
##                        34
```

## Normalization


```r
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.lawlor)
sce.lawlor <- computeSumFactors(sce.lawlor, clusters=clusters)
sce.lawlor <- logNormCounts(sce.lawlor)
```


```r
summary(sizeFactors(sce.lawlor))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.295   0.781   0.963   1.000   1.182   2.629
```


```r
plot(librarySizeFactors(sce.lawlor), sizeFactors(sce.lawlor), pch=16,
    xlab="Library size factors", ylab="Deconvolution factors", log="xy")
```

<div class="figure">
<img src="lawlor-pancreas_files/figure-html/unref-lawlor-norm-1.png" alt="Relationship between the library size factors and the deconvolution size factors in the Lawlor pancreas dataset." width="672" />
<p class="caption">(\#fig:unref-lawlor-norm)Relationship between the library size factors and the deconvolution size factors in the Lawlor pancreas dataset.</p>
</div>

## Variance modelling

Using age as a proxy for the donor.


```r
dec.lawlor <- modelGeneVar(sce.lawlor, block=sce.lawlor$`islet unos id`)
chosen.genes <- getTopHVGs(dec.lawlor, n=2000)
```


```r
par(mfrow=c(4,2))
blocked.stats <- dec.lawlor$per.block
for (i in colnames(blocked.stats)) {
    current <- blocked.stats[[i]]
    plot(current$mean, current$total, main=i, pch=16, cex=0.5,
        xlab="Mean of log-expression", ylab="Variance of log-expression")
    curfit <- metadata(current)
    curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}
```

<div class="figure">
<img src="lawlor-pancreas_files/figure-html/unnamed-chunk-4-1.png" alt="Per-gene variance as a function of the mean for the log-expression values in the Lawlor pancreas dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted separately for each donor." width="672" />
<p class="caption">(\#fig:unnamed-chunk-4)Per-gene variance as a function of the mean for the log-expression values in the Lawlor pancreas dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted separately for each donor.</p>
</div>

## Dimensionality reduction


```r
library(BiocSingular)
set.seed(101011001)
sce.lawlor <- runPCA(sce.lawlor, subset_row=chosen.genes, ncomponents=25)
sce.lawlor <- runTSNE(sce.lawlor, dimred="PCA")
```

## Clustering


```r
snn.gr <- buildSNNGraph(sce.lawlor, use.dimred="PCA")
colLabels(sce.lawlor) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
```


```r
table(colLabels(sce.lawlor), sce.lawlor$`cell type`)
```

```
##    
##     Acinar Alpha Beta Delta Ductal Gamma/PP None/Other Stellate
##   1      1     0    0    13      2       16          2        0
##   2      0     1   76     1      0        0          0        0
##   3      0   161    1     0      0        1          2        0
##   4      0     1    0     1      0        0          5       19
##   5      0     0  175     4      1        0          1        0
##   6     22     0    0     0      0        0          0        0
##   7      0    75    0     0      0        0          0        0
##   8      0     0    0     1     20        0          2        0
```


```r
table(colLabels(sce.lawlor), sce.lawlor$`islet unos id`)
```

```
##    
##     ACCG268 ACCR015A ACEK420A ACEL337 ACHY057 ACIB065 ACIW009 ACJV399
##   1       8        2        2       4       4       4       9       1
##   2      14        3        2      33       3       2       4      17
##   3      36       23       14      13      14      14      21      30
##   4       7        1        0       1       0       4       9       4
##   5      34       10        4      39       7      23      24      40
##   6       0        2       13       0       0       0       5       2
##   7      32       12        0       5       6       7       4       9
##   8       1        1        2       1       2       1      12       3
```


```r
gridExtra::grid.arrange(
    plotTSNE(sce.lawlor, colour_by="label"),
    plotTSNE(sce.lawlor, colour_by="islet unos id"),
    ncol=2
)
```

<div class="figure">
<img src="lawlor-pancreas_files/figure-html/unref-grun-tsne-1.png" alt="Obligatory $t$-SNE plots of the Lawlor pancreas dataset. Each point represents a cell that is colored by cluster (left) or batch (right)." width="672" />
<p class="caption">(\#fig:unref-grun-tsne)Obligatory $t$-SNE plots of the Lawlor pancreas dataset. Each point represents a cell that is colored by cluster (left) or batch (right).</p>
</div>

## Session Info {-}

<button class="aaron-collapse">View session info</button>
<div class="aaron-content">
```
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

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
 [1] BiocSingular_1.5.0          scran_1.17.3               
 [3] scater_1.17.2               ggplot2_3.3.2              
 [5] ensembldb_2.13.1            AnnotationFilter_1.13.0    
 [7] GenomicFeatures_1.41.0      AnnotationDbi_1.51.1       
 [9] AnnotationHub_2.21.1        BiocFileCache_1.13.0       
[11] dbplyr_1.4.4                scRNAseq_2.3.8             
[13] SingleCellExperiment_1.11.6 SummarizedExperiment_1.19.5
[15] DelayedArray_0.15.6         matrixStats_0.56.0         
[17] Matrix_1.2-18               Biobase_2.49.0             
[19] GenomicRanges_1.41.5        GenomeInfoDb_1.25.5        
[21] IRanges_2.23.10             S4Vectors_0.27.12          
[23] BiocGenerics_0.35.4         BiocStyle_2.17.0           
[25] simpleSingleCell_1.13.5    

loaded via a namespace (and not attached):
  [1] Rtsne_0.15                    ggbeeswarm_0.6.0             
  [3] colorspace_1.4-1              ellipsis_0.3.1               
  [5] scuttle_0.99.10               XVector_0.29.3               
  [7] BiocNeighbors_1.7.0           farver_2.0.3                 
  [9] bit64_0.9-7                   interactiveDisplayBase_1.27.5
 [11] codetools_0.2-16              knitr_1.29                   
 [13] Rsamtools_2.5.3               graph_1.67.1                 
 [15] shiny_1.5.0                   BiocManager_1.30.10          
 [17] compiler_4.0.2                httr_1.4.1                   
 [19] dqrng_0.2.1                   assertthat_0.2.1             
 [21] fastmap_1.0.1                 lazyeval_0.2.2               
 [23] limma_3.45.7                  later_1.1.0.1                
 [25] htmltools_0.5.0               prettyunits_1.1.1            
 [27] tools_4.0.2                   igraph_1.2.5                 
 [29] rsvd_1.0.3                    gtable_0.3.0                 
 [31] glue_1.4.1                    GenomeInfoDbData_1.2.3       
 [33] dplyr_1.0.0                   rappdirs_0.3.1               
 [35] Rcpp_1.0.4.6                  vctrs_0.3.1                  
 [37] Biostrings_2.57.2             ExperimentHub_1.15.0         
 [39] rtracklayer_1.49.3            DelayedMatrixStats_1.11.1    
 [41] xfun_0.15                     stringr_1.4.0                
 [43] ps_1.3.3                      mime_0.9                     
 [45] lifecycle_0.2.0               irlba_2.3.3                  
 [47] statmod_1.4.34                XML_3.99-0.3                 
 [49] edgeR_3.31.4                  zlibbioc_1.35.0              
 [51] scales_1.1.1                  hms_0.5.3                    
 [53] promises_1.1.1                ProtGenerics_1.21.0          
 [55] yaml_2.2.1                    curl_4.3                     
 [57] memoise_1.1.0                 gridExtra_2.3                
 [59] biomaRt_2.45.1                stringi_1.4.6                
 [61] RSQLite_2.2.0                 highr_0.8                    
 [63] BiocVersion_3.12.0            BiocParallel_1.23.0          
 [65] rlang_0.4.6                   pkgconfig_2.0.3              
 [67] bitops_1.0-6                  evaluate_0.14                
 [69] lattice_0.20-41               purrr_0.3.4                  
 [71] labeling_0.3                  GenomicAlignments_1.25.3     
 [73] CodeDepends_0.6.5             cowplot_1.0.0                
 [75] bit_1.1-15.2                  processx_3.4.2               
 [77] tidyselect_1.1.0              magrittr_1.5                 
 [79] bookdown_0.20                 R6_2.4.1                     
 [81] generics_0.0.2                DBI_1.1.0                    
 [83] pillar_1.4.4                  withr_2.2.0                  
 [85] RCurl_1.98-1.2                tibble_3.0.1                 
 [87] crayon_1.3.4                  rmarkdown_2.3                
 [89] viridis_0.5.1                 progress_1.2.2               
 [91] locfit_1.5-9.4                grid_4.0.2                   
 [93] blob_1.2.1                    callr_3.4.3                  
 [95] digest_0.6.25                 xtable_1.8-4                 
 [97] httpuv_1.5.4                  openssl_1.4.2                
 [99] munsell_0.5.0                 beeswarm_0.2.3               
[101] viridisLite_0.3.0             vipor_0.4.5                  
[103] askpass_1.1                  
```
</div>
