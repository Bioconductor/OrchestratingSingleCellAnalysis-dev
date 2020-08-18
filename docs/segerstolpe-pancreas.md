# Segerstolpe human pancreas (Smart-seq2)

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

This performs an analysis of the @segerstolpe2016singlecell dataset,
consisting of human pancreas cells from various donors.

## Data loading


```r
library(scRNAseq)
sce.seger <- SegerstolpePancreasData()
```


```r
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce.seger)$symbol
ens.id <- mapIds(edb, keys=symbols, keytype="SYMBOL", column="GENEID")
ens.id <- ifelse(is.na(ens.id), symbols, ens.id)

# Removing duplicated rows.
keep <- !duplicated(ens.id)
sce.seger <- sce.seger[keep,]
rownames(sce.seger) <- ens.id[keep]
```

We simplify the names of some of the relevant column metadata fields for ease of access.
Some editing of the cell type labels is necessary for consistency with other data sets.


```r
emtab.meta <- colData(sce.seger)[,c("cell type", 
    "individual", "single cell well quality")]
colnames(emtab.meta) <- c("CellType", "Donor", "Quality")
colData(sce.seger) <- emtab.meta

sce.seger$CellType <- gsub(" cell", "", sce.seger$CellType)
sce.seger$CellType <- paste0(
    toupper(substr(sce.seger$CellType, 1, 1)),
    substring(sce.seger$CellType, 2))
```

## Quality control


```r
unfiltered <- sce.seger
```

We remove low quality cells that were marked by the authors.
We then perform additional quality control as some of the remaining cells still have very low counts and numbers of detected features.
For some batches that seem to have a majority of low-quality cells (Figure \@ref(fig:unref-seger-qc-dist)), we use the other batches to define an appropriate threshold via `subset=`.


```r
low.qual <- sce.seger$Quality == "low quality cell"

library(scater)
stats <- perCellQCMetrics(sce.seger)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
    batch=sce.seger$Donor,
    subset=!sce.seger$Donor %in% c("HP1504901", "HP1509101"))

sce.seger <- sce.seger[,!(qc$discard | low.qual)]
```


```r
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, x="Donor", y="sum", colour_by="discard") +
        scale_y_log10() + ggtitle("Total count") +
        theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x="Donor", y="detected", colour_by="discard") +
        scale_y_log10() + ggtitle("Detected features") +
        theme(axis.text.x = element_text(angle = 90)),
    plotColData(unfiltered, x="Donor", y="altexps_ERCC_percent",
        colour_by="discard") + ggtitle("ERCC percent") +
        theme(axis.text.x = element_text(angle = 90)),
    ncol=2
)
```

<div class="figure">
<img src="segerstolpe-pancreas_files/figure-html/unref-seger-qc-dist-1.png" alt="Distribution of each QC metric across cells from each donor of the Segerstolpe pancreas dataset. Each point represents a cell and is colored according to whether that cell was discarded." width="960" />
<p class="caption">(\#fig:unref-seger-qc-dist)Distribution of each QC metric across cells from each donor of the Segerstolpe pancreas dataset. Each point represents a cell and is colored according to whether that cell was discarded.</p>
</div>


```r
colSums(as.matrix(qc))
```

```
##              low_lib_size            low_n_features high_altexps_ERCC_percent 
##                       788                      1056                      1031 
##                   discard 
##                      1246
```

## Normalization

We don't normalize the spike-ins as there are some cells with no spike-in counts.


```r
library(scran)
clusters <- quickCluster(sce.seger)
sce.seger <- computeSumFactors(sce.seger, clusters=clusters)
sce.seger <- logNormCounts(sce.seger) 
```


```r
summary(sizeFactors(sce.seger))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.014   0.390   0.708   1.000   1.332  11.182
```


```r
plot(librarySizeFactors(sce.seger), sizeFactors(sce.seger), pch=16,
    xlab="Library size factors", ylab="Deconvolution factors", log="xy")
```

<div class="figure">
<img src="segerstolpe-pancreas_files/figure-html/unref-seger-norm-1.png" alt="Relationship between the library size factors and the deconvolution size factors in the Segerstolpe pancreas dataset." width="672" />
<p class="caption">(\#fig:unref-seger-norm)Relationship between the library size factors and the deconvolution size factors in the Segerstolpe pancreas dataset.</p>
</div>

## Variance modelling

We do not use cells with no spike-ins for variance modelling.
Donor AZ also has very low spike-in counts and is subsequently ignored.


```r
for.hvg <- sce.seger[,librarySizeFactors(altExp(sce.seger)) > 0
    & sce.seger$Donor!="AZ"]
dec.seger <- modelGeneVarWithSpikes(for.hvg, "ERCC", block=for.hvg$Donor)
chosen.hvgs <- getTopHVGs(dec.seger, n=2000)
```


```r
par(mfrow=c(3,3))
blocked.stats <- dec.seger$per.block
for (i in colnames(blocked.stats)) {
    current <- blocked.stats[[i]]
    plot(current$mean, current$total, main=i, pch=16, cex=0.5,
        xlab="Mean of log-expression", ylab="Variance of log-expression")
    curfit <- metadata(current)
    points(curfit$mean, curfit$var, col="red", pch=16)
    curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}
```

<div class="figure">
<img src="segerstolpe-pancreas_files/figure-html/unref-seger-variance-1.png" alt="Per-gene variance as a function of the mean for the log-expression values in the Grun pancreas dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted to the spike-in transcripts (red) separately for each donor." width="672" />
<p class="caption">(\#fig:unref-seger-variance)Per-gene variance as a function of the mean for the log-expression values in the Grun pancreas dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted to the spike-in transcripts (red) separately for each donor.</p>
</div>

## Dimensionality reduction


```r
library(BiocSingular)
set.seed(101011001)
sce.seger <- runPCA(sce.seger, subset_row=chosen.hvgs, ncomponents=25)
sce.seger <- runTSNE(sce.seger, dimred="PCA")
```

## Clustering


```r
snn.gr <- buildSNNGraph(sce.seger, use.dimred="PCA")
colLabels(sce.seger) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
```

We see a strong donor effect, which suggests that should have called `fastMNN()` at some point.
(But hey, we already did that for the Muraro and Grun analyses, so where's the fun in doing that again?)


```r
tab <- table(Cluster=colLabels(sce.seger), Donor=sce.seger$Donor)
library(pheatmap)
pheatmap(log10(tab+10), color=viridis::viridis(100))
```

<div class="figure">
<img src="segerstolpe-pancreas_files/figure-html/unref-seger-heat-1-1.png" alt="Heatmap of the frequency of cells from each donor in each cluster." width="672" />
<p class="caption">(\#fig:unref-seger-heat-1)Heatmap of the frequency of cells from each donor in each cluster.</p>
</div>


```r
tab <- table(Cluster=colLabels(sce.seger), Donor=sce.seger$CellType)
pheatmap(log10(tab+10), color=viridis::viridis(100))
```

<div class="figure">
<img src="segerstolpe-pancreas_files/figure-html/unref-seger-heat-2-1.png" alt="Heatmap of the frequency of cells from each cell type label in each cluster." width="672" />
<p class="caption">(\#fig:unref-seger-heat-2)Heatmap of the frequency of cells from each cell type label in each cluster.</p>
</div>


```r
gridExtra::grid.arrange(
    plotTSNE(sce.seger, colour_by="label"),
    plotTSNE(sce.seger, colour_by="Donor"),
    ncol=2
)
```

<div class="figure">
<img src="segerstolpe-pancreas_files/figure-html/unref-grun-tsne-1.png" alt="Obligatory $t$-SNE plots of the Segerstolpe pancreas dataset. Each point represents a cell that is colored by cluster (left) or batch (right)." width="672" />
<p class="caption">(\#fig:unref-grun-tsne)Obligatory $t$-SNE plots of the Segerstolpe pancreas dataset. Each point represents a cell that is colored by cluster (left) or batch (right).</p>
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
 [1] pheatmap_1.0.12             BiocSingular_1.5.0         
 [3] scran_1.17.3                scater_1.17.2              
 [5] ggplot2_3.3.2               ensembldb_2.13.1           
 [7] AnnotationFilter_1.13.0     GenomicFeatures_1.41.0     
 [9] AnnotationDbi_1.51.1        AnnotationHub_2.21.1       
[11] BiocFileCache_1.13.0        dbplyr_1.4.4               
[13] scRNAseq_2.3.8              SingleCellExperiment_1.11.6
[15] SummarizedExperiment_1.19.5 DelayedArray_0.15.6        
[17] matrixStats_0.56.0          Matrix_1.2-18              
[19] Biobase_2.49.0              GenomicRanges_1.41.5       
[21] GenomeInfoDb_1.25.5         IRanges_2.23.10            
[23] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[25] BiocStyle_2.17.0            simpleSingleCell_1.13.5    

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
 [55] RColorBrewer_1.1-2            yaml_2.2.1                   
 [57] curl_4.3                      memoise_1.1.0                
 [59] gridExtra_2.3                 biomaRt_2.45.1               
 [61] stringi_1.4.6                 RSQLite_2.2.0                
 [63] highr_0.8                     BiocVersion_3.12.0           
 [65] BiocParallel_1.23.0           rlang_0.4.6                  
 [67] pkgconfig_2.0.3               bitops_1.0-6                 
 [69] evaluate_0.14                 lattice_0.20-41              
 [71] purrr_0.3.4                   labeling_0.3                 
 [73] GenomicAlignments_1.25.3      CodeDepends_0.6.5            
 [75] cowplot_1.0.0                 bit_1.1-15.2                 
 [77] processx_3.4.2                tidyselect_1.1.0             
 [79] magrittr_1.5                  bookdown_0.20                
 [81] R6_2.4.1                      generics_0.0.2               
 [83] DBI_1.1.0                     pillar_1.4.4                 
 [85] withr_2.2.0                   RCurl_1.98-1.2               
 [87] tibble_3.0.1                  crayon_1.3.4                 
 [89] rmarkdown_2.3                 viridis_0.5.1                
 [91] progress_1.2.2                locfit_1.5-9.4               
 [93] grid_4.0.2                    blob_1.2.1                   
 [95] callr_3.4.3                   digest_0.6.25                
 [97] xtable_1.8-4                  httpuv_1.5.4                 
 [99] openssl_1.4.2                 munsell_0.5.0                
[101] beeswarm_0.2.3                viridisLite_0.3.0            
[103] vipor_0.4.5                   askpass_1.1                  
```
</div>
