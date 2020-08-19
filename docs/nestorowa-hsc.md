# Nestorowa mouse HSC (Smart-seq2) 

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

This performs an analysis of the mouse haematopoietic stem cell (HSC) dataset generated with Smart-seq2 [@nestorowa2016singlecell].

## Data loading


```r
library(scRNAseq)
sce.nest <- NestorowaHSCData()
```


```r
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
anno <- select(ens.mm.v97, keys=rownames(sce.nest), 
    keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
rowData(sce.nest) <- anno[match(rownames(sce.nest), anno$GENEID),]
```

After loading and annotation, we inspect the resulting `SingleCellExperiment` object:


```r
sce.nest
```

```
## class: SingleCellExperiment 
## dim: 46078 1920 
## metadata(0):
## assays(1): counts
## rownames(46078): ENSMUSG00000000001 ENSMUSG00000000003 ...
##   ENSMUSG00000107391 ENSMUSG00000107392
## rowData names(3): GENEID SYMBOL SEQNAME
## colnames(1920): HSPC_007 HSPC_013 ... Prog_852 Prog_810
## colData names(2): cell.type FACS
## reducedDimNames(1): diffusion
## altExpNames(1): ERCC
```

## Quality control


```r
unfiltered <- sce.nest
```

For some reason, no mitochondrial transcripts are available, so we will perform quality control using the spike-in proportions only.


```r
library(scater)
stats <- perCellQCMetrics(sce.nest)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent")
sce.nest <- sce.nest[,!qc$discard]
```

We examine the number of cells discarded for each reason.


```r
colSums(as.matrix(qc))
```

```
##              low_lib_size            low_n_features high_altexps_ERCC_percent 
##                       146                        28                       241 
##                   discard 
##                       264
```

We create some diagnostic plots for each metric (Figure \@ref(fig:unref-nest-qc-dist)).


```r
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, y="sum", colour_by="discard") +
        scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, y="detected", colour_by="discard") +
        scale_y_log10() + ggtitle("Detected features"),
    plotColData(unfiltered, y="altexps_ERCC_percent",
        colour_by="discard") + ggtitle("ERCC percent"),
    ncol=2
)
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-nest-qc-dist-1.png" alt="Distribution of each QC metric across cells in the Nestorowa HSC dataset. Each point represents a cell and is colored according to whether that cell was discarded." width="672" />
<p class="caption">(\#fig:unref-nest-qc-dist)Distribution of each QC metric across cells in the Nestorowa HSC dataset. Each point represents a cell and is colored according to whether that cell was discarded.</p>
</div>

## Normalization


```r
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce.nest)
sce.nest <- computeSumFactors(sce.nest, clusters=clusters)
sce.nest <- logNormCounts(sce.nest)
```

We examine some key metrics for the distribution of size factors, and compare it to the library sizes as a sanity check (Figure \@ref(fig:unref-nest-norm)).


```r
summary(sizeFactors(sce.nest))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.044   0.422   0.748   1.000   1.249  15.927
```


```r
plot(librarySizeFactors(sce.nest), sizeFactors(sce.nest), pch=16,
    xlab="Library size factors", ylab="Deconvolution factors", log="xy")
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-nest-norm-1.png" alt="Relationship between the library size factors and the deconvolution size factors in the Nestorowa HSC dataset." width="672" />
<p class="caption">(\#fig:unref-nest-norm)Relationship between the library size factors and the deconvolution size factors in the Nestorowa HSC dataset.</p>
</div>

## Variance modelling

We use the spike-in transcripts to model the technical noise as a function of the mean (Figure \@ref(fig:unref-nest-var)).


```r
set.seed(00010101)
dec.nest <- modelGeneVarWithSpikes(sce.nest, "ERCC")
top.nest <- getTopHVGs(dec.nest, prop=0.1)
```


```r
plot(dec.nest$mean, dec.nest$total, pch=16, cex=0.5,
    xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec.nest)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
points(curfit$mean, curfit$var, col="red")
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-nest-var-1.png" alt="Per-gene variance as a function of the mean for the log-expression values in the Nestorowa HSC dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted to the spike-ins (red)." width="672" />
<p class="caption">(\#fig:unref-nest-var)Per-gene variance as a function of the mean for the log-expression values in the Nestorowa HSC dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted to the spike-ins (red).</p>
</div>

## Dimensionality reduction


```r
set.seed(101010011)
sce.nest <- denoisePCA(sce.nest, technical=dec.nest, subset.row=top.nest)
sce.nest <- runTSNE(sce.nest, dimred="PCA")
```

We check that the number of retained PCs is sensible.


```r
ncol(reducedDim(sce.nest, "PCA"))
```

```
## [1] 9
```

## Clustering


```r
snn.gr <- buildSNNGraph(sce.nest, use.dimred="PCA")
colLabels(sce.nest) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
```


```r
table(colLabels(sce.nest))
```

```
## 
##   1   2   3   4   5   6   7   8   9 
## 203 472 258 175 142 229  20  83  74
```


```r
plotTSNE(sce.nest, colour_by="label")
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-nest-tsne-1.png" alt="Obligatory $t$-SNE plot of the Nestorowa HSC dataset, where each point represents a cell and is colored according to the assigned cluster." width="672" />
<p class="caption">(\#fig:unref-nest-tsne)Obligatory $t$-SNE plot of the Nestorowa HSC dataset, where each point represents a cell and is colored according to the assigned cluster.</p>
</div>

## Marker gene detection


```r
markers <- findMarkers(sce.nest, colLabels(sce.nest), 
    test.type="wilcox", direction="up", lfc=0.5,
    row.data=rowData(sce.nest)[,"SYMBOL",drop=FALSE])
```



To illustrate the manual annotation process, we examine the marker genes for one of the clusters.
Upregulation of _Car2_, _Hebp1_ amd hemoglobins indicates that cluster 8 contains erythroid precursors.


```r
chosen <- markers[['8']]
best <- chosen[chosen$Top <= 10,]
aucs <- getMarkerEffects(best, prefix="AUC")
rownames(aucs) <- best$SYMBOL

library(pheatmap)
pheatmap(aucs, color=viridis::plasma(100))
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-heat-nest-markers-1.png" alt="Heatmap of the AUCs for the top marker genes in cluster 8 compared to all other clusters." width="672" />
<p class="caption">(\#fig:unref-heat-nest-markers)Heatmap of the AUCs for the top marker genes in cluster 8 compared to all other clusters.</p>
</div>



## Cell type annotation


```r
library(SingleR)
mm.ref <- MouseRNAseqData()

# Renaming to symbols to match with reference row names.
renamed <- sce.nest
rownames(renamed) <- uniquifyFeatureNames(rownames(renamed),
    rowData(sce.nest)$SYMBOL)
labels <- SingleR(renamed, mm.ref, labels=mm.ref$label.fine)
```

Most clusters are not assigned to any single lineage (Figure \@ref(fig:unref-assignments-nest)), which is perhaps unsurprising given that HSCs are quite different from their terminal fates.
Cluster 8 is considered to contain erythrocytes, which is roughly consistent with our conclusions from the marker gene analysis above.


```r
tab <- table(labels$labels, colLabels(sce.nest))
pheatmap(log10(tab+10), color=viridis::viridis(100))
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-assignments-nest-1.png" alt="Heatmap of the distribution of cells for each cluster in the Nestorowa HSC dataset, based on their assignment to each label in the mouse RNA-seq references from the _SingleR_ package." width="672" />
<p class="caption">(\#fig:unref-assignments-nest)Heatmap of the distribution of cells for each cluster in the Nestorowa HSC dataset, based on their assignment to each label in the mouse RNA-seq references from the _SingleR_ package.</p>
</div>



## Miscellaneous analyses

This dataset also contains information about the protein abundances in each cell from FACS.
There is barely any heterogeneity in the chosen markers across the clusters (Figure \@ref(fig:unref-nest-facs));
this is perhaps unsurprising given that all cells should be HSCs of some sort.


```r
Y <- colData(sce.nest)$FACS
keep <- rowSums(is.na(Y))==0 # Removing NA intensities.

se.averaged <- sumCountsAcrossCells(t(Y[keep,]), 
    colLabels(sce.nest)[keep], average=TRUE)
averaged <- assay(se.averaged)

log.intensities <- log2(averaged+1)
centered <- log.intensities - rowMeans(log.intensities)
pheatmap(centered, breaks=seq(-1, 1, length.out=101))
```

<div class="figure">
<img src="nestorowa-hsc_files/figure-html/unref-nest-facs-1.png" alt="Heatmap of the centered log-average intensity for each target protein quantified by FACS in the Nestorowa HSC dataset." width="672" />
<p class="caption">(\#fig:unref-nest-facs)Heatmap of the centered log-average intensity for each target protein quantified by FACS in the Nestorowa HSC dataset.</p>
</div>



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
 [1] celldex_0.99.1              SingleR_1.3.7              
 [3] pheatmap_1.0.12             scran_1.17.15              
 [5] scater_1.17.4               ggplot2_3.3.2              
 [7] AnnotationHub_2.21.2        BiocFileCache_1.13.1       
 [9] dbplyr_1.4.4                ensembldb_2.13.1           
[11] AnnotationFilter_1.13.0     GenomicFeatures_1.41.2     
[13] AnnotationDbi_1.51.3        scRNAseq_2.3.12            
[15] SingleCellExperiment_1.11.6 SummarizedExperiment_1.19.6
[17] DelayedArray_0.15.7         matrixStats_0.56.0         
[19] Matrix_1.2-18               Biobase_2.49.0             
[21] GenomicRanges_1.41.6        GenomeInfoDb_1.25.10       
[23] IRanges_2.23.10             S4Vectors_0.27.12          
[25] BiocGenerics_0.35.4         BiocStyle_2.17.0           
[27] simpleSingleCell_1.13.16   

loaded via a namespace (and not attached):
  [1] Rtsne_0.15                    ggbeeswarm_0.6.0             
  [3] colorspace_1.4-1              ellipsis_0.3.1               
  [5] scuttle_0.99.12               bluster_0.99.1               
  [7] XVector_0.29.3                BiocNeighbors_1.7.0          
  [9] farver_2.0.3                  bit64_4.0.2                  
 [11] interactiveDisplayBase_1.27.5 codetools_0.2-16             
 [13] knitr_1.29                    Rsamtools_2.5.3              
 [15] graph_1.67.1                  shiny_1.5.0                  
 [17] BiocManager_1.30.10           compiler_4.0.2               
 [19] httr_1.4.2                    dqrng_0.2.1                  
 [21] assertthat_0.2.1              fastmap_1.0.1                
 [23] lazyeval_0.2.2                limma_3.45.10                
 [25] later_1.1.0.1                 BiocSingular_1.5.0           
 [27] htmltools_0.5.0               prettyunits_1.1.1            
 [29] tools_4.0.2                   igraph_1.2.5                 
 [31] rsvd_1.0.3                    gtable_0.3.0                 
 [33] glue_1.4.1                    GenomeInfoDbData_1.2.3       
 [35] dplyr_1.0.1                   rappdirs_0.3.1               
 [37] Rcpp_1.0.5                    vctrs_0.3.2                  
 [39] Biostrings_2.57.2             ExperimentHub_1.15.1         
 [41] rtracklayer_1.49.5            DelayedMatrixStats_1.11.1    
 [43] xfun_0.16                     stringr_1.4.0                
 [45] ps_1.3.4                      mime_0.9                     
 [47] lifecycle_0.2.0               irlba_2.3.3                  
 [49] statmod_1.4.34                XML_3.99-0.5                 
 [51] edgeR_3.31.4                  zlibbioc_1.35.0              
 [53] scales_1.1.1                  hms_0.5.3                    
 [55] promises_1.1.1                ProtGenerics_1.21.0          
 [57] RColorBrewer_1.1-2            yaml_2.2.1                   
 [59] curl_4.3                      memoise_1.1.0                
 [61] gridExtra_2.3                 biomaRt_2.45.2               
 [63] stringi_1.4.6                 RSQLite_2.2.0                
 [65] highr_0.8                     BiocVersion_3.12.0           
 [67] BiocParallel_1.23.2           rlang_0.4.7                  
 [69] pkgconfig_2.0.3               bitops_1.0-6                 
 [71] evaluate_0.14                 lattice_0.20-41              
 [73] purrr_0.3.4                   labeling_0.3                 
 [75] GenomicAlignments_1.25.3      CodeDepends_0.6.5            
 [77] cowplot_1.0.0                 bit_4.0.4                    
 [79] processx_3.4.3                tidyselect_1.1.0             
 [81] magrittr_1.5                  bookdown_0.20                
 [83] R6_2.4.1                      generics_0.0.2               
 [85] DBI_1.1.0                     pillar_1.4.6                 
 [87] withr_2.2.0                   RCurl_1.98-1.2               
 [89] tibble_3.0.3                  crayon_1.3.4                 
 [91] rmarkdown_2.3                 viridis_0.5.1                
 [93] progress_1.2.2                locfit_1.5-9.4               
 [95] grid_4.0.2                    blob_1.2.1                   
 [97] callr_3.4.3                   digest_0.6.25                
 [99] xtable_1.8-4                  httpuv_1.5.4                 
[101] openssl_1.4.2                 munsell_0.5.0                
[103] beeswarm_0.2.3                viridisLite_0.3.0            
[105] vipor_0.4.5                   askpass_1.1                  
```
</div>
