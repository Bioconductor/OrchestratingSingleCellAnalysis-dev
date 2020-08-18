# Zeisel mouse brain (STRT-Seq)

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

## Introduction

Here, we examine a heterogeneous dataset from a study of cell types in the mouse brain [@zeisel2015brain].
This contains approximately 3000 cells of varying types such as oligodendrocytes, microglia and neurons.
Individual cells were isolated using the Fluidigm C1 microfluidics system [@pollen2014lowcoverage] and library preparation was performed on each cell using a UMI-based protocol.
After sequencing, expression was quantified by counting the number of unique molecular identifiers (UMIs) mapped to each gene.

## Data loading

We obtain a `SingleCellExperiment` object for this dataset using the relevant function from the *[scRNAseq](https://bioconductor.org/packages/3.12/scRNAseq)* package.
The idiosyncrasies of the published dataset means that we need to do some extra work to merge together redundant rows corresponding to alternative genomic locations for the same gene.


```r
library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

library(scater)
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
    id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))
```

We also fetch the Ensembl gene IDs, just in case we need them later.


```r
library(org.Mm.eg.db)
rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db, 
    keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")
```

## Quality control


```r
unfiltered <- sce.zeisel
```

The original authors of the study have already removed low-quality cells prior to data publication.
Nonetheless, we compute some quality control metrics to check whether the remaining cells are satisfactory.


```r
stats <- perCellQCMetrics(sce.zeisel, subsets=list(
    Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
    "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]
```


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
    plotColData(unfiltered, y="subsets_Mt_percent",
        colour_by="discard") + ggtitle("Mito percent"),
    ncol=2
)
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-qc-dist-1.png" alt="Distribution of each QC metric across cells in the Zeisel brain dataset. Each point represents a cell and is colored according to whether that cell was discarded." width="672" />
<p class="caption">(\#fig:unref-zeisel-qc-dist)Distribution of each QC metric across cells in the Zeisel brain dataset. Each point represents a cell and is colored according to whether that cell was discarded.</p>
</div>


```r
gridExtra::grid.arrange(
    plotColData(unfiltered, x="sum", y="subsets_Mt_percent",
        colour_by="discard") + scale_x_log10(),
    plotColData(unfiltered, x="altexps_ERCC_percent", y="subsets_Mt_percent",
        colour_by="discard"),
    ncol=2
)
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-qc-comp-1.png" alt="Percentage of mitochondrial reads in each cell in the Zeisel brain dataset, compared to the total count (left) or the percentage of spike-in reads (right). Each point represents a cell and is colored according to whether that cell was discarded." width="672" />
<p class="caption">(\#fig:unref-zeisel-qc-comp)Percentage of mitochondrial reads in each cell in the Zeisel brain dataset, compared to the total count (left) or the percentage of spike-in reads (right). Each point represents a cell and is colored according to whether that cell was discarded.</p>
</div>

We also examine the number of cells removed for each reason.


```r
colSums(as.matrix(qc))
```

```
##              low_lib_size            low_n_features high_altexps_ERCC_percent 
##                         0                         3                        65 
##   high_subsets_Mt_percent                   discard 
##                       128                       189
```

## Normalization


```r
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters) 
sce.zeisel <- logNormCounts(sce.zeisel)
```


```r
summary(sizeFactors(sce.zeisel))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.119   0.486   0.831   1.000   1.321   4.509
```


```r
plot(librarySizeFactors(sce.zeisel), sizeFactors(sce.zeisel), pch=16,
    xlab="Library size factors", ylab="Deconvolution factors", log="xy")
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-norm-1.png" alt="Relationship between the library size factors and the deconvolution size factors in the Zeisel brain dataset." width="672" />
<p class="caption">(\#fig:unref-zeisel-norm)Relationship between the library size factors and the deconvolution size factors in the Zeisel brain dataset.</p>
</div>

## Variance modelling

In theory, we should block on the plate of origin for each cell.
However, only 20-40 cells are available on each plate, and the population is also highly heterogeneous.
This means that we cannot assume that the distribution of sampled cell types on each plate is the same.
Thus, to avoid regressing out potential biology, we will not block on any factors in this analysis.


```r
dec.zeisel <- modelGeneVarWithSpikes(sce.zeisel, "ERCC")
top.hvgs <- getTopHVGs(dec.zeisel, prop=0.1)
```

We see from Figure \@ref(fig:unref-zeisel-var) that the technical and total variances are much smaller than those in the read-based datasets.
This is due to the use of UMIs, which reduces the noise caused by variable PCR amplification.
Furthermore, the spike-in trend is consistently lower than the variances of the endogenous gene, which reflects the heterogeneity in gene expression across cells of different types.


```r
plot(dec.zeisel$mean, dec.zeisel$total, pch=16, cex=0.5,
    xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec.zeisel)
points(curfit$mean, curfit$var, col="red", pch=16)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-var-1.png" alt="Per-gene variance as a function of the mean for the log-expression values in the Zeisel brain dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted to the spike-in transcripts (red)." width="672" />
<p class="caption">(\#fig:unref-zeisel-var)Per-gene variance as a function of the mean for the log-expression values in the Zeisel brain dataset. Each point represents a gene (black) with the mean-variance trend (blue) fitted to the spike-in transcripts (red).</p>
</div>

## Dimensionality reduction


```r
library(BiocSingular)
set.seed(101011001)
sce.zeisel <- denoisePCA(sce.zeisel, technical=dec.zeisel, subset.row=top.hvgs)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA")
```

We have a look at the number of PCs retained by `denoisePCA()`.


```r
ncol(reducedDim(sce.zeisel, "PCA"))
```

```
## [1] 50
```

## Clustering


```r
snn.gr <- buildSNNGraph(sce.zeisel, use.dimred="PCA")
colLabels(sce.zeisel) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
```


```r
table(colLabels(sce.zeisel))
```

```
## 
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
## 283 451 114 143 599 167 191 128 350  70 199  58  39  24
```


```r
plotTSNE(sce.zeisel, colour_by="label")
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-tsne-1.png" alt="Obligatory $t$-SNE plot of the Zeisel brain dataset, where each point represents a cell and is colored according to the assigned cluster." width="672" />
<p class="caption">(\#fig:unref-zeisel-tsne)Obligatory $t$-SNE plot of the Zeisel brain dataset, where each point represents a cell and is colored according to the assigned cluster.</p>
</div>

## Interpretation



We focus on upregulated marker genes as these can quickly provide positive identification of cell type in a heterogeneous population.
We examine the table for cluster 1, in which log-fold changes are reported between cluster 1 and every other cluster.
The same output is provided for each cluster in order to identify genes that discriminate between clusters.


```r
markers <- findMarkers(sce.zeisel, direction="up")
marker.set <- markers[["1"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity
```

```
## DataFrame with 10 rows and 8 columns
##                Top      p.value          FDR summary.logFC   logFC.2   logFC.3
##          <integer>    <numeric>    <numeric>     <numeric> <numeric> <numeric>
## Atp1a3           1 1.45982e-282 7.24035e-279       3.45669 0.0398568 0.0893943
## Celf4            1 2.27030e-246 4.50404e-243       3.10465 0.3886716 0.6145023
## Gad1             1 7.44925e-232 1.34351e-228       4.57719 4.5392751 4.3003280
## Gad2             1 2.88086e-207 3.57208e-204       4.25393 4.2322487 3.8884654
## Mllt11           1 1.72982e-249 3.81309e-246       2.88363 0.5782719 1.4933128
## Ndrg4            1  0.00000e+00  0.00000e+00       3.84337 0.8887239 1.0183408
## Slc32a1          1 2.38276e-110 4.04030e-108       1.92859 1.9196173 1.8252062
## Syngr3           1 3.68257e-143 1.30462e-140       2.55531 1.0981258 1.1994793
## Atp6v1g2         2 3.04451e-204 3.55295e-201       2.50875 0.0981706 0.5203760
## Napb             2 1.10402e-231 1.82522e-228       2.81533 0.1774508 0.3046901
##            logFC.4   logFC.5
##          <numeric> <numeric>
## Atp1a3    1.241388   3.45669
## Celf4     0.869334   3.10465
## Gad1      4.050305   4.47236
## Gad2      3.769556   4.16902
## Mllt11    0.951649   2.88363
## Ndrg4     1.140041   3.84337
## Slc32a1   1.804311   1.92426
## Syngr3    1.188856   2.47696
## Atp6v1g2  0.616391   2.50875
## Napb      0.673772   2.81533
```



Figure \@ref(fig:unref-zeisel-heat-cell) indicates that most of the top markers are strongly DE in cells of cluster 1 compared to some or all of the other clusters.
We can use these markers to identify cells from cluster 1 in validation studies with an independent population of cells.
A quick look at the markers suggest that cluster 1 represents interneurons based on expression of *Gad1* and *Slc6a1* [@zeng2012largescale].


```r
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce.zeisel, features=top.markers, order_columns_by="label")
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-heat-cell-1.png" alt="Heatmap of the log-expression of the top markers for cluster 1 compared to each other cluster. Cells are ordered by cluster and the color is scaled to the log-expression of each gene in each cell." width="672" />
<p class="caption">(\#fig:unref-zeisel-heat-cell)Heatmap of the log-expression of the top markers for cluster 1 compared to each other cluster. Cells are ordered by cluster and the color is scaled to the log-expression of each gene in each cell.</p>
</div>

An alternative visualization approach is to plot the log-fold changes to all other clusters directly (Figure \@ref(fig:unref-zeisel-heat-lfc)).
This is more concise and is useful in situations involving many clusters that contain different numbers of cells.


```r
library(pheatmap)
logFCs <- getMarkerEffects(marker.set[1:50,])
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<div class="figure">
<img src="zeisel-brain_files/figure-html/unref-zeisel-heat-lfc-1.png" alt="Heatmap of the log-fold changes of the top markers for cluster 1 compared to each other cluster." width="672" />
<p class="caption">(\#fig:unref-zeisel-heat-lfc)Heatmap of the log-fold changes of the top markers for cluster 1 compared to each other cluster.</p>
</div>

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
 [1] pheatmap_1.0.12             BiocSingular_1.5.0         
 [3] scran_1.17.15               org.Mm.eg.db_3.11.4        
 [5] AnnotationDbi_1.51.3        scater_1.17.4              
 [7] ggplot2_3.3.2               scRNAseq_2.3.12            
 [9] SingleCellExperiment_1.11.6 SummarizedExperiment_1.19.6
[11] DelayedArray_0.15.7         matrixStats_0.56.0         
[13] Matrix_1.2-18               Biobase_2.49.0             
[15] GenomicRanges_1.41.6        GenomeInfoDb_1.25.10       
[17] IRanges_2.23.10             S4Vectors_0.27.12          
[19] BiocGenerics_0.35.4         BiocStyle_2.17.0           
[21] rebook_0.99.4              

loaded via a namespace (and not attached):
 [1] Rtsne_0.15                    ggbeeswarm_0.6.0             
 [3] colorspace_1.4-1              ellipsis_0.3.1               
 [5] scuttle_0.99.13               bluster_0.99.1               
 [7] XVector_0.29.3                BiocNeighbors_1.7.0          
 [9] farver_2.0.3                  bit64_4.0.2                  
[11] interactiveDisplayBase_1.27.5 codetools_0.2-16             
[13] knitr_1.29                    dbplyr_1.4.4                 
[15] graph_1.67.1                  shiny_1.5.0                  
[17] BiocManager_1.30.10           compiler_4.0.0               
[19] httr_1.4.2                    dqrng_0.2.1                  
[21] assertthat_0.2.1              fastmap_1.0.1                
[23] limma_3.45.10                 later_1.1.0.1                
[25] htmltools_0.5.0               tools_4.0.0                  
[27] rsvd_1.0.3                    igraph_1.2.5                 
[29] gtable_0.3.0                  glue_1.4.1                   
[31] GenomeInfoDbData_1.2.3        dplyr_1.0.1                  
[33] rappdirs_0.3.1                Rcpp_1.0.5                   
[35] vctrs_0.3.2                   ExperimentHub_1.15.1         
[37] DelayedMatrixStats_1.11.1     xfun_0.16                    
[39] stringr_1.4.0                 ps_1.3.4                     
[41] mime_0.9                      lifecycle_0.2.0              
[43] irlba_2.3.3                   statmod_1.4.34               
[45] XML_3.99-0.5                  AnnotationHub_2.21.2         
[47] edgeR_3.31.4                  zlibbioc_1.35.0              
[49] scales_1.1.1                  promises_1.1.1               
[51] RColorBrewer_1.1-2            yaml_2.2.1                   
[53] curl_4.3                      memoise_1.1.0                
[55] gridExtra_2.3                 stringi_1.4.6                
[57] RSQLite_2.2.0                 BiocVersion_3.12.0           
[59] highr_0.8                     BiocParallel_1.23.2          
[61] rlang_0.4.7                   pkgconfig_2.0.3              
[63] bitops_1.0-6                  evaluate_0.14                
[65] lattice_0.20-41               purrr_0.3.4                  
[67] CodeDepends_0.6.5             labeling_0.3                 
[69] cowplot_1.0.0                 bit_4.0.4                    
[71] processx_3.4.3                tidyselect_1.1.0             
[73] magrittr_1.5                  bookdown_0.20                
[75] R6_2.4.1                      generics_0.0.2               
[77] DBI_1.1.0                     pillar_1.4.6                 
[79] withr_2.2.0                   RCurl_1.98-1.2               
[81] tibble_3.0.3                  crayon_1.3.4                 
[83] BiocFileCache_1.13.1          rmarkdown_2.3                
[85] viridis_0.5.1                 locfit_1.5-9.4               
[87] grid_4.0.0                    blob_1.2.1                   
[89] callr_3.4.3                   digest_0.6.25                
[91] xtable_1.8-4                  httpuv_1.5.4                 
[93] munsell_0.5.0                 beeswarm_0.2.3               
[95] viridisLite_0.3.0             vipor_0.4.5                  
```
</div>
