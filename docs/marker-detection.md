---
output:
  html_document
bibliography: ref.bib
---

# Marker gene detection {#marker-detection}

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

To interpret our clustering results from Chapter \@ref(clustering), we identify the genes that drive separation between clusters.
These marker genes allow us to assign biological meaning to each cluster based on their functional annotation.
In the most obvious case, the marker genes for each cluster are _a priori_ associated with particular cell types, allowing us to treat the clustering as a proxy for cell type identity.
The same principle can be applied to discover more subtle differences between clusters (e.g., changes in activation or differentiation state) based on the behavior of genes in the affected pathways.

Identification of marker genes is usually based around the retrospective detection of differential expression between clusters.
Genes that are more strongly DE are more likely to have caused separate clustering of cells in the first place.
Several different statistical tests are available to quantify the differences in expression profiles, and different approaches can be used to consolidate test results into a single ranking of genes for each cluster.
These choices parametrize the theoretical differences between the various marker detection strategies presented in this chapter.
We will demonstrate using the 10X PBMC dataset:

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

## Pairwise tests between clusters

### Motivation 

Our general strategy is to perform DE tests between pairs of clusters and then combine results into a single ranking of marker genes for each cluster.
We deliberately use pairwise comparisons rather than comparing each cluster to the average of all other cells; the latter approach is sensitive to the population composition, which introduces an element of unpredictability to the marker sets due to variation in cell type abundances.
(In the worst case, the presence of one subpopulation containing a majority of the cells will drive the selection of top markers for every other cluster, pushing out useful genes that can distinguish between the smaller subpopulations.)
Moreover, pairwise comparisons naturally provide more information to interpret of the utility of a marker, e.g., by providing log-fold changes to indicate which clusters are distinguished by each gene.

For this section, we will use the Welch $t$-test to perform our DE testing between clusters.
This is an easy choice as it is quickly computed and has good statistical properties for large numbers of cells [@soneson2018bias].
However, the same approach can also be applied with any pairwise statistical test, as discussed in Section \@ref(marker-tests).

### Combining pairwise statistics per cluster

#### Looking for any differences

We perform pairwise $t$-tests between clusters for each gene using the `findMarkers()` function, which returns a list of `DataFrame`s containing ranked candidate markers for each cluster.
The function will automatically retrieve the cluster identities from `sce.pbmc` using the `colLabels()` function, though we can easily specify other clustering schemes by explicitly supplying them via the `groups=` argument.


```r
library(scran)
markers.pbmc <- findMarkers(sce.pbmc)
markers.pbmc
```

```
## List of length 16
## names(16): 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
```



The default philosophy of `findMarkers()` is to identify a combination of marker genes that - together - uniquely define one cluster against the rest.
To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
We will demonstrate on cluster 7; the relevant `DataFrame` contains log~2~-fold changes of expression in cluster 7 over each other cluster, along with several statistics obtained by combining $p$-values [@simes1986improved] across the pairwise comparisons involving 7.


```r
chosen <- "7"
interesting <- markers.pbmc[[chosen]]
colnames(interesting)
```

```
##  [1] "Top"           "p.value"       "FDR"           "summary.logFC"
##  [5] "logFC.1"       "logFC.2"       "logFC.3"       "logFC.4"      
##  [9] "logFC.5"       "logFC.6"       "logFC.8"       "logFC.9"      
## [13] "logFC.10"      "logFC.11"      "logFC.12"      "logFC.13"     
## [17] "logFC.14"      "logFC.15"      "logFC.16"
```

Of particular interest is the `Top` field.
The set of genes with `Top` $\le X$ is the union of the top $X$ genes (ranked by $p$-value) from each pairwise comparison involving cluster 7.
For example, the set of all genes with `Top` values of 1 contains the gene with the lowest $p$-value from each comparison.
Similarly, the set of genes with `Top` values less than or equal to 10 contains the top 10 genes from each comparison.
The `Top` field represents `findMarkers()`'s approach to consolidating multiple pairwise comparisons into a single ranking for each cluster; each `DataFrame` produced by `findMarkers()` will order genes based on the `Top` value by default.


```r
interesting[1:10,1:4]
```

```
## DataFrame with 10 rows and 4 columns
##                Top      p.value          FDR summary.logFC
##          <integer>    <numeric>    <numeric>     <numeric>
## S100A4           1  2.59737e-38  1.27018e-36      -4.27560
## TAGLN2           1  8.65033e-28  2.44722e-26       5.07327
## FCGR3A           1  8.84356e-63  1.15048e-60      -3.07121
## GZMA             1 1.15392e-120 7.20000e-118      -1.92877
## HLA-DQA1         1  3.43640e-83  8.90663e-81      -3.54890
## TMSB4X           1  9.83227e-36  4.25820e-34       4.28970
## FCN1             1 1.74313e-239 9.78883e-236      -2.77594
## TRAC             1  0.00000e+00  0.00000e+00      -2.44793
## RPL17            1  2.95529e-71  5.18622e-69      -2.86310
## CD79A            1  0.00000e+00  0.00000e+00      -2.98030
```



We use the `Top` field to identify a set of genes that is guaranteed to distinguish cluster 7 from any other cluster.
Here, we examine the top 6 genes from each pairwise comparison (Figure \@ref(fig:heat-basic-pbmc)).
Some inspection of the most upregulated genes suggest that cluster 9 contains platelets or their precursors, based on the expression of platelet factor 4 (_PF4_) and pro-platelet basic protein (_PPBP_).


```r
best.set <- interesting[interesting$Top <= 6,]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<div class="figure">
<img src="marker-detection_files/figure-html/heat-basic-pbmc-1.png" alt="Heatmap of log-fold changes for cluster 7 over all other clusters. Colours are capped at -5 and 5 to preserve dynamic range." width="672" />
<p class="caption">(\#fig:heat-basic-pbmc)Heatmap of log-fold changes for cluster 7 over all other clusters. Colours are capped at -5 and 5 to preserve dynamic range.</p>
</div>



Each `DataFrame` also contains several other statistics that may be of interest.
The `summary.logFC` field provides a convenient summary of the direction and effect size for each gene, and is defined here as the log-fold change from the comparison with the lowest $p$-value.
The `p.value` field contains the combined $p$-value that is obtained by applying Simes' method to the pairwise $p$-values for each gene and represents the evidence against the joint null hypothesis, i.e., that the gene is not DE between cluster 7 and any other cluster.
Examination of these statistics permits a quick evaluation of the suitability of a candidate marker; if both of these metrics are poor (small log-fold change, large $p$-value), the gene can most likely be dismissed.

#### Finding cluster-specific markers 

By default, `findMarkers()` will give a high ranking to genes that are differentially expressed in any pairwise comparison.
This is because a gene only needs a very low $p$-value in a single pairwise comparison to achieve a low `Top` value.
A more stringent approach would only consider genes that are differentially expressed in all pairwise comparisons involving the cluster of interest.
To achieve this, we set `pval.type="all"` in `findMarkers()` to use an intersection-union test [@berger1996bioequivalence] where the combined $p$-value for each gene is the maximum of the $p$-values from all pairwise comparisons.
A gene will only achieve a low combined $p$-value if it is strongly DE in all comparisons to other clusters.


```r
# Set direction='up' to only consider upregulated genes as potential markers.
markers.pbmc.up3 <- findMarkers(sce.pbmc, pval.type="all", direction="up")
interesting.up3 <- markers.pbmc.up3[[chosen]]
interesting.up3[1:10,1:3]
```

```
## DataFrame with 10 rows and 3 columns
##               p.value         FDR summary.logFC
##             <numeric>   <numeric>     <numeric>
## SDPR      2.86451e-23 9.65166e-19       5.49695
## NRGN      5.91075e-23 9.95784e-19       4.71034
## TAGLN2    4.41617e-21 4.95995e-17       3.61864
## PPBP      1.36171e-20 1.14703e-16       6.34717
## GNG11     1.23155e-19 8.29918e-16       5.35904
## HIST1H2AC 3.94013e-19 2.21265e-15       5.46400
## TUBB1     9.64049e-19 4.64038e-15       4.92204
## PF4       1.87045e-14 7.87785e-11       6.49672
## CLU       8.75900e-13 3.27918e-09       3.90273
## RGS18     8.88042e-12 2.99217e-08       3.63236
```

This strategy will only report genes that are highly specific to the cluster of interest.
When it works, it can be highly effective as it generates a small focused set of candidate markers. 
However, any gene that is expressed at the same level in two or more clusters will simply not be detected. 
This is likely to discard many interesting genes, especially if the clusters are finely resolved with weak separation.
To give a concrete example, consider a mixed population of CD4^+^-only, CD8^+^-only, double-positive and double-negative T cells.
With `pval.type="all"`, neither _Cd4_ or _Cd8_ would be detected as subpopulation-specific markers because each gene is expressed in two subpopulations.
In comparison, `pval.type="any"` will detect both of these genes as they will be DE between at least one pair of subpopulations.

#### Balancing stringency and generality

If `pval.type="all"` is too stringent yet `pval.type="any"` is too generous, a compromise is to set `pval.type="some"`.
For each gene, we apply the Holm-Bonferroni correction across its $p$-values and take the middle-most value as the combined $p$-value.
This effectively tests the global null hypothesis that at least 50% of the individual pairwise comparisons exhibit no DE.
We then rank the genes by their combined $p$-values to obtain an ordered set of marker candidates.
The aim is to improve the conciseness of the top markers for defining a cluster while mitigating the risk of discarding useful genes that are not DE to all other clusters.
The downside is that taking this compromise position sacrifices the theoretical guarantees offered at the other two extremes.


```r
markers.pbmc.up4 <- findMarkers(sce.pbmc, pval.type="some", direction="up")
interesting.up4 <- markers.pbmc.up4[[chosen]]
interesting.up4[1:10,1:3]
```

```
## DataFrame with 10 rows and 3 columns
##               p.value         FDR summary.logFC
##             <numeric>   <numeric>     <numeric>
## PF4       5.23414e-32 1.76359e-27       6.86288
## TMSB4X    4.52854e-25 7.62923e-21       2.90202
## TAGLN2    2.31252e-24 2.59727e-20       4.88268
## NRGN      1.08964e-22 9.17861e-19       5.00827
## SDPR      2.47896e-22 1.67052e-18       5.60445
## PPBP      8.57360e-20 4.81465e-16       6.50103
## CCL5      5.66181e-19 2.72527e-15       5.30774
## GNG11     8.98381e-19 3.59373e-15       5.47403
## GPX1      9.59922e-19 3.59373e-15       4.86299
## HIST1H2AC 2.85071e-18 9.60519e-15       5.53275
```

In both cases, a different method is used to compute the summary effect size compared to `pval.type="any"`.
For `pval.type="all"`, the summary log-fold change is defined as that corresponding to the pairwise comparison with the largest $p$-value, while for `pval.type="some"`, it is defined as the log-fold change for the comparison with the middle-most $p$-value.
This reflects the calculation of the combined $p$-value and avoids focusing on genes with strong changes in only one comparison.

### Using the log-fold change 

The default `findMarkers()` call considers both up- and downregulated genes to be potential markers.
However, downregulated genes are less appealing as markers as it is more difficult to interpret and experimentally validate an absence of expression.
To focus on up-regulated markers, we can instead perform a one-sided $t$-test to identify genes that are upregulated in each cluster compared to the others.
This is achieved by setting `direction="up"` in the `findMarkers()` call.


```r
markers.pbmc.up <- findMarkers(sce.pbmc, direction="up")
interesting.up <- markers.pbmc.up[[chosen]]
interesting.up[1:10,1:4]
```

```
## DataFrame with 10 rows and 4 columns
##              Top     p.value         FDR summary.logFC
##        <integer>   <numeric>   <numeric>     <numeric>
## TAGLN2         1 4.32517e-28 4.85774e-24       5.07327
## PF4            1 4.78929e-35 8.06851e-31       6.71811
## TMSB4X         1 4.91613e-36 1.65644e-31       4.28970
## NRGN           2 1.35810e-23 9.15195e-20       4.86347
## B2M            2 8.10863e-25 6.83030e-21       2.40365
## SDPR           3 2.32759e-23 1.30710e-19       5.54225
## GPX1           4 4.74597e-21 2.28444e-17       5.71604
## PPBP           5 8.85410e-21 3.31478e-17       6.41411
## ACTB           6 6.22981e-21 2.62384e-17       3.79868
## GNG11          6 9.05522e-20 3.05107e-16       5.48735
```

The $t$-test also allows us to specify a non-zero log-fold change as the null hypothesis.
This allows us to consider the magnitude of the log-fold change in our $p$-value calculations, in a manner that is more rigorous than simply filtering directly on the log-fold changes [@mccarthy2009treat].
(Specifically, a simple threshold does not consider the variance and can enrich for genes that have both large log-fold changes and large variances.) 
We perform this by setting `lfc=` in our `findMarkers()` call - when combined with `direction=`, this tests for genes with log-fold changes that are significantly greater than 1:


```r
markers.pbmc.up2 <- findMarkers(sce.pbmc, direction="up", lfc=1)
interesting.up2 <- markers.pbmc.up2[[chosen]]
interesting.up2[1:10,1:4]
```

```
## DataFrame with 10 rows and 4 columns
##                 Top     p.value         FDR summary.logFC
##           <integer>   <numeric>   <numeric>     <numeric>
## TAGLN2            1 9.48392e-23 1.06517e-18       5.07327
## PF4               1 2.19317e-31 7.38966e-27       6.71811
## SDPR              2 5.42215e-20 4.56735e-16       5.54225
## TMSB4X            2 9.90003e-28 1.66786e-23       4.28970
## NRGN              3 9.24786e-20 6.23195e-16       4.95866
## GPX1              4 7.73653e-18 3.72392e-14       5.71604
## PPBP              4 5.53317e-18 3.10725e-14       6.51025
## GNG11             6 1.56486e-16 6.59079e-13       5.48735
## CCL5              6 2.72759e-16 1.02115e-12       5.39815
## HIST1H2AC         7 5.56042e-16 1.87353e-12       5.57765
```

These two settings yield a more focused set of candidate marker genes that are upregulated in cluster 7 (Figure \@ref(fig:heat-focused-pbmc)).


```r
best.set <- interesting.up2[interesting.up2$Top <= 5,]
logFCs <- getMarkerEffects(best.set)

library(pheatmap)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

<div class="figure">
<img src="marker-detection_files/figure-html/heat-focused-pbmc-1.png" alt="Heatmap of log-fold changes for cluster 7 over all other clusters. Colours are capped at -5 and 5 to preserve dynamic range." width="672" />
<p class="caption">(\#fig:heat-focused-pbmc)Heatmap of log-fold changes for cluster 7 over all other clusters. Colours are capped at -5 and 5 to preserve dynamic range.</p>
</div>

Of course, this increased stringency is not without cost.
If only upregulated genes are requested from `findMarkers()`, any cluster defined by downregulation of a marker gene will not contain that gene among the top set of features in its `DataFrame`.
This is occasionally relevant for subtypes or other states that are defined by low expression of particular genes^[Standard operating procedure is to (i) experience a brief but crushing bout of disappointment due to the poor quality of upregulated candidate markers, (ii) rage-quit, and (iii) remember to check the genes that are changing in the other direction.].
Similarly, setting an excessively high log-fold change threshold may discard otherwise useful genes.
For example, a gene upregulated in a small proportion of cells of a cluster will have a small log-fold change but can still be an effective marker if the focus is on specificity rather than sensitivity.

## Alternative testing regimes {#marker-tests}

### Using the Wilcoxon rank sum test

The Wilcoxon rank sum test (also known as the Wilcoxon-Mann-Whitney test, or WMW test) is another widely used method for pairwise comparisons between groups of observations.
Its strength lies in the fact that it directly assesses separation between the expression distributions of different clusters.
The WMW test statistic is proportional to the area-under-the-curve (AUC), i.e., the concordance probability, which is the probability of a random cell from one cluster having higher expression than a random cell from another cluster.
In a pairwise comparison, AUCs of 1 or 0 indicate that the two clusters have perfectly separated expression distributions.
Thus, the WMW test directly addresses the most desirable property of a candidate marker gene, while the $t$ test only does so indirectly via the difference in the means and the intra-group variance.

We perform WMW tests by again using the `findMarkers()` function, this time with `test="wilcox"`.
This returns a list of `DataFrame`s containing ranked candidate markers for each cluster.
The `direction=`, `lfc=` and `pval.type=` arguments can be specified and have the same interpretation as described for $t$-tests.
We demonstrate below by detecting upregulated genes in each cluster with `direction="up"`.


```r
markers.pbmc.wmw <- findMarkers(sce.pbmc, test="wilcox", direction="up")
names(markers.pbmc.wmw)
```

```
##  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
## [16] "16"
```

To explore the results in more detail, we focus on the `DataFrame` for cluster 7.
The interpretation of `Top` is the same as described for $t$-tests, and Simes' method is again used to combine $p$-values across pairwise comparisons.
If we want more focused sets, we can also change `pval.type=` as previously described.


```r
interesting.wmw <- markers.pbmc.wmw[[chosen]]
interesting.wmw[1:10,1:4]
```

```
## DataFrame with 10 rows and 4 columns
##                 Top      p.value          FDR summary.AUC
##           <integer>    <numeric>    <numeric>   <numeric>
## PF4               1 1.02312e-179 3.44731e-175    0.989080
## TMSB4X            1  3.14604e-29  1.20457e-26    0.998195
## SDPR              2 1.36598e-159 2.30126e-155    0.956221
## NRGN              2 2.77288e-142 1.86859e-138    0.966865
## TAGLN2            3  1.58373e-29  6.20491e-27    0.967680
## PPBP              3 1.35961e-147 1.52702e-143    0.934256
## GNG11             3 4.00798e-139 2.25075e-135    0.934030
## TUBB1             3 3.23282e-146 2.72317e-142    0.923386
## HIST1H2AC         5  2.49447e-97  5.60325e-94    0.932300
## B2M               5  3.15826e-25  9.85320e-23    0.968938
```

The `DataFrame` contains the AUCs from comparing cluster 7 to every other cluster (Figure \@ref(fig:heat-wmw-pbmc)).
A value greater than 0.5 indicates that the gene is upregulated in the current cluster compared to the other cluster,
while values less than 0.5 correspond to downregulation.
We would typically expect AUCs of 0.7-0.8 for a strongly upregulated candidate marker.


```r
best.set <- interesting.wmw[interesting.wmw$Top <= 5,]
AUCs <- getMarkerEffects(best.set, prefix="AUC")

library(pheatmap)
pheatmap(AUCs, breaks=seq(0, 1, length.out=21),
    color=viridis::viridis(21))
```

<div class="figure">
<img src="marker-detection_files/figure-html/heat-wmw-pbmc-1.png" alt="Heatmap of AUCs for cluster 7 compared to all other clusters." width="672" />
<p class="caption">(\#fig:heat-wmw-pbmc)Heatmap of AUCs for cluster 7 compared to all other clusters.</p>
</div>

One practical advantage of the WMW test over the Welch $t$-test is that it is symmetric with respect to differences in the size of the groups being compared.
This means that, all else being equal, the top-ranked genes on each side of a DE comparison will have similar expression profiles regardless of the number of cells in each group.
In contrast, the $t$-test will favor genes where the larger group has the higher relative variance as this increases the estimated degrees of freedom and decreases the resulting $p$-value.
This can lead to unappealing rankings when the aim is to identify genes upregulated in smaller groups.
The WMW test is not completely immune to variance effects - for example, it will slightly favor detection of DEGs at low average abundance where the greater number of ties at zero deflates the approximate variance of the rank sum statistic - but this is relatively benign as the selected genes are still fairly interesting.
We observe both of these effects in a comparison between alpha and gamma cells in the human pancreas data set from @lawlor2017singlecell (Figure \@ref(fig:comparative-markers-tw)).



<button class="aaron-collapse">View history</button>
<div class="aaron-content">
   
```r
#--- loading ---#
library(scRNAseq)
sce.lawlor <- LawlorPancreasData()

#--- gene-annotation ---#
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(edb, keys=rownames(sce.lawlor), keytype="GENEID", 
    columns=c("SYMBOL", "SEQNAME"))
rowData(sce.lawlor) <- anno[match(rownames(sce.lawlor), anno[,1]),-1]

#--- quality-control ---#
library(scater)
stats <- perCellQCMetrics(sce.lawlor, 
    subsets=list(Mito=which(rowData(sce.lawlor)$SEQNAME=="MT")))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
    batch=sce.lawlor$`islet unos id`)
sce.lawlor <- sce.lawlor[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.lawlor)
sce.lawlor <- computeSumFactors(sce.lawlor, clusters=clusters)
sce.lawlor <- logNormCounts(sce.lawlor)
```

</div>


```r
marker.lawlor.t <- findMarkers(sce.lawlor, groups=sce.lawlor$`cell type`, 
    direction="up", restrict=c("Alpha", "Gamma/PP"))
marker.lawlor.w <- findMarkers(sce.lawlor, groups=sce.lawlor$`cell type`, 
    direction="up", restrict=c("Alpha", "Gamma/PP"), test.type="wilcox")

# Upregulated in alpha:
marker.alpha.t <- marker.lawlor.t$Alpha
marker.alpha.w <- marker.lawlor.w$Alpha
chosen.alpha.t <- rownames(marker.alpha.t)[1:20]
chosen.alpha.w <- rownames(marker.alpha.w)[1:20]
u.alpha.t <- setdiff(chosen.alpha.t, chosen.alpha.w)
u.alpha.w <- setdiff(chosen.alpha.w, chosen.alpha.t)

# Upregulated in gamma:
marker.gamma.t <- marker.lawlor.t$`Gamma/PP`
marker.gamma.w <- marker.lawlor.w$`Gamma/PP`
chosen.gamma.t <- rownames(marker.gamma.t)[1:20]
chosen.gamma.w <- rownames(marker.gamma.w)[1:20]
u.gamma.t <- setdiff(chosen.gamma.t, chosen.gamma.w)
u.gamma.w <- setdiff(chosen.gamma.w, chosen.gamma.t)

# Examining all uniquely detected markers in each direction.
library(scater)
subset <- sce.lawlor[,sce.lawlor$`cell type` %in% c("Alpha", "Gamma/PP")]
gridExtra::grid.arrange(
    plotExpression(subset, x="cell type", features=u.alpha.t, ncol=2) +
        ggtitle("Upregulated in alpha, t-test-only"),
    plotExpression(subset, x="cell type", features=u.alpha.w, ncol=2) +
        ggtitle("Upregulated in alpha, WMW-test-only"),
    plotExpression(subset, x="cell type", features=u.gamma.t, ncol=2) +
        ggtitle("Upregulated in gamma, t-test-only"),
    plotExpression(subset, x="cell type", features=u.gamma.w, ncol=2) +
        ggtitle("Upregulated in gamma, WMW-test-only"),
    ncol=2
)
```

<div class="figure">
<img src="marker-detection_files/figure-html/comparative-markers-tw-1.png" alt="Distribution of expression values for alpha or gamma cell-specific markers in the GSE86469 human pancreas dataset. Each panel focuses on the genes that were uniquely ranked in the top 20 candidate markers by either the t-test or WMW test." width="672" />
<p class="caption">(\#fig:comparative-markers-tw)Distribution of expression values for alpha or gamma cell-specific markers in the GSE86469 human pancreas dataset. Each panel focuses on the genes that were uniquely ranked in the top 20 candidate markers by either the t-test or WMW test.</p>
</div>



The main disadvantage of the WMW test is that the AUCs are much slower to compute compared to $t$-statistics.
This may be inconvenient for interactive analyses involving multiple iterations of marker detection.
We can mitigate this to some extent by parallelizing these calculations using the `BPPARAM=` argument in `findMarkers()`.

### Using a binomial test

The binomial test identifies genes that differ in the proportion of expressing cells between clusters.
(For the purposes of this section, a cell is considered to express a gene simply if it has non-zero expression for that gene.)
This represents a much more stringent definition of marker genes compared to the other methods, as differences in expression between clusters are effectively ignored if both distributions of expression values are not near zero.
The premise is that genes are more likely to contribute to important biological decisions if they were active in one cluster and silent in another, compared to more subtle "tuning" effects from changing the expression of an active gene.
From a practical perspective, a binary measure of presence/absence is easier to validate.

We perform pairwise binomial tests between clusters using the `findMarkers()` function with `test="binom"`.
This returns a list of `DataFrame`s containing marker statistics for each cluster such as the `Top` rank and its $p$-value.
Here, the effect size is reported as the log-fold change in this proportion between each pair of clusters.
Large positive log-fold changes indicate that the gene is more frequently expressed in one cluster compared to the other.
We focus on genes that are upregulated in each cluster compared to the others by setting `direction="up"`.


```r
markers.pbmc.binom <- findMarkers(sce.pbmc, test="binom", direction="up")
names(markers.pbmc.binom)
```

```
##  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
## [16] "16"
```

```r
interesting.binom <- markers.pbmc.binom[[chosen]]
colnames(interesting.binom)
```

```
##  [1] "Top"           "p.value"       "FDR"           "summary.logFC"
##  [5] "logFC.1"       "logFC.2"       "logFC.3"       "logFC.4"      
##  [9] "logFC.5"       "logFC.6"       "logFC.8"       "logFC.9"      
## [13] "logFC.10"      "logFC.11"      "logFC.12"      "logFC.13"     
## [17] "logFC.14"      "logFC.15"      "logFC.16"
```

Figure \@ref(fig:viol-de-binom) confirms that the top genes exhibit strong differences in the proportion of expressing cells in cluster 7 compared to the others. 


```r
library(scater)
top.genes <- head(rownames(interesting.binom))
plotExpression(sce.pbmc, x="label", features=top.genes)
```

<div class="figure">
<img src="marker-detection_files/figure-html/viol-de-binom-1.png" alt="Distribution of log-normalized expression values for the top 10 DE genes involving cluster 7 with the binomial test, stratified by cluster assignment and coloured by the plate of origin for each cell." width="672" />
<p class="caption">(\#fig:viol-de-binom)Distribution of log-normalized expression values for the top 10 DE genes involving cluster 7 with the binomial test, stratified by cluster assignment and coloured by the plate of origin for each cell.</p>
</div>

The disadvantage of the binomial test is that its increased stringency can lead to the loss of good candidate markers.
For example, _GCG_ is a known marker for pancreatic alpha cells but is expressed in almost every other cell of the @lawlor2017singlecell pancreas data (Figure \@ref(fig:viol-gcg-lawlor)) and would not be highly ranked by the binomial test.


```r
plotExpression(sce.lawlor, x="cell type", features="ENSG00000115263")
```

<div class="figure">
<img src="marker-detection_files/figure-html/viol-gcg-lawlor-1.png" alt="Distribution of log-normalized expression values for _GCG_ across different pancreatic cell types in the Lawlor pancreas data." width="672" />
<p class="caption">(\#fig:viol-gcg-lawlor)Distribution of log-normalized expression values for _GCG_ across different pancreatic cell types in the Lawlor pancreas data.</p>
</div>

Another property of the binomial test is that it will not respond to scaling normalization.
Systematic differences in library size between clusters will not be considered when computing $p$-values or effect sizes.
This is not necessarily problematic for marker gene detection -
users can treat this as retaining information about the total RNA content, analogous to spike-in normalization in Section \@ref(spike-norm).

### Using custom DE methods

It is also possible to perform marker gene detection based on precomputed DE statistics, which allows us to take advantage of more sophisticated tests in dedicated DE analysis packages in the Bioconductor ecosystem.
To demonstrate, consider the `voom()` approach from the *[limma](https://bioconductor.org/packages/3.12/limma)* package [@law2014voom].
We first process our `SingleCellExperiment` to obtain a `fit` object as shown below.


```r
library(limma)
design <- model.matrix(~0 + label, data=colData(sce.pbmc))
colnames(design)
```

```
##  [1] "label1"  "label2"  "label3"  "label4"  "label5"  "label6"  "label7" 
##  [8] "label8"  "label9"  "label10" "label11" "label12" "label13" "label14"
## [15] "label15" "label16"
```

```r
# Removing very low-abundance genes.
keep <- calculateAverage(sce.pbmc) > 0.1 
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   29465    4229
```

```r
y <- convertTo(sce.pbmc, subset.row=keep)
v <- voom(y, design)
fit <- lmFit(v, design)
```

We then perform pairwise comparisons between clusters using the TREAT strategy [@mccarthy2009treat] to test for log-fold changes that are significantly greater than 0.5.
For each comparison, we store the corresponding data frame of statistics in `all.results`, along with the identities of the clusters involved in `all.pairs`.


```r
nclust <- length(unique(colLabels(sce.pbmc)))
all.results <- all.pairs <- list()
counter <- 1L

# Iterating across the first 'nclust' coefficients in design,
# and comparing them to each other in a pairwise manner.
for (x in seq_len(nclust)) {
    for (y in seq_len(x-1L)) {
        con <- integer(ncol(design))
        con[x] <- 1
        con[y] <- -1
        fit2 <- contrasts.fit(fit, con)
        fit2 <- treat(fit2, robust=TRUE, lfc=0.5)

        res <- topTreat(fit2, n=Inf, sort.by="none")
        all.results[[counter]] <- res
        all.pairs[[counter]] <- colnames(design)[c(x, y)]
        counter <- counter+1L

        # Also filling the reverse comparison.
        res$logFC <- -res$logFC
        all.results[[counter]] <- res
        all.pairs[[counter]] <- colnames(design)[c(y, x)]
        counter <- counter+1L
    }
}
```

These custom results are consolidated into a single marker list for each cluster with the `combineMarkers()` function.
This combines test statistics across all pairwise comparisons involving a single cluster,
yielding a per-cluster `DataFrame` that can be interpreted in the same manner as discussed previously.


```r
all.pairs <- do.call(rbind, all.pairs)
combined <- combineMarkers(all.results, all.pairs, pval.field="P.Value")

# Inspecting results for our cluster of interest again.
interesting.voom <- combined[[paste0("cluster", chosen)]] 
colnames(interesting.voom)
```

```
## NULL
```

```r
head(interesting.voom[,1:4])
```

```
## NULL
```

By default, we do not use custom DE methods to perform marker detection, for several reasons.
Many of these methods rely on empirical Bayes shrinkage to share information across genes in the presence of limited replication. 
However, this is unnecessary when there are large numbers of "replicate" cells in each group (Section \@ref(false-replicates)).
These methods also make stronger assumptions about the data (e.g., equal variances for linear models, the distribution of variances during empirical Bayes) that are more likely to be violated in noisy scRNA-seq contexts.
From a practical perspective, they require more work to set up and take more time to run.
Nonetheless, some custom methods (e.g., *[MAST](https://bioconductor.org/packages/3.12/MAST)*) may provide a useful point of difference from the simpler tests, in which case they can be converted into a marker detection scheme as described above.

### Combining multiple marker statistics

On occasion, we might want to combine marker statistics from several testing regimes into a single `DataFrame`.
This allows us to easily inspect multiple statistics at once to verify that a particular gene is a strong candidate marker.
For example, a large AUC from the WMW test indicates that the expression distributions are well-separated between clusters, while the log-fold change reported with the $t$-test provides a more interpretable measure of the magnitude of the change in expression.
We use the `multiMarkerStats()` to merge the results of separate `findMarkers()` calls into one `DataFrame` per cluster, with statistics interleaved to facilitate a direct comparison between different test regimes.


```r
combined <- multiMarkerStats(
    t=findMarkers(sce.pbmc, direction="up"),
    wilcox=findMarkers(sce.pbmc, test="wilcox", direction="up"),
    binom=findMarkers(sce.pbmc, test="binom", direction="up")
)

# Interleaved marker statistics from both tests for each cluster.
colnames(combined[["1"]])
```

```
##  [1] "Top"                 "p.value"             "FDR"                
##  [4] "t.Top"               "wilcox.Top"          "binom.Top"          
##  [7] "t.p.value"           "wilcox.p.value"      "binom.p.value"      
## [10] "t.FDR"               "wilcox.FDR"          "binom.FDR"          
## [13] "t.summary.logFC"     "wilcox.summary.AUC"  "binom.summary.logFC"
## [16] "t.logFC.2"           "wilcox.AUC.2"        "binom.logFC.2"      
## [19] "t.logFC.3"           "wilcox.AUC.3"        "binom.logFC.3"      
## [22] "t.logFC.4"           "wilcox.AUC.4"        "binom.logFC.4"      
## [25] "t.logFC.5"           "wilcox.AUC.5"        "binom.logFC.5"      
## [28] "t.logFC.6"           "wilcox.AUC.6"        "binom.logFC.6"      
## [31] "t.logFC.7"           "wilcox.AUC.7"        "binom.logFC.7"      
## [34] "t.logFC.8"           "wilcox.AUC.8"        "binom.logFC.8"      
## [37] "t.logFC.9"           "wilcox.AUC.9"        "binom.logFC.9"      
## [40] "t.logFC.10"          "wilcox.AUC.10"       "binom.logFC.10"     
## [43] "t.logFC.11"          "wilcox.AUC.11"       "binom.logFC.11"     
## [46] "t.logFC.12"          "wilcox.AUC.12"       "binom.logFC.12"     
## [49] "t.logFC.13"          "wilcox.AUC.13"       "binom.logFC.13"     
## [52] "t.logFC.14"          "wilcox.AUC.14"       "binom.logFC.14"     
## [55] "t.logFC.15"          "wilcox.AUC.15"       "binom.logFC.15"     
## [58] "t.logFC.16"          "wilcox.AUC.16"       "binom.logFC.16"
```

```r
head(combined[["1"]][,1:9])
```

```
## DataFrame with 6 rows and 9 columns
##              Top     p.value         FDR     t.Top wilcox.Top binom.Top
##        <integer>   <numeric>   <numeric> <integer>  <integer> <integer>
## TYROBP         1 1.36219e-37 1.31136e-34         1          1         1
## FCER1G         2 5.54939e-48 8.90386e-45         1          1         2
## GZMA           2 7.10783e-83 2.39491e-78         1          2         1
## HOPX           2 1.25041e-79 2.10656e-75         2          1         1
## CTSW           3 2.51098e-71 1.20864e-67         1          1         3
## KLRF1          3 5.69193e-66 2.39730e-62         3          1         1
##           t.p.value wilcox.p.value binom.p.value
##           <numeric>      <numeric>     <numeric>
## TYROBP 3.93768e-112   2.99215e-124   1.36219e-37
## FCER1G  4.67496e-82   1.73332e-116   5.54939e-48
## GZMA    1.06381e-88   1.00829e-165   7.10783e-83
## HOPX    1.25041e-79   3.23816e-190  2.40034e-111
## CTSW   1.13373e-107   7.90522e-131   2.51098e-71
## KLRF1   5.69193e-66   2.73030e-184  2.77818e-113
```

In addition, `multiMarkerStats()` will compute a number of new statistics by combining the per-regime statistics.
The combined `Top` value is obtained by simply taking the largest `Top` value across all tests for a given gene, while the reported `p.value` is obtained by taking the largest $p$-value.
Ranking on either metric focuses on genes with robust differences that are highly ranked and detected by each of the individual testing regimes.
Of course, this might be considered an overly conservative approach in practice, so it is entirely permissible to re-rank the `DataFrame` according to the `Top` or `p.value` for an individual regime (effectively limiting the use of the other regimes' statistics to diagnostics only).

## Handling blocking factors {#marker-batch}

### Using the `block=` argument

Large studies may contain factors of variation that are known and not interesting (e.g., batch effects, sex differences).
If these are not modelled, they can interfere with marker gene detection - most obviously by inflating the variance within each cluster, but also by distorting the log-fold changes if the cluster composition varies across levels of the blocking factor.
To avoid these issues, we set the `block=` argument in the `findMarkers()` call, as demonstrated below for the 416B data set.

<button class="aaron-collapse">View history</button>
<div class="aaron-content">
   
```r
#--- loading ---#
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b") 
sce.416b$block <- factor(sce.416b$block)

#--- gene-annotation ---#
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
    keytype="GENEID", column="SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
    keytype="GENEID", column="SEQNAME")

library(scater)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, 
    rowData(sce.416b)$SYMBOL)

#--- quality-control ---#
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
    "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

#--- normalization ---#
library(scran)
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)

#--- variance-modelling ---#
dec.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", block=sce.416b$block)
chosen.hvgs <- getTopHVGs(dec.416b, prop=0.1)

#--- batch-correction ---#
library(limma)
assay(sce.416b, "corrected") <- removeBatchEffect(logcounts(sce.416b), 
    design=model.matrix(~sce.416b$phenotype), batch=sce.416b$block)

#--- dimensionality-reduction ---#
sce.416b <- runPCA(sce.416b, ncomponents=10, subset_row=chosen.hvgs,
    exprs_values="corrected", BSPARAM=BiocSingular::ExactParam())

set.seed(1010)
sce.416b <- runTSNE(sce.416b, dimred="PCA", perplexity=10)

#--- clustering ---#
my.dist <- dist(reducedDim(sce.416b, "PCA"))
my.tree <- hclust(my.dist, method="ward.D2")

library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist),
    minClusterSize=10, verbose=0))
colLabels(sce.416b) <- factor(my.clusters)
```

</div>


```r
m.out <- findMarkers(sce.416b, block=sce.416b$block, direction="up") 
```

For each gene, each pairwise comparison between clusters is performed separately in each level of the blocking factor - in this case, the plate of origin.
The function will then combine $p$-values from different plates using Stouffer's Z method to obtain a single $p$-value per pairwise comparison.
(These $p$-values are further combined across comparisons to obtain a single $p$-value per gene, using either Simes' method or an intersection-union test depending on the value of `pval.type=`.)
This approach favours genes that exhibit consistent DE in the same direction in each plate.


```r
demo <- m.out[["1"]] 
demo[demo$Top <= 5,1:4]
```

```
## DataFrame with 13 rows and 4 columns
##                          Top     p.value         FDR summary.logFC
##                    <integer>   <numeric>   <numeric>     <numeric>
## Foxs1                      1 1.37387e-12 4.35563e-10       3.07058
## Pirb                       1 2.08277e-33 1.21332e-29       5.87820
## Myh11                      1 6.44327e-47 3.00282e-42       4.38182
## Tmsb4x                     2 3.22944e-44 7.52525e-40       1.47689
## Ctsd                       2 6.78109e-38 7.90065e-34       2.89152
## ...                      ...         ...         ...           ...
## Tob1                       4 6.63870e-09 1.18088e-06       2.74161
## Pi16                       4 1.69247e-32 7.88758e-29       5.76914
## Cd53                       5 1.08574e-27 2.97646e-24       5.75200
## Alox5ap                    5 1.33791e-28 4.15679e-25       1.36676
## CBFB-MYH11-mcherry         5 3.75556e-35 3.50049e-31       3.01677
```

The `block=` argument works with all tests shown above and is robust to difference in the log-fold changes or variance between batches.
However, it assumes that each pair of clusters is present in at least one batch.
In scenarios where cells from two clusters never co-occur in the same batch, the comparison will be impossible and `NA`s will be reported in the output.

### Using the `design=` argument

Another approach is to define a design matrix containing the batch of origin as the sole factor.
`findMarkers()` will then fit a linear model to the log-expression values, similar to the use of *[limma](https://bioconductor.org/packages/3.12/limma)* for bulk RNA sequencing data [@ritchie2015limma].
This handles situations where multiple batches contain unique clusters, as comparisons can be implicitly performed via shared cell types in each batch.
There is also a slight increase in power when information is shared across clusters for variance estimation.


```r
# Setting up the design matrix (we remove intercept for full rank
# in the final design matrix with the cluster-specific terms).
design <- model.matrix(~sce.416b$block)
design <- design[,-1,drop=FALSE]

m.alt <- findMarkers(sce.416b, design=design, direction="up")
demo <- m.alt[["1"]]
demo[demo$Top <= 5,1:4]
```

```
## DataFrame with 12 rows and 4 columns
##                          Top     p.value         FDR summary.logFC
##                    <integer>   <numeric>   <numeric>     <numeric>
## Gm6977                     1 7.15187e-24 8.77120e-21      0.810553
## Myh11                      1 4.56882e-64 2.12925e-59      4.381806
## Tmsb4x                     2 9.48997e-46 2.21135e-41      1.478213
## Cd63                       2 1.80446e-15 7.85933e-13      0.813016
## Cd200r3                    2 2.40861e-45 3.74170e-41      6.684003
## ...                      ...         ...         ...           ...
## Actb                       4 5.61751e-36 2.90887e-32      0.961762
## Ctsd                       4 2.08646e-42 2.43094e-38      2.893014
## Fth1                       4 1.83949e-23 2.14319e-20      0.797407
## Ccl9                       5 1.75378e-30 3.71514e-27      5.396347
## CBFB-MYH11-mcherry         5 9.09026e-39 8.47285e-35      3.017758
```

The use of a linear model makes some strong assumptions, necessitating some caution when interpreting the results.
If the batch effect is not consistent across clusters, the variance will be inflated and the log-fold change estimates will be distorted.
Variances are also assumed to be equal across groups, which is not true in general.
In particular, the presence of clusters in which a gene is silent will shrink the residual variance towards zero, preventing the model from penalizing genes with high variance in other clusters.
Thus, we generally recommend the use of `block=` where possible.

## Invalidity of $p$-values {#p-value-invalidity}

### From data snooping

All of our DE strategies for detecting marker genes between clusters are statistically flawed to some extent.
The DE analysis is performed on the same data used to obtain the clusters, which represents "data dredging" (also known as fishing or data snooping).
The hypothesis of interest - are there differences between clusters? - is formulated from the data, so we are more likely to get a positive result when we re-use the data set to test that hypothesis.

The practical effect of data dredging is best illustrated with a simple simulation.
We simulate i.i.d. normal values, perform $k$-means clustering and test for DE between clusters of cells with `findMarkers()`.
The resulting distribution of $p$-values is heavily skewed towards low values (Figure \@ref(fig:pval-dist)).
Thus, we can detect "significant" differences between clusters even in the absence of any real substructure in the data.
This effect arises from the fact that clustering, by definition, yields groups of cells that are separated in expression space.
Testing for DE genes between clusters will inevitably yield some significant results as that is how the clusters were defined.


```r
library(scran)
set.seed(0)
y <- matrix(rnorm(100000), ncol=200)
clusters <- kmeans(t(y), centers=2)$cluster
out <- findMarkers(y, clusters)
hist(out[[1]]$p.value, col="grey80", xlab="p-value")
```

<div class="figure">
<img src="marker-detection_files/figure-html/pval-dist-1.png" alt="Distribution of $p$-values from a DE analysis between two clusters in a simulation with no true subpopulation structure." width="672" />
<p class="caption">(\#fig:pval-dist)Distribution of $p$-values from a DE analysis between two clusters in a simulation with no true subpopulation structure.</p>
</div>

For marker gene detection, this effect is largely harmless as the $p$-values are used only for ranking.
However, it becomes an issue when the $p$-values are used to define "significant differences" between clusters with respect to an error rate threshold.
Meaningful interpretation of error rates require consideration of the long-run behavior, i.e., the rate of incorrect rejections if the experiment were repeated many times.
The concept of statistical significance for differences between clusters is not applicable if clusters and their interpretations are not stably reproducible across (hypothetical) replicate experiments.

### Nature of replication {#false-replicates}

The naive application of DE analysis methods will treat counts from the same cluster of cells as replicate observations.
This is not the most relevant level of replication when cells are derived from the same biological sample (i.e., cell culture, animal or patient).
DE analyses that treat cells as replicates fail to properly model the sample-to-sample variability [@lun2017overcoming].
The latter is arguably the more important level of replication as different samples will necessarily be generated if the experiment is to be replicated.
Indeed, the use of cells as replicates only masks the fact that the sample size is actually one in an experiment involving a single biological sample.
This reinforces the inappropriateness of using the marker gene $p$-values to perform statistical inference.

We strongly recommend selecting some markers for use in validation studies with an independent replicate population of cells.
A typical strategy is to identify a corresponding subset of cells that express the upregulated markers and do not express the downregulated markers.
Ideally, a different technique for quantifying expression would also be used during validation, e.g., fluorescent _in situ_ hybridisation or quantitative PCR.
This confirms that the subpopulation genuinely exists and is not an artifact of the scRNA-seq protocol or the computational analysis.

## Further comments

One consequence of the DE analysis strategy is that markers are defined relative to subpopulations in the same dataset.
Biologically meaningful genes will not be detected if they are expressed uniformly throughout the population, e.g., T cell markers will not be detected if only T cells are present in the dataset.
In practice, this is usually only a problem when the experimental data are provided without any biological context - certainly, we would hope to have some _a priori_ idea about what cells have been captured.
For most applications, it is actually desirable to avoid detecting such genes as we are interested in characterizing heterogeneity  within the context of a known cell population.
Continuing from the example above, the failure to detect T cell markers is of little consequence if we already know we are working with T cells.
Nonetheless, if "absolute" identification of cell types is necessary, we discuss some strategies for doing so in Chapter \@ref(cell-type-annotation).

Alternatively, marker detection can be performed by treating gene expression as a predictor variable for cluster assignment.
For a pair of clusters, we can find genes that discriminate between them by performing inference with a logistic model where the outcome for each cell is whether it was assigned to the first cluster and the lone predictor is the expression of each gene.
Treating the cluster assignment as the dependent variable is more philosophically pleasing in some sense, as the clusters are indeed defined from the expression data rather than being known in advance.
(Note that this does not solve the data snooping problem.)
In practice, this approach effectively does the same task as a Wilcoxon rank sum test in terms of quantifying separation between clusters.
Logistic models have the advantage in that they can easily be extended to block on multiple nuisance variables, though this is not typically necessary in most use cases.
Even more complex strategies use machine learning methods to determine which features contribute most to successful cluster classification, but this is probably unnecessary for routine analyses.

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
 [1] limma_3.45.10               scater_1.17.4              
 [3] ggplot2_3.3.2               pheatmap_1.0.12            
 [5] scran_1.17.15               SingleCellExperiment_1.11.6
 [7] SummarizedExperiment_1.19.6 DelayedArray_0.15.7        
 [9] matrixStats_0.56.0          Matrix_1.2-18              
[11] Biobase_2.49.0              GenomicRanges_1.41.6       
[13] GenomeInfoDb_1.25.10        IRanges_2.23.10            
[15] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[17] BiocStyle_2.17.0            simpleSingleCell_1.13.16   

loaded via a namespace (and not attached):
 [1] viridis_0.5.1             edgeR_3.31.4             
 [3] BiocSingular_1.5.0        viridisLite_0.3.0        
 [5] DelayedMatrixStats_1.11.1 scuttle_0.99.12          
 [7] statmod_1.4.34            BiocManager_1.30.10      
 [9] highr_0.8                 dqrng_0.2.1              
[11] vipor_0.4.5               GenomeInfoDbData_1.2.3   
[13] yaml_2.2.1                pillar_1.4.6             
[15] lattice_0.20-41           glue_1.4.1               
[17] digest_0.6.25             RColorBrewer_1.1-2       
[19] XVector_0.29.3            colorspace_1.4-1         
[21] cowplot_1.0.0             htmltools_0.5.0          
[23] XML_3.99-0.5              pkgconfig_2.0.3          
[25] bookdown_0.20             zlibbioc_1.35.0          
[27] purrr_0.3.4               scales_1.1.1             
[29] processx_3.4.3            BiocParallel_1.23.2      
[31] tibble_3.0.3              farver_2.0.3             
[33] generics_0.0.2            ellipsis_0.3.1           
[35] withr_2.2.0               magrittr_1.5             
[37] crayon_1.3.4              CodeDepends_0.6.5        
[39] evaluate_0.14             ps_1.3.4                 
[41] bluster_0.99.1            beeswarm_0.2.3           
[43] graph_1.67.1              tools_4.0.2              
[45] lifecycle_0.2.0           stringr_1.4.0            
[47] munsell_0.5.0             locfit_1.5-9.4           
[49] irlba_2.3.3               callr_3.4.3              
[51] compiler_4.0.2            rsvd_1.0.3               
[53] rlang_0.4.7               grid_4.0.2               
[55] RCurl_1.98-1.2            BiocNeighbors_1.7.0      
[57] igraph_1.2.5              labeling_0.3             
[59] bitops_1.0-6              rmarkdown_2.3            
[61] gtable_0.3.0              codetools_0.2-16         
[63] R6_2.4.1                  gridExtra_2.3            
[65] knitr_1.29                dplyr_1.0.1              
[67] ggbeeswarm_0.6.0          stringi_1.4.6            
[69] Rcpp_1.0.5                vctrs_0.3.2              
[71] tidyselect_1.1.0          xfun_0.16                
```
</div>
