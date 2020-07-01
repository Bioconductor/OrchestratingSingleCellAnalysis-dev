---
output: html_document
bibliography: ref.bib
---

# Cell cycle assignment

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

On occasion, it can be desirable to determine cell cycle activity from scRNA-seq data.
In and of itself, the distribution of cells across phases of the cell cycle is not usually informative, but we can use this to determine if there are differences in proliferation between subpopulations or across treatment conditions.
Many of the key events in the cell cycle (e.g., passage through checkpoints) are post-translational and thus not directly visible in transcriptomic data; nonetheless, there are enough changes in expression that can be exploited to determine cell cycle phase.
We demonstrate using the 416B dataset, which is known to contain actively cycling cells after oncogene induction.

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
sce.416b
```

```
## class: SingleCellExperiment 
## dim: 46604 185 
## metadata(0):
## assays(3): counts logcounts corrected
## rownames(46604): 4933401J01Rik Gm26206 ... CAAA01147332.1
##   CBFB-MYH11-mcherry
## rowData names(4): Length ENSEMBL SYMBOL SEQNAME
## colnames(185): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
##   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ...
##   SLX-11312.N712_S507.H5H5YBBXX.s_8.r_1
##   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
## colData names(11): Source Name cell line ... sizeFactor label
## reducedDimNames(2): PCA TSNE
## altExpNames(2): ERCC SIRV
```

## Using the cyclins

The cyclins control progression through the cell cycle and have well-characterized patterns of expression across cell cycle phases.
Cyclin D is expressed throughout but peaks at G1; cyclin E is expressed highest in the G1/S transition; cyclin A is expressed across S and G2; and cyclin B is expressed highest in late G2 and mitosis.
Inspection of the relative expression of cyclins across the population can often be sufficient to determine the relative cell cycle activity in each cluster (Figure \@ref(fig:heat-cyclin)).
For example, cluster 1 is likely to be in G1 while the other clusters are scattered across the later phases.


```r
cyclin.genes <- grep("^Ccn[abde][0-9]$", rowData(sce.416b)$SYMBOL)
cyclin.genes <- rownames(sce.416b)[cyclin.genes]
cyclin.genes
```

```
##  [1] "Ccnb3" "Ccna2" "Ccna1" "Ccne2" "Ccnd2" "Ccne1" "Ccnd1" "Ccnb2" "Ccnb1"
## [10] "Ccnd3"
```

```r
library(scater)
plotHeatmap(sce.416b, order_columns_by="label", 
    cluster_rows=FALSE, features=sort(cyclin.genes))
```

<div class="figure">
<img src="cell-cycle_files/figure-html/heat-cyclin-1.png" alt="Heatmap of the log-normalized expression values of the cyclin genes in the 416B dataset. Each column represents a cell that is sorted by the cluster of origin." width="672" />
<p class="caption">(\#fig:heat-cyclin)Heatmap of the log-normalized expression values of the cyclin genes in the 416B dataset. Each column represents a cell that is sorted by the cluster of origin.</p>
</div>



We can use this approach to make statements about the relative cell cycle activity across clusters.
In this case, we apply standard DE methods (Chapter \@ref(marker-detection)) to look for upregulation of each cyclin between clusters, which would imply that a subpopulation contains more cells in the corresponding cell cycle phase.
The same logic applies to comparisons between treatment conditions as described in Chapter \@ref(multi-sample-comparisons).


```r
library(scran)
markers <- findMarkers(sce.416b, subset.row=cyclin.genes, 
    test.type="wilcox", direction="up")

# We can infer that cluster 4 has more cells in G2/M than the other clusters,
# based on higher expression of the cyclin B's.
markers[[4]]
```

```
## DataFrame with 10 rows and 7 columns
##             Top     p.value         FDR summary.AUC     AUC.1     AUC.2
##       <integer>   <numeric>   <numeric>   <numeric> <numeric> <numeric>
## Ccna2         1 4.47082e-09 4.47082e-08    0.996337  0.996337  0.641822
## Ccnd1         1 2.27713e-04 5.69283e-04    0.822981  0.368132  0.822981
## Ccnb1         1 1.19027e-07 5.95137e-07    0.949634  0.949634  0.519669
## Ccnb2         2 3.87799e-07 1.29266e-06    0.934066  0.934066  0.781573
## Ccna1         4 2.96992e-02 5.93985e-02    0.535714  0.535714  0.495342
## Ccne2         5 6.56983e-02 1.09497e-01    0.641941  0.641941  0.447205
## Ccne1         6 5.85979e-01 8.37113e-01    0.564103  0.564103  0.366460
## Ccnd3         7 9.94578e-01 1.00000e+00    0.402930  0.402930  0.283644
## Ccnd2         8 9.99993e-01 1.00000e+00    0.306548  0.134615  0.327122
## Ccnb3        10 1.00000e+00 1.00000e+00    0.500000  0.500000  0.500000
##           AUC.3
##       <numeric>
## Ccna2  0.925595
## Ccnd1  0.776786
## Ccnb1  0.934524
## Ccnb2  0.898810
## Ccna1  0.535714
## Ccne2  0.455357
## Ccne1  0.473214
## Ccnd3  0.273810
## Ccnd2  0.306548
## Ccnb3  0.500000
```



This approach assumes that cyclin expression is not affected by biological processes other than the cell cycle.
This is a strong assumption in highly heterogeneous populations where cyclins may perform cell-type-specific roles.
For example, using the Grun HSC dataset [@grun2016denovo], we see an upregulation of cyclin D2 in sorted HSCs (Figure \@ref(fig:heat-cyclin-grun)) that is consistent with a particular reliance on D-type cyclins in these cells [@steinman2002cell;@kozar2004mouse].


```r
extractCached("grun-hsc.Rmd", "clustering", "sce.grun.hsc")
```

<button class="aaron-collapse">View history</button>
<div class="aaron-content">
   
```r
#--- data-loading ---#
library(scRNAseq)
sce.grun.hsc <- GrunHSCData(ensembl=TRUE)

#--- gene-annotation ---#
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
anno <- select(ens.mm.v97, keys=rownames(sce.grun.hsc), 
    keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
rowData(sce.grun.hsc) <- anno[match(rownames(sce.grun.hsc), anno$GENEID),]

#--- quality-control ---#
library(scuttle)
stats <- perCellQCMetrics(sce.grun.hsc)
qc <- quickPerCellQC(stats, batch=sce.grun.hsc$protocol,
    subset=grepl("sorted", sce.grun.hsc$protocol))
sce.grun.hsc <- sce.grun.hsc[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce.grun.hsc)
sce.grun.hsc <- computeSumFactors(sce.grun.hsc, clusters=clusters)
sce.grun.hsc <- logNormCounts(sce.grun.hsc)

#--- variance-modelling ---#
set.seed(00010101)
dec.grun.hsc <- modelGeneVarByPoisson(sce.grun.hsc) 
top.grun.hsc <- getTopHVGs(dec.grun.hsc, prop=0.1)

#--- dimensionality-reduction ---#
set.seed(101010011)
sce.grun.hsc <- denoisePCA(sce.grun.hsc, technical=dec.grun.hsc, subset.row=top.grun.hsc)
sce.grun.hsc <- runTSNE(sce.grun.hsc, dimred="PCA")

#--- clustering ---#
snn.gr <- buildSNNGraph(sce.grun.hsc, use.dimred="PCA")
colLabels(sce.grun.hsc) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
```

</div>


```r
# Switching the row names for a nicer plot.
rownames(sce.grun.hsc) <- uniquifyFeatureNames(rownames(sce.grun.hsc),
    rowData(sce.grun.hsc)$SYMBOL)

cyclin.genes <- grep("^Ccn[abde][0-9]$", rowData(sce.grun.hsc)$SYMBOL)
cyclin.genes <- rownames(sce.grun.hsc)[cyclin.genes]

plotHeatmap(sce.grun.hsc, order_columns_by="label",
    cluster_rows=FALSE, features=sort(cyclin.genes),
    colour_columns_by="protocol")
```

<div class="figure">
<img src="cell-cycle_files/figure-html/heat-cyclin-grun-1.png" alt="Heatmap of the log-normalized expression values of the cyclin genes in the Grun HSC dataset. Each column represents a cell that is sorted by the cluster of origin and extraction protocol." width="672" />
<p class="caption">(\#fig:heat-cyclin-grun)Heatmap of the log-normalized expression values of the cyclin genes in the Grun HSC dataset. Each column represents a cell that is sorted by the cluster of origin and extraction protocol.</p>
</div>



Admittedly, this is merely a symptom of a more fundamental issue -
that the cell cycle is not independent of the other processes that are occurring in a cell.
This will be a recurring theme throughout the chapter, which suggests that cell cycle inferences are best used in comparisons between closely related cell types where there are fewer changes elsewhere that might interfere with interpretation.

## Using reference profiles

Cell cycle assignment can be considered a specialized case of cell annotation, which suggests that the strategies described in Chapter \@ref(cell-type-annotation) can also be applied here.
Given a reference dataset containing cells of known cell cycle phase, we could use methods like *[SingleR](https://bioconductor.org/packages/3.12/SingleR)* to determine the phase of each cell in a test dataset.
We demonstrate on a reference of mouse ESCs from @buettner2015computational that were sorted by cell cycle phase prior to scRNA-seq.


```r
library(scRNAseq)
sce.ref <- BuettnerESCData()
sce.ref <- logNormCounts(sce.ref)
sce.ref
```

```
## class: SingleCellExperiment 
## dim: 38293 288 
## metadata(0):
## assays(2): counts logcounts
## rownames(38293): ENSMUSG00000000001 ENSMUSG00000000003 ...
##   ENSMUSG00000097934 ENSMUSG00000097935
## rowData names(3): EnsemblTranscriptID AssociatedGeneName GeneLength
## colnames(288): G1_cell1_count G1_cell2_count ... G2M_cell95_count
##   G2M_cell96_count
## colData names(2): phase sizeFactor
## reducedDimNames(0):
## altExpNames(1): ERCC
```

We will restrict the annotation process to a subset of genes with _a priori_ known roles in cell cycle.
This aims to avoid detecting markers for other biological processes that happen to be correlated with the cell cycle in the reference dataset, which would reduce classification performance if those processes are absent or uncorrelated in the test dataset.


```r
# Find genes that are cell cycle-related.
library(org.Mm.eg.db)
cycle.anno <- select(org.Mm.eg.db, keytype="GOALL", keys="GO:0007049", 
    columns="ENSEMBL")[,"ENSEMBL"]

# Find the genes that are present in both datasets as well.
candidates <- intersect(cycle.anno, rownames(sce.ref))
candidates <- intersect(candidates, rowData(sce.416b)$ENSEMBL)
str(candidates)
```

```
##  chr [1:1606] "ENSMUSG00000026842" "ENSMUSG00000029580" ...
```

We use the `SingleR()` function to assign labels to the 416B data based on the cell cycle phases in the ESC reference.
Cluster 1 mostly consists of G1 cells while the other clusters have more cells in the other phases, which is broadly consistent with our conclusions from the cyclin-based analysis.
Unlike the cyclin-based analysis, this approach yields "absolute" assignments of cell cycle phase that do not need to be interpreted relative to other cells in the same dataset.


```r
# Switching row names back to Ensembl to match the reference.
test.data <- logcounts(sce.416b)
rownames(test.data) <- rowData(sce.416b)$ENSEMBL

library(SingleR)
assignments <- SingleR(test.data[candidates,], ref=sce.ref[candidates,],
    de.method="wilcox", label=sce.ref$phase)

tab <- table(assignments$labels, colLabels(sce.416b))
tab
```

```
##      
##        1  2  3  4
##   G1  71  7 19  1
##   G2M  2 60  1 13
##   S    5  2  4  0
```



The key assumption here is that the cell cycle is orthogonal to cell type and other aspects of cell behavior.
This justifies the use of a reference involving cell types that are quite different from the cells in the test dataset, provided that the cell cycle transcriptional program is conserved across datasets [@bertoli2013control;@conboy2007cell].
However, it is not difficult to find routine violations of this assumption - for example, _Lef1_ is detected as one of the top markers to distinguish between G1 from G2/M in the reference but has no detectable expression in the 416B dataset (Figure \@ref(fig:dist-lef1)).


```r
gridExtra::grid.arrange(
    plotExpression(sce.ref, features="ENSMUSG00000027985", x="phase"),
    plotExpression(sce.416b, features="Lef1", x="label"),
    ncol=2)
```

<div class="figure">
<img src="cell-cycle_files/figure-html/dist-lef1-1.png" alt="Distribution of log-normalized expression values for _Lef1_ in the reference dataset (left) and in the 416B dataset (right)." width="672" />
<p class="caption">(\#fig:dist-lef1)Distribution of log-normalized expression values for _Lef1_ in the reference dataset (left) and in the 416B dataset (right).</p>
</div>



Thus, a healthy dose of skepticism is required when interpreting these assignments.
Our hope is that any systematic assignment error is consistent across clusters and conditions such that they cancel out in comparisons of phase frequencies, which is the more interesting analysis anyway. 
Indeed, while the availability of absolute phase calls may be more appealing, it may not make much practical difference to the conclusions if the frequencies are ultimately interpreted in a relative sense (e.g., using a chi-squared test). 


```r
# Test for differences in phase distributions between clusters 1 and 2.
chisq.test(tab[,1:2])
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  tab[, 1:2]
## X-squared = 108, df = 2, p-value <2e-16
```

## Using the `cyclone()` classifier

The method described by @scialdone2015computational is yet another approach for classifying cells into cell cycle phases.
Using a reference dataset, we first compute the sign of the difference in expression between each pair of genes.
Pairs with changes in the sign across cell cycle phases are chosen as markers.
Cells in a test dataset can then be classified into the appropriate phase, based on whether the observed sign for each marker pair is consistent with one phase or another.
This approach is implemented in the `cyclone()` function from the *[scran](https://bioconductor.org/packages/3.12/scran)* package, which also contains pre-trained set of marker pairs for mouse and human data.


```r
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
    package="scran"))

# Using Ensembl IDs to match up with the annotation in 'mm.pairs'.
assignments <- cyclone(sce.416b, mm.pairs, gene.names=rowData(sce.416b)$ENSEMBL)
```

The phase assignment result for each cell in the 416B dataset is shown in Figure \@ref(fig:phaseplot416b).
For each cell, a higher score for a phase corresponds to a higher probability that the cell is in that phase.
We focus on the G1 and G2/M scores as these are the most informative for classification.


```r
plot(assignments$score$G1, assignments$score$G2M,
    xlab="G1 score", ylab="G2/M score", pch=16)
```

<div class="figure">
<img src="cell-cycle_files/figure-html/phaseplot416b-1.png" alt="Cell cycle phase scores from applying the pair-based classifier on the 416B dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases." width="672" />
<p class="caption">(\#fig:phaseplot416b)Cell cycle phase scores from applying the pair-based classifier on the 416B dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases.</p>
</div>

Cells are classified as being in G1 phase if the G1 score is above 0.5 and greater than the G2/M score;
    in G2/M phase if the G2/M score is above 0.5 and greater than the G1 score;
    and in S phase if neither score is above 0.5.
We see that the results are quite similar to those from `SingleR()`, which is reassuring.


```r
table(assignments$phases, colLabels(sce.416b))
```

```
##      
##        1  2  3  4
##   G1  74  8 20  0
##   G2M  1 48  0 13
##   S    3 13  4  1
```



The same considerations and caveats described for the *[SingleR](https://bioconductor.org/packages/3.12/SingleR)*-based approach are also applicable here.
From a practical perspective, `cyclone()` takes much longer but does not require an explicit reference as the marker pairs are already computed.

## Regressing out cell cycle phase

For some time, it was popular to regress out the cell cycle phase prior to downstream analyses.
The aim was to remove uninteresting variation due to cell cycle, thus improving resolution of other biological processes of interest.
We could implement this by performing cell cycle phase assignment as described above, treating each phase as a separate batch and applying any of the batch correction strategies described in Chapter \@ref(data-integration).
The most common approach is to use a linear model to simply regress out the phase effect, e.g., via `regressBatches()`.


```r
library(batchelor)
sce.nocycle <- regressBatches(sce.416b, batch=assignments$phases)

# The corrected matrix can then be used for downstream analyses:
sce.nocycle <- runPCA(sce.nocycle, exprs_values="corrected")
```

Similarly, for functions that support blocking, we can use the phase assignments as a blocking factor.


```r
# Similar use in related functions that support blocking:
dec.nocycle <- modelGeneVarWithSpikes(sce.416b, "ERCC", 
    block=assignments$phases)
marker.nocycle <- findMarkers(sce.416b, block=assignments$phases)
```

That said, we do not consider cell cycle adjustment to be necessary for routine scRNA-seq analyses.

- In most applications, the cell cycle is a minor factor of variation, secondary to differences between cell types.
It will often have no effect on many analyses that focus on broader aspects of heterogeneity.
More subtle heterogeneity may be masked by cell cycle variation but this should be demonstrated rather than assumed by default.
- Any attempt at removal assumes that the cell cycle effect is orthogonal to other biological processes.
Regression will remove interesting signal if cell cycle activity varies across clusters or conditions. 
This is not an uncommon occurence with, e.g., increased proliferation of T cells upon activation [@richard2018tcell] and changes in cell cycle phase progression across developmental stages [@roccio2013predicting].
Violations of this assumption may also introduce spurious signal within clusters, interfering with any interpretation of subtle variation.
- If adjustment is truly necessary, it should be applied separately to the subset of cells in each cluster.
This avoids the worst violations of the orthogonality assumption due to differences in cell cycle behavior across clusters.
Similarly, gene-based analyses should use the uncorrected data with blocking where possible (Section \@ref(using-corrected-values)), which provides a sanity check that protects against distortions introduced by the adjustment.

It can also be an informative exercise to repeat the analysis after removing all known cell cycle-related genes.
This allows us to explore other factors of variation that are correlated with but distinct from the cell cycle, such as cell fate decisions [@soufi2016cycling] that would otherwise have been eliminated by regression.
We demonstrate below with the @leng2015oscope dataset containing phase-sorted ESCs, where removal of a variety of cell cycle-related genes does not eliminate the separation between G1 and S populations (Figure \@ref(fig:leng-nocycle)).
The persistence of this separation is driven by differential expression in genes without any direct role in cell cyle progression, possibly indicative of a correlated biological process (or a sorting artifact).


```r
library(org.Hs.eg.db)
go.genes <- select(org.Hs.eg.db, keys="GO:0007049", # cell cycle
    keytype="GOALL", column="ENSEMBL")[,"ENSEMBL"]

library(reactome.db)
rct.genes <- select(reactome.db, keys="R-HSA-1640170", # cell cycle
    keytype="PATHID", column="ENTREZID")[,"ENTREZID"]
rct.genes <- select(org.Hs.eg.db, keys=as.character(rct.genes), 
    keytype="ENTREZID", column="ENSEMBL")[,"ENSEMBL"]

combined <- union(rct.genes, go.genes)
length(combined)
```

```
## [1] 2244
```

```r
# Performing an analysis without the cell cycle-related genes.
library(scRNAseq)
sce.leng <- LengESCData(ensembl=TRUE)
leftovers <- setdiff(rownames(sce.leng), combined)
sce.nocycle <- sce.leng[leftovers,]

sce.nocycle <- logNormCounts(sce.nocycle, assay.type="normcounts")
dec.nocycle <- modelGeneVar(sce.nocycle)
sce.nocycle <- runPCA(sce.nocycle, subset_row=getTopHVGs(dec.nocycle, n=1000))
plotPCA(sce.nocycle, colour_by="Phase")
```

<div class="figure">
<img src="cell-cycle_files/figure-html/leng-nocycle-1.png" alt="PCA plot of the Leng ESC dataset, generated after comprehensive removal of cell cycle-related genes. Each point corresponds to a cell that is colored by the sorted cell cycle phase." width="672" />
<p class="caption">(\#fig:leng-nocycle)PCA plot of the Leng ESC dataset, generated after comprehensive removal of cell cycle-related genes. Each point corresponds to a cell that is colored by the sorted cell cycle phase.</p>
</div>

```r
diff <- findMarkers(sce.nocycle, sce.nocycle$Phase, direction="up", 
    row.data=rowData(sce.nocycle)[,"originalName",drop=FALSE]) 
as.data.frame(diff$S[1:20,])
```

```
##                 originalName Top   p.value       FDR summary.logFC logFC.G1
## ENSG00000168298     HIST1H1E   1 1.151e-40 1.808e-36        3.9911   3.9911
## ENSG00000184357     HIST1H1B   2 6.582e-38 5.168e-34        4.1561   4.1561
## ENSG00000172006       ZNF554   3 1.083e-28 5.670e-25        1.3533   1.3533
## ENSG00000105426        PTPRS   4 3.425e-27 1.345e-23        1.3334   1.3334
## ENSG00000124575     HIST1H1D   5 6.914e-24 2.172e-20        2.9226   2.9226
## ENSG00000124610     HIST1H1A   6 4.335e-22 1.135e-18        3.8796   3.8796
## ENSG00000122787       AKR1D1   7 1.944e-20 4.362e-17        1.7696   1.7696
## ENSG00000187837     HIST1H1C   8 2.367e-20 4.645e-17        2.6925   2.6925
## ENSG00000112599       GUCA1B   9 5.005e-20 8.733e-17        1.6067   1.6067
## ENSG00000196912     ANKRD36B  10 9.344e-20 1.467e-16        1.9969   1.9969
## ENSG00000140505       CYP1A2  11 5.226e-18 7.460e-15        1.0224   1.0224
## ENSG00000244694       PTCHD4  12 6.086e-18 7.964e-15        2.2297   2.2297
## ENSG00000135976      ANKRD36  13 1.822e-17 2.175e-14        3.3136   3.3136
## ENSG00000105392          CRX  14 1.939e-17 2.175e-14        1.1622   1.1622
## ENSG00000134757         DSG3  15 2.701e-17 2.655e-14        1.4882   1.4882
## ENSG00000113946       CLDN16  16 2.705e-17 2.655e-14        1.3234   1.3234
## ENSG00000147697        GSDMC  17 3.850e-17 3.455e-14        1.2055   1.2055
## ENSG00000142408       CACNG8  18 3.961e-17 3.455e-14        1.2537   1.2537
## ENSG00000180616        SSTR2  19 5.355e-17 4.425e-14        1.6697   1.6697
## ENSG00000175544        CABP4  20 1.297e-16 1.019e-13        0.9409   0.9409
```

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
 [1] reactome.db_1.70.0          org.Hs.eg.db_3.11.4        
 [3] batchelor_1.5.1             SingleR_1.3.6              
 [5] org.Mm.eg.db_3.11.4         ensembldb_2.13.1           
 [7] AnnotationFilter_1.13.0     GenomicFeatures_1.41.0     
 [9] AnnotationDbi_1.51.1        scRNAseq_2.3.8             
[11] scran_1.17.3                scater_1.17.2              
[13] ggplot2_3.3.2               SingleCellExperiment_1.11.6
[15] SummarizedExperiment_1.19.5 DelayedArray_0.15.6        
[17] matrixStats_0.56.0          Matrix_1.2-18              
[19] Biobase_2.49.0              GenomicRanges_1.41.5       
[21] GenomeInfoDb_1.25.5         IRanges_2.23.10            
[23] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[25] BiocStyle_2.17.0            simpleSingleCell_1.13.5    

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0              colorspace_1.4-1             
  [3] ellipsis_0.3.1                scuttle_0.99.10              
  [5] XVector_0.29.3                BiocNeighbors_1.7.0          
  [7] farver_2.0.3                  bit64_0.9-7                  
  [9] interactiveDisplayBase_1.27.5 codetools_0.2-16             
 [11] knitr_1.29                    Rsamtools_2.5.3              
 [13] dbplyr_1.4.4                  pheatmap_1.0.12              
 [15] graph_1.67.1                  shiny_1.5.0                  
 [17] BiocManager_1.30.10           compiler_4.0.2               
 [19] httr_1.4.1                    dqrng_0.2.1                  
 [21] lazyeval_0.2.2                assertthat_0.2.1             
 [23] fastmap_1.0.1                 limma_3.45.7                 
 [25] later_1.1.0.1                 BiocSingular_1.5.0           
 [27] prettyunits_1.1.1             htmltools_0.5.0              
 [29] tools_4.0.2                   rsvd_1.0.3                   
 [31] igraph_1.2.5                  gtable_0.3.0                 
 [33] glue_1.4.1                    GenomeInfoDbData_1.2.3       
 [35] dplyr_1.0.0                   rappdirs_0.3.1               
 [37] Rcpp_1.0.4.6                  vctrs_0.3.1                  
 [39] Biostrings_2.57.2             rtracklayer_1.49.3           
 [41] ExperimentHub_1.15.0          DelayedMatrixStats_1.11.1    
 [43] xfun_0.15                     stringr_1.4.0                
 [45] ps_1.3.3                      mime_0.9                     
 [47] lifecycle_0.2.0               irlba_2.3.3                  
 [49] statmod_1.4.34                XML_3.99-0.3                 
 [51] AnnotationHub_2.21.1          edgeR_3.31.4                 
 [53] zlibbioc_1.35.0               scales_1.1.1                 
 [55] ProtGenerics_1.21.0           hms_0.5.3                    
 [57] promises_1.1.1                RColorBrewer_1.1-2           
 [59] yaml_2.2.1                    curl_4.3                     
 [61] memoise_1.1.0                 gridExtra_2.3                
 [63] biomaRt_2.45.1                stringi_1.4.6                
 [65] RSQLite_2.2.0                 BiocVersion_3.12.0           
 [67] highr_0.8                     BiocParallel_1.23.0          
 [69] rlang_0.4.6                   pkgconfig_2.0.3              
 [71] bitops_1.0-6                  evaluate_0.14                
 [73] lattice_0.20-41               purrr_0.3.4                  
 [75] labeling_0.3                  GenomicAlignments_1.25.3     
 [77] CodeDepends_0.6.5             cowplot_1.0.0                
 [79] bit_1.1-15.2                  processx_3.4.2               
 [81] tidyselect_1.1.0              magrittr_1.5                 
 [83] bookdown_0.20                 R6_2.4.1                     
 [85] generics_0.0.2                DBI_1.1.0                    
 [87] pillar_1.4.4                  withr_2.2.0                  
 [89] RCurl_1.98-1.2                tibble_3.0.1                 
 [91] crayon_1.3.4                  BiocFileCache_1.13.0         
 [93] rmarkdown_2.3                 progress_1.2.2               
 [95] viridis_0.5.1                 locfit_1.5-9.4               
 [97] grid_4.0.2                    blob_1.2.1                   
 [99] callr_3.4.3                   digest_0.6.25                
[101] xtable_1.8-4                  httpuv_1.5.4                 
[103] openssl_1.4.2                 munsell_0.5.0                
[105] beeswarm_0.2.3                viridisLite_0.3.0            
[107] vipor_0.4.5                   askpass_1.1                  
```
</div>
