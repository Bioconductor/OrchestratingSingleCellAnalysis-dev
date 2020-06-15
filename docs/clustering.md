---
bibliography: ref.bib
---

# Clustering

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

Clustering is an unsupervised learning procedure that is used in scRNA-seq data analysis to empirically define groups of cells with similar expression profiles.
Its primary purpose is to summarize the data in a digestible format for human interpretation. 
This allows us to describe population heterogeneity in terms of discrete labels that are easily understood, rather than attempting to comprehend the high-dimensional manifold on which the cells truly reside.
After annotation based on marker genes, the clusters can be treated as proxies for more abstract biological concepts such as cell types or states.
Clustering is thus a critical step for extracting biological insights from scRNA-seq data.
Here, we demonstrate the application of several commonly used methods with the 10X PBMC dataset. 

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
```

</div>


```r
sce.pbmc
```

```
## class: SingleCellExperiment 
## dim: 33694 3922 
## metadata(1): Samples
## assays(2): counts logcounts
## rownames(33694): RP11-34P13.3 FAM138A ... AC213203.1 FAM231B
## rowData names(2): ID Symbol
## colnames(3922): AAACCTGAGAAGGCCT-1 AAACCTGAGACAGACC-1 ...
##   TTTGTCACAGGTCCAC-1 TTTGTCATCCCAAGAT-1
## colData names(3): Sample Barcode sizeFactor
## reducedDimNames(3): PCA TSNE UMAP
## altExpNames(0):
```

## What is the "true clustering"?

At this point, it is worth stressing the distinction between clusters and cell types.
The former is an empirical construct while the latter is a biological truth (albeit a vaguely defined one).
For this reason, questions like "what is the true number of clusters?" are usually meaningless.
We can define as many clusters as we like, with whatever algorithm we like - each clustering will represent its own partitioning of the high-dimensional expression space, and is as "real" as any other clustering.

A more relevant question is "how well do the clusters approximate the cell types?"
Unfortunately, this is difficult to answer given the context-dependent interpretation of biological truth.
Some analysts will be satisfied with resolution of the major cell types; other analysts may want resolution of subtypes; and others still may require resolution of different states (e.g., metabolic activity, stress) within those subtypes.
Moreover, two clusterings can be highly inconsistent yet both valid, simply partitioning the cells based on different aspects of biology.
Indeed, asking for an unqualified "best" clustering is akin to asking for the best magnification on a microscope without any context.

It is helpful to realize that clustering, like a microscope, is simply a tool to explore the data.
We can zoom in and out by changing the resolution of the clustering parameters, and we can experiment with different clustering algorithms to obtain alternative perspectives of the data.
This iterative approach is entirely permissible for data exploration, which constitutes the majority of all scRNA-seq data analyses.
   
## Graph-based clustering {#clustering-graph}

### Background

Popularized by its use in *[Seurat](https://CRAN.R-project.org/package=Seurat)*, graph-based clustering is a flexible and scalable technique for clustering large scRNA-seq datasets.
We first build a graph where each node is a cell that is connected to its nearest neighbors in the high-dimensional space.
Edges are weighted based on the similarity between the cells involved, with higher weight given to cells that are more closely related.
We then apply algorithms to identify "communities" of cells that are more connected to cells in the same community than they are to cells of different communities.
Each community represents a cluster that we can use for downstream interpretation.

The major advantage of graph-based clustering lies in its scalability.
It only requires a $k$-nearest neighbor search that can be done in log-linear time on average, in contrast to hierachical clustering methods with runtimes that are quadratic with respect to the number of cells.
Graph construction avoids making strong assumptions about the shape of the clusters or the distribution of cells within each cluster, compared to other methods like $k$-means (that favor spherical clusters) or Gaussian mixture models (that require normality).
From a practical perspective, each cell is forcibly connected to a minimum number of neighboring cells, which reduces the risk of generating many uninformative clusters consisting of one or two outlier cells.

The main drawback of graph-based methods is that, after graph construction, no information is retained about relationships beyond the neighboring cells^[Sten Linarrsson talked about this in SCG2018, but I don't know where that work ended up. So this is what passes as a reference for the time being.].
This has some practical consequences in datasets that exhibit differences in cell density, as more steps through the graph are required to move the same distance through a region of higher cell density.
From the perspective of community detection algorithms, this effect "inflates" the high-density regions such that any internal substructure or noise is more likely to cause formation of subclusters.
The resolution of clustering thus becomes dependent on the density of cells, which can occasionally be misleading if it overstates the heterogeneity in the data.



### Implementation

There are several considerations in the practical execution of a graph-based clustering method:

- How many neighbors are considered when constructing the graph.
- What scheme is used to weight the edges.
- Which community detection algorithm is used to define the clusters.

For example, the following code uses the 10 nearest neighbors of each cell to construct a shared nearest neighbor graph.
Two cells are connected by an edge if any of their nearest neighbors are shared,
with the edge weight defined from the highest average rank of the shared neighbors [@xu2015identification].
The Walktrap method from the *[igraph](https://CRAN.R-project.org/package=igraph)* package is then used to identify communities.
All calculations are performed using the top PCs to take advantage of data compression and denoising.


```r
library(scran)
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
```

```
## clust
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
## 585 518 364 458 170 791 295 107  45  46 152  84  40  60 142  16  28  21
```

We assign the cluster assignments back into our `SingleCellExperiment` object as a factor in the column metadata.
This allows us to conveniently visualize the distribution of clusters in a $t$-SNE plot (Figure \@ref(fig:tsne-clust-graph)).


```r
library(scater)
colLabels(sce.pbmc) <- factor(clust)
plotReducedDim(sce.pbmc, "TSNE", colour_by="label")
```

<div class="figure">
<img src="clustering_files/figure-html/tsne-clust-graph-1.png" alt="$t$-SNE plot of the 10X PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from graph-based clustering." width="672" />
<p class="caption">(\#fig:tsne-clust-graph)$t$-SNE plot of the 10X PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from graph-based clustering.</p>
</div>

One of the most important parameters is `k`, the number of nearest neighbors used to construct the graph.
This controls the resolution of the clustering where higher `k` yields a more inter-connected graph and broader clusters.
Users can exploit this by experimenting with different values of `k` to obtain a satisfactory resolution.


```r
# More resolved.
g.5 <- buildSNNGraph(sce.pbmc, k=5, use.dimred = 'PCA')
clust.5 <- igraph::cluster_walktrap(g.5)$membership
table(clust.5)
```

```
## clust.5
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
##  81  45 457 296 168 350  79 432 428 897  64 171  68 135  79  26  18  30  21  16 
##  21  22  23 
##  36   9  16
```

```r
# Less resolved.
g.50 <- buildSNNGraph(sce.pbmc, k=50, use.dimred = 'PCA')
clust.50 <- igraph::cluster_walktrap(g.50)$membership
table(clust.50)
```

```
## clust.50
##   1   2   3   4   5   6   7   8 
## 307 729 789 187 516 524 825  45
```



The graph itself can be visualized using a force-directed layout (Figure \@ref(fig:pbmc-force)).
This yields a dimensionality reduction result that is closely related to $t$-SNE and UMAP, though which of these is the most aesthetically pleasing is left to the eye of the beholder.


```r
set.seed(2000)
reducedDim(sce.pbmc, "force") <- igraph::layout_with_fr(g)
plotReducedDim(sce.pbmc, colour_by="label", dimred="force")
```

<div class="figure">
<img src="clustering_files/figure-html/pbmc-force-1.png" alt="Force-directed layout for the shared nearest-neighbor graph of the PBMC dataset. Each point represents a cell and is coloured according to its assigned cluster identity." width="672" />
<p class="caption">(\#fig:pbmc-force)Force-directed layout for the shared nearest-neighbor graph of the PBMC dataset. Each point represents a cell and is coloured according to its assigned cluster identity.</p>
</div>

### Other parameters

Further tweaking can be performed by changing the edge weighting scheme during graph construction.
Setting `type="number"` will weight edges based on the number of nearest neighbors that are shared between two cells.
Similarly, `type="jaccard"` will weight edges according to the Jaccard index of the two sets of neighbors.
We can also disable weighting altogether by using `buildKNNGraph()`, which is occasionally useful for downstream graph operations that do not support weights.


```r
g.num <- buildSNNGraph(sce.pbmc, use.dimred="PCA", type="number")
g.jaccard <- buildSNNGraph(sce.pbmc, use.dimred="PCA", type="jaccard")
g.none <- buildKNNGraph(sce.pbmc, use.dimred="PCA")
```

All of these `g` variables are `graph` objects from the *[igraph](https://CRAN.R-project.org/package=igraph)* package and can be used with any of the community detection algorithms provided by *[igraph](https://CRAN.R-project.org/package=igraph)*.
We have already mentioned the Walktrap approach, but many others are available to choose from:


```r
clust.louvain <- igraph::cluster_louvain(g)$membership
clust.infomap <- igraph::cluster_infomap(g)$membership
clust.fast <- igraph::cluster_fast_greedy(g)$membership
clust.labprop <- igraph::cluster_label_prop(g)$membership
clust.eigen <- igraph::cluster_leading_eigen(g)$membership
```

It is then straightforward to compare two clustering strategies to see how they differ.
For example, Figure \@ref(fig:walktrap-v-others) suggests that Infomap yields finer clusters than Walktrap while fast-greedy yields coarser clusters.


```r
library(pheatmap)

# Using a large pseudo-count for a smoother color transition
# between 0 and 1 cell in each 'tab'.
tab <- table(paste("Infomap", clust.infomap), 
    paste("Walktrap", clust))
ivw <- pheatmap(log10(tab+10), main="Infomap vs Walktrap",
    color=viridis::viridis(100), silent=TRUE)

tab <- table(paste("Fast", clust.fast), 
    paste("Walktrap", clust))
fvw <- pheatmap(log10(tab+10), main="Fast-greedy vs Walktrap",
    color=viridis::viridis(100), silent=TRUE)

gridExtra::grid.arrange(ivw[[4]], fvw[[4]])
```

<div class="figure">
<img src="clustering_files/figure-html/walktrap-v-others-1.png" alt="Number of cells assigned to combinations of cluster labels with different community detection algorithms in the PBMC dataset. Each entry of each heatmap represents a pair of labels, coloured proportionally to the log-number of cells with those labels." width="672" />
<p class="caption">(\#fig:walktrap-v-others)Number of cells assigned to combinations of cluster labels with different community detection algorithms in the PBMC dataset. Each entry of each heatmap represents a pair of labels, coloured proportionally to the log-number of cells with those labels.</p>
</div>



Some community detection algorithms operate by agglomeration and thus can be used to construct a hierarchical dendrogram based on the pattern of merges between clusters.
The dendrogram itself is not particularly informative as it simply describes the order of merge steps performed by the algorithm; unlike the dendrograms produced by hierarchical clustering (Section \@ref(hierarchical-clustering)), it does not capture the magnitude of differences between subpopulations.
However, it does provide a convenient avenue for manually tuning the clustering resolution by generating nested clusterings using the `cut_at()` function, as shown below.


```r
community.walktrap <- igraph::cluster_walktrap(g)
table(igraph::cut_at(community.walktrap, n=5))
```

```
## 
##    1    2    3    4    5 
## 3612  198   45   46   21
```

```r
table(igraph::cut_at(community.walktrap, n=20))
```

```
## 
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
## 533 364 458 442 170 791 295 107  45  46 152  84  40  60 142  76  52  16  28  21
```

Pipelines involving *[scran](https://bioconductor.org/packages/3.12/scran)* default to rank-based weights followed by Walktrap clustering.
In contrast, *[Seurat](https://CRAN.R-project.org/package=Seurat)* uses Jaccard-based weights followed by Louvain clustering.
Both of these strategies work well, and it is likely that the same could be said for many other combinations of weighting schemes and community detection algorithms.

### Assessing cluster separation

When dealing with graphs, the modularity is a natural metric for evaluating the separation between communities/clusters.
This is defined as the (scaled) difference between the observed total weight of edges between nodes in the same cluster and the expected total weight if edge weights were randomly distributed across all pairs of nodes.
Larger modularity values indicate that there most edges occur within clusters, suggesting that the clusters are sufficiently well separated to avoid edges forming between neighboring cells in different clusters.

The standard approach is to report a single modularity value for a clustering on a given graph.
This is useful for comparing different clusterings on the same graph - and indeed, some community detection algorithms are designed with the aim of maximizing the modularity - but it is less helpful for interpreting a given clustering.
Rather, we use the `clusterModularity()` function with `as.ratio=TRUE`, which returns the ratio of the observed to expected sum of weights between each pair of clusters.
We use the ratio instead of the difference as the former is less dependent on the number of cells in each cluster.


```r
ratio <- clusterModularity(g, clust, as.ratio=TRUE)
dim(ratio)
```

```
## [1] 18 18
```

In this matrix, each row/column corresponds to a cluster and each entry contains the ratio of the observed to total weight of edges between cells in the respective clusters.
A dataset containing well-separated clusters should contain most of the observed total weight on the diagonal entries, i.e., most edges occur between cells in the same cluster.
Indeed, concentration of the weight on the diagonal of (Figure \@ref(fig:cluster-mod)) indicates that most of the clusters are well-separated, while some modest off-diagonal entries represent closely related clusters with more inter-connecting edges.


```r
library(pheatmap)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100))
```

<div class="figure">
<img src="clustering_files/figure-html/cluster-mod-1.png" alt="Heatmap of the log~2~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links." width="672" />
<p class="caption">(\#fig:cluster-mod)Heatmap of the log~2~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links.</p>
</div>

One useful approach is to use the `ratio` matrix to form another graph where the nodes are clusters rather than cells.
Edges between nodes are weighted according to the ratio of observed to expected edge weights between cells in those clusters.
We can then repeat our graph operations on this new cluster-level graph to explore the relationships between clusters.
For example, we could obtain clusters of clusters, or we could simply create a new cluster-based layout for visualization (Figure \@ref(fig:cluster-graph)).
This is analogous to the "graph abstraction" approach described by @wolf2017graph, which can be used to identify trajectories in the data based on high-weight paths between clusters.


```r
cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1), 
    mode="upper", weighted=TRUE, diag=FALSE)

# Increasing the weight to increase the visibility of the lines.
set.seed(11001010)
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*5,
    layout=igraph::layout_with_lgl)
```

<div class="figure">
<img src="clustering_files/figure-html/cluster-graph-1.png" alt="Force-based layout showing the relationships between clusters based on the log-ratio of observed to expected total weights between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding log-ratio." width="672" />
<p class="caption">(\#fig:cluster-graph)Force-based layout showing the relationships between clusters based on the log-ratio of observed to expected total weights between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding log-ratio.</p>
</div>

Incidentally, some readers may have noticed that all *[igraph](https://CRAN.R-project.org/package=igraph)* commands were prefixed with `igraph::`.
We have done this deliberately to avoid bringing `igraph::normalize` into the global namespace.
Rather unfortunately, this `normalize` function accepts any argument and returns `NULL`, which causes difficult-to-diagnose bugs when it overwrites `normalize` from *[BiocGenerics](https://bioconductor.org/packages/3.12/BiocGenerics)*.

## $k$-means clustering 

### Background

$k$-means clustering is a classic technique that aims to partition cells into $k$ clusters.
Each cell is assigned to the cluster with the closest centroid, which is done by minimizing the within-cluster sum of squares using a random starting configuration for the $k$ centroids.
The main advantage of this approach lies in its speed, given the simplicity and ease of implementation of the algorithm.
However, it suffers from a number of serious shortcomings that reduce its appeal for obtaining interpretable clusters:

- It implicitly favors spherical clusters of equal radius.
This can lead to unintuitive partitionings on real datasets that contain groupings with irregular sizes and shapes.
- The number of clusters $k$ must be specified beforehand and represents a hard cap on the resolution of the clustering..
For example, setting $k$ to be below the number of cell types will always lead to co-clustering of two cell types, regardless of how well separated they are.
In contrast, other methods like graph-based clustering will respect strong separation even if the relevant resolution parameter is set to a low value.
- It is dependent on the randomly chosen initial coordinates.
This requires multiple runs to verify that the clustering is stable.

That said, $k$-means clustering is still one of the best approaches for sample-based data compression. 
In this application, we set $k$ to a large value such as the square root of the number of cells to obtain fine-grained clusters.
These are not meant to be interpreted directly, but rather, the centroids are treated as "samples" for further analyses.
The idea here is to obtain a single representative of each region of the expression space, reducing the number of samples and computational work in later steps like, e.g., trajectory reconstruction [@ji2016tscan].
This approach will also eliminate differences in cell density across the expression space, ensuring that the most abundant cell type does not dominate downstream results. 

### Base implementation 

Base R provides the `kmeans()` function that does as its name suggests.
We call this on our top PCs to obtain a clustering for a specified number of clusters in the `centers=` argument, after setting the random seed to ensure that the results are reproducible.
In general, the $k$-means clusters correspond to the visual clusters on the $t$-SNE plot in Figure \@ref(fig:tsne-clust-kmeans), though there are some divergences that are not observed in, say, Figure \@ref(fig:tsne-clust-graph).
(This is at least partially due to the fact that $t$-SNE is itself graph-based and so will naturally agree more with a graph-based clustering strategy.)


```r
set.seed(100)
clust.kmeans <- kmeans(reducedDim(sce.pbmc, "PCA"), centers=10)
table(clust.kmeans$cluster)
```

```
## 
##   1   2   3   4   5   6   7   8   9  10 
## 472 200  46 515  90 320 241 865 735 438
```

```r
colLabels(sce.pbmc) <- factor(clust.kmeans$cluster)
plotReducedDim(sce.pbmc, "TSNE", colour_by="label")
```

<div class="figure">
<img src="clustering_files/figure-html/tsne-clust-kmeans-1.png" alt="$t$-SNE plot of the 10X PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from $k$-means clustering." width="672" />
<p class="caption">(\#fig:tsne-clust-kmeans)$t$-SNE plot of the 10X PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from $k$-means clustering.</p>
</div>

If we were so inclined, we could obtain a "reasonable" choice of $k$ by computing the gap statistic using methods from the *[cluster](https://CRAN.R-project.org/package=cluster)* package.
This is the log-ratio of the expected to observed within-cluster sum of squares, where the expected value is computed by randomly distributing cells within the minimum bounding box of the original data.
A larger gap statistic represents a lower observed sum of squares - and thus better clustering - compared to a population with no structure.
Ideally, we would choose the $k$ that maximizes the gap statistic, but this is often unhelpful as the tendency of $k$-means to favor spherical clusters drives a large $k$ to capture different cluster shapes.
Instead, we choose the most parsimonious $k$ beyond which the increases in the gap statistic are considered insignificant (Figure \@ref(fig:kmeans-gap)).


```r
library(cluster)
set.seed(110010101)
gaps <- clusGap(reducedDim(sce.pbmc, "PCA"), kmeans, K.max=20)
best.k <- maxSE(gaps$Tab[,"gap"], gaps$Tab[,"SE.sim"])
best.k
```

```
## [1] 9
```

```r
plot(gaps$Tab[,"gap"], xlab="Number of clusters", ylab="Gap statistic")
abline(v=best.k, col="red")
```

<div class="figure">
<img src="clustering_files/figure-html/kmeans-gap-1.png" alt="Gap statistic with respect to increasing number of $k$-means clusters in the 10X PBMC dataset. The red line represents the chosen $k$." width="672" />
<p class="caption">(\#fig:kmeans-gap)Gap statistic with respect to increasing number of $k$-means clusters in the 10X PBMC dataset. The red line represents the chosen $k$.</p>
</div>



A more practical use of $k$-means is to deliberately set $k$ to a large value to achieve overclustering.
This will forcibly partition cells inside broad clusters that do not have well-defined internal structure.
For example, we might be interested in the change in expression from one "side" of a cluster to the other, but the lack of any clear separation within the cluster makes it difficult to separate with graph-based methods, even at the highest resolution.
$k$-means has no such problems and will readily split these broad clusters for greater resolution, though obviously one must be prepared for the additional work involved in interpreting a greater number of clusters.


```r
set.seed(100)
clust.kmeans2 <- kmeans(reducedDim(sce.pbmc, "PCA"), centers=20)
table(clust.kmeans2$cluster)
```

```
## 
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
## 153 172  47 254 125 207 160 334 204 442 163  68 192 271 113 168 124 420  45 260
```

```r
colLabels(sce.pbmc) <- factor(clust.kmeans2$cluster)
plotTSNE(sce.pbmc, colour_by="label", text_by="label")
```

<div class="figure">
<img src="clustering_files/figure-html/tsne-clust-kmeans-best-1.png" alt="$t$-SNE plot of the 10X PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from $k$-means clustering with $k=20$." width="672" />
<p class="caption">(\#fig:tsne-clust-kmeans-best)$t$-SNE plot of the 10X PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from $k$-means clustering with $k=20$.</p>
</div>

### Assessing cluster separation

The within-cluster sum of squares (WCSS) for each cluster is the most relevant diagnostic for $k$-means, given that the algorithm aims to find a clustering that minimizes the WCSS.
Specifically, we use the WCSS to compute the root-mean-squared deviation (RMSD) that represents the spread of cells within each cluster.
A cluster is more likely to have a low RMSD if it has no internal structure and is separated from other clusters (such that there are not many cells on the boundaries between clusters, which would result in a higher sum of squares from the centroid).


```r
ncells <- tabulate(clust.kmeans2$cluster)
tab <- data.frame(wcss=clust.kmeans2$withinss, ncells=ncells)
tab$rms <- sqrt(tab$wcss/tab$ncells)
tab
```

```
##     wcss ncells   rms
## 1   2872    153 4.333
## 2   4204    172 4.944
## 3   1443     47 5.542
## 4   4275    254 4.103
## 5   1711    125 3.700
## 6   3055    207 3.842
## 7   1633    160 3.194
## 8   2160    334 2.543
## 9   7027    204 5.869
## 10  2858    442 2.543
## 11  4259    163 5.112
## 12  2289     68 5.802
## 13  2112    192 3.317
## 14  6391    271 4.856
## 15  1603    113 3.767
## 16 12400    168 8.591
## 17  1823    124 3.834
## 18  7026    420 4.090
## 19  2591     45 7.587
## 20  3640    260 3.742
```

(As an aside, the RMSDs of the clusters are poorly correlated with their sizes in Figure \@ref(fig:tsne-clust-kmeans-best).
This highlights the risks of attempting to quantitatively interpret the sizes of visual clusters in $t$-SNE plots.)

To explore the relationships between $k$-means clusters, a natural approach is to compute distances between their centroids.
This directly lends itself to visualization as a tree after hierarchical clustering (Figure \@ref(fig:kmeans-tree)). 


```r
cent.tree <- hclust(dist(clust.kmeans2$centers), "ward.D2")
plot(cent.tree)
```

<div class="figure">
<img src="clustering_files/figure-html/kmeans-tree-1.png" alt="Hierarchy of $k$-means cluster centroids, using Ward's minimum variance method." width="672" />
<p class="caption">(\#fig:kmeans-tree)Hierarchy of $k$-means cluster centroids, using Ward's minimum variance method.</p>
</div>

### In two-step procedures

As previously mentioned, $k$-means is most effective in its role of vector quantization, 
i.e., compressing adjacent cells into a single representative point.
This allows $k$-means to be used as a prelude to more sophisticated and interpretable (but expensive) clustering algorithms.
For example, the `clusterSNNGraph()` function from *[scran](https://bioconductor.org/packages/3.12/scran)* can use $k$-means as an initial step to obtain representative centroids that are then subjected to graph-based clustering (Figure \@ref(fig:tsne-kmeans-graph-pbmc)).
Each cell is then placed in the same graph-based cluster that its $k$-means centroid was assigned to.


```r
# Setting the seed due to the randomness of k-means.
set.seed(0101010)
kgraph.clusters <- clusterSNNGraph(sce.pbmc, use.dimred="PCA", 
    use.kmeans=TRUE, kmeans.centers=1000, k=5)
table(kgraph.clusters)
```

```
## kgraph.clusters
##   1   2   3   4   5   6   7   8   9  10  11 
## 840 137 550 517 220 528 829  46 127  83  45
```

```r
plotTSNE(sce.pbmc, colour_by=I(kgraph.clusters))
```

<div class="figure">
<img src="clustering_files/figure-html/tsne-kmeans-graph-pbmc-1.png" alt="$t$-SNE plot of the PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from combined $k$-means/graph-based clustering." width="672" />
<p class="caption">(\#fig:tsne-kmeans-graph-pbmc)$t$-SNE plot of the PBMC dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from combined $k$-means/graph-based clustering.</p>
</div>

The obvious benefit of this approach over direct graph-based clustering is the speed improvement.
We avoid the need to identifying nearest neighbors for each cell and the construction of a large intermediate graph,
while benefiting from the relative interpretability of graph-based clusters compared to those from $k$-means.
This approach also mitigates the "inflation" effect discussed in Section \@ref(clustering-graph).
Each centroid serves as a representative of a region of space that is roughly similar in volume,
ameliorating differences in cell density that can cause (potentially undesirable) differences in resolution.

The choice of the number of $k$-means clusters (defined here by the `kmeans.clusters=` argument) determines the trade-off between speed and fidelity.
Larger values provide a more faithful representation of the underlying distribution of cells,
at the cost of requiring more computational work by the second-stage clustering procedure.
Note that the second step operates on the centroids, so increasing `kmeans.clusters=` may have further implications if the second-stage procedure is sensitive to the total number of input observations. 
For example, increasing the number of centroids would require an concomitant increase in `k=` (the number of neighbors in graph construction) to maintain the same level of resolution in the final output.

## Hierarchical clustering

### Background

Hierarchical clustering is an ancient technique that aims to generate a dendrogram containing a hierarchy of samples.
This is most commonly done by greedily agglomerating samples into clusters, then agglomerating those clusters into larger clusters, and so on until all samples belong to a single cluster. 
Variants of hierarchical clustering methods primarily differ in how they choose to perform the agglomerations.
For example, complete linkage aims to merge clusters with the smallest maximum distance between their elements, while Ward's method aims to minimize the increase in within-cluster variance.

In the context of scRNA-seq, the main advantage of hierarchical clustering lies in the production of the dendrogram.
This is a rich summary that describes the relationships between cells and subpopulations at various resolutions and in a quantitative manner based on the branch lengths.
Users can easily "cut" the tree at different heights to define clusters with different granularity, where clusters defined at high resolution are guaranteed to be nested within those defined at a lower resolution.
(Guaranteed nesting can be helpful for interpretation, as discussed in Section \@ref(subclustering).)
The dendrogram is also a natural representation of the data in situations where cells have descended from a relatively recent common ancestor.

In practice, hierachical clustering is too slow to be used for anything but the smallest scRNA-seq datasets.
Most variants require a cell-cell distance matrix that is prohibitively expensive to compute for many cells.
Greedy agglomeration is also likely to result in a quantitatively suboptimal partitioning (as defined by the agglomeration measure) at higher levels of the dendrogram when the number of cells and merge steps is high.
Nonetheless, we will still demonstrate the application of hierarchical clustering here, as it can occasionally be useful for squeezing more information out of datasets with very few cells.

### Implementation

As the PBMC dataset is too large, we will demonstrate on the 416B dataset instead.

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
## colData names(10): Source Name cell line ... block sizeFactor
## reducedDimNames(2): PCA TSNE
## altExpNames(2): ERCC SIRV
```

We compute a cell-cell distance matrix using the top PCs and we apply hierarchical clustering with Ward's method.
The resulting tree in Figure \@ref(fig:dend-416b) shows a clear split in the population caused by oncogene induction. 
While both Ward's method and complete linkage (`hclust()`'s default) yield compact clusters, we prefer the former it is less affected by differences in variance between clusters.


```r
dist.416b <- dist(reducedDim(sce.416b, "PCA"))
tree.416b <- hclust(dist.416b, "ward.D2")

# Making a prettier dendrogram.
library(dendextend)
tree.416b$labels <- seq_along(tree.416b$labels)
dend <- as.dendrogram(tree.416b, hang=0.1)

combined.fac <- paste0(sce.416b$block, ".", 
    sub(" .*", "", sce.416b$phenotype))
labels_colors(dend) <- c(
    `20160113.wild`="blue",
    `20160113.induced`="red",
    `20160325.wild`="dodgerblue",
    `20160325.induced`="salmon"
)[combined.fac][order.dendrogram(dend)]

plot(dend)
```

<div class="figure">
<img src="clustering_files/figure-html/dend-416b-1.png" alt="Hierarchy of cells in the 416B data set after hierarchical clustering, where each leaf node is a cell that is coloured according to its oncogene induction status (red is induced, blue is control) and plate of origin (light or dark)." width="672" />
<p class="caption">(\#fig:dend-416b)Hierarchy of cells in the 416B data set after hierarchical clustering, where each leaf node is a cell that is coloured according to its oncogene induction status (red is induced, blue is control) and plate of origin (light or dark).</p>
</div>



To obtain explicit clusters, we "cut" the tree by removing internal branches such that every subtree represents a distinct cluster.
This is most simply done by removing internal branches above a certain height of the tree, as performed by the `cutree()` function.
A more sophisticated variant of this approach is implemented in the *[dynamicTreeCut](https://CRAN.R-project.org/package=dynamicTreeCut)* package,
which uses the shape of the branches to obtain a better partitioning for complex dendrograms (Figure \@ref(fig:dend-cluster)).


```r
library(dynamicTreeCut)

# minClusterSize needs to be turned down for small datasets.
# deepSplit controls the resolution of the partitioning.
clust.416b <- cutreeDynamic(tree.416b, distM=as.matrix(dist.416b),
    minClusterSize=10, deepSplit=1)
```

```
##  ..cutHeight not given, setting it to 783  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
table(clust.416b)
```

```
## clust.416b
##  1  2  3  4 
## 78 69 24 14
```

```r
labels_colors(dend) <- clust.416b[order.dendrogram(dend)]
plot(dend)
```

<div class="figure">
<img src="clustering_files/figure-html/dend-cluster-1.png" alt="Hierarchy of cells in the 416B data set after hierarchical clustering, where each leaf node is a cell that is coloured according to its assigned cluster identity from a dynamic tree cut." width="672" />
<p class="caption">(\#fig:dend-cluster)Hierarchy of cells in the 416B data set after hierarchical clustering, where each leaf node is a cell that is coloured according to its assigned cluster identity from a dynamic tree cut.</p>
</div>

This generally corresponds well to the grouping of cells on a $t$-SNE plot (Figure \@ref(fig:tsne-416b)).
The exception is cluster 2, which is split across two visual clusters in the plot.
We attribute this to a distortion introduced by $t$-SNE rather than inappropriate behavior of the clustering algorithm, based on the examination of some later diagnostics.


```r
colLabels(sce.416b) <- factor(clust.416b)
plotReducedDim(sce.416b, "TSNE", colour_by="label")
```

<div class="figure">
<img src="clustering_files/figure-html/tsne-416b-1.png" alt="$t$-SNE plot of the 416B dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from hierarchical clustering." width="672" />
<p class="caption">(\#fig:tsne-416b)$t$-SNE plot of the 416B dataset, where each point represents a cell and is coloured according to the identity of the assigned cluster from hierarchical clustering.</p>
</div>



### Assessing cluster separation {#silhouette-width}

We check the separation of the clusters using the silhouette width (Figure \@ref(fig:silhouette416b)).
For each cell, we compute the average distance to all cells in the same cluster.
We also compute the average distance to all cells in another cluster, taking the minimum of the averages across all other clusters.
The silhouette width for each cell is defined as the difference between these two values divided by their maximum.
Cells with large positive silhouette widths are closer to other cells in the same cluster than to cells in different clusters.

Each cluster would ideally contain large positive silhouette widths, indicating that it is well-separated from other clusters.
This is indeed the case in Figure \@ref(fig:silhouette416b) - and in fact, cluster 2 has the largest width of all, indicating that it is a more coherent cluster than portrayed in Figure \@ref(fig:tsne-416b).
Smaller widths can arise from the presence of internal subclusters, which inflates the within-cluster distance; or overclustering, where cells at the boundary of a partition are closer to the neighboring cluster than their own cluster.


```r
sil <- silhouette(clust.416b, dist = dist.416b)
plot(sil)
```

<div class="figure">
<img src="clustering_files/figure-html/silhouette416b-1.png" alt="Silhouette widths for cells in each cluster in the 416B dataset. Each bar represents a cell, grouped by the cluster to which it is assigned." width="672" />
<p class="caption">(\#fig:silhouette416b)Silhouette widths for cells in each cluster in the 416B dataset. Each bar represents a cell, grouped by the cluster to which it is assigned.</p>
</div>



For a more detailed examination, we identify the closest neighboring cluster for cells with negative widths.
This provides a perspective on the relationships between clusters that is closer to the raw data than the dendrogram in Figure \@ref(fig:dend-cluster).


```r
neg.widths <- sil[,3] < 0
table(Cluster=sil[neg.widths,1], Neighbor=sil[neg.widths,2])
```

```
##        Neighbor
## Cluster 1 2 3
##       2 0 0 3
##       3 1 3 0
```

The average silhouette width across all cells can also be used to choose clustering parameters.
The aim is to maximize the average silhouette width in order to obtain well-separated clusters.
This can be helpful to automatically obtain a "reasonable" clustering, though in practice, the clustering that yields the strongest separation often does not provide the most biological insight.

## General-purpose cluster diagnostics

### Cluster separation, redux

We previously introduced the silhouette width in the context of hierarchical clustering (Section \@ref(silhouette-width)).
While this can be applied with other clustering algorithms, 
it requires calculation of all pairwise distances between cells and is not scalable for larger cdatasets.
In such cases, we instead use an approximate approach that replaces the average of the distances with the distance to the average (i.e., centroid) of each cluster, with some tweaks to account for the distance due to the within-cluster variance.
This is implemented in the `clusterSilhouette()` function from *[scran](https://bioconductor.org/packages/3.12/scran)*,
allowing us to quickly identify poorly separate clusters with mostly negative widths (Figure \@ref(fig:pbmc-silhouette)).


```r
sil.approx <- clusterSilhouette(sce.pbmc, clusters=clust, use.dimred="PCA") 

sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust)

ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")
```

<div class="figure">
<img src="clustering_files/figure-html/pbmc-silhouette-1.png" alt="Distribution of the approximate silhouette width across cells in each cluster of the PBMC dataset. Each point represents a cell and colored with the identity of its own cluster if its silhouette width is positive and that of the closest other cluster if the width is negative." width="672" />
<p class="caption">(\#fig:pbmc-silhouette)Distribution of the approximate silhouette width across cells in each cluster of the PBMC dataset. Each point represents a cell and colored with the identity of its own cluster if its silhouette width is positive and that of the closest other cluster if the width is negative.</p>
</div>

Alternatively, we can quantify the degree to which cells from multiple clusters intermingle in expression space.
The "clustering purity" is defined for each cell as the proportion of neighboring cells that are assigned to the same cluster.
Well-separated clusters should exhibit little intermingling and thus high purity values for all member cells,
as demonstrated below in Figure \@ref(fig:pbmc-box-purities).
Median purity values are consistently greater than 0.9,
indicating that most cells in each cluster are primarily surrounded by other cells from the same cluster.


```r
pure.pbmc <- clusterPurity(sce.pbmc, clusters=clust, use.dimred="PCA") 

pure.data <- as.data.frame(pure.pbmc)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- factor(clust)

ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley")
```

<div class="figure">
<img src="clustering_files/figure-html/pbmc-box-purities-1.png" alt="Distribution of cluster purities across cells in each cluster of the PBMC dataset. Each point represents a cell and colored with the identity of the cluster contributing the largest proportion of its neighbors." width="672" />
<p class="caption">(\#fig:pbmc-box-purities)Distribution of cluster purities across cells in each cluster of the PBMC dataset. Each point represents a cell and colored with the identity of the cluster contributing the largest proportion of its neighbors.</p>
</div>



The main difference between these two methods is that the purity is ignorant of the intra-cluster variance.
This may or may not be desirable depending on what level of heterogeneity is of interest.
In addition, the purity will - on average - only decrease with increasing cluster number/resolution,
making it less effective for choosing between different clusterings.
However, regardless of the chosen method, it is worth keeping in mind that poor separation is not synonymous with poor quality. 
In fact, poorly separated clusters will often be observed in non-trivial analyses of scRNA-seq data 
where the aim is to characterize closely related subtypes or states.
These diagnostics are best used to guide interpretation by highlighting clusters that require more investigation
rather than to rule out poorly separated clusters altogether. 

### Comparing different clusterings

As previously mentioned, clustering's main purpose is to obtain a discrete summary of the data for further interpretation.
The diversity of available methods (and the subsequent variation in the clustering results)
reflects the many different "perspectives" that can be derived from a high-dimensional scRNA-seq dataset.
It is helpful to determine how these perspectives relate to each other by comparing the clustering results.
More concretely, we want to know which clusters map to each other across algorithms;
inconsistencies may be indicative of complex variation that is summarized differently by each clustering procedure.

A simple yet effective approach for comparing two clusterings of the same dataset
is to create a 2-dimensional table of label frequencies (Figure \@ref(fig:walktrap-v-others)).
We can further improve the interpretability of this table by computing the proportions of cell assignments,
which avoids difficulties with dynamic range when visualizing clusters of differing abundances,
For example, we may be interested in how our Walktrap clusters from Section \@ref(clustering-graph)
are redistributed when we switch to using Louvain community detection (Figure \@ref(fig:walktrap-v-louvain-prop)).
Note that this heatmap is best interpreted on a row-by-row basis as the proportions are computed per row
and there is no guarantee that cell abundances are comparable across columns.

<!--
No point interpreting it by column, as you'd get multiple 100% values.
An entry with a larger proportion may have a lower number of cells along a column.
-->


```r
tab <- table(Walktrap=clust, Louvain=clust.louvain)
tab <- tab/rowSums(tab)
pheatmap(tab, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)
```

<div class="figure">
<img src="clustering_files/figure-html/walktrap-v-louvain-prop-1.png" alt="Heatmap of the proportions of cells from each Walktrap cluster (rows) across the Louvain clusters (columns) in the PBMC dataset." width="672" />
<p class="caption">(\#fig:walktrap-v-louvain-prop)Heatmap of the proportions of cells from each Walktrap cluster (rows) across the Louvain clusters (columns) in the PBMC dataset.</p>
</div>

For clusterings that differ primarily in resolution (usually from different parameterizations of the same algorithm),
we can use the *[clustree](https://CRAN.R-project.org/package=clustree)* package to visualize the relationships between them.
Here, the aim is to capture the redistribution of cells from one clustering to another at progressively higher resolution, 
providing a convenient depiction of how clusters split apart (Figure \@ref(fig:walktrap-res-clustree)).
This approach is most effective with multiple nested mappings but is less useful when the mappings are "crossed",
e.g., for comparisons involving theoretically distinct clustering procedures.


```r
library(clustree)
combined <- cbind(k.50=clust.50, k.10=clust, k.5=clust.5)
clustree(combined, prefix="k.", edge_arrow=FALSE)
```

<div class="figure">
<img src="clustering_files/figure-html/walktrap-res-clustree-1.png" alt="Graph of the relationships between the Walktrap clusterings of the PBMC dataset, generated with varying $k$ during the nearest-neighbor graph construction. (A higher $k$ generally corresponds to a lower resolution clustering.) The size of the nodes is proportional to the number of cells in each cluster, and the edges depict cells in one cluster that are reassigned to another cluster at a different resolution. The color of the edges is defined according to the number of reassigned cells and the opacity is defined from the corresponding proportion relative to the size of the lower-resolution cluster." width="960" />
<p class="caption">(\#fig:walktrap-res-clustree)Graph of the relationships between the Walktrap clusterings of the PBMC dataset, generated with varying $k$ during the nearest-neighbor graph construction. (A higher $k$ generally corresponds to a lower resolution clustering.) The size of the nodes is proportional to the number of cells in each cluster, and the edges depict cells in one cluster that are reassigned to another cluster at a different resolution. The color of the edges is defined according to the number of reassigned cells and the opacity is defined from the corresponding proportion relative to the size of the lower-resolution cluster.</p>
</div>

We can quantify the agreement between two clusterings by computing the Rand index with *[scran](https://bioconductor.org/packages/3.12/scran)*'s `clusterRand()`.
This is defined as the proportion of pairs of cells that retain the same status
(i.e., both cells in the same cluster, or each cell in different clusters) in both clusterings.
A larger Rand index indicates that the clusters are preserved though this tends to be driven by the most abundant clusters.


```r
clusterRand(clust, clust.5, mode="index")
```

```
## [1] 0.9555
```

A more granular perspective can be generated by `clusterRand()` with `mode="ratio"`,
where the Rand index is broken down into its contributions from each cluster pair.
Low values on the diagonal in Figure \@ref(fig:pbmc-rand-breakdown) indicate that 
cells from the corresponding cluster in `clust` are redistributed to multiple other clusters in `clust.50`. 
Conversely, low off-diagonal values indicate that the corresponding pair of clusters in `clust` are merged together in `clust.50`.


```r
breakdown <- clusterRand(ref=clust, alt=clust.5, mode="ratio")
pheatmap(breakdown, color=viridis::magma(100), 
    cluster_rows=FALSE, cluster_cols=FALSE)
```

<div class="figure">
<img src="clustering_files/figure-html/pbmc-rand-breakdown-1.png" alt="Breakdown of the Rand index into its contributions from each pair of clusters in the reference Walktrap clustering compared to a higher-resolution alternative clustering for the PBMC dataset. Rows and columns of the heatmap represent clusters in the reference clustering. Each entry represents the proportion of pairs of cells involving the row/column clusters that retain the same status in the alternative clustering." width="672" />
<p class="caption">(\#fig:pbmc-rand-breakdown)Breakdown of the Rand index into its contributions from each pair of clusters in the reference Walktrap clustering compared to a higher-resolution alternative clustering for the PBMC dataset. Rows and columns of the heatmap represent clusters in the reference clustering. Each entry represents the proportion of pairs of cells involving the row/column clusters that retain the same status in the alternative clustering.</p>
</div>

Alternatively, we can compute coassignment probabilities between every pair of original clusters in `originals`.
For clusters $X$ and $Y$ in a reference clustering, the coassignment probability is defined as the probability that a randomly chosen cell from $X$ and a randomly chosen cell from $Y$ are assigned to the same cluster in the alternative clustering.
High co-assignment probabilities indicate that $X$ is not stable with respect to its separation from $Y$, given that their cells are liable to cluster together in the replicates.
This metric is closely related to the Rand index when the latter is broken down by cluster pairs 
though the interpretation of off-diagonal elements is reversed, 
i.e., low off-diagonal values are indicative of continued separation between the corresponding pair of clusters
(Figure \@ref(fig:pbmc-coassignment)).


```r
coassign <- coassignProb(clust, clust.5)
pheatmap(coassign, color=viridis::magma(100), 
    cluster_rows=FALSE, cluster_cols=FALSE)
```

<div class="figure">
<img src="clustering_files/figure-html/pbmc-coassignment-1.png" alt="Coassignment probabilities between each pair of clusters in the reference Walktrap clustering compared to a higher-resolution alternative clustering for the PBMC dataset. Rows and columns of the heatmap represent clusters in the reference clustering." width="672" />
<p class="caption">(\#fig:pbmc-coassignment)Coassignment probabilities between each pair of clusters in the reference Walktrap clustering compared to a higher-resolution alternative clustering for the PBMC dataset. Rows and columns of the heatmap represent clusters in the reference clustering.</p>
</div>

### Evaluating cluster stability {#cluster-bootstrapping}

A desirable property of a given clustering is that it is stable to perturbations to the input data [@luxburg2010clustering].
Stable clusters are logistically convenient as small changes to upstream processing will not change the conclusions;
greater stability also increases the likelihood that those conclusions can be reproduced in an independent replicate study.
*[scran](https://bioconductor.org/packages/3.12/scran)* uses bootstrapping to evaluate the stability of a clustering algorithm on a given dataset - that is, cells are sampled with replacement to create a "bootstrap replicate" dataset, and clustering is repeated on this replicate to see if the same clusters can be reproduced.
We demonstrate below for graph-based clustering on the PBMC dataset.


```r
myClusterFUN <- function(x) {
    g <- buildSNNGraph(x, use.dimred="PCA", type="jaccard")
    igraph::cluster_louvain(g)$membership
}

originals <- myClusterFUN(sce.pbmc)

set.seed(0010010100)
coassign <- bootstrapCluster(sce.pbmc, FUN=myClusterFUN, clusters=originals)
dim(coassign)
```

```
## [1] 19 19
```

The function returns a matrix of coassignment probabilities between every pair of original clusters in `originals` (Figure \@ref(fig:bootstrap-matrix)).
Ideally, we would hope for high coassignment probabilities on the diagonal (i.e., $X$ cells cluster with themselves) and low probabilities off the diagonal.
This would indicate that the clustering in the bootstrap replicates are highly consistent with that of the original dataset.


```r
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE,
    color=rev(viridis::magma(100)))
```

<div class="figure">
<img src="clustering_files/figure-html/bootstrap-matrix-1.png" alt="Heatmap of coassignment probabilities from bootstrapping of graph-based clustering in the PBMC dataset. Each row and column represents an original cluster and each entry is colored according to the value of the coassignment probability between that pair of clusters." width="672" />
<p class="caption">(\#fig:bootstrap-matrix)Heatmap of coassignment probabilities from bootstrapping of graph-based clustering in the PBMC dataset. Each row and column represents an original cluster and each entry is colored according to the value of the coassignment probability between that pair of clusters.</p>
</div>

Bootstrapping is a general approach for evaluating cluster stability that is compatible with any clustering algorithm.
The coassignment probability is also more informative than a single per-cluster stability measure as the former considers the relationships between clusters, e.g., unstable separation between $X$ and $Y$ does not penalize the stability of separation between $X$ and another cluster $Z$.
Of course, one should take these probabilities with a grain of salt, as bootstrapping only considers the effect of sampling noise and ignores other factors that affect reproducibility in an independent study (e.g., batch effects, donor variation).
In addition, it is possible for a poor separation to be highly stable, so highly stable cluster may not necessarily represent some distinct subpopulation.

## Subclustering {#subclustering}

Another simple approach to improving resolution is to repeat the feature selection and clustering _within_ a single cluster.
This aims to select HVGs and PCs that are more relevant to internal structure, improving resolution by avoiding noise from unnecessary features.
Subsetting also encourages clustering methods to separate cells according to more modest heterogeneity in the absence of distinct subpopulations.
We demonstrate with a cluster of putative memory T cells from the PBMC dataset, identified according to several markers (Figure \@ref(fig:ccr7-dist-memory)).


```r
g.full <- buildSNNGraph(sce.pbmc, use.dimred = 'PCA')
clust.full <- igraph::cluster_walktrap(g.full)$membership
plotExpression(sce.pbmc, features=c("CD3E", "CCR7", "CD69", "CD44"),
    x=I(factor(clust.full)), colour_by=I(factor(clust.full)))
```

<div class="figure">
<img src="clustering_files/figure-html/ccr7-dist-memory-1.png" alt="Distribution of log-normalized expression values for several T cell markers within each cluster in the 10X PBMC dataset. Each cluster is color-coded for convenience." width="672" />
<p class="caption">(\#fig:ccr7-dist-memory)Distribution of log-normalized expression values for several T cell markers within each cluster in the 10X PBMC dataset. Each cluster is color-coded for convenience.</p>
</div>

```r
# Repeating modelling and PCA on the subset.
memory <- 6L
sce.memory <- sce.pbmc[,clust.full==memory]
dec.memory <- modelGeneVar(sce.memory)
sce.memory <- denoisePCA(sce.memory, technical=dec.memory,
    subset.row=getTopHVGs(dec.memory, prop=0.1))
```



We apply graph-based clustering within this memory subset to obtain _CD4_^+^ and _CD8_^+^ subclusters (Figure \@ref(fig:cd48-memory-expression)).
Admittedly, the expression of _CD4_ is so low that the change is rather modest, but the interpretation is clear enough.


```r
g.memory <- buildSNNGraph(sce.memory, use.dimred="PCA")
clust.memory <- igraph::cluster_walktrap(g.memory)$membership
plotExpression(sce.memory, features=c("CD8A", "CD4"),
    x=I(factor(clust.memory)))
```

<div class="figure">
<img src="clustering_files/figure-html/cd48-memory-expression-1.png" alt="Distribution of _CD4_ and _CD8A_ log-normalized expression values within each cluster in the memory T cell subset of the 10X PBMC dataset." width="672" />
<p class="caption">(\#fig:cd48-memory-expression)Distribution of _CD4_ and _CD8A_ log-normalized expression values within each cluster in the memory T cell subset of the 10X PBMC dataset.</p>
</div>



For subclustering analyses, it is helpful to define a customized function that calls our desired algorithms to obtain a clustering from a given `SingleCellExperiment`.
This function can then be applied multiple times on different subsets without having to repeatedly copy and modify the code for each subset.
For example, `quickSubCluster()` loops over all subsets and executes this user-specified function to generate a list of `SingleCellExperiment` objects containing the subclustering results.
(Of course, the downside is that this assumes that a similar analysis is appropriate for each subset.
If different subsets require extensive reparametrization, copying the code may actually be more straightforward.)


```r
set.seed(1000010)
subcluster.out <- quickSubCluster(sce.pbmc, groups=clust.full,
    prepFUN=function(x) { # Preparing the subsetted SCE for clustering.
        dec <- modelGeneVar(x)
        input <- denoisePCA(x, technical=dec,
            subset.row=getTopHVGs(dec, prop=0.1),
            BSPARAM=BiocSingular::IrlbaParam())
    },
    clusterFUN=function(x) { # Performing the subclustering in the subset.
        g <- buildSNNGraph(x, use.dimred="PCA", k=20)
        igraph::cluster_walktrap(g)$membership
    }
)

# One SingleCellExperiment object per parent cluster:
names(subcluster.out)
```

```
##  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
## [16] "16" "17" "18"
```

```r
# Looking at the subclustering for one example:
table(subcluster.out[[1]]$subcluster)
```

```
## 
## 1.1 1.2 1.3 1.4 1.5 
## 179 113  44 117 132
```

Subclustering is a general and conceptually straightforward procedure for increasing resolution. 
It can also simplify the interpretation of the subclusters, which only need to be considered in the context of the parent cluster's identity - for example, we did not have to re-identify the cells in cluster 6 as T cells.
However, this is a double-edged sword as it is difficult for practitioners to consider the uncertainty of identification for parent clusters when working with deep nesting.
If cell types or states span cluster boundaries, conditioning on the putative cell type identity of the parent cluster can encourage the construction of a "house of cards" of cell type assignments, e.g., where a subcluster of one parent cluster is actually contamination from a cell type in a separate parent cluster.

## Session Info {-}

<button class="aaron-collapse">View session info</button>
<div class="aaron-content">
```
R version 4.0.0 Patched (2020-05-01 r78341)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

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
 [1] clustree_0.4.3              ggraph_2.0.3               
 [3] dynamicTreeCut_1.63-1       dendextend_1.13.4          
 [5] cluster_2.1.0               pheatmap_1.0.12            
 [7] scater_1.17.3               ggplot2_3.3.1              
 [9] scran_1.17.2                SingleCellExperiment_1.11.4
[11] SummarizedExperiment_1.19.5 DelayedArray_0.15.4        
[13] matrixStats_0.56.0          Matrix_1.2-18              
[15] Biobase_2.49.0              GenomicRanges_1.41.5       
[17] GenomeInfoDb_1.25.2         IRanges_2.23.10            
[19] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[21] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] bitops_1.0-6              RColorBrewer_1.1-2       
 [3] backports_1.1.7           tools_4.0.0              
 [5] R6_2.4.1                  irlba_2.3.3              
 [7] vipor_0.4.5               colorspace_1.4-1         
 [9] withr_2.2.0               tidyselect_1.1.0         
[11] gridExtra_2.3             processx_3.4.2           
[13] compiler_4.0.0            graph_1.67.1             
[15] BiocNeighbors_1.7.0       labeling_0.3             
[17] bookdown_0.19             checkmate_2.0.0          
[19] scales_1.1.1              callr_3.4.3              
[21] stringr_1.4.0             digest_0.6.25            
[23] rmarkdown_2.2             XVector_0.29.2           
[25] pkgconfig_2.0.3           htmltools_0.4.0          
[27] limma_3.45.7              highr_0.8                
[29] rlang_0.4.6               DelayedMatrixStats_1.11.0
[31] generics_0.0.2            farver_2.0.3             
[33] BiocParallel_1.23.0       dplyr_1.0.0              
[35] RCurl_1.98-1.2            magrittr_1.5             
[37] BiocSingular_1.5.0        GenomeInfoDbData_1.2.3   
[39] scuttle_0.99.9            Rcpp_1.0.4.6             
[41] ggbeeswarm_0.6.0          munsell_0.5.0            
[43] viridis_0.5.1             lifecycle_0.2.0          
[45] stringi_1.4.6             yaml_2.2.1               
[47] edgeR_3.31.4              MASS_7.3-51.6            
[49] zlibbioc_1.35.0           grid_4.0.0               
[51] ggrepel_0.8.2             dqrng_0.2.1              
[53] crayon_1.3.4              lattice_0.20-41          
[55] graphlayouts_0.7.0        cowplot_1.0.0            
[57] locfit_1.5-9.4            CodeDepends_0.6.5        
[59] knitr_1.28                ps_1.3.3                 
[61] pillar_1.4.4              igraph_1.2.5             
[63] codetools_0.2-16          XML_3.99-0.3             
[65] glue_1.4.1                evaluate_0.14            
[67] BiocManager_1.30.10       tweenr_1.0.1             
[69] vctrs_0.3.1               polyclip_1.10-0          
[71] tidyr_1.1.0               gtable_0.3.0             
[73] purrr_0.3.4               ggforce_0.3.1            
[75] xfun_0.14                 rsvd_1.0.3               
[77] tidygraph_1.2.0           viridisLite_0.3.0        
[79] tibble_3.0.1              beeswarm_0.2.3           
[81] statmod_1.4.34            ellipsis_0.3.1           
```
</div>
