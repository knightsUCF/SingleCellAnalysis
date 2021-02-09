# Overview

Single cell analysis aims to analyze the frequency of gene groups in a cell. In the case of a cell injury, we can detect which genes seem significant in the regeneration process.

The following data and methods are based on the research paper:

https://www.biorxiv.org/content/10.1101/2020.05.13.094854v1


<h3>Processing Data</h3>

```R
library(dplyr)
library(Seurat)
library(patchwork)


path = '../Data/Original/'
pbmc.data <- Read10X(data.dir = path, gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "SC", min.cells = 3, min.features = 200)


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
# plot(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot(VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2))
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/QCMetrics.png)

```R
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot(plot1 + plot2)

plot(plot2)
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/ScatterCountvsFeature.png)


<h3>Identifying Highly Variable Features</h3>

```R
# remove unwanted cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# normalize the rest
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot(plot3 + plot4)
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/HighlyVariableFeatures.png)



<h3>Scaling the Data</h3>

```R
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

<h3>Dimensional reduction</h3>

```R
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

output:

PC_ 1 
Positive:  Crip2, Tm4sf1, Hspb1, Cavin3, Pcp4l1 
Negative:  Ctss, Fcgr3, Cotl1, Cd68, Gm49339 
PC_ 2 
Positive:  Esam, Ly6c1, Egfl7, Flt1, Eng 
Negative:  Mt3, Clu, Gfap, Chchd10, Rsph1 
PC_ 3 
Positive:  C1qa, C1qc, Cx3cr1, Ly86, Cd81 
Negative:  S100a6, Thbs1, S100a10, Lgals3, Crip1 
PC_ 4 
Positive:  Cspg5, Dlgap1, Gm3764, C1ql1, Ptprz1 
Negative:  Ccdc153, Ccdc113, 1110017D15Rik, Rsph9, Tmem212 
PC_ 5 
Positive:  Cspg5, C1ql1, Olig1, Tmod2, Dlgap1 
Negative:  Pcolce, Bgn, Smoc2, Gpc3, Fbln1


plot(DimPlot(pbmc, reduction = "pca"))
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/pca.png)

<h3>Dimensionality</h3>

```R
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

plot(JackStrawPlot(pbmc, dims = 1:15))
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/Dimensionality.png)


<h3>Elbow Plot</h3>

```R
plot(ElbowPlot(pbmc))
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/ElbowPlot.png)


<h3>Clustering</h3>

```R
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

output:

Computing nearest neighbor graph
Computing SNN
Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

Number of nodes: 344
Number of edges: 8121

Running Louvain algorithm...
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Maximum modularity in 10 random starts: 0.8660
Number of communities: 5
```

<h3>Non-linear Dimensional Reduction (UMAP)</h3>

```R
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that we can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot(DimPlot(pbmc, reduction = "umap"))
```

![](https://github.com/knightsUCF/SingleCellAnalysis/blob/main/images/UMAP.png)

<h3>Finding differentially expressed features (cluster biomarkers)</h3>

```R
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
print(head(cluster1.markers, n = 5))

output

               p_val avg_log2FC pct.1 pct.2    p_val_adj
Esam    1.279957e-59   3.595577 0.933 0.042 1.611209e-55
Adgrl4  6.386108e-59   2.659595 0.867 0.018 8.038833e-55
Tmem252 8.195863e-58   3.703981 0.883 0.028 1.031695e-53
Vwa1    5.332391e-56   3.521792 0.900 0.046 6.712414e-52
Egfl7   5.689981e-56   3.957891 0.983 0.092 7.162548e-52
```


```R
# find all markers distinguishing cluster 3 from clusters 0 and 3
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, ident.2 = c(0, 3), min.pct = 0.25) # last 3 here does not denote clusters but the dimensions of data
print(head(cluster5.markers, n = 3))

output

               p_val avg_log2FC pct.1 pct.2     p_val_adj
FCGR3A 7.583625e-209   4.274913 0.975 0.037 1.040018e-204
IFITM3 2.500844e-199   3.892661 0.975 0.046 3.429657e-195
CFD    1.763722e-195   3.408196 0.938 0.037 2.418768e-191
```

