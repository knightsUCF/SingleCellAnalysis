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




