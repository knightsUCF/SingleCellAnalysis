

<h3>What NicheNet Can Predict</h3>

1) which ligands from one cell population (“sender/niche”) are most likely to affect target gene expression in an interacting cell population (“receiver/target”)

2) which specific target genes are affected by which of these predicted ligands

<h3>Dependencies</h3>

```R
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
```

<h3>Cell Types Post Cluster Analysis</h3>

pbmc <- readRDS("pbmc.rds")


# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
# print(head(cluster1.markers, n = 5))

# find all markers distinguishing cluster 3 from clusters 0 and 3
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, ident.2 = c(0, 3), min.pct = 0.25) # last 3 here does not denote clusters but the dimensions of data
# print(head(cluster5.markers, n = 3))


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print(pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))

# Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# plot(VlnPlot(pbmc, features = c("Esam", "Tm4sf1")))

# we can plot raw counts as well
# plot(VlnPlot(pbmc, features = c("Esam", "Tm4sf1"), slot = "counts", log = TRUE))

# plot(FeaturePlot(pbmc, features = c("Esam", "Vwa1", "Egfl7", "Flt1", "Cd81", "C1ql1", "Cspg5", "Dlgap1", "Smoc2")))

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# plot(DoHeatmap(pbmc, features = top10$gene) + NoLegend())


new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)


hn1 <- RunTSNE(pbmc, dims.use = 1:20,reduction.use = "pca", dim_embed = 2)

plot(DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())




