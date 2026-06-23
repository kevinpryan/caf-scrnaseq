library(Seurat)
library(dplyr)
install.packages("clustree")
install.packages("devtools")
library(clustree)


#install.packages("devtools")
devtools::install_github('immunogenomics/presto')
# read in data
integrated_seurat <- readRDS("../results/02_integration/integrated_seurat_object.Rds")
# from the elbow plot it looks like the elbow is at 11
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:11)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:11)
VizDimLoadings(integrated_seurat, dims = 1:2, reduction = "pca")
DimPlot(integrated_seurat, reduction = "pca", group.by = "orig.ident") 
DimHeatmap(integrated_seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(integrated_seurat, dims = 1:11, cells = 500, balanced = TRUE)
DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated_seurat, reduction = "umap")
integrated_seurat <- FindClusters(integrated_seurat, resolution = seq(from = 0.1, to = 2, by = 0.1))
p <- clustree(integrated_seurat, prefix = "integrated_snn_res.", exprs = "scale.data")
plot(p)
# value of 0.4 gives same number of clusters as Cords breast, start there
Idents(integrated_seurat) <- integrated_seurat$integrated_snn_res.0.4
DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.4")
DimPlot(integrated_seurat, reduction = "pca", group.by = "integrated_snn_res.0.4")

p2 <- DimPlot(
  integrated_seurat,
  reduction = "umap",
  group.by = c("integrated_snn_res.0.4", "orig.ident"),
  combine = TRUE, label.size = 2
)
p2
#wrap_plots(c(p1, p2), ncol = 2, byrow = F)
saveRDS(integrated_seurat, file = "../03_processed_data/03-clustering/integreated_seurat_object_clustered.Rds")
integrated_seurat$integrated_snn_res.
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
integrated_seurat_markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE)
sigs <- integrated_seurat_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
