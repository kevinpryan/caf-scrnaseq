library(SingleR)
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(vegan)
library(tidyr)
library(ggthemes)
library(forcats)
integrated_seurat <- readRDS("../../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat_redo-integration-20251128.Rds")
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
raw_counts <- GetAssayData(integrated_seurat, assay = "RNA", layer = "counts")
cell_metadata <- integrated_seurat@meta.data
norm_integrated_data <- GetAssayData(integrated_seurat, assay = "integrated", layer = "data")
cell_metadata <- integrated_seurat@meta.data
common_features <- rownames(norm_integrated_data)
common_cells <- colnames(norm_integrated_data)
raw_counts_matched <- raw_counts[common_features, common_cells]
cat("Verifying matrix dimensions:\n")
cat("Raw Counts (Matched): "); print(dim(raw_counts_matched))
cat("Integrated Data: "); print(dim(norm_integrated_data))
cell_metadata <- integrated_seurat@meta.data[common_cells, ] # Ensure metadata also matches

sce_object <- SingleCellExperiment(
  assays = list(
    counts = raw_counts_matched,      # Use the matched matrix
    logcounts = norm_integrated_data
  ),
  colData = cell_metadata
)
breast <- as.SingleCellExperiment(readRDS("/data/scRNA-seq-cords/BREAST_fibro_tumour.rds"))
rm(integrated_seurat, raw_counts, cell_metadata)
gc()
# SingleR() expects reference datasets to be normalized and log-transformed.
breast <- logNormCounts(breast)
pred.grun <- SingleR(test=sce_object, ref=breast, labels=breast$CAFtype, de.method="wilcox")
saveRDS(object = pred.grun, file = "../../04_results/singler/singler-default-params-labels-breast-ref-use-variable-genes-integrated.Rds")
