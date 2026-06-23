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
de_genes <- read.csv("../../00_raw_data/reference-data/41467_2023_39762_MOESM6_ESM-cordsetal-supp3.csv") %>% dplyr::select(gene) %>% pull()
de_genes <- unique(de_genes)
integrated_seurat <- readRDS("../../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat_redo-integration-20251128.Rds")
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
raw_counts <- GetAssayData(integrated_seurat, assay = "RNA", layer = "counts")
cell_metadata <- integrated_seurat@meta.data
sce_object_all_genes <- SingleCellExperiment(
  assays = list(
    counts = raw_counts#,      # Use the matched matrix
    #logcounts = norm_integrated_data
  ),
  colData = cell_metadata
)
sce_object_all_genes <- logNormCounts(sce_object_all_genes)
saveRDS(sce_object_all_genes, file = "../../04_results/singler/sce_object_for_singleR_all_genes.Rds")
breast <- as.SingleCellExperiment(readRDS("/data/scRNA-seq-cords/BREAST_fibro_tumour.rds"))
rm(integrated_seurat, raw_counts, cell_metadata)
gc()
# SingleR() expects reference datasets to be normalized and log-transformed.
breast <- logNormCounts(breast)
breast_genes <- rownames(assays(breast)[[2]])
sce_genes <- rownames(sce_object_all_genes)
de_genes_overlap <- intersect(de_genes, breast_genes)
de_genes_overlap <- intersect(de_genes_overlap, sce_genes)
breast_subset <- breast[de_genes_overlap, ]
sce_object_subset <- sce_object_all_genes[de_genes_overlap, ]
pred.grun <- SingleR(test=sce_object_subset, 
                     ref=breast_subset, 
                     labels=breast$CAFtype, 
                     de.method="wilcox"#, 
                     #genes = de_genes_overlap
                     )
saveRDS(object = pred.grun, file = "../../04_results/singler/singler-default-params-labels-breast-ref-use-de-genes-supp3-cords.Rds")
saveRDS(object = sce_object_subset, file = "../../04_results/singler/sce_object_for_singleR_subset_de_genes_cords_supp3.Rds")
