library(SingleR)
library(Seurat)
library(SingleCellExperiment)
# read in single cell data
integrated_seurat <- as.SingleCellExperiment(readRDS("../03_processed_data/03-clustering/integreated_seurat_object_clustered.Rds"))
# read in reference data
breast <- readRDS("~/Document")