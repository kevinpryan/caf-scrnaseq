#!/usr/bin/env Rscript

library(argparse)
library(Seurat)
library(dplyr)
library(ggplot2)

# --- Argument Parsing ---
parser <- ArgumentParser(description='Integrate multiple normalized Seurat objects')
parser$add_argument('--input_files', type='character', nargs='+', required=TRUE, help='List of paths to the input normalized RDS files')
parser$add_argument('--output_rds', type='character', required=TRUE, help='Path for the final integrated RDS file')
parser$add_argument('--output_plot', type='character', required=TRUE, help='Path for the output UMAP plot')
args <- parser$parse_args()

# --- Load and Prepare Data ---
# Read the list of RDS files into a list of Seurat objects
seurat_object_list <- lapply(args$input_files, readRDS)

cat("All samples loaded. Preparing for integration...\n")

# --- Integration using SCTransform workflow ---
# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_object_list, nfeatures = 3000)

# Prepare the objects for integration
seurat_object_list <- PrepSCTIntegration(object.list = seurat_object_list, anchor.features = features)

# Find integration anchors
# Using rpca for speed and memory efficiency, recommended for large datasets
anchors <- FindIntegrationAnchors(object.list = seurat_object_list, normalization.method = "SCT", anchor.features = features, reduction = "rpca")

# Integrate data
integrated_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

cat("Integration complete. Running downstream analysis...\n")

# --- Downstream Analysis on the Integrated Object ---
# IMPORTANT: Perform analysis on the "integrated" assay
DefaultAssay(integrated_seurat) <- "integrated"

# Run PCA, FindNeighbors, FindClusters, and UMAP
integrated_seurat <- RunPCA(integrated_seurat, verbose = FALSE)
png("elbow_plot.png")
ElbowPlot(integrated_seurat, ndims = 50)
dev.off()

#integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30)
#integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)
#integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30)

# --- Save Final Integrated Object ---
saveRDS(integrated_seurat, file = args$output_rds)

# --- Create and Save a Visualization ---
# UMAP plot colored by sample
#p <- DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident") +
#    ggtitle("UMAP of Integrated Samples")

#ggsave(args$output_plot, plot = p, width = 8, height = 6)

cat("Workflow complete!\n")
