#!/usr/bin/env Rscript

library(argparse)
library(Seurat)
library(dplyr)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)

# --- Argument Parsing ---
parser <- ArgumentParser(description='Run QC, filtering, doublet detection, and normalization for a single 10x sample')
parser$add_argument('--sample_id', type='character', required=TRUE, help='The sample identifier')
parser$add_argument('--input_file', type='character', required=TRUE, help='Path to the raw data H5 file')
parser$add_argument('--output_rds', type='character', required=TRUE, help='Path for the output normalized RDS file')
parser$add_argument('--output_qc_metrics', type='character', required=TRUE, help='Path for the output QC metrics CSV file') # <-- NEW ARGUMENT
args <- parser$parse_args()

# --- Load Data ---
counts <- Read10X_h5(filename = args$input_file)
obj <- CreateSeuratObject(counts = counts, project = args$sample_id)
initial_cell_count <- ncol(obj)

# --- QC Metric Calculation ---
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

# --- Doublet Detection ---
cat("Running Doublet Detection...\n")
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)
obj$scDblFinder.class <- sce$scDblFinder.class
obj$scDblFinder.score <- sce$scDblFinder.score

# ==============================================================================
# --- NEW: GENERATE QC REPORT ---
# Calculate the number of cells that would be removed by each filter independently.
# Note: A single cell can fail multiple filters.
# ==============================================================================
cat("Generating QC metrics report...\n")
cells_fail_nFeature <- sum(obj@meta.data$nFeature_RNA < 250)
cells_fail_nCount <- sum(obj@meta.data$nCount_RNA < 500)
cells_fail_log10Genes <- sum(obj@meta.data$log10GenesPerUMI <= 0.80)
cells_fail_mito <- sum(obj@meta.data$percent.mt >= 5)
cells_fail_doublet <- sum(obj@meta.data$scDblFinder.class == "doublet")

# --- Perform the actual filtering ---
cat("Filtering cells based on QC metrics and doublet calls...\n")
filtered_seurat <- subset(x = obj,
                         subset = (nFeature_RNA >= 250) &
                                  (nCount_RNA >= 500) &
                                  (log10GenesPerUMI > 0.80) &
                                  (percent.mt < 5) &
                                  (scDblFinder.class == "singlet"))
final_cell_count <- ncol(filtered_seurat)
total_cells_removed <- initial_cell_count - final_cell_count

# Create a data frame with the QC summary
qc_summary <- data.frame(
  Metric = c(
    "Initial number of cells",
    "Cells failing nFeature_RNA (< 250)",
    "Cells failing nCount_RNA (< 500)",
    "Cells failing log10GenesPerUMI (<= 0.80)",
    "Cells failing percent.mt (>= 5%)",
    "Cells identified as doublets",
    "Total cells removed (combination of all filters)",
    "Final number of cells remaining"
  ),
  Count = c(
    initial_cell_count,
    cells_fail_nFeature,
    cells_fail_nCount,
    cells_fail_log10Genes,
    cells_fail_mito,
    cells_fail_doublet,
    total_cells_removed,
    final_cell_count
  )
)

# Calculate and format percentages based on the initial cell count
qc_summary$Percentage_of_Initial <- sprintf("%.2f%%", (qc_summary$Count / initial_cell_count) * 100)
qc_summary$Percentage_of_Initial[1] <- "100.00%" # Set initial to 100%

# Write the summary to the output file
write.csv(qc_summary, file = args$output_qc_metrics, row.names = FALSE)

# --- Gene Filtering ---
counts <- GetAssayData(object = filtered_seurat, assay = "RNA", layer = "counts")
genes_to_keep <- which(Matrix::rowSums(counts > 0) >= 10)
filtered_seurat <- filtered_seurat[genes_to_keep, ]

# --- Normalization ---
cat("Running SCTransform...\n")
filtered_seurat <- SCTransform(filtered_seurat, vars.to.regress = "percent.mt", verbose = TRUE)

# --- RunPCA ---
cat("Running PCA...\n")
filtered_seurat <- RunPCA(filtered_seurat, verbose = FALSE)
# --- Save Output ---
cat(paste("Saving normalized object for sample:", args$sample_id, "\n"))
saveRDS(filtered_seurat, file = args$output_rds)
