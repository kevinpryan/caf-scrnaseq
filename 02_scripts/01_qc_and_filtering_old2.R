library(argparse)
library(Seurat)
library(dplyr)
library(Matrix)

# --- Argument Parsing ---
parser <- ArgumentParser(description='Run QC, filtering, and normalization for a single 10x sample')
parser$add_argument('--sample_id', type='character', required=TRUE, help='The sample identifier')
parser$add_argument('--input_file', type='character', required=TRUE, help='Path to the raw data H5 file')
parser$add_argument('--output_rds', type='character', required=TRUE, help='Path for the output normalized RDS file')
args <- parser$parse_args()

# --- Load Data ---
counts <- Read10X_h5(filename = args$input_file)
obj <- CreateSeuratObject(counts = counts, project = args$sample_id)

# --- QC Metric Calculation ---
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

# --- Cell Filtering ---
filtered_seurat <- subset(x = obj, 
                         subset = (nFeature_RNA >= 250) & 
                                  (nCount_RNA >= 500) & 
                                  (log10GenesPerUMI > 0.80) & 
                                  (percent.mt < 5))

# --- Gene Filtering ---
# Remove genes expressed in fewer than 10 cells
counts <- GetAssayData(object = filtered_seurat, assay = "RNA", layer = "counts")
genes_to_keep <- which(Matrix::rowSums(counts > 0) >= 10)
filtered_seurat <- filtered_seurat[genes_to_keep, ]

# --- Normalization ---
# Run SCTransform. This replaces NormalizeData, ScaleData, and FindVariableFeatures.
# It's more memory-efficient and is the recommended workflow.
# We also regress out mitochondrial mapping percentage as a confounding factor.
cat("Running SCTransform...\n")
filtered_seurat <- SCTransform(filtered_seurat, vars.to.regress = "percent.mt", verbose = FALSE)

# --- Save Output ---
cat(paste("Saving normalized object for sample:", args$sample_id, "\n"))
saveRDS(filtered_seurat, file = args$output_rds)
