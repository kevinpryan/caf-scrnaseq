library(argparse)
library(Seurat)
library(dplyr)
library(ggplot2)
parser <- ArgumentParser(description='Run QC and filtering for a single sample')
#parser$add_argument('--samplesheet', type='character', required=TRUE, help='The sample identifier')
parser$add_argument('--sample_id', type='character', required=TRUE, help='The sample identifier')
parser$add_argument('--input_file', type='character', required=TRUE, help='Path to the raw data directory (matrix, features, barcodes)')
parser$add_argument('--output_rds', type='character', required=TRUE, help='Path for the output filtered RDS file')
args <- parser$parse_args()
sample_id <- args$sample_id
input_file <- args$input_file
output_rds <- args$output_rds
#samplesheet <- read.csv(args$samplesheet, header = TRUE)
#samplesheet$path <- paste("/home/rstudio/caf-scrnaseq", samplesheet$path, sep = "/")
#print("samplesheet...")
#print(samplesheet)
counts <- Read10X_h5(filename = input_file)
obj <- CreateSeuratObject(counts = counts, project = sample_id)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
obj$mitoRatio <- PercentageFeatureSet(object = obj, pattern = "^MT-")
obj$mitoRatio <- obj@meta.data$mitoRatio / 100
metadata <- obj@meta.data
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
obj@meta.data <- metadata
filtered_seurat <- subset(x = obj, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.05))
filtered_seurat[["RNA_joined"]] <- JoinLayers(filtered_seurat[["RNA"]])
DefaultAssay(filtered_seurat) <- "RNA_joined"
counts <- filtered_seurat[["RNA_joined"]]$counts 
#LayerData(object = filtered_seurat, layer = "counts")
n_cells_per_gene <- Matrix::rowSums(counts > 0)

# Output a logical vector for every gene on whether the more than zero counts per cell
genes_to_keep <- names(n_cells_per_gene[n_cells_per_gene >= 10])
seu_obj_filtered_genes <- filtered_seurat[genes_to_keep, ]
DefaultAssay(seu_obj_filtered_genes) <- "RNA"

print(paste("Original number of genes:", nrow(filtered_seurat)))
print(paste("Number of genes to keep after filtering:", length(genes_to_keep)))
metadata_clean <- seu_obj_filtered_genes@meta.data
saveRDS(seu_obj_filtered_genes, file = output_rds)


#saveRDS(scrnaseq_data, file = args$output_rds)
  
