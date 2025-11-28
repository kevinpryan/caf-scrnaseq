# File: join_layers_script.R

library(Seurat)

# 1. Load the object
cat("Loading Seurat object...\n")
integrated_seurat <- readRDS("../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat.Rds")

# 2. Perform the memory-intensive operation
cat("Joining layers in the RNA assay...\n")
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")

# 3. Save the NEW, modified object to a different file
cat("Saving the object with joined layers...\n")
saveRDS(integrated_seurat, file = "../03_processed_data/03-clustering/integrated_seurat_object_LAYERS_JOINED.Rds")

cat("Script finished successfully!\n")