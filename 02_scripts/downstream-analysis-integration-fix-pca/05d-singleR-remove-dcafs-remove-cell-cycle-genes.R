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
# ==============================================================================
# 1. SETUP & LOAD DATA
# ==============================================================================
integrated_seurat <- readRDS("../../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat_redo-integration-20251128.Rds")
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
# Load Reference
breast <- as.SingleCellExperiment(readRDS("/data/scRNA-seq-cords/BREAST_fibro_tumour.rds"))
breast <- logNormCounts(breast)
# ==============================================================================
# 2. CALCULATE CELL CYCLE STATE (The "Dividing" Label)
# ==============================================================================
# A. Get standard cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# B. Score the cells (Seurat adds 'S.Score', 'G2M.Score', and 'Phase' to metadata)
# Note: Run this on the RNA assay (raw counts normalized), not integrated residuals
DefaultAssay(integrated_seurat) <- "RNA"
integrated_seurat <- NormalizeData(integrated_seurat)
integrated_seurat <- CellCycleScoring(integrated_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
# C. Create a custom "Status" column
# "Phase" gives G1, S, or G2M. Let's simplify to "Dividing" vs "Resting"
integrated_seurat$Cycle_Status <- ifelse(integrated_seurat$Phase == "G1", "Resting", "Dividing")
# ==============================================================================
# 3. PREPARE FOR SINGLER (The "Lineage" Label)
# ==============================================================================
# A. Filter the Reference: Remove dCAF cells entirely
breast_clean <- breast[, breast$CAFtype != "dCAF"]
breast_clean$CAFtype <- droplevels(breast_clean$CAFtype) # Clean up empty factor level
# B. Filter the Genes: Remove cell cycle genes from the comparison
# Even though we removed dCAF cells, removing the genes ensures the "Dividing"
# signal doesn't confuse the mCAF vs iCAF decision.
common_genes <- intersect(rownames(integrated_seurat), rownames(breast_clean))
genes_no_cc <- setdiff(common_genes, c(s.genes, g2m.genes))
# C. Prepare Query SCE (Raw counts -> LogNorm)
sce_query <- SingleCellExperiment(
  assays = list(counts = GetAssayData(integrated_seurat, assay = "RNA", layer = "counts")),
  colData = integrated_seurat@meta.data
)
sce_query <- logNormCounts(sce_query)
# ==============================================================================
# 4. RUN SINGLER
# ==============================================================================
rm(breast)
gc()
pred <- SingleR(test = sce_query[genes_no_cc, ], 
                ref = breast_clean[genes_no_cc, ], 
                labels = breast_clean$CAFtype, 
                de.method = "wilcox")

# Add predictions to Seurat
integrated_seurat$SingleR_Lineage <- pred$labels
# ==============================================================================
# 5. CREATE COMPOSITE LABELS
# ==============================================================================
# Combine the two columns
integrated_seurat$Final_Annotation <- paste(integrated_seurat$Cycle_Status, 
                                            integrated_seurat$SingleR_Lineage, 
                                            sep = " ")
saveRDS(object = pred, file = "../../04_results/singler/singler-remove-dcafs-annotate-dividing-output.Rds")
saveRDS(object = sce_query, file = "../../04_results/singler/sce_object_for_singleR_remove-dcafs-annotate-dividing.Rds")
saveRDS(object = integrated_seurat, file = "../../04_results/singler/integrated_seurat_annotated_singleR-20251202.Rds")
