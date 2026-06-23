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

# --- Load Data ---
integrated_seurat <- readRDS("../../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat_redo-integration-20251128.Rds")
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
# Load Reference
breast <- as.SingleCellExperiment(readRDS("/data/scRNA-seq-cords/BREAST_fibro_tumour.rds"))
breast <- logNormCounts(breast)

# ==============================================================================
# DEFINING THE "BLACKLIST" (Artifact Removal)
# ==============================================================================

# 1. Cell Cycle Genes (Seurat Standard)
# Removes the "dCAF" proliferation signal
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

# 2. Ribosomal Genes
# Removes the "High Protein Synthesis" culture artifact
# Matches RPL..., RPS..., MRPL...
ribo_genes <- grep("^(RP[LS]|MRP[LS])", rownames(integrated_seurat), value = TRUE)

# 3. Mitochondrial Genes (Optional, but good practice)
mt_genes <- grep("^MT-", rownames(integrated_seurat), value = TRUE)

# 4. Hypoxia & Glycolysis (The "tpCAF" Artifact)
# Based on your list + Buffa et al. signature
hypoxia_genes <- c(
  "ENO1", "ENO2", "GAPDH", "PGK1", "ALDOA", "LDHA", "CA9", "VEGFA", 
  "HK2", "PKM", "BNIP3", "DDIT4", "TPI1", "P4HA1", "ADM", "NDRG1"
)

# 5. Stress / Heat Shock / Dissociation (The "hsp_tpCAF" Artifact)
# A. Heat Shock Proteins (Regex for HSPs and DNAJ co-chaperones)
hsp_genes <- grep("^(HSP|DNAJ|CRYAB)", rownames(integrated_seurat), value = TRUE)

# B. Immediate Early Genes (Dissociation stress - van den Brink et al 2017)
ieg_genes <- c(
  "FOS", "FOSB", "JUN", "JUND", "JUNB", "EGR1", "EGR3", 
  "NR4A1", "NR4A2", "NR4A3", "ATF3", "DUSP1", "IER2", "IER3", "BTG2"
)

# --- COMBINE ALL EXCLUSIONS ---
blacklist <- unique(c(s_genes, g2m_genes, ribo_genes, mt_genes, hypoxia_genes, hsp_genes, ieg_genes))

#====================================================================
# GENERATE SCE_QUERY OBJECT AND CALCULATE CELL CYCLE SCORES
#====================================================================

DefaultAssay(integrated_seurat) <- "RNA"
integrated_seurat <- NormalizeData(integrated_seurat)
integrated_seurat <- CellCycleScoring(integrated_seurat, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
integrated_seurat$Cycle_Status <- ifelse(integrated_seurat$Phase == "G1", "Resting", "Dividing")
# C. Prepare Query SCE (Raw counts -> LogNorm)
sce_query <- SingleCellExperiment(
  assays = list(counts = GetAssayData(integrated_seurat, assay = "RNA", layer = "counts")),
  colData = integrated_seurat@meta.data
)
sce_query <- logNormCounts(sce_query)

# ==============================================================================
# FILTER & RUN SINGLER
# ==============================================================================

# 1. Intersect genes
common_genes <- intersect(rownames(sce_query), rownames(breast))

# 2. Remove the Blacklist
genes_to_use <- setdiff(common_genes, blacklist)

cat("Total genes used for classification:", length(genes_to_use), "\n")

# 3. Clean Reference
# Remove dCAF (dividing) and hsp_tpCAF/tpCAF (if you consider them purely pathological states)
# However, filtering the GENES is usually sufficient. 
# Removing the reference labels is an extra safety step.
types_to_remove <- c("dCAF", "hsp_tpCAF", "tpCAF") 
breast_clean <- breast[, !breast$CAFtype %in% types_to_remove]
breast_clean$CAFtype <- droplevels(breast_clean$CAFtype)

# 4. Run SingleR
rm(breast)
gc()
pred <- SingleR(test = sce_query[genes_to_use, ], 
                ref = breast_clean[genes_to_use, ], 
                labels = breast_clean$CAFtype, 
                de.method = "wilcox")

# View Result
table(pred$labels)

integrated_seurat$SingleR_Lineage <- pred$labels

# ==============================================================================
# 5. CREATE COMPOSITE LABELS
# ==============================================================================
# Combine the two columns
integrated_seurat$Final_Annotation <- paste(integrated_seurat$Cycle_Status, 
                                            integrated_seurat$SingleR_Lineage, 
                                            sep = " ")

saveRDS(object = pred, file = "../../04_results/singler/singler-remove-culture-artifacts-output.Rds")
saveRDS(object = sce_query, file = "../../04_results/singler/sce_object_for_singleR_remove-artifacts.Rds")
saveRDS(object = integrated_seurat, file = "../../04_results/singler/integrated_seurat_annotated_singleR-remove-artifacts-20251202.Rds")
