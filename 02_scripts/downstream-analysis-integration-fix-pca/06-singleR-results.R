library(SingleR)
library(Seurat)
library(ggplot2)
# all genes log normalised un-integrated counts
pred.grun <- readRDS("../../04_results/singler/singler-default-params-labels-breast-ref-use-all-genes-norm-reference.Rds")
scoreheatmap <- plotScoreHeatmap(pred.grun)
scoreheatmap
pred.grun
table(pred.grun$labels)
plotDeltaDistribution(pred.grun, ncol = 3)

# highly variable genes used in data integration
pred.grun2 <- readRDS("../../04_results/singler/singler-default-params-labels-breast-ref-use-variable-genes-integrated.Rd")
sce_object <- readRDS("../../04_results/singler/sce_object_for_singleR_integrated_data.Rds")
table(pred.grun2$labels)
scoreheatmap2 <- plotScoreHeatmap(pred.grun2)
scoreheatmap2
plotDeltaDistribution(pred.grun2, ncol = 3)
plotMarkerHeatmap(pred.grun2, sce_object, label="rCAF")



# DE genes Cords et al used in data integration
# allowed SingleR to find DE genes
pred.grun3 <- readRDS("../../04_results/singler/singler-default-params-labels-breast-ref-use-de-genes-supp3-cords.Rds")
sce_object <- readRDS("../../04_results/singler/sce_object_for_singleR_subset_de_genes_cords_supp3.Rds")
table(pred.grun3$labels)
scoreheatmap3 <- plotScoreHeatmap(pred.grun3)
plotMarkerHeatmap(pred.grun3, sce_object, label="apCAF")
plotDeltaDistribution(pred.grun3, ncol = 3)

# DE genes Cords et al used in data integration
# specified which CAF subset DE genes belonged to
pred.grun4 <- readRDS("../../04_results/singler/singler-default-params-labels-breast-ref-use-de-genes-supp3-cords-specify-subpops.Rds")
sce_object <- readRDS("../../04_results/singler/sce_object_for_singleR_all_genes.Rds")
table(pred.grun4$labels)
scoreheatmap4 <- plotScoreHeatmap(pred.grun4)
plotMarkerHeatmap(pred.grun4, sce_object, label="apCAF")
plotDeltaDistribution(pred.grun4, ncol = 3)

# remove dCAFs from reference and remove cell division-related genes when doing annotation
# label cells as dividing or not dividing
pred <- readRDS("../../04_results/singler/singler-remove-dcafs-annotate-dividing-output.Rds")
table(pred$labels)
seurat_obj <- readRDS("../../04_results/singler/integrated_seurat_annotated_singleR-20251202.Rds")
table(seurat_obj@meta.data$Final_Annotation)
scoreheatmap5 <- plotScoreHeatmap(pred)

# remove dCAFs, tpCAFs, hsp_tpCAFs from reference
pred2 <- readRDS("../../04_results/singler/singler-remove-culture-artifacts-output.Rds")
table(pred2$labels)
seurat_obj <- readRDS("../../04_results/singler/integrated_seurat_annotated_singleR-remove-artifacts-20251202.Rds")
sce_object <- readRDS("../../04_results/singler/sce_object_for_singleR_remove-artifacts.Rds")
table(seurat_obj@meta.data$Final_Annotation)
scoreheatmap6 <- plotScoreHeatmap(pred2)
plotMarkerHeatmap(pred2, sce_object, label="apCAF")

# Check expression of the defining apCAF markers in your query object
# (Assuming your Seurat object is named 'seurat_obj')

markers <- c("CD74", "HLA-DRA", "HLA-DRB1", "MIF")

# DotPlot to see who is expressing what
DotPlot(seurat_obj, features = markers, group.by = "SingleR_Lineage") + 
  RotatedAxis()

# Check if your cells are actually fibroblasts (they should be)
# and which direction they lean
features <- c(
  "COL1A1", "PDGFRA", "FAP",  # Pan-CAF markers (Should be High)
  "POSTN", "MMP11", "LRRC15", # mCAF markers
  "CFD", "C3", "CXCL12",       # iCAF markers
  "PDGFRB", "RGS5", "CSPG4" # pericyte markers
)

DotPlot(seurat_obj, features = features, group.by = "SingleR_Lineage") + 
  RotatedAxis()

# Define Gene Signatures (using lists effectively averages out dropout noise)
mCAF_sig <- list(c("MMP11", "POSTN", "COL1A1", "COL1A2", "LRRC15", "FN1"))
iCAF_sig <- list(c("CXCL12", "CFD", "C3", "IL6", "CXCL14", "DCN"))

# Score the cells
seurat_obj <- AddModuleScore(seurat_obj, features = mCAF_sig, name = "mCAF_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = iCAF_sig, name = "iCAF_Score")

# Visualize
FeaturePlot(seurat_obj, features = c("mCAF_Score1", "iCAF_Score1"))
VlnPlot(seurat_obj, features = c("mCAF_Score1", "iCAF_Score1"), group.by = "SingleR_Lineage")

features_check <- c(
  # 1. Pan-Fibroblast markers (Should be present in iCAF, even if lower than mCAF)
  "PDGFRA", "VIM", "FAP", 
  
  # 2. iCAF specific markers (Should be High in iCAF)
  "CXCL12", "CFD", "C3", "IL6",
  
  # 3. Contaminants (Should be ZERO)
  "PTPRC",  # CD45 (Immune cells)
  "EPCAM",  # Epithelial/Tumor cells
  "PECAM1"  # CD31 (Endothelial cells)
)

p <- DotPlot(seurat_obj, 
             features = c("COL1A1", "PDGFRA", "CXCL12", "IL6"), 
             group.by = "SingleR_Lineage")

# 2. Extract the underlying data frame
plot_data <- p$data

# 3. View the first few rows
head(plot_data)
percent_table <- plot_data %>%
  select(features.plot, id, pct.exp) %>%
  pivot_wider(names_from = id, values_from = pct.exp)


# remove dCAFs, apCAFs, hsp_tpCAFs, tpCAFs

pred3 <- readRDS("../../04_results/singler/singler-remove-all-culture-artifacts-output-remove-apCAF.Rds")
table(pred3$labels)
seurat_obj <- readRDS("../../04_results/singler/integrated_seurat_annotated_singleR-remove-all-artifacts-remove-apCAFs-20251202.Rds")
sce_object <- readRDS("../../04_results/singler/sce_object_for_singleR_remove-all-artifacts-remove-apCAF.Rds")
# Visualize CXCL12 on your UMAP
FeaturePlot(seurat_obj, features = "CXCL12", order = TRUE) 
# Visualize it specifically within the mCAF population
# This plot will likely show a "tail" or sub-blob of mCAFs that are CXCL12+
p_cxcl12 <- VlnPlot(seurat_obj, features = "CXCL12", group.by = "SingleR_Lineage")
VlnPlot(seurat_obj, features = "CXCL12", group.by = "SingleR_Lineage")
p_cxcl12_data <- p_cxcl12$data
library(dplyr)
p_cxcl12_data_mCAF <- p_cxcl12_data %>% dplyr::filter(ident == "mCAF")
p_cxcl12_data_mCAF$is_expressed <- ifelse(p_cxcl12_data_mCAF$CXCL12 > 0, "yes", "no")
table(p_cxcl12_data_mCAF$is_expressed)
mCAF_cells <- WhichCells(seurat_obj, expression = SingleR_Lineage == "mCAF")
cxcl12_pos_cells <- WhichCells(seurat_obj, expression = CXCL12 > 0, cells = mCAF_cells, slot = "counts")
seurat_obj$Refined_Annotation <- seurat_obj$SingleR_Lineage
seurat_obj$Refined_Annotation[cxcl12_pos_cells] <- "iCAF"
table(seurat_obj$Refined_Annotation)

# use module score method to refine annotation
# Define a "Culture-Adapted" iCAF signature
# Focus on the chemokines that persist in culture (CXCL12, IL6, CXCL14)
# Exclude C3/CFD
culture_iCAF_genes <- list(c("CXCL12", "IL6", "CXCL14", "LIF", "IL11"))

# Score the cells
seurat_obj <- AddModuleScore(seurat_obj, features = culture_iCAF_genes, name = "Culture_iCAF_Score")

# Find mCAFs with high iCAF scores
# We pick a threshold based on the distribution (e.g., top 20% or > 0)
scores <- seurat_obj$Culture_iCAF_Score1
hist(scores)
mCAF_indices <- which(seurat_obj$SingleR_Lineage == "mCAF")
# from violin plot, there is a separation at zero, call anything <= 0 mCCAF, anything > 0 iCAF
VlnPlot(seurat_obj, features = "Culture_iCAF_Score1", group.by = "SingleR_Lineage") & geom_hline(yintercept = 0)
icaf_score_data <- VlnPlot(seurat_obj, features = "Culture_iCAF_Score1", group.by = "SingleR_Lineage")
icaf_score_data <- icaf_score_data$data
icaf_score_data_mcaf <- icaf_score_data %>% dplyr::filter(ident == "mCAF")

high_score_indices <- intersect(mCAF_indices, which(scores > 0)) # based on VlnPlot

# Re-label
cells_to_rescue <- colnames(seurat_obj)[high_score_indices]
seurat_obj$Refined_Annotation[cells_to_rescue] <- "iCAF"
VlnPlot(seurat_obj, features = "Culture_iCAF_Score1", group.by = "Refined_Annotation")
DimPlot(seurat_obj, group.by = "Refined_Annotation", label = TRUE) 


# Plot 1: Original SingleR calls (mostly mCAF)
p1 <- DimPlot(seurat_obj, 
              group.by = "SingleR_Lineage", 
              label = TRUE
              ) + 
  ggtitle("Original SingleR Labels") #+ NoLegend()

# Plot 2: Your Refined Annotation (iCAF rescued)
p2 <- DimPlot(seurat_obj, group.by = "Refined_Annotation", label = TRUE) + 
  ggtitle("Refined Annotation (CXCL12+ Gated)") # + NoLegend()

# Combine them side-by-side
p1 + p2
"CXCL12" %in% VariableFeatures(seurat_obj)
# This looks much cleaner than a UMAP for shallow data
VlnPlot(seurat_obj, 
        features = c("CXCL12", "IL6", "MMP11", "COL1A1"), 
        group.by = "Refined_Annotation", 
        pt.size = 0, # Hide the messy dots
        ncol = 2)
# possibly missing IL6 and MMP11 due to lack of depth
# try alternative markers mCAF
VlnPlot(seurat_obj, 
        features = c("FN1", "POSTN", "COL1A2", "LRRC15"), 
        group.by = "Refined_Annotation", 
        pt.size = 0, # Hide the messy dots
        ncol = 2)
# try alternative markers iCAF
VlnPlot(seurat_obj, 
        features = c("CXCL14", "DCN", "LIF", "CCL2"), 
        group.by = "Refined_Annotation", 
        pt.size = 0, # Hide the messy dots
        ncol = 2)
# The "Backup" Panel
backup_markers <- c(
  # The Splitter (Confirmed)
  "CXCL12", 
  
  # Backup iCAF (Try to find one that lights up)
  "CXCL14", "DCN", "LIF", "CCL2",
  
  # Backup mCAF (Look for higher intensity here vs iCAF)
  "FN1", "POSTN", "LRRC15", "ACTA2"
)

DotPlot(seurat_obj, features = backup_markers, group.by = "Refined_Annotation") + 
  RotatedAxis()

# expression of DCN seems higher in iCAF vs mCAF, check for significance:
# 1. Set the identity to your new annotations
Idents(seurat_obj) <- "Refined_Annotation"

# 2. Run the test specifically for DCN
# ident.1 is the group you are testing FOR (iCAF)
# ident.2 is the group you are comparing AGAINST (mCAF)
# assay = "RNA" ensures you use the normalized data, not integrated residuals
stat_result <- FindMarkers(seurat_obj, 
                           ident.1 = "iCAF", 
                           ident.2 = "mCAF", 
                           features = c("DCN"),
                           logfc.threshold = 0,
                           min.pct = 0,
                           assay = "RNA") # Use RNA assay for DE testing

# 3. View the result
print(stat_result)
#            p_val avg_log2FC pct.1 pct.2    p_val_adj
# DCN 1.145057e-83  0.2739568 0.709 0.621 2.307176e-79
FeaturePlot(seurat_obj, features = c("DCN"))

cells_per_celltype <- seurat_obj@meta.data %>% group_by(Refined_Annotation, orig.ident) %>% summarise(ncells = n()) %>% ungroup()
proportions_df <- cells_per_celltype %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells_in_sample = sum(ncells),
    proportion = ncells / total_cells_in_sample
  ) %>%
  ungroup()
library(tidyr)
counts_wide <- proportions_df %>%
  dplyr::select(orig.ident, Refined_Annotation, ncells) %>%
  spread(key = Refined_Annotation, value = ncells, fill = 0) %>%
  # The first column is the sample ID, so we need to move it to rownames
  tibble::column_to_rownames("orig.ident")
library(vegan)
shannon_diversity <- diversity(counts_wide, index = "shannon")
simpson_diversity <- diversity(counts_wide, index = "invsimpson") # Use Inverse Simpson
diversity_summary <- data.frame(
  orig.ident = names(shannon_diversity),
  shannon = shannon_diversity,
  inv_simpson = simpson_diversity
)




desired_order <- c("mCAF", "iCAF", "vCAF", "Pericyte", "rCAF")
proportions_df$Refined_Annotation <- as.factor(proportions_df$Refined_Annotation)

library(forcats)
proportions_df <- 
  proportions_df %>% 
  mutate(
    Refined_Annotation = fct_relevel(
      Refined_Annotation, desired_order
    )
  )

library(ggthemes)
plt <- ggplot(proportions_df, aes(x = orig.ident, y = proportion, fill = `Refined_Annotation`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Cluster",
    title = "Cellular Composition by Sample"
  ) +
  theme_classic() +
  scale_fill_tableau() +  
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_manual(values = my_colors) + # Use the generated palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for 
plt

p_acta2 <- VlnPlot(seurat_obj, features = "ACTA2", group.by = "Refined_Annotation")
p_acta2_data <- p_acta2$data
p_acta2_data$is_acta2_expressed <- ifelse(p_acta2_data$ACTA2 > 0, "yes", "no")
p_acta2_data %>% group_by(ident, is_acta2_expressed) %>% summarise(n = n())

stat_result_acta2 <- FindMarkers(seurat_obj, 
                           ident.1 = "iCAF", 
                           ident.2 = "mCAF", 
                           features = c("ACTA2"),
                           logfc.threshold = 0,
                           min.pct = 0,
                           assay = "RNA") # Use RNA assay for DE testing
