library(Seurat)
library(dplyr)
integrated_seurat <- readRDS("../../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat_redo-integration-20251128.Rds")

orig.ident <- c("4027", "4299", "DF_CAF", "PC_CAF", "RM_CAF")
cancer.type <- c("breast", "breast", "lung", "lung", "lung")
recovery <- c("high", "high", "high", "low", "low")
metadata <- data.frame(orig.ident, cancer.type, recovery)
metadata_full <- left_join(integrated_seurat@meta.data, metadata)
integrated_seurat <- AddMetaData(integrated_seurat, metadata = metadata_full)

p2 <- DimPlot(
  integrated_seurat,
  reduction = "umap",
  group.by = c("integrated_snn_res.0.4", "orig.ident", "cancer.type", "recovery"),
  combine = TRUE, label.size = 2
)
p2

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
integrated_seurat_markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE)
all_markers <- FindAllMarkers(
  integrated_seurat,
  assay = "integrated",         # Use the integrated assay
  only.pos = TRUE,
  # only test genes expressed in this percent of cells
  min.pct = 0.1,
  # stick to default lfc threshold
  logfc.threshold = 0.1
)

sigs <- integrated_seurat_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

sigs_padj <- integrated_seurat_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.1)
#fwrite(sigs_padj, file = "../04_results/sigs-clustering-integrated_snn_res.0.5_lfc_1_padj_0.1.csv", sep = ",", quote = F, row.names = F)
table(Idents(integrated_seurat))

# look at marker expression across clusters
gene <- "FAP"
expr_values <- FetchData(integrated_seurat, vars = c(gene, "integrated_snn_res.0.4"))
expr_values <- expr_values %>%
  mutate(expressed = !!sym(gene) > 0)
proportions <- expr_values %>%
  group_by(integrated_snn_res.0.4) %>%
  summarise(
    n_cells = n(),
    n_expressing = sum(expressed),
    prop_expressing = n_expressing / n_cells
  ) %>%
  arrange(desc(prop_expressing))
ggplot(proportions, aes(x = reorder(integrated_snn_res.0.4, prop_expressing), y = prop_expressing)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = paste0("Proportion expressing FAP"),
    title = paste0("Expression of FAP across cell types")
  ) +
  theme_minimal()



