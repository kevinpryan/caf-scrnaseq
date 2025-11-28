library(Seurat)
library(dplyr)
install.packages("clustree")
install.packages("devtools")
library(clustree)
library(future)  # For parallelization
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2) # Set max size to 8GB
library(presto)
#install.packages("devtools")
devtools::install_github('immunogenomics/presto')
# read in data
integrated_seurat <- readRDS("../results/02_integration/integrated_seurat_object.Rds")
# from the elbow plot it looks like the elbow is at 11
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:11)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:11)
VizDimLoadings(integrated_seurat, dims = 1:2, reduction = "pca")
DimPlot(integrated_seurat, reduction = "pca", group.by = "orig.ident") 
DimHeatmap(integrated_seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(integrated_seurat, dims = 1:11, cells = 500, balanced = TRUE)
DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated_seurat, reduction = "umap")
integrated_seurat <- FindClusters(integrated_seurat, resolution = seq(from = 0.1, to = 2, by = 0.1))
p <- clustree(integrated_seurat, prefix = "integrated_snn_res.", exprs = "scale.data")
plot(p)
ggsave(p, file = "../04_results/clustree-integrated_snn_res.0.5.png")
ggsave(p, file = "../04_results/clustree-integrated_snn_res.0.5.svg")
ggsave(p, file = "../04_results/clustree-integrated_snn_res.0.5.pdf")

# value of 0.5 gives same number of clusters as Cords breast, start there
Idents(integrated_seurat) <- integrated_seurat$integrated_snn_res.0.5
DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.5")
DimPlot(integrated_seurat, reduction = "pca", group.by = "integrated_snn_res.0.5")

p2 <- DimPlot(
  integrated_seurat,
  reduction = "umap",
  group.by = c("integrated_snn_res.0.5", "orig.ident"),
  combine = TRUE, label.size = 2
)
p2
#wrap_plots(c(p1, p2), ncol = 2, byrow = F)
#saveRDS(integrated_seurat, file = "../03_processed_data/03-clustering/integreated_seurat_object_clustered.Rds")
#integrated_seurat <- readRDS("../03_processed_data/03-clustering/integrated_seurat_object_clustered.Rds")
integrated_seurat <- DietSeurat(
  integrated_seurat,
  assays = c("RNA", "integrated"), # Keep raw counts and integrated data
  dimreducs = c("pca", "umap")     # Keep your reductions
)
saveRDS(integrated_seurat, file = "../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat.Rds")

integrated_seurat$integrated_snn_res.
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

sigs <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

sigs_padj <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.1)
fwrite(sigs_padj, file = "../04_results/sigs-clustering-integrated_snn_res.0.5_lfc_1_padj_0.1.csv", sep = ",", quote = F, row.names = F)
table(Idents(integrated_seurat))

table(sigs_padj$gene)

head(integrated_seurat@meta.data)

metadata <- integrated_seurat@meta.data
metadata$barcode <- rownames(metadata)
cells_per_cluster <- metadata %>% group_by(integrated_snn_res.0.5, orig.ident) %>% summarise(ncells = n())
head(cells_per_cluster)
proportions_df <- cells_per_cluster %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells_in_sample = sum(ncells),
    proportion = ncells / total_cells_in_sample
  ) %>%
  ungroup() # It's good practice to ungroup after the operation
library(data.table)
#fwrite(proportions_df, file = "../04_results/cellular-composition-by-sample-integrated_snn_res.0.5_table.txt", sep = "\t", row.names = F, quote = F)
#fwrite(cells_per_cluster, file = "../04_results/cells-per-cluster-by-sample-integrated_snn_res.0.5_table.txt", sep = "\t", row.names = F, quote = F)

install.packages("ggthemes")
library(ggthemes)
library(RColorBrewer)
n_clusters <- 12

#my_colors <- brewer.pal(n = n_clusters, name = "Set2")
#color_generator <- colorRampPalette(brewer.pal(8, "Set2") + brewer.pal(8, "Dark2"))
color_generator <- colorRampPalette(c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2")))

# Generate exactly 12 unique colors from this function
my_colors <- color_generator(n_clusters)

plt <- ggplot(proportions_df, aes(x = orig.ident, y = proportion, fill = `integrated_snn_res.0.5`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Cluster",
    title = "Cellular Composition by Sample"
  ) +
  theme_classic() +
  #scale_fill_tableau() +  
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = my_colors) + # Use the generated palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for 
ggsave(plt, file = "../04_results/cellular-composition-by-sample-integrated_snn_res.0.5.png")
install.packages("svglite")
ggsave(plt, file = "../04_results/cellular-composition-by-sample-integrated_snn_res.0.5.svg")
ggsave(plt, file = "../04_results/cellular-composition-by-sample-integrated_snn_res.0.5.pdf")
