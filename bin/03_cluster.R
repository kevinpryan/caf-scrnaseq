#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
install.packages("clustree")
install.packages("devtools")
library(clustree)
library(future)  # For parallelization
library(data.table)
library(ggthemes)
library(RColorBrewer)
library(SingleCellExperiment)
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2) # Set max size to 8GB

parser$add_argument('--input_rds', type='character', required=TRUE, help='Path for the input integrated RDS file')

integrated_seurat <- readRDS(parser$input_rds)
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:11)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:11)
integrated_seurat <- FindClusters(integrated_seurat, resolution = seq(from = 0.1, to = 2, by = 0.1))
p <- clustree(integrated_seurat, prefix = "integrated_snn_res.", exprs = "scale.data")
pdf("clustree_plot.pdf")
plot(p)
dev.off()
Idents(integrated_seurat) <- integrated_seurat$integrated_snn_res.0.5
integrated_seurat <- DietSeurat(
  integrated_seurat,
  assays = c("RNA", "integrated"), # Keep raw counts and integrated data
  dimreducs = c("pca", "umap")     # Keep your reductions
)
saveRDS(integrated_seurat, file = "clustered_integrated_seurat_dietseurat.rds")

all_markers <- FindAllMarkers(
  integrated_seurat,
  assay = "integrated",         # Use the integrated assay
  only.pos = TRUE,
  # only test genes expressed in this percent of cells
  min.pct = 0.1,
  # stick to default lfc threshold
  logfc.threshold = 0.1
)
# use fwrite to write the markers to a file
fwrite(all_markers, file = "all_markers.csv", sep = ",", quote = F, row.names = F)
sigs <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# use fwrite to write the significant markers to a file
fwrite(sigs, file = "significant_markers.csv", sep = ",", quote = F, row.names = F)
sigs_padj <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.1)
# use fwrite to write the significant markers with adjusted p-value to a file
fwrite(sigs_padj, file = "significant_markers_padj.csv", sep = ",", quote = F, row.names = F)

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
  ungroup()

n_clusters <- 12

color_generator <- colorRampPalette(c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2")))

# Generate exactly 12 unique colors from this function
my_colors <- color_generator(n_clusters)
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
ggsave(plt, file = "cellular-composition-by-sample-integrated_snn_res.0.5.png")
ggsave(plt, file = "cellular-composition-by-sample-integrated_snn_res.0.5.svg")
ggsave(plt, file = "cellular-composition-by-sample-integrated_snn_res.0.5.pdf")

# convert integrated_seurate to SingleCellExperiment
integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
raw_counts <- GetAssayData(integrated_seurat, assay = "RNA", layer = "counts")
norm_integrated_data <- GetAssayData(integrated_seurat, assay = "integrated", layer = "data")
cell_metadata <- integrated_seurat@meta.data

common_features <- rownames(norm_integrated_data)
common_cells <- colnames(norm_integrated_data)
raw_counts_matched <- raw_counts[common_features, common_cells]
cat("Verifying matrix dimensions:\n")
cat("Raw Counts (Matched): "); print(dim(raw_counts_matched))
cat("Integrated Data: "); print(dim(norm_integrated_data))
cell_metadata <- integrated_seurat@meta.data[common_cells, ] # Ensure metadata also matches

sce_object <- SingleCellExperiment(
  assays = list(
    counts = raw_counts_matched,      # Use the matched matrix
    logcounts = norm_integrated_data
  ),
  colData = cell_metadata
)

saveRDS(sce_object, file = "clustered_integrated_sce_dietseurat.rds")
