BiocManager::install(c("sp", "scran", "scrapper"))
#BiocManager::install("scran")
#BiocManager::install("scrapper")
#install.packages("vegan")
#install.packages("spam")
install.packages(c("spam", "vegan"))
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
# read in single cell data
integrated_seurat <- readRDS("../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat.Rds")
cat("Class of the loaded object is:", class(integrated_seurat), "\n")

integrated_seurat <- JoinLayers(integrated_seurat, assay = "RNA")
integrated_seurat <- readRDS("../03_processed_data/03-clustering/integrated_seurat_object_LAYERS_JOINED.Rds")

# seu_v4_compat <- JoinLayers(
#   integrated_seurat,
#   assays = c("RNA", "integrated")
# )
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
#seu_v4_compat <- JoinLayers(integrated_seurat)
#ce <- as.SingleCellExperiment(seu_v4_compat)

# read in reference data
breast <- as.SingleCellExperiment(readRDS("/data/scRNA-seq-cords/BREAST_fibro_tumour.rds"))
rm(integrated_seurat, raw_counts, raw_counts_matched, norm_integrated_data, cell_metadata, common_features, common_cells)
gc()
pred.grun <- SingleR(test=sce_object, ref=breast, labels=breast$CAFtype, de.method="wilcox")
#pred.grun
saveRDS(object = pred.grun, file = "../04_results/singler/singler-default-params-labels-breast-ref.Rds")
pred.grun <- readRDS("../04_results/singler/singler-default-params-labels-breast-ref.Rds")
#cells_per_cluster <- metadata %>% group_by(integrated_snn_res.0.5, orig.ident) %>% summarise(ncells = n())

table(pred.grun$labels)
pred.grun$delta.next
pred.grun@rownames
pred.grun$labels

scoreheatmap <- plotScoreHeatmap(pred.grun)
ggsave(scoreheatmap, file = "../04_results/singler/plot-score-heatmap.pdf")
ggsave(scoreheatmap, file = "../04_results/singler/plot-score-heatmap.svg")

plotDeltaDistribution(pred.grun, ncol = 3)
summary(is.na(pred.grun$pruned.labels))

# log-transformed counts??
plotMarkerHeatmap(pred.grun, 
                  logcounts(sce_object), 
                  label = "Pericyte"
                  )

all.markers <- metadata(pred.grun)$de.genes
sce_object$labels <- pred.grun$labels
library(scater)

plotHeatmap(sce_object, order_columns_by="labels",
            features=unique(unlist(all.markers$apCAF))) 
plotHeatmap(sce_object, order_columns_by="labels",
            features=unique(unlist(all.markers$dCAF))) 

to.remove.default <- pruneScores(pred.grun)
table(Label=pred.grun$labels, Removed=to.remove.default)

to.remove2 <- pruneScores(pred.grun, min.diff.next=0.1)
table(Label=pred.grun$labels, Removed=to.remove2)
barcode_label_prune <- data.frame(Barcode = rownames(pred.grun), CAFType_SingleR = pred.grun$labels, ToRemove = to.remove2)

pred.grun@metadata

#log_transformed_counts <- GetAssayData(sce_object,)
#stopifnot(length(pred.grun$labels) == length(pred.grun@rownames))
output_celltypes <- data.frame(Barcode = rownames(pred.grun), CAFType_SingleR = pred.grun$labels)
original_coldata <- as.data.frame(sce_object@colData)
original_coldata$Barcode <- rownames(original_coldata)

output_celltypes_df <- left_join(output_celltypes, as.data.frame(original_coldata))
output_celltypes_df_alldata <- left_join(original_coldata, as.data.frame(barcode_label_prune))
original_coldata
cells_per_celltype <- output_celltypes_df_alldata %>% group_by(CAFType_SingleR, orig.ident) %>% summarise(ncells = n())

prune.per.sample.df <- output_celltypes_df_alldata %>% 
                       group_by(orig.ident, CAFType_SingleR, ToRemove) %>%
                       summarise(N_to_remove_per_cell_type_sample = n()) %>% 
                      dplyr::filter(ToRemove == TRUE) %>% 
                       ungroup() 
              
                              #,
                              #total_cells_in_sample = sum(ncells),
                              #proportion_to_remove = N_To_Remove/total_cells_in_sample) 
                              #%>% 
                       #summarise(to_remove = n()) %>% 

prune.per.sample.df.proprotion <- prune.per.sample.df %>% 
                                  full_join(cells_per_celltype) %>% 
                                  mutate(proprotion_pruned = N_to_remove_per_cell_type_sample/ncells)
    
proportions_df <- cells_per_celltype %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells_in_sample = sum(ncells),
    proportion = ncells / total_cells_in_sample
  ) %>%
  ungroup()

counts_wide <- proportions_df %>%
  dplyr::select(orig.ident, CAFType_SingleR, ncells) %>%
  spread(key = CAFType_SingleR, value = ncells, fill = 0) %>%
  # The first column is the sample ID, so we need to move it to rownames
  tibble::column_to_rownames("orig.ident")

shannon_diversity <- diversity(counts_wide, index = "shannon")
simpson_diversity <- diversity(counts_wide, index = "invsimpson") # Use Inverse Simpson

diversity_summary <- data.frame(
  orig.ident = names(shannon_diversity),
  shannon = shannon_diversity,
  inv_simpson = simpson_diversity
)
fwrite(x = diversity_summary, file = "../04_results/singler/diversity-summary.csv", sep = ",", quote = F, row.names = F)

diversity_summary <- read.csv("../04_results/singler/diversity-summary.csv")
diversity_summary_long <- diversity_summary %>% pivot_longer(names_to = "diversity_metric", cols = c("shannon", "inv_simpson"))
diversity <- ggplot(diversity_summary_long, aes(fill=diversity_metric, y=value, x=orig.ident)) + 
  geom_bar(position="dodge", stat="identity") +
  ylab("Diversity metric value") +
  xlab("Sample") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("CAF subpopulation diversity of samples") +
  theme_classic() 
plot(diversity)
ggsave(plot = diversity, filename = "../04_results/singler/cell-type-diversity-per-sample.png")  
ggsave(plot = diversity, filename = "../04_results/singler/cell-type-diversity-per-sample.pdf")  
ggsave(plot = diversity, filename = "../04_results/singler/cell-type-diversity-per-sample.svg") 
n_clusters <- 10

#my_colors <- brewer.pal(n = n_clusters, name = "Set2")
#color_generator <- colorRampPalette(brewer.pal(8, "Set2") + brewer.pal(8, "Dark2"))
color_generator <- colorRampPalette(c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2")))

# Generate exactly 12 unique colors from this function
my_colors <- color_generator(n_clusters)

proportions_df$CAFType_SingleR <- gsub("tpCAF", "tCAF", proportions_df$CAFType_SingleR)
proportions_df$CAFType_SingleR <- gsub("IDO_CAF", "ifnCAF", proportions_df$CAFType_SingleR)
proportions_df$CAFType_SingleR <- gsub("tpCAF", "tCAF", proportions_df$CAFType_SingleR)

desired_order <- c("mCAF", "iCAF", "vCAF", "Pericyte", "tCAF", "hsp_tCAF", "ifnCAF", "apCAF", "rCAF", "dCAF")
proportions_df$CAFType_SingleR <- as.factor(proportions_df$CAFType_SingleR)

proportions_df <- 
  proportions_df %>% 
  mutate(
    CAFType_SingleR = fct_relevel(
      CAFType_SingleR, desired_order
    )
  )

plt <- ggplot(proportions_df, aes(x = orig.ident, y = proportion, fill = `CAFType_SingleR`)) +
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
ggsave(plt, file = "../04_results/singler/cellular-composition-by-sample-labelled.png")
ggsave(plt, file = "../04_results/singler/cellular-composition-by-sample-labelled.svg")
ggsave(plt, file = "../04_results/singler/cellular-composition-by-sample-labelled.pdf")

pred.grun

#read in integrated_seurat again
integrated_seurat <- readRDS("../03_processed_data/03-clustering/integrated_seurat_object_clustered_dietseurat.Rds")
integrated_seurat@meta.data
output_celltypes_df_to_add <- output_celltypes_df
rownames(output_celltypes_df_to_add) <- output_celltypes_df_to_add$Barcode
output_celltypes_df_to_add <- output_celltypes_df_to_add %>% dplyr::select(CAFType_SingleR)
integrated_seurat <- AddMetaData(integrated_seurat, output_celltypes_df_to_add)
gene <- "PBK"
expr_values <- FetchData(integrated_seurat, vars = c(gene, "CAFType_SingleR"))
expr_values <- expr_values %>%
  mutate(expressed = !!sym(gene) > 0)
proportions <- expr_values %>%
  group_by(CAFType_SingleR) %>%
  summarise(
    n_cells = n(),
    n_expressing = sum(expressed),
    prop_expressing = n_expressing / n_cells
  ) %>%
  arrange(desc(prop_expressing))
ggplot(proportions, aes(x = reorder(CAFType_SingleR, prop_expressing), y = prop_expressing)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = paste0("Proportion expressing PBK"),
    title = paste0("Expression of PBK across cell types")
  ) +
  theme_minimal()

DimPlot(integrated_seurat, reduction = "umap", group.by = "CAFType_SingleR")
DimPlot(integrated_seurat, reduction = "pca", group.by = "CAFType_SingleR") 

