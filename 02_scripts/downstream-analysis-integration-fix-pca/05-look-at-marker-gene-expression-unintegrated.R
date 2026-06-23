samplesheet <- read.csv("../../01_metadata/samplesheet.csv")
files <- list.files(path = "../../results-updated-20251128/01_normalized_rds", full.names = T)
names(files) <- c("4027", "4299", "DF_CAF", "PC_CAF", "RM_CAF")
files
seurat_list <- lapply(seq_along(files), function(i) {
  obj <- readRDS(files[i])
  return(obj)
})
merged_obj <- merge(x = seurat_list[[1]], y = list(seurat_list[[2]], seurat_list[[3]], seurat_list[[4]], seurat_list[[5]]))

plot_expression_across_samples <- function(gene, obj) {
  expr_values <- FetchData(obj, vars = c(gene, "integrated_snn_res.0.4"), layer = "RNA")
  ylabel <- paste("Proportion expressing", gene, sep = " ")
  title_plot <- paste("Expression of", gene, "across clusters", sep = " ")
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
      y = paste0(ylabel),
      title = paste0(title_plot)
    ) +
    theme_minimal()
  
}
plot_expression_across_clusters("ACTA2", integrated_seurat)
plot_expression_across_clusters("VIM", integrated_seurat)
plot_expression_across_clusters("PDPN", integrated_seurat)