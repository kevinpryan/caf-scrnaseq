# read in cells per cluster data
data_in <- read.table("../04_results/cells-per-cluster-by-sample-integrated_snn_res.0.5_table.txt", header = T)
counts_wide <- data_in %>%
  dplyr::select(orig.ident, integrated_snn_res.0.5, ncells) %>%
  spread(key = integrated_snn_res.0.5, value = ncells, fill = 0) %>%
  # The first column is the sample ID, so we need to move it to rownames
  tibble::column_to_rownames("orig.ident")
shannon_diversity <- diversity(counts_wide, index = "shannon")
simpson_diversity <- diversity(counts_wide, index = "invsimpson") # Use Inverse Simpson

diversity_summary <- data.frame(
  orig.ident = names(shannon_diversity),
  shannon = shannon_diversity,
  inv_simpson = simpson_diversity
)
diversity_summary_long <- diversity_summary %>% pivot_longer(names_to = "diversity_metric", cols = c("shannon", "inv_simpson"))
diversity <- ggplot(diversity_summary_long, aes(fill=diversity_metric, y=value, x=orig.ident)) + 
  geom_bar(position="dodge", stat="identity") +
  ylab("Diversity metric value") +
  xlab("Sample") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Cluster diversity of samples") +
  theme_classic() 
plot(diversity)
ggsave(plot = diversity, filename = "../04_results/cluster-diversity-per-sample.png")  
ggsave(plot = diversity, filename = "../04_results/cluster-diversity-per-sample.pdf")  
ggsave(plot = diversity, filename = "../04_results/cluster-diversity-per-sample.svg")  


fwrite(x = diversity_summary, file = "../04_results/diversity-summary-clustering.csv", sep = ",", quote = F, row.names = F)
