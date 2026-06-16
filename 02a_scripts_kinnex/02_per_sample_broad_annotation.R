library(Seurat)
library(SingleCellExperiment)
library(dplyr)
#install.packages("clustree")
library(clustree)
library(msigdbr)
library(fgsea)

library(UCell)
samplesheet <- read.csv("../01_metadata/samplesheet-kinnex-20260204.csv")
samplesheet$path <- paste("/home/rstudio/caf-scrnaseq", samplesheet$path, sep = "/")
# define signatures
xu_data <- readxl::read_xlsx(path = "../00_raw_data/reference-data/mmc3-xu-et-al-canonical-markers-2024.xlsx")
xu_data_endothelial_up <- xu_data %>% dplyr::filter(`cell type` == "Endothelial cells" & `up/down` == "upregulated") %>% dplyr::select(gene) %>%  pull() %>% paste(., "+", sep = "")
xu_data_epithelial_up <- xu_data %>% dplyr::filter(`cell type` == "Epithelial cells" & `up/down` == "upregulated") %>% dplyr::select(gene) %>%  pull() %>% paste(., "+", sep = "")
xu_data_fibroblast_up <- xu_data %>% dplyr::filter(`cell type` == "Fibroblasts" & `up/down` == "upregulated") %>% dplyr::select(gene) %>%  pull() %>% paste(., "+", sep = "")
xu_data_fibroblast_down <- xu_data %>% dplyr::filter(`cell type` == "Fibroblasts" & `up/down` == "downregulated") %>% dplyr::select(gene) %>%  pull() %>% paste(., "-", sep = "")
xu_data_fibroblast <- c(xu_data_fibroblast_up, xu_data_fibroblast_down)
#xu_data_epithelial_up_exclude_krts <- c("EGFR",  "FZR1", "ITGA6", "TP63", "MME", "FOXA1", "GATA3", "MUC1",  "CD24", "GABRP", "EPCAM")
# immune_genes <- c(
#   "PTPRC",
#   "LST1",
#   "TYROBP",
#   "FCER1G",
#   "HLA-DRA",
#   "CD74"
# )
immune_genes <- c(
  "PTPRC",
  "LYZ",
  "MS4A1",
  "CD3D",
  "CD3E"
)
gene_sets <- list(
  #Fibroblast = xu_data_fibroblast_up,
  Fibroblast = xu_data_fibroblast,
  Epithelial = xu_data_epithelial_up,
  
  Endothelial = xu_data_endothelial_up,

  Immune = immune_genes
)


# read in normalised data
seurat_qc_list <- readRDS(file = "../03_processed_data/kinnex/01_qc_and_filtering/caf_kinnex_seurat_qc_list_sctransform.Rds")
seurat_qc_list
gc()

medians <- lapply(seurat_qc_list, function(x) {
  med <- median(x$nFeature_RNA)
  med
})
medians <- unlist(medians)
mean(medians)
# mean median number of genes is 1737, can stick to default of 1500 - should be the same across samples (see supplementary section 5 UCell and pyUCell: single-cell gene signature scoring for R and Python)
seurat_qc_list <- lapply(seurat_qc_list, function(x) {
  DefaultAssay(x) <- "SCT"
  AddModuleScore_UCell(
    x,
    features = gene_sets,
    name = NULL
  )
})

add_caf_score <- function(obj) {
  
  obj$CAF_score <-
    obj$Fibroblast -
    pmax(
      obj$Epithelial,
      obj$Immune,
      obj$Endothelial
    )
  
  obj$CAF_confident <- 
    obj$CAF_score > 0.1 &
    obj$Fibroblast > obj$Epithelial &
    obj$Fibroblast > obj$Immune &
    obj$Fibroblast > obj$Endothelial 
    # obj$Fibroblast > obj$Epithelial &
    # obj$Fibroblast > obj$Immune &
    # obj$Fibroblast > obj$Endothelial
  
  obj$CAF_low_score <- 
    obj$CAF_score <= 0.1 &
    obj$Fibroblast > obj$Epithelial &
    obj$Fibroblast > obj$Immune &
    obj$Fibroblast > obj$Endothelial 
  
  obj$max_celltype <- apply(
    obj@meta.data[, c(
      "Fibroblast",
      "Epithelial",
      "Immune",
      "Endothelial"
    )],
    1,
    function(x) names(which.max(x))
  )
  
  print(table(obj$max_celltype))
  print(table(obj$CAF_confident))
  return(obj)
}

seurat_qc_list <- lapply(seurat_qc_list, add_caf_score)
names(seurat_qc_list) <- samplesheet$samples
table(seurat_qc_list$`4027`$CAF_low_score)
seurat_qc_list <- readRDS(file = "../03_processed_data/kinnex/01_qc_and_filtering/caf_kinnex_seurat_qc_list_sctransform_caf_scores.Rds")
#seurat_qc_list <- lapply(seurat_qc_list, add_caf_score)
# cluster and look at UMAP per sample based on CAF_confident
gc()
seurat_qc_list <- lapply(seurat_qc_list, function(x) {
  # x <- NormalizeData(x)
  # x <- FindVariableFeatures(x, 
  #                           selection.method = "vst", nfeatures = 500)
  # x <- ScaleData(x)
  RunPCA(x, npcs = 20, features=VariableFeatures(x))
})
  
ElbowPlot(seurat_qc_list$`4027`) 
# 11
ElbowPlot(seurat_qc_list$`4299`) 
# 11
ElbowPlot(seurat_qc_list$DF_CAF) 
# 18
ElbowPlot(seurat_qc_list$RM_CAF) 
# 13
ElbowPlot(seurat_qc_list$PC_CAF) 
# 14

seurat_qc_list$`4027` <- RunUMAP(seurat_qc_list$`4027`, reduction = "pca", 
                          dims = 1:11, seed.use=123)
seurat_qc_list$`4299` <- RunUMAP(seurat_qc_list$`4299`, reduction = "pca", 
                                 dims = 1:11, seed.use=123)
seurat_qc_list$DF_CAF <- RunUMAP(seurat_qc_list$DF_CAF, reduction = "pca", 
                                 dims = 1:18, seed.use=123)
seurat_qc_list$RM_CAF <- RunUMAP(seurat_qc_list$RM_CAF, reduction = "pca", 
                                 dims = 1:13, seed.use=123)
seurat_qc_list$PC_CAF <- RunUMAP(seurat_qc_list$PC_CAF, reduction = "pca", 
                                 dims = 1:14, seed.use=123)

plots <- lapply(names(seurat_qc_list), function(nm) {
  DimPlot(
    seurat_qc_list[[nm]],
    reduction = "umap",
    group.by = "CAF_confident",
    cols = c("firebrick", "grey90")
    # cols = c("grey90", "firebrick")
  ) +
    ggtitle(nm)
})

names(plots) <- names(seurat_qc_list)

plots_max_cell_type <- lapply(names(seurat_qc_list), function(nm) {
  #if (names(nm) == "DF_CAF") next
  DimPlot(
    seurat_qc_list[[nm]],
    reduction = "umap",
    group.by = "max_celltype",
    cols = c("firebrick", "grey90", "green")
    # cols = c("grey90", "firebrick")
  ) +
    ggtitle(nm)
})

names(plots_max_cell_type) <- names(seurat_qc_list)

VlnPlot(
  seurat_qc_list$DF_CAF,
  features = c("COL1A1","KRT8","KRT18","KRT7"),
  group.by = "max_celltype"
)

VlnPlot(
  seurat_qc_list$DF_CAF,
  features = c(
    "COL1A1","COL1A2","DCN","LUM",
    "EPCAM","KRT18","KRT19",
    "PECAM1","VWF",
    "PTPRC","LYZ"
  ),
  group.by = "max_celltype"
)

# in plots for PC_CAF, there is a cluster of cells with mostly CAF_confident FALSE
# perform clustering for this sample - is there something about the markers we find in this cluster that can tell us what these cells are?
seurat_qc_PC <- FindNeighbors(seurat_qc_list$PC_CAF, dims = 1:14)
#DefaultAssay(seurat_qc_PC) <- "SCT"
seurat_qc_PC <- FindClusters(seurat_qc_PC, resolution = seq(from = 0.1, to = 2, by = 0.1))
p_PC <- clustree(seurat_qc_PC, prefix = "SCT_snn_res.", exprs = "scale.data")
p_PC
# try res = 0.3
Idents(seurat_qc_PC) <- seurat_qc_PC$SCT_snn_res.0.3
DimPlot(seurat_qc_PC, reduction = "umap", group.by = "SCT_snn_res.0.3")

# cluster 5 has our CAF_confident FALSE cells - this could be two subclusters, do initial DE using 0.3
markers_PC_0.3 <- FindAllMarkers(
  seurat_qc_PC,
  only.pos = TRUE
)
markers_PC_0.3_sig <-  markers_PC_0.3 %>%
  dplyr::filter(avg_log2FC > 0.25 & 
                  p_val_adj < 0.05 &
                pct.1 > 0.25 &
                pct.1 - pct.2 > 0.2
                )
markers_PC_0.3_sig

# try more fine-grained clustering resolution to get this as a single cluster
Idents(seurat_qc_PC) <- seurat_qc_PC$SCT_snn_res.0.5
DimPlot(seurat_qc_PC, reduction = "umap", group.by = "SCT_snn_res.0.5")

# now we have cluster 7 with our cells of interest, do DE analysis
markers_PC <- FindAllMarkers(
  seurat_qc_PC,
  only.pos = TRUE
)
  
markers_PC_sig <-  markers_PC %>%
  dplyr::filter(avg_log2FC > 0.25 & 
                  p_val_adj < 0.05 &
                  pct.1 > 0.25 &
                  pct.1 - pct.2 > 0.2
  )
markers_PC_sig
# no significant DE genes in cluster 7

#seurat_qc_4027 <- seurat_qc_list$`4027`
hist(seurat_qc_4027$Fibroblast, breaks = 100)

hist(seurat_qc_4027$CAF_score, breaks = 100)
DefaultAssay(seurat_qc_4027) <- "SCT"
seurat_qc_4027 <- NormalizeData(seurat_qc_4027)
seurat_qc_4027 <- FindVariableFeatures(seurat_qc_4027, 
                                      selection.method = "vst", nfeatures = 500)

seurat_qc_4027 <- ScaleData(seurat_qc_4027)
seurat_qc_4027 <- RunPCA(seurat_qc_4027, npcs = 20, 
                        features=VariableFeatures(seurat_qc_list_4027))
seurat_qc_4027 <- RunUMAP(seurat_qc_4027, reduction = "pca", 
                         dims = 1:20, seed.use=123)
FeaturePlot(seurat_qc_4027, reduction = "umap", features = names(gene_sets)) &
  theme(aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

VlnPlot(
  seurat_qc_4027,
  features = c(
    "Fibroblast",
    "Epithelial",
    "Immune",
    "Endothelial"
  ),
  pt.size = 0
)

seurat_qc_list_4027_scores <- data.frame(
  seurat_qc_4027$Fibroblast,
  seurat_qc_4027$Endothelial,
  seurat_qc_4027$Epithelial#,
  # seurat_qc_list_4027$Fibroblast_endothelial_down,
  # seurat_qc_list_4027$Epithelial_down
)

ElbowPlot(seurat_qc_4027)
seurat_qc_4027 <- FindNeighbors(seurat_qc_4027, dims = 1:12)
seurat_qc_4027 <- FindClusters(seurat_qc_4027, resolution = seq(from = 0.1, to = 2, by = 0.1))
p_4027 <- clustree(seurat_qc_4027, prefix = "SCT_snn_res.", exprs = "scale.data")
p_4027
Idents(seurat_qc_4027) <- seurat_qc_4027$SCT_snn_res.0.2
seurat_qc_4027 <- RunUMAP(seurat_qc_4027, reduction = "pca", 
                          dims = 1:12, seed.use=123)

VlnPlot(
  seurat_qc_4027,
  features = c(
    "Fibroblast",
    "Epithelial",
    "Endothelial"
  ),
  group.by = "SCT_snn_res.0.2"
)


#lung_fibroblast_clean <- unique(c("THY1","PDGFRB","PDGFRA","VIM","ACTA2","DES","CD36","ADGRE5","VCAN","DKK3","PDGFRA","DDR2","TCF21","TWIST2","P4HA1","S100A4","CCN2","PRRX1","SNAI1","VIM","KMT2B","OSR1","MIF","COL1A2","COL1A1","GAS6","GSN","CDKN1A","MYC","PLIN2","TOPORS","FBN1","CD34","VIT","FBLN2","PDGFRB","COL1A1","ACTA2","COL3A1"))
samplesheet <- read.csv("../01_metadata/samplesheet-kinnex-20260204.csv")
samplesheet$path <- paste("/home/rstudio/caf-scrnaseq", samplesheet$path, sep = "/")
names(seurat_qc_list) <- samplesheet$samples

# 4027
seurat_qc_4027 <- seurat_qc_list$`4027`
DefaultAssay(seurat_qc_4027) <- "SCT"
seurat_qc_4027 <- RunPCA(seurat_qc_4027)
ElbowPlot(seurat_qc_4027)
seurat_qc_4027 <- FindNeighbors(seurat_qc_4027, dims = 1:11)
seurat_qc_4027 <- FindClusters(seurat_qc_4027, resolution = seq(from = 0.1, to = 2, by = 0.1))
p_4027 <- clustree(seurat_qc_4027, prefix = "SCT_snn_res.", exprs = "scale.data")
p_4027
Idents(seurat_qc_4027) <- seurat_qc_4027$SCT_snn_res.0.1
seurat_qc_4027 <- RunUMAP(seurat_qc_4027, dims = 1:11)
DimPlot(seurat_qc_4027, reduction = "umap", group.by = "SCT_snn_res.0.1")
DimPlot(seurat_qc_4027, reduction = "pca", group.by = "SCT_snn_res.0.1")

Idents(seurat_qc_4027) <- seurat_qc_4027$SCT_snn_res.0.2
seurat_qc_4027 <- RunUMAP(seurat_qc_4027, dims = 1:11)
DimPlot(seurat_qc_4027, reduction = "umap", group.by = "SCT_snn_res.0.2")
DimPlot(seurat_qc_4027, reduction = "pca", group.by = "SCT_snn_res.0.2")
markers_4027 <- FindAllMarkers(
  seurat_qc_4027,
  only.pos = TRUE
)

markers_4027_sig <- markers_4027 %>%
  dplyr::filter(avg_log2FC > 0.25 & p_val_adj < 0.05)
markers_4027_sig

DotPlot(
  seurat_qc_4027,
  features = c(
    "COL1A1","COL1A2","DCN","LUM",
    "EPCAM","KRT18","KRT19",
    "PECAM1","VWF",
    "PTPRC","LYZ"
  )
)

# cluster zero has highest KRT18, still seem to be expressing COL1A1 -> is there a subcluster here that is expressing COL1A1 and not KRT18
# going from 0.2 to 0.3 splits cluster zero in two - see if they have different expression of these markers
Idents(seurat_qc_4027) <- seurat_qc_4027$SCT_snn_res.0.3
seurat_qc_4027 <- RunUMAP(seurat_qc_4027, dims = 1:11)
DimPlot(seurat_qc_4027, reduction = "umap", group.by = "SCT_snn_res.0.3")
DimPlot(seurat_qc_4027, reduction = "pca", group.by = "SCT_snn_res.0.3")
markers_4027_SCT_snn_res.0.3 <- FindAllMarkers(
  seurat_qc_4027,
  only.pos = TRUE
)
markers_4027_SCT_snn_res.0.3.sig <- markers_4027_SCT_snn_res.0.3 %>% dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05)
DotPlot(
  seurat_qc_4027,
  features = c(
    "COL1A1","COL1A2","DCN","LUM",
    "EPCAM","KRT18","KRT19",
    "PECAM1","VWF",
    "PTPRC","LYZ"
  )
)

FeaturePlot(
  seurat_qc_4027,
  features = c(
    "COL1A1","DCN","LUM",
    "EPCAM","KRT18","KRT19",
    "PECAM1","VWF",
    "PTPRC"
  )
)

VlnPlot(
  seurat_qc_4027,
  features = c("COL1A1","KRT8","KRT18","KRT7"),
  idents = "1"
)

FeatureScatter(seurat_qc_4027, "COL1A1", "KRT18")

msig_hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)
hallmark_list <- msig_hallmark %>%
  split(x = .$gene_symbol, f = .$gs_name)
markers_cluster1 <- markers_4027_SCT_snn_res.0.3 %>%
  filter(cluster == 1)
ranks <- markers_cluster1$avg_log2FC
names(ranks) <- markers_cluster1$gene

ranks <- sort(ranks, decreasing = TRUE)
fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats = ranks,
  minSize = 15,
  maxSize = 500
)
fgsea_res %>%
  arrange(padj) %>%
  select(pathway, NES, padj) %>%
  head(20)
library(ggplot2)

fgsea_res %>%
  arrange(padj) %>%
  head(15) %>%
  ggplot(aes(
    x = reorder(pathway, NES),
    y = NES
  )) +
  geom_col() +
  coord_flip()

stress_genes <- c(
  "FOS","JUN","JUNB",
  "HSPA1A","HSPB1","DNAJB1",
  "ATF3","DDIT3",
  "GADD45A"
)

epithelial_genes <- c(
  "EPCAM","CDH1","MUC1"
)

keratin_genes <- c(
  "KRT7","KRT8","KRT18","KRT19"
)

seurat_qc_4027 <- AddModuleScore(
  seurat_qc_4027,
  features = list(
    stress = stress_genes,
    epithelial = epithelial_genes,
    keratin = keratin_genes
  )
)
valid_stress <- stress_genes[stress_genes %in% rownames(seurat_qc_4027)]
valid_keratin <- keratin_genes[keratin_genes %in% rownames(seurat_qc_4027)]
seurat_qc_4027 <- AddModuleScore(
  seurat_qc_4027,
  features = list(
    stress = valid_stress
    #keratin = valid_keratin
  ),
  name = "Stress_module_score"
)

seurat_qc_4027 <- AddModuleScore(
  seurat_qc_4027,
  features = list(
    keratin = valid_keratin
  ),
  name = "Keratin_module_score"
)

seurat_qc_4027 <- AddModuleScore(
  seurat_qc_4027,
  features = list(
    keratin = epithelial_genes
  ),
  name = "Epithelial_module_score"
)
VlnPlot(seurat_qc_4027, features = c("Stress_module_score1","Keratin_module_score1", "Epithelial_module_score1"), group.by = "SCT_snn_res.0.3")
FeaturePlot(seurat_qc_4027, features = c("Stress_module_score1", "Keratin_module_score1", "COL1A1"))


rm(seurat_qc_4027)
gc()


seurat_qc_4299 <- seurat_qc_list$`4299`
seurat_qc_4299 <- RunPCA(seurat_qc_4299)
ElbowPlot(seurat_qc_4299)
seurat_qc_4299 <- FindNeighbors(seurat_qc_4299, dims = 1:15)
seurat_qc_4299 <- FindClusters(seurat_qc_4299, resolution = seq(from = 0.1, to = 2, by = 0.1))
p_4299 <- clustree(seurat_qc_4299, prefix = "SCT_snn_res.", exprs = "scale.data")
p_4299
Idents(seurat_qc_4299) <- seurat_qc_4299$SCT_snn_res.0.3

add_lineage_scores <- function(obj, gene_sets, assay = "SCT") {
  
  DefaultAssay(obj) <- assay
  
  # ensure genes exist
  gene_sets <- lapply(gene_sets, function(gs) {
    intersect(gs, rownames(obj))
  })
  print("gene_sets...")
  print(gene_sets)
  obj <- AddModuleScore(
    obj,
    features = gene_sets,
    name = names(gene_sets)
  )
  
  return(obj)
}
gene_sets <- list(
  Fibroblast = c("COL1A1","COL1A2","DCN","LUM","PDGFRA","PDGFRB","VIM"),
  Epithelial = c("EPCAM","KRT7","KRT8","KRT18","KRT19","CDH1"),
  Immune = c("PTPRC","LYZ","MS4A1","CD3D","CD3E"),
  Endothelial = c("PECAM1","VWF","KDR","RGS5"),
  Stress = c("FOS","JUN","HSP90AA1","HSPA1A","ATF3")
)

seurat_qc_list <- lapply(
  seurat_qc_list,
  function(x) {
    add_lineage_scores(
      obj = x,
      gene_sets = gene_sets
    )
  }
)
gc()

obj$fib_score <- obj$Fibroblast1
obj$epi_score <- obj$Epithelial1
obj$imm_score <- obj$Immune1
obj$endo_score <- obj$Endothelial1

xu_data <- readxl::read_xlsx(path = "caf-scrnaseq/00_raw_data/reference-data/mmc3-xu-et-al-canonical-markers-2024.xlsx")
xu_data_endothelial_up <- xu_data %>% dplyr::filter(`cell type` == "Endothelial cells" & `up/down` == "upregulated") %>% pull()
