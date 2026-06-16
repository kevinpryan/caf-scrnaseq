library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
BiocManager::install('glmGamPoi')
setwd("~/caf-scrnaseq/02a_scripts_kinnex")
samplesheet <- read.csv("../01_metadata/samplesheet-kinnex-20260204.csv")
samplesheet$path <- paste("/home/rstudio/caf-scrnaseq", samplesheet$path, sep = "/")
process_sample <- function(path, sample_name) {
  
  # ---------------------------
  # 1. Read + create Seurat
  # ---------------------------
  counts <- Read10X(data.dir = path)
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name
  )
  
  # ---------------------------
  # 2. Basic QC metrics (NO filtering)
  # ---------------------------
  obj$barcode_raw <- Cells(obj)
  obj$sample <- sample_name
  
  obj$mt_percent <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$nCount_RNA <- obj$nCount_RNA
  obj$nFeature_RNA <- obj$nFeature_RNA
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA + 1) / log10(obj$nCount_RNA + 1)
  
  # ---------------------------
  # 3. QC flags 
  # ---------------------------
  obj$qc_pass <- (
    obj$mt_percent < 5 &
      obj$nFeature_RNA > 250 &
      obj$nCount_RNA > 500 &
      obj$log10GenesPerUMI > 0.80
  )
  
  # ---------------------------
  # 4. Doublet detection 
  # ---------------------------
  sce <- SingleCellExperiment(
    assays = list(counts = GetAssayData(obj, layer = "counts"))
  )
  
  sce <- scDblFinder(sce)
  
  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score
  
  obj$qc_pass <- obj$qc_pass & (obj$doublet_class != "doublet")
  
  # ---------------------------
  # 5. Return full annotated object
  # ---------------------------
  return(obj)
}

generate_qc_table <- function(obj){
  cat("Generating QC metrics report...\n")
  cells_fail_nFeature <- sum(obj@meta.data$nFeature_RNA < 250)
  cells_fail_nCount <- sum(obj@meta.data$nCount_RNA < 500)
  cells_fail_log10Genes <- sum(obj@meta.data$log10GenesPerUMI <= 0.80)
  cells_fail_mito <- sum(obj@meta.data$mt_percent >= 5)
  cells_fail_doublet <- sum(obj@meta.data$doublet_class == "doublet")
  meta <- obj@meta.data
  final_cell_filter <- meta %>% dplyr::filter(
      nFeature_RNA > 250 &
      nCount_RNA > 500 &
      log10GenesPerUMI > 0.80 &
        mt_percent < 5 &
      doublet_class != "doublet"
  ) 
  final_cell_count <- nrow(final_cell_filter)
  print(final_cell_count)
  initial_cell_count <- nrow(obj@meta.data)
  total_cells_removed <- initial_cell_count - final_cell_count
  
  # Create a data frame with the QC summary
  qc_summary <- data.frame(
    Metric = c(
      "Initial number of cells",
      "Cells failing nFeature_RNA (< 250)",
      "Cells failing nCount_RNA (< 500)",
      "Cells failing log10GenesPerUMI (<= 0.80)",
      "Cells failing percent.mt (>= 5%)",
      "Cells identified as doublets",
      "Total cells removed (combination of all filters)",
      "Final number of cells remaining"
    ),
    Count = c(
      initial_cell_count,
      cells_fail_nFeature,
      cells_fail_nCount,
      cells_fail_log10Genes,
      cells_fail_mito,
      cells_fail_doublet,
      total_cells_removed,
      final_cell_count
    )
  )
}
seurat_list <- mapply(
  process_sample,
  path = samplesheet$path,
  sample_name = samplesheet$samples,
  SIMPLIFY = FALSE
)
names(seurat_list) <- samplesheet$samples

qc_list <- mapply(
  generate_qc_table,
  obj = seurat_list,
  SIMPLIFY = FALSE
)

seurat_qc_list <- lapply(
  seurat_list,
  function(x) subset(x, subset = qc_pass)
)

gc()

library(future)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)

seurat_qc_list <- lapply(
  seurat_qc_list,
  function(x) {
    SCTransform(
      x,
      vars.to.regress = "mt_percent",
      method = "glmGamPoi",
      vst.flavor = "v2",
      verbose = FALSE
    )
  }
)
gc()

seurat_qc_list
saveRDS(seurat_qc_list, file = "../03_processed_data/kinnex/01_qc_and_filtering/caf_kinnex_seurat_qc_list_sctransform.Rds")
# seurat_merged <- merge(
#   x = seurat_list[[1]],
#   y = seurat_list[2:5]
# )
# 
# seurat_merged$cell_id <- paste(
#   seurat_merged$sample,
#   seurat_merged$barcode_raw,
#   sep = ":"
# )
# 
# seurat_qc <- subset(seurat_merged, subset = qc_pass)

# celltypedb
# breast -> fibroblast
# POSTN,PDGFRB,PDGFRA,CK18,CK19,EPCAM,CK8
# breast -> CAF
# PDGFRA, PDGFRB, LGALS1, CAV1, ACTA2
# lung - tissue lung - fibroblast 
lung_fibroblast_original <- c("CD90","PDGFRB","PDGFR-alpha","Vimentin","alpha-SMA","Desmin","CD36","CD97","VCAN","DKK3","PDGFRA","DDR2","TCF21","TWIST2","P4HA1","S100A4","CTGF","PRRX1","SNAI1","VIM","KMT2B","OSR1","MIF","COL1A2","COL1A1","GAS6","GSN","CDKN1A","MYC","PLIN2","LUN","FBN1","CD34","VIT","FBLN2","PDGFRB","COL1A1","ACTA2","COL3A1")
lung_fibroblast_clean <- unique(c("THY1","PDGFRB","PDGFRA","VIM","ACTA2","DES","CD36","ADGRE5","VCAN","DKK3","PDGFRA","DDR2","TCF21","TWIST2","P4HA1","S100A4","CCN2","PRRX1","SNAI1","VIM","KMT2B","OSR1","MIF","COL1A2","COL1A1","GAS6","GSN","CDKN1A","MYC","PLIN2","TOPORS","FBN1","CD34","VIT","FBLN2","PDGFRB","COL1A1","ACTA2","COL3A1"))
length(lung_fibroblast_clean)
# ? FAP-1

# normalise 
# seurat_qc <- NormalizeData(seurat_qc, normalization.method = "LogNormalize", scale.factor = 10000)
# seurat_qc <- FindVariableFeatures(seurat_qc, selection.method = "vst", nfeatures = 2000)
# seurat_qc <- ScaleData(seurat_qc)#, features = all.genes)
gc()
seurat_qc <- SCTransform(seurat_qc, vars.to.regress = "mt_percent")
top10 <- head(VariableFeatures(seurat_qc), 10)



seurat_list_merged <- merge(seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]], seurat_list[[4]], seurat_list[[5]]))
seurat_list_merged$log10GenesPerUMI <- log10(seurat_list_merged$nFeature_RNA) / log10(seurat_list_merged$nCount_RNA)
seurat_list_merged$mitoRatio <- PercentageFeatureSet(object = seurat_list_merged, pattern = "^MT-")
seurat_list_merged$mitoRatio <- seurat_list_merged@meta.data$mitoRatio / 100
