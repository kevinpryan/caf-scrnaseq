#!/usr/bin/env Rscript

library(SingleR)
library(Seurat)
library(SingleCellExperiment)
library(scran)

parser$add_argument('--input_rds', type='character', required=TRUE, help='Path for the input integrated RDS file')
parser$add_argument('--input_breast_rds', type='character', required=TRUE, help='Path for the breast RDS file')
parser$add_argument('--output_rds', type='character', required=TRUE, help='Path for the output annotated RDS file')
integrated_seurat_sce <- readRDS(parser$input_rds)


saveRDS(sce_object, file = "clustered_integrated_sce_dietseurat.rds")
