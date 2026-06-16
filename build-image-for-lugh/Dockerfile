FROM rocker/rstudio:4.4.3

LABEL description="R 4.4.3 + Seurat v5 single-cell analysis environment"

ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

# System dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libhdf5-dev \
    libv8-dev \
    libsodium-dev \
    libglpk40 \
    libgdal-dev \
    libboost-dev \
    libomp-dev \
    libbamtools-dev \
    libboost-iostreams-dev \
    libboost-log-dev \
    libboost-system-dev \
    libboost-test-dev \
    libarmadillo-dev \
    libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('remotes','devtools','BiocManager'), repos='https://cloud.r-project.org')"

RUN R -e "install.packages(c( \
    'hdf5r', \
    'dplyr', \
    'tidyverse', \
    'cowplot', \
    'knitr', \
    'msigdbr', \
    'metap', \
    'R.utils', \
    'ggalt', \
    'ggpubr', \
    'argparse', \
    'clustree', \
    'devtools', \
    'patchwork' \
), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
    'SingleR', \
    'scRNAseq', \
    'celldex', \
    'fgsea', \
    'multtest', \
    'scuttle', \
    'BiocGenerics', \
    'DelayedArray', \
    'DelayedMatrixStats', \
    'limma', \
    'S4Vectors', \
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'batchelor', \
    'org.Mm.eg.db', \
    'AnnotationHub', \
    'scater', \
    'edgeR', \
    'apeglm', \
    'DESeq2', \
    'clusterProfiler', \
    'PCAtools', \
    'slingshot', \
    'scDblFinder', \
    'SpatialExperiment', \
    'GSEABase', \
    'GSVA', \
    'singscore', \
    'UCell', \
    'glmGamPoi' \
), ask=FALSE)"

# Seurat v5
#RUN R -e "remotes::install_github('satijalab/seurat', ref='seurat5')"

RUN R -e "install.packages('Seurat', version = '5.5.0', repos='https://cloud.r-project.org')"

#RUN R -e "update.packages(ask=FALSE, checkBuilt=TRUE)"

# GitHub-only packages
RUN R -e "remotes::install_github(c( \
    'satijalab/seurat-wrappers', \
    'velocyto-team/velocyto.R', \
    'immunogenomics/presto', \
    'immunogenomics/harmony' \
), upgrade='never')"

#RUN R -e "devtools::install_github('immunogenomics/presto')"
