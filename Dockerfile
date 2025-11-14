FROM rocker/rstudio:4.3.2
LABEL description="Docker image for Seurat5"

RUN apt-get update && apt-get install -y \
    libhdf5-dev build-essential libxml2-dev \
    libssl-dev libv8-dev libsodium-dev libglpk40 \
    libgdal-dev libboost-dev libomp-dev \
    libbamtools-dev libboost-iostreams-dev \
    libboost-log-dev libboost-system-dev \
    libboost-test-dev libcurl4-openssl-dev libz-dev \
    libarmadillo-dev libhdf5-cpp-103

RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/')"

RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5', quiet = TRUE)"

RUN R -e "install.packages(c('hdf5r', 'dplyr', 'tidyverse', 'cowplot', 'knitr', 'slingshot', 'msigdbr', 'remotes', 'metap', 'devtools', 'R.utils', 'ggalt', 'ggpubr', 'BiocManager'), repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install(c('SingleR', 'slingshot', 'scRNAseq', 'celldex', 'fgsea', 'multtest', 'scuttle', 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'org.Mm.eg.db', 'AnnotationHub', 'scater', 'edgeR', 'apeglm', 'DESeq2', 'pcaMethods', 'clusterProfiler'))"

RUN R -e "remotes::install_github(c('satijalab/seurat-wrappers', 'kevinblighe/PCAtools', 'chris-mcginnis-ucsf/DoubletFinder', 'velocyto-team/velocyto.R'))"

RUN R -e "install.packages('argparse')"
