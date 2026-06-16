docker run --rm -it \
    -p 8787:8787 \
    -e PASSWORD=beachballs \
    -v "$(pwd)/../":/home/rstudio/caf-scrnaseq \
    -v /home/kevin/Documents/PhD/Projects/caf-bc/data/scRNA-seq-cords/:/data/scRNA-seq-cords/ \
    kevinr9525/rocker-rstudio-4.4.3-seurat-v5:dev
    #kevinr9525/rocker-rstudio-seurat-5:dev

#docker run --rm -it \
#    -p 8787:8787 \
#    -e PASSWORD=beachballs \
#    -v "$(pwd)":/home/rstudio/caf-scrnaseq \
#    -v /home/kevin/Documents/PhD/Projects/caf-bc/data/scRNA-seq-cords/:/data/scRNA-seq-cords/ \
#    kevinr9525/rocker-rstudio-4.3.2-seurat-v5:dev
    #kevinr9525/rocker-tidyverse-4.5-bioconductor-seurat:dev
