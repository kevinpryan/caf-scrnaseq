docker run --rm -it \
    -p 8787:8787 \
    -e PASSWORD=beachballs \
    -v "$(pwd)":/home/rstudio/caf-scrnaseq \
    kevinr9525/rocker-rstudio-seurat-5:dev
