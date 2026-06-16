docker build --build-arg GITHUB_PAT=$(cat ~/.config/github/pat) \
             -t kevinr9525/rocker-rstudio-4.4.3-seurat-v5:dev .
