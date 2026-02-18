options(repos = c(CRAN = "https://cloud.r-project.org"))

if(!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

options(repos=BiocManager::repositories())

install.packages(c("NMF", "dotenv", "rprojroot", "Seurat"), dependencies=TRUE)
BiocManager::install("rhdf5")
BiocManager::install("anndataR")
BiocManager::install("SingleCellExperiment")

devtools::install_github("immunogenomics/presto")
devtools::install_github("jokergoo/ComplexHeatmap", upgrade = "never")
devtools::install_github("jinworks/CellChat", upgrade = "never")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz")
devtools::install_github("cole-trapnell-lab/monocle3")
