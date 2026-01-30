options(repos = c(CRAN = "https://cloud.r-project.org"))

if(!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

options(repos=BiocManager::repositories())

install.packages(c("NMF", "dotenv", "rprojroot"), dependencies=TRUE)
BiocManager::install("rhdf5")
BiocManager::install("anndataR")

devtools::install_github("immunogenomics/presto")
devtools::install_github("jokergoo/ComplexHeatmap", upgrade = "never")
devtools::install_github("jinworks/CellChat", upgrade = "never")
