options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(
  c(
    "BiocManager",
    "ashr", "NMF", "dotenv", "rprojroot", "Seurat"
  ),
  dependencies = TRUE
)

options(repos = BiocManager::repositories())
BiocManager::install(
  c(
    "rhdf5", "anndataR", "SingleCellExperiment",
    "scDblFinder", "DropletUtils"
  ),
  dependencies = TRUE,
  update = FALSE,
  ask = FALSE
)

devtools::install_github("immunogenomics/presto")
devtools::install_github("jokergoo/ComplexHeatmap", upgrade = "never")
devtools::install_github("jinworks/CellChat", upgrade = "never")
remotes::install_url(
  "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz"
)
devtools::install_github("cole-trapnell-lab/monocle3")
devtools::install_github("myles-lewis/glmmSeq")
