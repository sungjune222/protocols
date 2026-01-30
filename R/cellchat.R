invisible(
  suppressPackageStartupMessages({
    library(anndataR)
    library(CellChat)
    library(Matrix)
  })
)

series_name <- "SCP1038"
clustered_data_dir <- get_env_dir("CLUSTERED_DATA_LOCATION")
cellchat_dir <- get_env_dir("CELLCHAT_DATA_LOCATION")
cellchat_dir <- file.path(cellchat_dir, series_name)
dir.create(cellchat_dir, recursive = TRUE, showWarnings = FALSE)

adata <- anndataR::read_h5ad(
  file.path(
    clustered_data_dir,
    paste0(series_name, "_clustered.h5ad")
  ),
  as = "HDF5AnnData"
)

group_key <- "leiden"
stopifnot(group_key %in% colnames(adata$obs))

provided_genes <- CellChat::extractGene(CellChatDB.mouse)
genes_use <- intersect(provided_genes, adata$var_names)
compact_adata <- adata[, genes_use]

raw_counts <- compact_adata$X
stopifnot(inherits(raw_counts, "dgRMatrix"))

cells <- compact_adata$obs_names
genes <- compact_adata$var_names

data_raw <- Matrix::t(raw_counts)
rownames(data_raw) <- genes
colnames(data_raw) <- cells

meta_data <- data.frame(
  group = paste0("cluster_", compact_adata$obs[[group_key]]),
  row.names = cells,
  samples = compact_adata$obs[["sample"]]
)

data_input <- CellChat::normalizeData(
  data_raw,
  scale.factor = 10000, do.log = TRUE
)

cellchat <- CellChat::createCellChat(
  object = data_input, meta = meta_data,
  group.by = "group"
)

# Set the CellChatDB
cellchat@DB <- CellChatDB.mouse
# Keep only signaling genes present in the CellChat database
cellchat <- subsetData(cellchat)

# For parallel processing
# This option must be inactivated if you have insufficient memory
# future::plan("multisession", workers = 4)
# options(future.globals.maxSize = Inf)
# Identify signaling genes that are overexpressed in each cell group
cellchat <- identifyOverExpressedGenes(cellchat)
# Select overexpressed ligand–receptor interactions
# based on the overexpressed genes and the CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute ligand–receptor communication probabilities between all interactions
cellchat <- computeCommunProb(cellchat)
# Remove communications involving groups with too few cells
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Aggregate L-R probabilities into pathway-level communication probabilities
cellchat <- computeCommunProbPathway(cellchat)
# Summarize results into 2D group-by-group networks
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = file.path(cellchat_dir, "cellchat.rds"))
