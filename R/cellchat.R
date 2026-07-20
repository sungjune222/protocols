run_cellchat <- function(
  h5ad_file,
  species,
  data_dir = "CLUSTERED_DATA",
  group_key = "leiden",
  sample_key = "sample",
  min_cells = 100
) {
  if (missing(species)) {
    stop("Error: 'species' argument is required. Please specify 'mouse' or 'human'")
  }

  species <- match.arg(species, choices = c("mouse", "human"))

  invisible(
    suppressPackageStartupMessages({
      library(anndataR)
      library(CellChat)
      library(Matrix)
    })
  )

  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else {
    CellChatDB <- CellChatDB.human
  }

  base_name <- tools::file_path_sans_ext(h5ad_file)

  data_dir <- get_env_dir(data_dir)
  cellchat_dir <- get_env_dir("CELLCHAT")
  cellchat_dir <- file.path(cellchat_dir, base_name)
  dir.create(cellchat_dir, recursive = TRUE, showWarnings = FALSE)

  adata <- anndataR::read_h5ad(
    file.path(data_dir, h5ad_file),
    as = "HDF5AnnData"
  )

  stopifnot(group_key %in% colnames(adata$obs))
  stopifnot(sample_key %in% colnames(adata$obs))

  raw_counts <- adata$X
  stopifnot(inherits(raw_counts, "dgRMatrix"))

  cells <- as.character(adata$obs_names)
  genes <- as.character(adata$var_names)

  if (anyDuplicated(cells)) {
    stop("Duplicated cell names detected.")
  }

  if (anyDuplicated(genes)) {
    genes <- make.unique(genes)
    warning("Duplicated gene names detected; make.unique() was applied.")
  }

  data_raw <- Matrix::t(raw_counts)
  rownames(data_raw) <- genes
  colnames(data_raw) <- cells

  group <- paste0("cluster_", as.character(adata$obs[[group_key]]))
  samples <- as.character(adata$obs[[sample_key]])

  meta_data <- data.frame(
    group = group,
    samples = samples,
    row.names = cells,
    stringsAsFactors = FALSE
  )

  meta_data$group <- droplevels(factor(meta_data$group))

  data_input <- CellChat::normalizeData(
    data_raw,
    scale.factor = 10000,
    do.log = TRUE
  )

  cellchat <- CellChat::createCellChat(
    object = data_input,
    meta = meta_data,
    group.by = "group"
  )

  cellchat@DB <- CellChatDB

  cellchat <- subsetData(cellchat)

  message("Genes in CellChat input: ", nrow(cellchat@data))
  message("Genes in data.signaling: ", nrow(cellchat@data.signaling))

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  n_lrsig <- ifelse(is.null(cellchat@LR$LRsig), 0, nrow(cellchat@LR$LRsig))
  message("LRsig after identifyOverExpressedInteractions: ", n_lrsig)

  if (n_lrsig == 0) {
    stop(
      "No ligand-receptor pairs were selected by identifyOverExpressedInteractions(). ",
      "This is the real failure point. Check gene symbols, species, group granularity, ",
      "and expression sparsity. Do not run computeCommunProb() with LRsig = 0."
    )
  }

  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  saveRDS(
    cellchat,
    file = file.path(cellchat_dir, paste0(base_name, ".rds"))
  )

  cellchat
}
