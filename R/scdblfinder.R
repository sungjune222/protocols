suppressPackageStartupMessages({
  library(scDblFinder)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    return(default)
  }
  args[i + 1]
}

manifest <- get_arg("--manifest")
outdir <- get_arg("--outdir")
threads <- as.integer(get_arg("--threads", "4"))
dbr_per1k <- as.numeric(get_arg("--dbr_per1k", "0.008"))
seed <- as.integer(get_arg("--seed", "1"))

if (is.null(manifest) || is.null(outdir)) {
  stop(
    "Usage: --manifest <success_libraries.csv> --outdir <dir>"
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

tab <- read.csv(manifest)
if (!all(c("LibraryID", "FilteredH5") %in% colnames(tab))) {
  stop("Manifest must contain LibraryID and FilteredH5 columns")
}

sce_list <- vector("list", nrow(tab))
removed_zero <- setNames(integer(nrow(tab)), tab$LibraryID)

for (i in seq_len(nrow(tab))) {
  library_id <- tab$LibraryID[i]
  message(sprintf("[%d/%d] Loading library: %s ...", i, nrow(tab), library_id))
  sce <- read10xCounts(tab$FilteredH5[i], col.names = TRUE, type = "HDF5")

  rownames(sce) <- make.unique(rownames(sce))
  orig_bc <- make.unique(colnames(sce))
  colnames(sce) <- paste0(orig_bc, "-", library_id)
  sce$orig_barcode <- orig_bc

  keep <- Matrix::colSums(counts(sce)) > 0
  removed_zero[library_id] <- sum(!keep)
  sce <- sce[, keep]

  if (ncol(sce) == 0) next
  sce$library_id <- library_id
  sce_list[[i]] <- sce
}

sce_list <- Filter(Negate(is.null), sce_list)
if (length(sce_list) == 0) stop("No cells remain after zero-count filtering.")
all_genes <- Reduce(union, lapply(sce_list, rownames))

sce_list <- lapply(sce_list, function(x) {
  missing_genes <- setdiff(all_genes, rownames(x))

  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(
      0,
      nrow = length(missing_genes), ncol = ncol(x), sparse = TRUE
    )
    rownames(zero_mat) <- missing_genes
    new_counts <- rbind(counts(x), zero_mat)
  } else {
    new_counts <- counts(x)
  }

  new_counts <- new_counts[all_genes, , drop = FALSE]
  SingleCellExperiment(list(counts = new_counts), colData = colData(x))
})
sce <- do.call(cbind, sce_list)

bp <- if (threads > 1) {
  MulticoreParam(workers = threads, progressbar = TRUE, RNGseed = seed)
} else {
  SerialParam(progressbar = TRUE, RNGseed = seed)
}

sce <- scDblFinder(
  sce,
  samples = sce$library_id,
  clusters = TRUE,
  dbr = NULL,
  dbr.per1k = dbr_per1k,
  threshold = TRUE,
  verbose = TRUE,
  BPPARAM = bp
)

calls <- data.frame(
  cell = colnames(sce),
  orig_barcode = sce$orig_barcode,
  library_id = sce$library_id,
  score = sce$scDblFinder.score,
  class = sce$scDblFinder.class
)

for (library_id in unique(calls$library_id)) {
  sdir <- file.path(outdir, library_id)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)

  sub <- calls[calls$library_id == library_id, , drop = FALSE]
  singlet <- sub$class == "singlet"
  doublet <- sub$class == "doublet"

  write.csv(
    data.frame(barcode = sub$orig_barcode[singlet]),
    file.path(sdir, sprintf("%s_singlet_barcodes.csv", library_id)),
    row.names = FALSE
  )

  expected_doublet_fraction <- min(1, (nrow(sub) / 1000) * dbr_per1k)
  write.csv(
    data.frame(
      LibraryID = library_id,
      TotalCells = nrow(sub),
      RemovedZeroCountCells = unname(removed_zero[library_id]),
      PredictedDoublets = sum(doublet, na.rm = TRUE),
      PredictedSinglets = sum(singlet, na.rm = TRUE),
      ObservedDoubletFraction = mean(doublet, na.rm = TRUE),
      ExpectedDoubletFraction = expected_doublet_fraction,
      dbr_per1k = dbr_per1k,
      scDblFinderVersion = as.character(packageVersion("scDblFinder")),
      Status = "SUCCESS",
      Note = "OK"
    ),
    file.path(sdir, sprintf("%s_summary.csv", library_id)),
    row.names = FALSE
  )
}
