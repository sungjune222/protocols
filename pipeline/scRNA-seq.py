import gc
import os
from pipeline.config.constants import SINGLE_CELL_VAE_BATCH_SIZE, CPU_CORE_COUNT

os.environ["NUMBA_NUM_THREADS"] = str(CPU_CORE_COUNT)
os.environ["OMP_NUM_THREADS"] = str(CPU_CORE_COUNT)
os.environ["OPENBLAS_NUM_THREADS"] = str(CPU_CORE_COUNT)
os.environ["MKL_NUM_THREADS"] = str(CPU_CORE_COUNT)

import cupy as cp
import numpy as np
import pandas as pd
import scvi
import scanpy as sc
import rapids_singlecell as rsc
from anndata import AnnData
from scipy.sparse import csr_matrix
from pipeline.utils import plot
from pipeline.utils.env import find_env_dir
from pipeline.utils.extract_mt_gene import extract_mt_genes
from pipeline.utils.fix_nullable_strings import fix_nullable_strings
from pipeline.config.machine_learning import DataLoader

sc.settings.n_jobs = CPU_CORE_COUNT

if __name__ == "__main__":
    pre_h5ad_dir = find_env_dir("PRE_H5AD")
    MAIN_SERIES = "pbmc_human"
    SUB_SERIES = ""
    SPECIES = "human"
    leiden_resolution = 3.5
    only_cpu = False

    if SUB_SERIES:
        SERIES = MAIN_SERIES + "_" + SUB_SERIES
    else:
        SERIES = MAIN_SERIES
    file = os.path.join(pre_h5ad_dir, SERIES + "_raw.h5ad")

    # Preprocessing each sample
    print("Loading data: " + SERIES + "...")
    loaded_adata = sc.read_h5ad(file)

    if not only_cpu:
        assert isinstance(loaded_adata.X, csr_matrix)
        loaded_adata.X = loaded_adata.X.astype(np.float32)
        rsc.get.anndata_to_GPU(loaded_adata)

    # Cell and gene filtering
    MIN_GENES = 200
    MIN_CELLS = 20
    MIN_COUNTS = 1000
    if SUB_SERIES == "":
        # Quality assessment by calculating QC metrics
        def quality_assess(adata: AnnData) -> AnnData:
            # Marking mitochondrial and ribosomal genes
            mt_genes = extract_mt_genes(species=SPECIES)
            mt_mask = adata.var.index.isin(mt_genes)
            ribo_mask = adata.var.index.str.upper().str.startswith(("RPS", "RPL"))

            adata.var["mt"] = np.asarray(mt_mask, dtype=bool)
            adata.var["ribo"] = np.asarray(ribo_mask, dtype=bool)

            # qc_vars: List of categories that you want to make as a QC metrics (It must be set as a boolean list in AnnData.obs)
            if not only_cpu:
                rsc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["mt", "ribo"],
                    log1p=True,
                )
            else:
                sc.pp.calculate_qc_metrics(
                    adata, qc_vars=["mt", "ribo"], log1p=True, inplace=True
                )
            return adata

        print("Assessing quality...")
        quality_assessed_adata = quality_assess(loaded_adata)

        if not only_cpu:
            rsc.pp.filter_cells(quality_assessed_adata, min_genes=MIN_GENES)
            rsc.pp.filter_genes(quality_assessed_adata, min_cells=MIN_CELLS)
            rsc.pp.filter_cells(quality_assessed_adata, min_counts=MIN_COUNTS)
        else:
            sc.pp.filter_cells(quality_assessed_adata, min_genes=MIN_GENES)
            sc.pp.filter_genes(quality_assessed_adata, min_cells=MIN_CELLS)
            sc.pp.filter_cells(quality_assessed_adata, min_counts=MIN_COUNTS)

        # This plot will take a lot of time
        # plot.plot_qc(quality_assessed_adata, SERIES)

        # Cytoplasmic RNA in dead cells leaks out, resulting in a higher proportion of remaining mitochondrial RNA
        filtered_adata = quality_assessed_adata[
            (quality_assessed_adata.obs["pct_counts_mt"] < 15)
        ].copy()
    else:
        quality_assessed_adata = loaded_adata
        filtered_adata = quality_assessed_adata

    del loaded_adata
    del quality_assessed_adata
    gc.collect()

    if only_cpu:
        filtered_adata = fix_nullable_strings(filtered_adata)
        gc.collect()

    filtered_h5ad_dir = find_env_dir("FILTERED_H5AD")
    filtered_adata.write_h5ad(
        os.path.join(filtered_h5ad_dir, SERIES + "_filtered.h5ad"), compression="gzip"
    )

    # %% Producing latent representation of cells using scVI
    print("Producing latent representation of cells using scVI...")
    N_TOP_GENES = 4000
    if not only_cpu:
        rsc.pp.highly_variable_genes(
            filtered_adata,
            n_top_genes=N_TOP_GENES,
            flavor="seurat_v3",
        )
    else:
        sc.pp.highly_variable_genes(
            filtered_adata, n_top_genes=N_TOP_GENES, flavor="seurat_v3"
        )
        filtered_adata = fix_nullable_strings(filtered_adata)

    if not only_cpu:
        rsc.get.anndata_to_CPU(filtered_adata, convert_all=True)
        gc.collect()
        cp.get_default_memory_pool().free_all_blocks()
    scvi_adata = filtered_adata[:, filtered_adata.var["highly_variable"]].copy()

    # Setting up AnnData for scVI model with batch effects
    scvi.model.SCVI.setup_anndata(
        scvi_adata,
        # The primary batch info
        batch_key="sample",
    )
    model = scvi.model.SCVI(scvi_adata, n_latent=50)
    model.train(
        accelerator="gpu",
        precision="bf16-mixed",
        batch_size=SINGLE_CELL_VAE_BATCH_SIZE,
        datasplitter_kwargs=DataLoader,
        train_size=0.9,
        check_val_every_n_epoch=5,
        max_epochs=400,
        early_stopping=True,
    )
    plot.plot_validation_loss(model, SERIES, file_info="model_validation_loss")

    # latent_representation: (cell, latent_space_dimension)
    if filtered_adata.obs_names.equals(scvi_adata.obs_names):
        filtered_adata.obsm["X_scvi"] = model.get_latent_representation()  # type: ignore
    else:
        raise ValueError(
            "Cell names do not match between filtered_adata and scvi_adata"
        )

    scvi_model_dir = find_env_dir("SCVI_MODEL")
    model.save(
        os.path.join(scvi_model_dir, SERIES),
        overwrite=True,
        save_anndata=True,
    )
    del scvi_adata
    gc.collect()

    dimension_reduced_h5ad_dir = find_env_dir("DIMENSION_REDUCED_H5AD")
    filtered_adata.write_h5ad(
        os.path.join(
            dimension_reduced_h5ad_dir,
            SERIES + "_dimension_reduced.h5ad",
        )
    )

    # %% Clustering
    if not only_cpu:
        rsc.get.anndata_to_GPU(filtered_adata)

    # Uses pre-calculated scVI latent representation for calculating similarity score and constructing neighborhood graph
    # Store settings in .uns["neighbors"] and connectivity matrices in .obsp
    print("Constructing neighborhood graph...")
    if not only_cpu:
        rsc.pp.neighbors(
            filtered_adata, n_neighbors=15, use_rep="X_scvi", metric="cosine"
        )
    else:
        sc.pp.neighbors(
            filtered_adata, n_neighbors=15, use_rep="X_scvi", metric="cosine"
        )

    # Embeds the neighborhood graph into 2D space using UMAP algorithm (optimized via SGD, Stochastic Gradient Descent)
    # You can change n_components to 3 for 3D UMAP
    print("Calculating UMAP...")
    if not only_cpu:
        rsc.tl.umap(filtered_adata, n_components=2, min_dist=0.3)
    else:
        sc.tl.umap(
            filtered_adata,
            n_components=2,
            min_dist=0.3,
            init_pos="random",
            random_state=0,
        )
    # If the clusters appear too clumped or merged, try decreasing n_neighbors and min_dist

    # Clustering cells using leiden algorithm, maximizes modularity which is defined based on its intergroup connectivity and expected (random) connectivity
    # High resolution value results in more clusters
    print("Clustering with Leiden algorithm...")
    if not only_cpu:
        rsc.tl.leiden(filtered_adata, resolution=leiden_resolution)
    else:
        sc.tl.leiden(
            filtered_adata,
            resolution=leiden_resolution,
            flavor="igraph",
            n_iterations=-1,
            directed=False,
        )

    plot.plot_umap(filtered_adata, SERIES)

    clustered_h5ad_dir = find_env_dir("CLUSTERED_H5AD")
    filtered_adata.write_h5ad(
        os.path.join(
            clustered_h5ad_dir,
            SERIES + ".h5ad",
        ),
        compression="gzip",
    )

    print("Clustering completed")
