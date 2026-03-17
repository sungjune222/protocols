# %% Environment Setup
import gc
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import rapids_singlecell as rsc
from anndata import AnnData
from pipeline.utils import plot
from pipeline.utils.env import find_env_dir
from pipeline.config.constants import (
    SINGLE_CELL_VAE_BATCH_SIZE,
)
from pipeline.config.machine_learning import DataLoader

if __name__ == "__main__":
    pre_h5ad_dir = find_env_dir("PRE_H5AD")
    series_name = "macnair"
    leiden_resolution = 3.5
    file = os.path.join(pre_h5ad_dir, series_name + "_raw.h5ad")

    # Preprocessing each sample
    print("Loading data...")
    loaded_adata = sc.read_h5ad(file)
    rsc.get.anndata_to_GPU(loaded_adata)

    # Quality assessment by calculating QC metrics
    def quality_assess(adata: AnnData, ribo_genes: pd.DataFrame) -> AnnData:
        # Marking mitochondrial and ribosomal genes
        mt_mask = adata.var.index.str.lower().str.startswith("mt-")
        ribo_mask = adata.var_names.isin(ribo_genes[0].values)
        adata.var["mt"] = np.asarray(mt_mask, dtype=bool)
        adata.var["ribo"] = np.asarray(ribo_mask, dtype=bool)

        # Generates QC metrics
        # qc_vars: List of categories that you want to make as a QC metrics (It must be set as a boolean list in AnnData.obs)
        rsc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo"], log1p=True,
        )
        return adata

    ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)

    print("Assessing quality...")
    quality_assessed_adata = quality_assess(loaded_adata, ribo_genes)

    # Cell and gene filtering
    rsc.pp.filter_cells(quality_assessed_adata, min_genes=200)
    rsc.pp.filter_genes(quality_assessed_adata, min_cells=10)
    # Filtered cells and genes necessitate re-calculation of QC metrics
    quality_assessed_adata = quality_assess(quality_assessed_adata, ribo_genes)

    plot.plot_qc(quality_assessed_adata, series_name)

    # Cytoplasmic RNA in dead cells leaks out, resulting in a higher proportion of remaining mitochondrial RNA
    filtered_adata = quality_assessed_adata[
        (quality_assessed_adata.obs["pct_counts_mt"] < 15)
    ]
    gc.collect()

    filtered_h5ad_dir = find_env_dir("FILTERED_H5AD")
    filtered_adata.write_h5ad(
        os.path.join(filtered_h5ad_dir, series_name + "_filtered.h5ad")
    )

    # %% Producing latent representation of cells using scVI
    print("Producing latent representation of cells using scVI...")
    
    rsc.pp.highly_variable_genes(
        filtered_adata,
        n_top_genes=4000,
        flavor="seurat_v3",
    )
    scvi_adata = filtered_adata[:, filtered_adata.var["highly_variable"]].copy()
    rsc.get.anndata_to_CPU(scvi_adata)

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
    plot.plot_validation_loss(model, series_name, file_info="model_validation_loss")

    # latent_representation: (cell, latent_space_dimension)
    if filtered_adata.obs_names.equals(scvi_adata.obs_names):
        filtered_adata.obsm["X_scvi"] = model.get_latent_representation()
    else:
        raise ValueError(
            "Cell names do not match between filtered_adata and scvi_adata"
        )

    scvi_model_dir = find_env_dir("SCVI_MODEL")
    model.save(
        os.path.join(scvi_model_dir, series_name),
        overwrite=True,
        save_anndata=True,
    )
    del scvi_adata
    gc.collect()

    dimension_reduced_h5ad_dir = find_env_dir("DIMENSION_REDUCED_H5AD")
    filtered_adata.write_h5ad(
        os.path.join(
            dimension_reduced_h5ad_dir,
            series_name + "_dimension_reduced.h5ad",
        )
    )

    # %% Clustering
    rsc.get.anndata_to_GPU(filtered_adata)

    # Uses pre-calculated scVI latent representation for calculating similarity score and constructing neighborhood graph
    # Store settings in .uns["neighbors"] and connectivity matrices in .obsp
    print("Constructing neighborhood graph...")
    rsc.pp.neighbors(filtered_adata, n_neighbors=15, use_rep="X_scvi", metric="cosine")
    # Embeds the neighborhood graph into 2D space using UMAP algorithm (optimized via SGD, Stochastic Gradient Descent)
    # You can change n_components to 3 for 3D UMAP
    print("Calculating UMAP...")
    rsc.tl.umap(filtered_adata, n_components=2, min_dist=0.3)
    # If the clusters appear too clumped or merged, try decreasing n_neighbors and min_dist

    # Clustering cells using leiden algorithm, maximizes modularity which is defined based on its intergroup connectivity and expected (random) connectivity
    # High resolution value results in more clusters
    print("Clustering with Leiden algorithm...")
    rsc.tl.leiden(filtered_adata, resolution=leiden_resolution)
    plot.plot_umap(filtered_adata, series_name)

    clustered_h5ad_dir = find_env_dir("CLUSTERED_H5AD")

    filtered_adata.write_h5ad(
        os.path.join(
            clustered_h5ad_dir,
            series_name + ".h5ad",
        )
    )

    print("Clustering completed")
