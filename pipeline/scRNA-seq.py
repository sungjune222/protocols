# %% Environment Setup
import anndata
import gc
import gseapy as gp
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from anndata import AnnData
from pipeline.utils import plot
from pipeline.utils.env import find_env_dir
from pipeline.config.constants import (
    SINGLE_CELL_VAE_BATCH_SIZE,
    CPU_CORE_COUNT,
)
from pipeline.config.machine_learning import DataLoader
from pipeline.utils.deseq2 import pseudobulk_deseq2

anndata.settings.allow_write_nullable_strings = True

if __name__ == "__main__":
    # Loading .env
    h5ad_matrix_location = find_env_dir("H5AD_MATRIX")
    series_name = "Zheng_opc"
    covariates = ["series"]
    group_keys = ["condition"]
    file = os.path.join(h5ad_matrix_location, series_name + ".h5ad")

    # %% Preprocessing functions

    # Quality assessment by calculating QC metrics
    def quality_assess(adata: AnnData, ribo_genes: pd.DataFrame) -> AnnData:
        # Marking mitochondrial and ribosomal genes
        mt_mask = adata.var.index.str.lower().str.startswith("mt-")
        ribo_mask = adata.var_names.isin(ribo_genes[0].values)
        adata.var["mt"] = np.asarray(mt_mask, dtype=bool)
        adata.var["ribo"] = np.asarray(ribo_mask, dtype=bool)

        # Generates QC metrics
        # qc_vars: List of categories that you want to make as a QC metrics (It must be set as a boolean list in AnnData.obs)
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo"], percent_top=[20], log1p=True, inplace=True
        )

        return adata

    ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)

    # Preprocessing each sample
    print("Loading data...")
    loaded_adata = sc.read_h5ad(file)

    print("Assessing quality...")
    quality_assessed_adata = quality_assess(loaded_adata, ribo_genes)

    # Cell and gene filtering
    sc.pp.filter_cells(quality_assessed_adata, min_genes=200)
    sc.pp.filter_genes(quality_assessed_adata, min_cells=10)
    # Filtered cells and genes necessitate re-calculation of QC metrics
    quality_assessed_adata = quality_assess(quality_assessed_adata, ribo_genes)

    plot.plot_qc(quality_assessed_adata, series_name)

    # Cytoplasmic RNA in dead cells leaks out, resulting in a higher proportion of remaining mitochondrial RNA
    quality_assessed_adata = quality_assessed_adata[
        (quality_assessed_adata.obs["pct_counts_mt"] < 15)
    ].copy()

    filtered_adata = quality_assessed_adata
    gc.collect()

    filtered_count_matrix_location = find_env_dir("FILTERED_COUNT_MATRIX")
    filtered_adata.write_h5ad(
        os.path.join(filtered_count_matrix_location, series_name + ".h5ad")
    )

    # %% Producing latent representation of cells using scVI
    print("Producing latent representation of cells using scVI...")

    sc.pp.highly_variable_genes(
        filtered_adata,
        n_top_genes=4000,
        subset=False,
        flavor="seurat_v3",
    )
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
        batch_size=SINGLE_CELL_VAE_BATCH_SIZE,
        datasplitter_kwargs=DataLoader,
        train_size=0.9,
        check_val_every_n_epoch=1,
        max_epochs=500,
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

    scvi_model_location = find_env_dir("SCVI_MODEL")
    model.save(
        os.path.join(scvi_model_location, series_name),
        overwrite=True,
        save_anndata=True,
    )
    del scvi_adata
    gc.collect()

    compressed_count_matrix_location = find_env_dir("COMPRESSED_COUNT_MATRIX")
    filtered_adata.write_h5ad(
        os.path.join(
            compressed_count_matrix_location,
            series_name + "_compressed.h5ad",
        )
    )

    # %% Clustering
    # Uses pre-calculated scVI latent representation for calculating similarity score and constructing neighborhood graph
    # Store settings in .uns["neighbors"] and connectivity matrices in .obsp
    print("Constructing neighborhood graph...")
    sc.pp.neighbors(filtered_adata, n_neighbors=15, use_rep="X_scvi", metric="cosine")
    # Embeds the neighborhood graph into 2D space using UMAP algorithm (optimized via SGD, Stochastic Gradient Descent)
    # You can change n_components to 3 for 3D UMAP
    print("Calculating UMAP...")
    sc.tl.umap(filtered_adata, n_components=2, min_dist=0.3)
    # If the clusters appear too clumped or merged, try decreasing n_neighbors and min_dist

    # Clustering cells using leiden algorithm, maximizes modularity which is defined based on its intergroup connectivity and expected (random) connectivity
    # High resolution value results in more clusters
    print("Clustering with Leiden algorithm...")
    sc.tl.leiden(
        filtered_adata, resolution=1.5, flavor="igraph", n_iterations=-1, directed=False
    )
    plot.plot_umap(filtered_adata, series_name)

    # Differential expression analysis is not performed on scVI denoised data because of imputation artifacts
    # Instead, pseudobulk differential expression analysis using DESeq2 is performed
    print("Differential expression analysis...")

    clustered_data_location = find_env_dir("CLUSTERED_DATA")

    filtered_adata.write_h5ad(
        os.path.join(
            clustered_data_location,
            series_name + "_clustered.h5ad",
        )
    )

    # %% Differential Expression Analysis using Pseudobulk DESeq2
    de_analysis_location = find_env_dir("DESEQ")
    de_result = pseudobulk_deseq2(filtered_adata, group_keys=group_keys, covariates=covariates)
    de_result.to_csv(os.path.join(de_analysis_location, series_name + "_deseq2.csv"))

    # %% Gene Set Enrichment Analysis (GSEA)
    def gsea(de_result: pd.DataFrame):
        enrichment_analysis_location = find_env_dir("ENRICHMENT_ANALYSIS")
        enrichment_analysis_location = os.path.join(
            enrichment_analysis_location, series_name
        )
        os.makedirs(enrichment_analysis_location, exist_ok=True)

        if "gene" in de_result.columns:
            de_result = de_result.set_index("gene")

        clean_de = de_result.replace([np.inf, -np.inf], np.nan).dropna(subset=["stat"])
        clean_de.index = clean_de.index.astype(str).str.upper()

        libraries = [
            "GO_Biological_Process_2025",
            "GO_Cellular_Component_2025",
            "GO_Molecular_Function_2025",
            "KEGG_2026",
        ]

        for (cluster, contrast), group_de in clean_de.groupby(["group", "contrast"]):
            rank = group_de["stat"].sort_values(ascending=False).rename("score")

            for lib in libraries:
                prerank = gp.prerank(
                    rnk=rank,
                    gene_sets=lib,
                    permutation_num=2000,
                    min_size=10,
                    max_size=500,
                    seed=0,
                    threads=CPU_CORE_COUNT,
                    verbose=True,
                )

                assert prerank.res2d is not None
                prerank.res2d.sort_values("FDR q-val").to_csv(
                    os.path.join(
                        enrichment_analysis_location,
                        f"{series_name}_{cluster}vs{contrast}_{lib}.csv",
                    )
                )

    gsea(de_result)
