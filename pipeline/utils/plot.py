import matplotlib

# Supports file saving only; GUI rendering is not available
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.axes import Axes
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import warnings
from anndata import AnnData
from scipy.sparse import spmatrix
from typing import List, Dict, Any, Optional
from pipeline.config.directory import (
    DOTPLOTS_DIR,
    QC_RIDGEPLOTS_DIR,
    UMAP_PLOTS_DIR,
    VALIDATION_LOSS_PLOTS_DIR,
)

# Setting seaborn global theme
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

warnings.filterwarnings("ignore", message="Tight layout not applied")


# Plots validation loss over training epochs
def plot_validation_loss(
    model: scvi.model.SCVI | scvi.external.SOLO, filename: str
) -> None:
    os.makedirs(VALIDATION_LOSS_PLOTS_DIR, exist_ok=True)

    if model.history is None:
        raise ValueError("Model history is not available.")

    plt.figure(figsize=(10, 6))
    plt.xlabel("Epoch", fontweight="bold")

    if isinstance(model, scvi.model.SCVI):
        elbo_history = model.history["elbo_validation"]
        plt.plot(elbo_history, label="ELBO Validation Loss", linewidth=2)
        plt.ylabel("ELBO Loss", fontweight="bold")
        plt.title("ELBO Validation Loss Over Epochs", fontweight="bold")
    elif isinstance(model, scvi.external.SOLO):
        elbo_history = model.history["validation_loss"]
        plt.plot(elbo_history, label="Validation Loss", linewidth=2)
        plt.ylabel("Validation Loss", fontweight="bold")
        plt.title("Validation Loss Over Epochs", fontweight="bold")

    plt.legend(frameon=False)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(
        os.path.join(VALIDATION_LOSS_PLOTS_DIR, filename),
        format="svg",
        bbox_inches="tight",
    )
    plt.close()


# %% Visualizes sample quality metrics across multiple samples
def plot_qc(adata: AnnData) -> None:
    os.makedirs(QC_RIDGEPLOTS_DIR, exist_ok=True)

    assert isinstance(adata.obs, pd.DataFrame)
    sample_quality = adata.obs.sort_values("sample")

    variables = [
        "singlet_probability",
        "pct_counts_mt",
        "n_genes_by_counts",
        "pct_counts_in_top_20_genes",
        "log1p_total_counts",
    ]
    pretty_names = {
        "singlet_probability": "Singlet Probability",
        "pct_counts_mt": "Mitochondrial Fraction (%)",
        "n_genes_by_counts": "Detected Genes",
        "pct_counts_in_top_20_genes": "Library Complexity (Top 20%)",
        "log1p_total_counts": "Sequencing Depth (Log1p UMI)",
    }
    series_name = adata[0].obs["series"].iloc[0]

    for variable in variables:
        # Initializes seaborn FacetGrid
        sns_grid = sns.FacetGrid(
            sample_quality,
            row="sample",
            hue="sample",
            aspect=15,
            height=0.6,
            palette="tab20",
            sharex=True,
        )
        # Draw KDE (Kernel Density Estimation) plot
        sns_grid.map(
            sns.kdeplot, variable, clip_on=False, fill=True, alpha=0.8, linewidth=1.5
        )
        # Outlining the KDE plot with white line
        sns_grid.map(sns.kdeplot, variable, clip_on=False, color="w", linewidth=2)
        # Depicting a y axis
        sns_grid.map(plt.axhline, y=0, linewidth=2, clip_on=False)

        # Write sample names on the left side of each KDE plot
        def label(_, color, label):
            kde_plot = plt.gca()
            text = kde_plot.text(
                1,
                0.2,
                label,
                fontweight="bold",
                color=color,
                ha="right",
                va="center",
                transform=kde_plot.transAxes,
            )
            text.set_path_effects([patheffects.withStroke(linewidth=3, foreground="w")])

        sns_grid.map(label, variable)

        # Allow subplots to overlap
        sns_grid.figure.subplots_adjust(hspace=-0.4)
        # Remove unnecessary subplot details
        sns_grid.set_titles("")
        # Remove y-axis ticks and labels
        sns_grid.set(yticks=[], ylabel="")
        # Remove spines (Outer box of each subplot)
        sns_grid.despine(bottom=True, left=True)

        median_val = sample_quality[variable].median()
        for i, kde_plot in enumerate(sns_grid.axes.flat):
            kde_plot: Axes = kde_plot
            # Draw median line (Red)
            kde_plot.axvline(
                x=median_val,
                color="#d62728",
                linestyle="-",
                alpha=1.0,
                linewidth=1.5,
                zorder=0,
            )
            if i == 0:
                text_obj = kde_plot.text(
                    median_val,
                    0.8,
                    f"{median_val:.2f}",
                    ha="left",
                    color="#d62728",
                    transform=kde_plot.get_xaxis_transform(),
                    fontsize=10,
                    fontweight="bold",
                    clip_on=False,
                )
                text_obj.set_path_effects(
                    [patheffects.withStroke(linewidth=2, foreground="black")]
                )

            # Special logic: Draw threshold line only for Singlet Probability
            if variable == "singlet_probability":
                kde_plot.axvline(x=0.6, color="#1f77b4", linestyle="-", linewidth=1.5)
                if i == 0:
                    cutoff_text = kde_plot.text(
                        0.6,
                        1.1,
                        " Cutoff (0.6)",
                        color="#1f77b4",
                        transform=kde_plot.get_xaxis_transform(),
                        fontsize=7,
                        fontweight="bold",
                        clip_on=False,
                    )
                    cutoff_text.set_path_effects(
                        [patheffects.withStroke(linewidth=2, foreground="black")]
                    )

            # Set X-axis label
            kde_plot.set_xlabel(pretty_names.get(variable, variable))

        filename = os.path.join(
            QC_RIDGEPLOTS_DIR,
            f"{series_name}_{variable}.svg",
        )
        plt.savefig(filename, format="svg", bbox_inches="tight")
        plt.close()


# If celltype information is already annotated in adata.obs["celltype"], set has_celltype=True to visualize
# If you want to highlight specific cell types, provide their names in highlight_cells list
# highlight_cells only work when has_celltype=True
# additional_config can be used to provide extra plotting configurations
def plot_umap(
    adata: AnnData,
    has_celltype: bool = False,
    highlight_cells: Optional[List[str]] = None,
    additional_config: Optional[List[Dict[str, Any]]] = None,
    dot_size: int = 7,
) -> None:
    os.makedirs(UMAP_PLOTS_DIR, exist_ok=True)

    assert isinstance(adata.obs, pd.DataFrame)

    idx = np.random.permutation(adata.n_obs)
    plot_adata = AnnData(obs=adata.obs.iloc[idx].copy())
    plot_adata.obsm["X_umap"] = adata.obsm["X_umap"][idx].copy()

    series_name = plot_adata.obs["series"].iloc[0]

    plot_configs = [
        {
            "color": "leiden",
            "title": "Leiden Clustering",
            "suffix": "_Leiden_Clustering",
            "legend_loc": "on data",
            "palette": None,
        },
        {
            "color": "sample",
            "title": "Sample Distribution",
            "suffix": "_Sample_Distribution",
            "legend_loc": "best",
            "palette": None,
        },
    ]

    if has_celltype:
        plot_configs.append(
            {
                "color": "celltype",
                "title": "Cell Type Distribution",
                "suffix": "_Cell_Type_Distribution",
                "legend_loc": "best",
                "palette": None,
            }
        )
    if has_celltype and highlight_cells is not None:
        cell_types = plot_adata.obs["celltype"].unique()

        for cell in highlight_cells:
            custom_palette = {ct: "lightgray" for ct in cell_types}
            if cell in custom_palette:
                custom_palette[cell] = "#FF0000"
            else:
                print(f"Warning: '{cell}' not found in cell types.")

            plot_configs.append(
                {
                    "color": "celltype",
                    "title": f"{cell} Distribution",
                    "suffix": f"_{cell}_Highlight",
                    "legend_loc": "best",
                    "palette": custom_palette,
                },
            )

    if additional_config is not None:
        plot_configs.extend(additional_config)

    for config in plot_configs:
        fig, ax = plt.subplots(figsize=(16, 10))

        sc.pl.umap(
            plot_adata,
            color=config["color"],
            projection="2d",
            ax=ax,
            palette=config["palette"],
            legend_loc=config["legend_loc"],
            legend_fontoutline=2,
            size=dot_size,
            frameon=False,
            show=False,
        )
        ax.set_title(config["title"], fontweight="bold", fontsize=24)

        legend = ax.get_legend()
        if legend:
            legend.set_frame_on(False)

        filename = f"{series_name}_umap_{config['suffix']}.svg"
        save_path = os.path.join(UMAP_PLOTS_DIR, filename)
        plt.savefig(save_path, format="svg", bbox_inches="tight")

        plt.close(fig)


def plot_dotplot(adata: AnnData, target_genes_dict: Dict[str, List[str]]) -> None:
    os.makedirs(DOTPLOTS_DIR, exist_ok=True)
    series_name = adata.obs["series"].iloc[0]

    all_target_genes = []
    for genes in target_genes_dict.values():
        all_target_genes.extend(genes)

    if len(all_target_genes) != len(set(all_target_genes)):
        seen = set()
        duplicates = set()
        for gene in all_target_genes:
            if gene in seen:
                duplicates.add(gene)
            seen.add(gene)
        raise ValueError(
            f"Duplicate genes found in target_genes_dict: {list(duplicates)}. Please ensure each gene appears only once."
        )

    missing_genes = [gene for gene in all_target_genes if gene not in adata.var_names]
    if missing_genes:
        raise ValueError(f"Such genes are missing in the data: {missing_genes}")

    if "total_counts" in adata.obs:
        library_size = adata.obs["total_counts"].to_numpy()
    else:
        raise ValueError(
            "total_counts not found in adata.obs, quality check did not run properly"
        )

    target_gene_expression = adata[:, all_target_genes].X
    if isinstance(target_gene_expression, spmatrix):
        target_gene_expression = target_gene_expression.toarray()  # type: ignore
    else:
        target_gene_expression = np.asarray(target_gene_expression)

    assert isinstance(adata.obs, pd.DataFrame)
    core_adata = AnnData(
        X=target_gene_expression,
        obs=adata.obs[["leiden"]].copy(),
        var=pd.DataFrame(index=all_target_genes),
    )

    library_size = np.where(library_size == 0, 1, library_size)
    core_adata.X = (core_adata.X / library_size[:, np.newaxis]) * 1e4
    sc.pp.log1p(core_adata)

    sc.pl.dotplot(
        core_adata,
        var_names=target_genes_dict,
        groupby="leiden",
        standard_scale="var",  # None: Absolute expression values, "var": Relative to gene
        show=False,
        use_raw=False,
    )

    suffix = "_".join(target_genes_dict.keys())
    filename = f"{series_name}_{suffix}_dotplot.svg"
    plt.savefig(os.path.join(DOTPLOTS_DIR, filename), format="svg", bbox_inches="tight")
    plt.close()
