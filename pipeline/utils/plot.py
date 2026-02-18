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
from scipy.sparse import issparse, csr_matrix
from typing import List, Dict, Any, Optional
from pipeline.utils.env import find_env_dir

warnings.filterwarnings("ignore", message="Tight layout not applied")


# Plots validation loss over training epochs
def plot_validation_loss(
    model: scvi.model.SCVI | scvi.external.SOLO, series_name: str, file_info: str
) -> None:
    validation_loss_plots_dir = find_env_dir("VALIDATION_LOSS_PLOTS")

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
        os.path.join(validation_loss_plots_dir, f"{series_name}_{file_info}.svg"),
        format="svg",
        bbox_inches="tight",
    )
    plt.close()


# %% Visualizes sample quality metrics across multiple samples
def plot_qc(adata: AnnData, series_name: str) -> None:
    qc_ridgeplots_dir = find_env_dir("QC_RIDGEPLOTS")
    qc_ridgeplots_dir = os.path.join(qc_ridgeplots_dir, series_name)
    os.makedirs(
        qc_ridgeplots_dir,
        exist_ok=True,
    )

    assert isinstance(adata.obs, pd.DataFrame)
    sample_quality = adata.obs.sort_values("sample")

    variables = [
        "pct_counts_mt",
        "n_genes_by_counts",
        "pct_counts_in_top_20_genes",
        "log1p_total_counts",
    ]
    pretty_names = {
        "pct_counts_mt": "Mitochondrial Fraction (%)",
        "n_genes_by_counts": "Detected Genes",
        "pct_counts_in_top_20_genes": "Library Complexity (Top 20%)",
        "log1p_total_counts": "Sequencing Depth (Log1p UMI)",
    }

    # Setting seaborn global theme
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

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
            sns.kdeplot,
            variable,
            clip_on=False,
            fill=True,
            alpha=0.8,
            linewidth=1.5,
            warn_singular=False,
        )
        # Outlining the KDE plot with white line
        sns_grid.map(
            sns.kdeplot,
            variable,
            clip_on=False,
            color="w",
            linewidth=2,
            warn_singular=False,
        )
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
            qc_ridgeplots_dir,
            f"{variable}.svg",
        )
        plt.savefig(filename, format="svg", bbox_inches="tight")
        plt.close()
    sns.reset_defaults()


# If celltype information is already annotated in adata.obs["celltype"], set has_celltype=True to visualize
# If you want to highlight specific cell types, provide their names in highlight_cells list
# highlight_cells only work when has_celltype=True
# additional_config can be used to provide extra plotting configurations
def plot_umap(
    adata: AnnData,
    series_name: str,
    has_celltype: bool = False,
    highlight_cells: Optional[List[str]] = None,
    additional_config: Optional[List[Dict[str, Any]]] = None,
    dot_size: int = 7,
) -> None:
    umap_plots_dir = find_env_dir("UMAP_PLOTS")
    umap_plots_dir = os.path.join(umap_plots_dir, series_name)
    os.makedirs(
        umap_plots_dir,
        exist_ok=True,
    )
    assert isinstance(adata.obs, pd.DataFrame)

    idx = np.random.permutation(adata.n_obs)
    plot_adata = AnnData(obs=adata.obs.iloc[idx].copy())
    plot_adata.obsm["X_umap"] = adata.obsm["X_umap"][idx].copy()

    n_samples = len(adata.obs['sample'].unique())
    sample_palette = sns.color_palette("husl", n_samples)

    plot_configs = [
        {
            "color": "leiden",
            "title": "Leiden Clustering",
            "legend_loc": "on data",
            "palette": None,
        },
        {
            "color": "sample",
            "title": "Sample Distribution",
            "legend_loc": "best",
            "palette": sample_palette,
        },
    ]

    if has_celltype:
        plot_configs.append(
            {
                "color": "celltype",
                "title": "Cell Type Distribution",
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
                return

            plot_configs.append(
                {
                    "color": "celltype",
                    "title": f"{cell} Distribution",
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

        filename = f"umap_{config['title']}.svg"
        save_path = os.path.join(umap_plots_dir, filename)
        plt.savefig(save_path, format="svg", bbox_inches="tight")

        plt.close(fig)


def plot_dotplot(
    adata: AnnData, series_name: str, target_genes_dict: Dict[str, List[str]], group: str
) -> None:
    dotplots_dir = find_env_dir("DOTPLOTS")
    suffix = "_".join(target_genes_dict.keys())
    dotplots_dir = os.path.join(dotplots_dir, series_name + "_" + suffix)
    os.makedirs(
        dotplots_dir,
        exist_ok=True,
    )

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

    assert isinstance(adata.X, csr_matrix)
    if issparse(adata.X):
        library_size = np.array(adata.X.sum(axis=1)).ravel()
    else:
        library_size = adata.X.sum(axis=1)

    target_gene_expression = adata[:, all_target_genes].X
    if isinstance(target_gene_expression, csr_matrix):
        target_gene_expression = target_gene_expression.toarray()
    else:
        target_gene_expression = np.asarray(target_gene_expression)

    assert isinstance(adata.obs, pd.DataFrame)
    core_adata = AnnData(
        X=target_gene_expression,
        obs=adata.obs[[group]].copy(),
        var=pd.DataFrame(index=all_target_genes),
    )

    library_size = np.where(library_size == 0, 1, library_size)
    core_adata.X = (core_adata.X / library_size[:, np.newaxis]) * 1e4
    sc.pp.log1p(core_adata)

    sc.pl.dotplot(
        core_adata,
        var_names=target_genes_dict,
        groupby=group,
        standard_scale="var",  # None: Absolute expression values, "var": Relative to gene
        show=False,
        use_raw=False,
        var_group_rotation=0,
    )

    filename = "Var_dotplot.svg"
    plt.savefig(os.path.join(dotplots_dir, filename), format="svg", bbox_inches="tight")

    sc.pl.dotplot(
        core_adata,
        var_names=target_genes_dict,
        groupby=group,
        standard_scale=None,  # None: Absolute expression values, "var": Relative to gene
        show=False,
        use_raw=False,
        var_group_rotation=0,
    )

    suffix = "_".join(target_genes_dict.keys())
    filename = "None_dotplot.svg"
    plt.savefig(os.path.join(dotplots_dir, filename), format="svg", bbox_inches="tight")
    plt.close()


def plot_violin(adata: AnnData, gene: str) -> None:
    violin_plots_dir = find_env_dir("VIOLIN_PLOTS")
    # series_name = adata.obs["series"].iloc[0]
    series_name = "SCP1038"

    num_categories = len(adata.obs["leiden"].unique())
    width_per_category = 0.6
    fig_width = max(5, num_categories * width_per_category)
    fig_height = 6

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_yscale("symlog", linthresh=1)

    sc.pl.violin(
        adata,
        keys=gene,
        groupby="leiden",
        rotation=90,
        ax=ax,
        show=False,
        use_raw=False,
        stripplot=False,
        density_norm="count",
    )
    ax.set_ylim((0, 100))

    filename = f"{series_name}_{gene}_violin.svg"
    fig.savefig(
        os.path.join(violin_plots_dir, filename), format="svg", bbox_inches="tight"
    )
    plt.close()

def plot_proportions(
        adata: AnnData,
        series_name: str,
        group_key: str,
        sample_key: str
    ):
    proportions_plots_dir = find_env_dir("PROPORTION_PLOTS")
    proportions_plots_dir = os.path.join(proportions_plots_dir, series_name)
    os.makedirs(proportions_plots_dir, exist_ok=True)

    prop_df = pd.crosstab(
        adata.obs[group_key].to_numpy(), 
        adata.obs[sample_key].to_numpy(), 
        normalize='index'
    )

    fig, ax = plt.subplots(figsize=(10, 6))

    prop_df.plot(
        kind='bar', 
        stacked=True, 
        ax=ax, 
        colormap='tab20',
        edgecolor='none'
    )

    plt.title(f"Proportion of {sample_key} by {group_key}", fontsize=14)
    plt.xlabel(group_key.capitalize(), fontsize=12)
    plt.ylabel("Proportion", fontsize=12)
    plt.xticks(rotation=45, ha='right') 
    
    plt.legend(
        title=sample_key, 
        bbox_to_anchor=(1.05, 1), 
        loc='upper left', 
        borderaxespad=0.
    )

    plt.tight_layout()
    filename = f"proportion_{group_key}_by_{sample_key}.svg"
    fig.savefig(
        os.path.join(proportions_plots_dir, filename), format="svg", bbox_inches="tight"
    )
    plt.close()
