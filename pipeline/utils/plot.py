import matplotlib

# Supports file saving only; GUI rendering is not available
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import numpy as np
import os
import pandas as pd
import re
import scanpy as sc
import scvi
import seaborn as sns
import textwrap
import warnings
from adjustText import adjust_text
from anndata import AnnData
from matplotlib.axes import Axes
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr, spearmanr
from typing import List, Dict, Any, Optional
from pipeline.utils.env import find_env_dir
from pipeline.utils.pseudobulk import pseudobulk
from pipeline.config.constants import FIG_FORMAT

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
        os.path.join(validation_loss_plots_dir, f"{series_name}_{file_info}.{FIG_FORMAT}"),
        format=FIG_FORMAT,
        bbox_inches="tight",
    )
    plt.close("all")


# %% Visualizes sample quality metrics across multiple samples
def plot_qc(adata: AnnData, series_name: str, max_cells_per_sample: int = 5000) -> None:
    qc_ridgeplots_dir = find_env_dir("QC_RIDGEPLOTS")
    qc_ridgeplots_dir = os.path.join(qc_ridgeplots_dir, series_name)
    os.makedirs(
        qc_ridgeplots_dir,
        exist_ok=True,
    )

    assert isinstance(adata.obs, pd.DataFrame)
    rng = np.random.default_rng(0)
    sample_quality = adata.obs.sort_values("sample").copy()
    sample_quality["rand"] = rng.random(len(sample_quality))

    sample_quality = (
        sample_quality
        .sort_values(["sample", "rand"])
        .groupby("sample", group_keys=False, observed=True)
        .head(max_cells_per_sample)
        .drop(columns="rand")
        .reset_index(drop=True)
    )

    variables = [
        "pct_counts_mt",
        "pct_counts_ribo",
        "n_genes_by_counts",
        "log1p_total_counts",
    ]
    pretty_names = {
        "pct_counts_mt": "Mitochondrial Fraction (%)",
        "pct_counts_ribo": "Ribosomal Fraction (%)",
        "n_genes_by_counts": "Detected Genes",
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
        sns_grid.map_dataframe(
            sns.kdeplot,
            x=variable,
            clip_on=False,
            fill=True,
            alpha=0.8,
            linewidth=1.5,
            warn_singular=False,
            bw_adjust=0.8,
            cut=0,
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
            f"{variable}.{FIG_FORMAT}",
        )
        plt.savefig(filename, format=FIG_FORMAT, bbox_inches="tight")
        plt.close("all")
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

        filename = f"umap_{config['title']}.{FIG_FORMAT}"
        save_path = os.path.join(umap_plots_dir, filename)
        plt.savefig(save_path, format=FIG_FORMAT, bbox_inches="tight")

        plt.close(fig)


def plot_dotplot(
    adata: AnnData, series_name: str, target_genes_dict: Dict[str, List[str]], group: str, 
    filter_dict: Optional[dict] = None,
    project: Optional[str] = None,
    is_pseudobulk: bool = False,
) -> None:
    if filter_dict is None:
        filter_dict = {}

    dotplots_dir = find_env_dir("DOTPLOTS")
    if is_pseudobulk:
        series_name = series_name + "_pseudobulk"

    suffix = "_".join(target_genes_dict.keys())
    if filter_dict:
        flattened = [item for sublist in filter_dict.values() for item in sublist]
        suffix = "_".join(flattened) + "_" + suffix

    if project is not None:
        dotplots_dir = os.path.join(dotplots_dir, project)

    dotplots_dir = os.path.join(dotplots_dir, "_".join([series_name, suffix]))
    os.makedirs(
        dotplots_dir,
        exist_ok=True,
    )

    mask = pd.Series(True, index=adata.obs.index)
    for col, values in filter_dict.items():
        if col not in adata.obs.columns:
            raise ValueError(f"Column '{col}' not found in adata.obs")
        
        if isinstance(values, str):
            values = [values]
        mask &= adata.obs[col].isin(values)
    adata = adata[mask]

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

    if is_pseudobulk:
        group_keys = ["sample", group]
        pb_counts_df, pb_meta = pseudobulk(adata, group_keys=group_keys)
        pb_counts_target = pb_counts_df[all_target_genes]

        core_adata = AnnData(
            X=pb_counts_target.to_numpy().astype(np.float32),
            obs=pb_meta.copy(),
            var=pd.DataFrame(index=all_target_genes)
        )

        library_size = pb_meta["library_size"].to_numpy().astype(np.float32)

        assert isinstance(core_adata.obs, pd.DataFrame)
        core_adata.obs[group] = core_adata.obs[group].astype(str).astype('category')
    
    else:
        target_gene_expression = adata[:, all_target_genes].X
        if isinstance(target_gene_expression, csr_matrix):
            target_gene_expression = target_gene_expression.toarray()
        else:
            target_gene_expression = np.asarray(target_gene_expression)
            
        assert isinstance(adata.X, csr_matrix)
        library_size = np.array(adata.X.sum(axis=1)).ravel()
        
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
    plt.savefig(os.path.join(dotplots_dir, f"Var_dotplot.{FIG_FORMAT}"), format=FIG_FORMAT, bbox_inches="tight")

    sc.pl.dotplot(
        core_adata,
        var_names=target_genes_dict,
        groupby=group,
        standard_scale=None,  # None: Absolute expression values, "var": Relative to gene
        show=False,
        use_raw=False,
        var_group_rotation=0,
    )
    plt.savefig(os.path.join(dotplots_dir, f"None_dotplot.{FIG_FORMAT}"), format=FIG_FORMAT, bbox_inches="tight")
    plt.close("all")


def plot_violin(
    adata: AnnData,
    series: str,
    gene: list[str],
    group: str = "leiden",
    group_order: list | None = None, 
) -> None:
    out_dir = os.path.join(find_env_dir("VIOLIN_PLOTS"), series)
    os.makedirs(out_dir, exist_ok=True)

    data = adata.copy()
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)

    df = sc.get.obs_df(data, keys=[group, *gene], use_raw=False)
    df[group] = df[group].astype("category").cat.remove_unused_categories()
    groups = df[group].cat.categories

    if group_order is None:
        try:
            group_order = sorted(groups, key=lambda x: int(x))
        except (ValueError, TypeError):
            group_order = sorted(groups)
    else: 
        missing = [g for g in groups if g not in group_order]
        extra = [g for g in group_order if g not in groups]
        order_index = pd.Index(group_order)
        duplicated = order_index[order_index.duplicated()].tolist()

        if missing or extra or duplicated:
            raise ValueError(
                f"Invalid group_order: missing={missing}, "
                f"extra={extra}, duplicated={duplicated}"
            )
    ymax = max(1, int(np.ceil(df[gene].to_numpy().max())))

    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
    }), sns.axes_style("ticks"):
        fig, axes = plt.subplots(
            len(gene), 1,
            sharex=True, sharey=True, squeeze=False,
            figsize=(
                max(5, df[group].nunique() * 0.8 + 1.5),
                max(2.5, len(gene) * 0.7 + 1),
            ),
            layout="constrained",
        )

        colors = sns.husl_palette(len(gene), s=0.55, l=0.60)

        for ax, g, color in zip(axes[:, 0], gene, colors):
            sns.violinplot(
                data=df, x=group, y=g, ax=ax,
                order=group_order, 
                color=color, linecolor="#333333",
                inner=None, cut=0,
                density_norm="width", linewidth=0.6,
            )
            ax.set(
                xlabel="", ylabel="",
                ylim=(0, ymax * 1.05), yticks=[0, ymax],
            )
            ax.text(
                1.02, 0.5, g, transform=ax.transAxes,
                ha="left", va="center", fontstyle="italic",
            )
            ax.tick_params(axis="both", labelsize=9, length=3, width=0.6)
            sns.despine(ax=ax)

        axes[-1, 0].tick_params(axis="x", labelrotation=0)
        fig.supylabel("Log-normalized expression", fontsize=10)

        fig.savefig(
            os.path.join(out_dir, f"violin_by_{group}.{FIG_FORMAT}"),
            format=FIG_FORMAT, dpi=300,
            bbox_inches="tight", facecolor="white",
        )
        plt.close(fig)

def plot_proportions(
    adata: AnnData,
    series: str,
    group_key: str,
    sample_key: str,
    exclude_group: list[str] | None = None,
    group_order: list | None = None,
    sample_order: list | None = None,
    connect: bool = False,
) -> None:
    out_dir = os.path.join(find_env_dir("PROPORTION_PLOTS"), series)
    os.makedirs(out_dir, exist_ok=True)

    if exclude_group is None:
        exclude_group = ["Doublet", "LowCount", "LowQuality"]

    obs_df = adata.obs[~adata.obs[group_key].isin(exclude_group)]
    prop_df = pd.crosstab(
        obs_df[group_key].to_numpy(),
        obs_df[sample_key].to_numpy(),
        normalize="index",
    )

    if prop_df.empty:
        raise ValueError("No valid data available for plotting")

    if group_order is not None:
        missing = [g for g in prop_df.index if g not in group_order]
        extra = [g for g in group_order if g not in prop_df.index]
        order_index = pd.Index(group_order)
        duplicated = order_index[order_index.duplicated()].tolist()

        if missing or extra or duplicated:
            raise ValueError(
                f"group_order missing: missing={missing}, "
                f"extra groups={extra}, duplicated={duplicated}"
            )
        prop_df = prop_df.loc[group_order]
    else:
        try:
            prop_df = prop_df.loc[sorted(prop_df.index, key=lambda x: int(x))]
        except (ValueError, TypeError):
            prop_df = prop_df.sort_index()

    if sample_order is not None:
        missing = [s for s in prop_df.columns if s not in sample_order]
        extra = [s for s in sample_order if s not in prop_df.columns]
        order_index = pd.Index(sample_order)
        duplicated = order_index[order_index.duplicated()].tolist()

        if missing or extra or duplicated:
            raise ValueError(
                f"Invalid sample_order: missing={missing}, "
                f"extra={extra}, duplicated={duplicated}"
            )

        prop_df = prop_df.loc[:, sample_order]

    palette = [
        "#587FA2", "#D5A15B", "#77A69A", "#B98191",
        "#9487B2", "#9AA6AF", "#C5BE89", "#C48D70",
    ]
    n_colors = prop_df.shape[1]
    colors = (
        palette[:n_colors]
        if n_colors <= len(palette)
        else sns.husl_palette(n_colors, s=0.55, l=0.65)
    )

    width = 0.7
    fig, ax = plt.subplots(figsize=(10, 6))

    prop_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        width=width,
        edgecolor="none",
    )

    if connect:
        bounds = np.column_stack([
            np.zeros(len(prop_df)),
            prop_df.cumsum(axis=1).to_numpy(),
        ])
        for i in range(len(prop_df) - 1):
            for j, color in enumerate(colors):
                ax.fill_between(
                    [i + width / 2, i + 1 - width / 2],
                    bounds[i:i + 2, j],
                    bounds[i:i + 2, j + 1],
                    color=color,
                    alpha=0.5,
                    linewidth=0,
                    zorder=0,
                )

    ax.set_title(f"Proportion of {sample_key} by {group_key}", fontsize=14)
    ax.set_xlabel(group_key.capitalize(), fontsize=12)
    ax.set_ylabel("Proportion", fontsize=12)
    ax.set_ylim(0, 1)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.spines[["top", "right"]].set_visible(False)

    ax.legend(
        title=sample_key,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
    )

    fig.tight_layout()
    suffix = "_connected" if connect else ""
    filename = f"proportion_{group_key}_by_{sample_key}{suffix}.{FIG_FORMAT}"
    fig.savefig(
        os.path.join(out_dir, filename),
        format=FIG_FORMAT,
        bbox_inches="tight",
    )
    plt.close(fig)

def plot_volcano(
        deg: pd.DataFrame,
        series: str,
        name: str,
        genes: List[str],
        lfc_threshold: float = 0.5,
        dot_size: int = 3,
        text_size: int = 16,
        xlim=(-5, 5),
        ylim=(0, 50),
    ):
    volcano_plots_dir = find_env_dir("VOLCANO_PLOTS")
    out_dir = os.path.join(volcano_plots_dir, series)
    os.makedirs(out_dir, exist_ok=True)

    df = deg.dropna(subset=["log2FoldChange_shrunk", "padj"]).copy()
    df["logp"] = -np.log10(df["padj"].clip(lower=1e-300))

    LFC = "log2FoldChange_shrunk"
    up = (df["padj"] < 0.05) & (df[LFC] > lfc_threshold)
    down = (df["padj"] < 0.05) & (df[LFC] < -lfc_threshold)

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(df.loc[~(up | down), LFC], df.loc[~(up | down), "logp"],
               s=dot_size, c="lightgray", label="Not Sig")

    ax.scatter(df.loc[down, LFC], df.loc[down, "logp"],
               s=dot_size, c="cornflowerblue", label="Down")

    ax.scatter(df.loc[up, LFC], df.loc[up, "logp"],
               s=dot_size, c="firebrick", label="Up")

    ax.axhline(-np.log10(0.05), c="gray", ls="--", lw=1)
    ax.axvline(-lfc_threshold, c="gray", ls="--", lw=1)
    ax.axvline(lfc_threshold, c="gray", ls="--", lw=1)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    visible = df[df[LFC].between(*xlim) & df["logp"].between(*ylim)]
    texts, xs, ys = [], [], []

    dx = (xlim[1] - xlim[0]) * 0.12
    dy = (ylim[1] - ylim[0]) * 0.10

    for gene in genes:
        hit = df[df["gene"].str.upper() == gene.upper()]  
        if hit.empty:
            print(f"Gene not found: {gene}")
            continue

        row = hit.iloc[0]
        x, y = row[LFC], row["logp"]

        if not (xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]):
            print(f"Gene outside plot limits: {gene} (x={x:.2f}, y={y:.2f})")
            continue
        xs.append(x)
        ys.append(y)

        ax.scatter(x, y, s=30, c="black", zorder=5)
        texts.append(ax.text(
            x + (dx if x >= 0 else -dx),
            min(y + dy, ylim[1] - 0.05 * (ylim[1] - ylim[0])), 
            row["gene"],
            fontsize=text_size,
            fontweight="normal",
            ha="center", va="center", zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.3),
        ))

    ax.set_xlabel("log2 Fold Change Shrunk")
    ax.set_ylabel("-log10(padj)")
    legend = ax.legend(frameon=False)
    fig.tight_layout()

    if texts:
        adjust_text(
            texts,
            x=xs, y=ys, 
            target_x=xs, target_y=ys,
            objects=[legend],
            ax=ax,
            expand=(1.3, 1.6), 
            force_text=(0.8, 1.0), 
            force_pull=0, 
            ensure_inside_axes=True,
            min_arrow_len=0,
            iter_lim=1000, 
            arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8),
        )
    
    fig.savefig(
        os.path.join(out_dir, f"{name}.{FIG_FORMAT}"),
        format=FIG_FORMAT,
        bbox_inches="tight"
    )
    plt.close(fig)

def plot_gsea(
    result: pd.DataFrame,
    series: str,
    name: str,
    library: str | None = None,
    keywords: list[str] | None = None,
    fdr: float = 0.05,
    direction: str | None = None, 
    top_n: int | None = 10, 
) -> pd.DataFrame:
    if direction not in (None, "positive", "negative"):
        raise ValueError("direction must be None, 'positive', or 'negative'")
    if top_n is not None and top_n < 1:
        raise ValueError("top_n must be a positive integer or None.")

    out_dir = os.path.join(find_env_dir("ENRICHMENT"), series)
    os.makedirs(out_dir, exist_ok=True)

    df = result.copy()
    if library is not None:
        df = df[df["Library"] == library].copy()
    if df["Library"].nunique() > 1:
        raise ValueError("Select one library with library=...")

    df[["NES", "FDR q-val"]] = df[["NES", "FDR q-val"]].astype(float)
    df = df.dropna(subset=["NES", "FDR q-val"])

    if keywords is None:
        keywords = ["cell"]

    names = df["Term"].str.replace("_", " ", regex=False)
    related = names.str.contains(
        "|".join(map(re.escape, keywords)), case=False, na=False,
    )

    df = df[related & (df["FDR q-val"] < fdr)] 
    if direction == "positive":
        df = df[df["NES"] > 0]
    elif direction == "negative":
        df = df[df["NES"] < 0]

    df = df.sort_values("NES", key=lambda s: s.abs(), ascending=False)
    if top_n is not None:
        df = df.head(top_n)

    if df.empty:
        raise ValueError("No pathways pass the keyword, direction, and FDR filters.")

    labels = (
        df["Term"]
        .str.replace(r"^GOBP_", "", regex=True)
        .str.replace("_", " ", regex=False)
        .str.replace(r"\s*\(GO:\d+\)$", "", regex=True)
        .str.capitalize()
        .map(lambda s: textwrap.fill(s, width=40) if isinstance(s, str) else "")
    )

    fig, ax = plt.subplots(
        figsize=(8, max(3, len(df) * 0.6)),
        layout="constrained",
    )
    y = np.arange(len(df))
    colors = np.where(df["NES"] < 0, "#4C72B0", "#C44E52") 
    ax.barh(y, df["NES"], color=colors, height=0.7)
    ax.set_yticks(y, labels=labels)
    ax.invert_yaxis() 
    ax.axvline(0, color="gray", lw=0.8) 

    ax.set_xlabel("Normalized enrichment score (NES)")
    ax.set_title(f"{name.replace('_', ' ')} | FDR < {fdr}")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=10)

    suffix = f"_{direction}" if direction else "" 
    filename = f"gsea_{name}{suffix}"

    fig.savefig(
        os.path.join(out_dir, f"{filename}.{FIG_FORMAT}"),
        format=FIG_FORMAT, dpi=300, bbox_inches="tight",
    )
    plt.close(fig)
    return df

def plot_species_scatter(
    human_deg: pd.DataFrame,
    mouse_deg: pd.DataFrame,
    series: str,
    name: str,
    orthologs: pd.DataFrame | None = None,
    genes: list[str] | None = None,
    lfc_threshold: float = 0.5,
    fdr: float = 0.05,
    text_size: int = 8,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
) -> pd.DataFrame:
    out_dir = os.path.join(find_env_dir("SCATTER_PLOTS"), series)
    os.makedirs(out_dir, exist_ok=True)

    cols = ["gene", "log2FoldChange_shrunk", "padj"]
    human = human_deg[cols].set_axis(
        ["human_gene", "human_lfc", "human_padj"], axis=1,
    )
    mouse = mouse_deg[cols].set_axis(
        ["mouse_gene", "mouse_lfc", "mouse_padj"], axis=1,
    )

    if orthologs is None: 
        human["_match"] = human["human_gene"].str.upper()
        mouse["_match"] = mouse["mouse_gene"].str.upper()

        df = (
            human.dropna(subset=["_match"])
            .merge(
                mouse.dropna(subset=["_match"]),
                on="_match",
                how="inner", 
                validate="one_to_one",
            )
            .drop(columns="_match")
        )
    else: 
        mapping = orthologs[["human_gene", "mouse_gene"]].dropna().drop_duplicates()

        if mapping["human_gene"].duplicated().any() or mapping["mouse_gene"].duplicated().any():
            raise ValueError("orthologs must contain one-to-one gene pairs.")

        df = (
            mapping
            .merge(human, on="human_gene", how="inner", validate="one_to_one")
            .merge(mouse, on="mouse_gene", how="inner", validate="one_to_one")
        )

    numeric = ["human_lfc", "mouse_lfc", "human_padj", "mouse_padj"]
    df[numeric] = df[numeric].astype(float).replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["human_lfc", "mouse_lfc"]).copy()

    x, y = df["human_lfc"], df["mouse_lfc"]
    if len(df) < 3 or x.nunique() < 2 or y.nunique() < 2:
        raise ValueError("At least 3 matched genes with nonconstant LFC values are required.")

    sig = (df["human_padj"] < fdr) & (df["mouse_padj"] < fdr)
    up = sig & (x > lfc_threshold) & (y > lfc_threshold)
    down = sig & (x < -lfc_threshold) & (y < -lfc_threshold)
    opposite = (
        sig & (x * y < 0)
        & (x.abs() > lfc_threshold) & (y.abs() > lfc_threshold)
    )
    df["status"] = np.select(
        [up, down, opposite, sig],
        ["Up in both", "Down in both", "Opposite", "LFC cutoff not met"],
        default="padj cutoff not met",
    )

    r = float(pearsonr(x, y).statistic) #type: ignore
    rho = float(spearmanr(x, y).statistic) #type: ignore

    fig, ax = plt.subplots(figsize=(7, 7))
    palette = {
        "padj cutoff not met": "#D0D0D0",
        "LFC cutoff not met": "#929BA5",
        "Down in both": "#4C72B0",
        "Up in both": "#C44E52",
        "Opposite": "#DD9853",
    }
    for label, color in palette.items():
        d = df[df["status"] == label]
        ax.scatter(
            d["human_lfc"], d["mouse_lfc"],
            s=8, color=color, alpha=0.7, edgecolors="none",
            rasterized=True, label=f"{label} (n={len(d):,})",
        )

    lim = max(1, float(np.abs(df[["human_lfc", "mouse_lfc"]].to_numpy()).max()) * 1.15)
    xlim = xlim if xlim is not None else (-lim, lim)  # 추가
    ylim = ylim if ylim is not None else (-lim, lim)  # 추가

    if xlim[0] >= xlim[1] or ylim[0] >= ylim[1]:
        raise ValueError("Axis limits must satisfy min < max.")

    ax.axline((0, 0), slope=1, ls=":", color="#BBBBBB", lw=0.8, zorder=0)
    ax.axhline(0, color="#BBBBBB", lw=0.6, zorder=0)
    ax.axvline(0, color="#BBBBBB", lw=0.6, zorder=0)

    for cutoff in (-lfc_threshold, lfc_threshold):
        ax.axhline(cutoff, color="gray", ls="--", lw=0.8, zorder=0)
        ax.axvline(cutoff, color="gray", ls="--", lw=0.8, zorder=0)

    ax.set(
        xlim=xlim, ylim=ylim,
        xlabel="Human log2 fold change (shrunk)",
        ylabel="Mouse log2 fold change (shrunk)",
        title=name.replace("_", " "),
    )
    ax.set_aspect("equal")

    info = ax.text(
        0.03, 0.97,
        f"n = {len(df):,}\nPearson r = {r:.2f}\nSpearman ρ = {rho:.2f}",
        transform=ax.transAxes, va="top", fontsize=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9),
    )

    wanted = {g.upper() for g in genes or []}
    marked = df[
        df["human_gene"].str.upper().isin(wanted)
        | df["mouse_gene"].str.upper().isin(wanted)
    ]
    found = set(marked["human_gene"].str.upper()) | set(marked["mouse_gene"].str.upper())
    for g in sorted(wanted - found):
        print(f"Gene not found among matched pairs: {g}")

    visible = (
        marked["human_lfc"].between(*xlim)
        & marked["mouse_lfc"].between(*ylim)
    )
    for g in marked.loc[~visible, "human_gene"]:
        print(f"Gene outside plot limits: {g}")
    marked = marked.loc[visible]

    ax.scatter(
        marked["human_lfc"], marked["mouse_lfc"],
        s=45,
        c=marked["status"].map(palette).tolist(),
        edgecolors="black", linewidths=0.8, zorder=5,
    )

    texts = [
        ax.text(
            row.human_lfc, row.mouse_lfc, str(row.human_gene), #type: ignore
            fontsize=text_size, fontweight="normal", zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.3),
        )
        for row in marked.itertuples()
    ]

    ax.legend(
        bbox_to_anchor=(1.02, 1), loc="upper left",
        frameon=False, markerscale=2,
    )
    fig.tight_layout()

    if texts:
        adjust_text(
            texts,
            x=marked["human_lfc"].to_numpy(),
            y=marked["mouse_lfc"].to_numpy(),
            objects=[info], ax=ax,
            expand=(1.3, 1.6), force_text=(0.8, 1.0),
            force_pull=0, min_arrow_len=0, iter_lim=1000,
            arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8),
        )

    fig.savefig(
        os.path.join(out_dir, f"scatter_{name}.{FIG_FORMAT}"),
        format=FIG_FORMAT, dpi=300, bbox_inches="tight",
    )
    plt.close(fig)

    df.to_csv(os.path.join(out_dir, f"scatter_{name}.csv"), index=False)
    return df