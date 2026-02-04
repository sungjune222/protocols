import numpy as np
import os
import scanpy as sc
from anndata import AnnData
from cellrank import pl as cr_pl
from cellrank.estimators import GPCCA
from cellrank.kernels import PseudotimeKernel, ConnectivityKernel
from cellrank.models import GAM
from pipeline.utils.env import find_env_dir

 # %% CellRank-based Trajectory Inference and Gene Trend Analysis
def cellrank_analysis(adata: AnnData, series_name: str, start_cluster: str, target_gene: str, group_key: str = "leiden"):
    cellrank_analysis_location = find_env_dir("CELLRANK")
    cellrank_analysis_location = os.path.join(cellrank_analysis_location, series_name)
    os.makedirs(cellrank_analysis_location, exist_ok=True)

    # Calculating connectivity graph between clusters using PAGA
    sc.tl.paga(adata, groups=group_key)
    # Calculating Diffusion Pseudotime (DPT) for ordering cells along developmental trajectories
    root_cell_index = np.where(adata.obs[group_key] == start_cluster)[0][0]
    adata.uns["iroot"] = root_cell_index
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    # CellRank-based Trajectory Inference and Gene Trend Analysis
    pseudotime_kernel = PseudotimeKernel(adata, time_key="dpt_pseudotime")
    pseudotime_kernel.compute_transition_matrix()

    connectivity_kernel = ConnectivityKernel(adata)
    connectivity_kernel.compute_transition_matrix()

    kernel = 0.8 * pseudotime_kernel + 0.2 * connectivity_kernel
    kernel.compute_transition_matrix()
    kernel.plot_projection(
        basis="umap",
        color=group_key,
        title="Model-based Differentiation Flow",
        save=os.path.join(cellrank_analysis_location, "differentiation_flow_umap.svg"),
    )

    #  Generalized Perron Cluster Cluster Analysis (GPCCA)
    gpcca = GPCCA(kernel)
    gpcca.fit(cluster_key=group_key, n_states=None)
    gpcca.predict_terminal_states()
    gpcca.predict_initial_states(allow_overlap=True)

    gpcca.compute_fate_probabilities()  # type: ignore
    gpcca.plot_fate_probabilities(  # type: ignore
        same_plot=True,
        basis="umap",
        save=os.path.join(cellrank_analysis_location, "fate_probabilities_umap.svg"),
    )

    assert(gpcca.initial_states is not None)
    assert (gpcca.terminal_states is not None)
    init = gpcca.initial_states.copy() 
    term = gpcca.terminal_states

    trend_model = GAM(adata)
    lineages = gpcca.terminal_states.cat.categories

    cr_pl.gene_trends(
        adata,
        model=trend_model,
        genes=[target_gene],
        time_key="dpt_pseudotime",
        lineages=lineages,  # type: ignore
        data_key="X",
        same_plot=True,
        hide_cells=True,
        figsize=(10, 6),
        save=os.path.join(cellrank_analysis_location, f"gene_trends_{target_gene}.svg"),
    )

    drivers = gpcca.compute_lineage_drivers()  # type: ignore
    target_lineage = lineages[0]
    print(drivers.head(10))