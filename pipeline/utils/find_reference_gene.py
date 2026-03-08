import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy.stats import gmean

def find_pseudobulk_reference_genes(
    adata: AnnData,  
    min_sample_fraction: float = 0.7, 
    min_mean_expr: float = 10.0,
    min_cells_per_sample: int = 100
) -> pd.DataFrame:
    sample_counts = adata.obs["sample"].value_counts()
    valid_samples = sample_counts[sample_counts >= min_cells_per_sample].index
    
    if len(valid_samples) == 0:
        raise ValueError(f"No samples with at least {min_cells_per_sample} cells")
        
    adata = adata[adata.obs["sample"].isin(valid_samples)]
    samples = adata.obs["sample"].astype("category")
    samples = samples.cat.remove_unused_categories()
    
    row_idx = samples.cat.codes.values
    col_idx = np.arange(adata.n_obs)
    val = np.ones(adata.n_obs, dtype=np.float32)
    n_samples = len(samples.cat.categories)

    M = csr_matrix((val, (row_idx, col_idx)), shape=(n_samples, adata.n_obs))
    pb_X = M.dot(adata.X) 
    pb_dense = np.asarray(pb_X.todense())

    assert isinstance(adata.X, csr_matrix)
    X_binary = adata.X.copy()
    X_binary.data = np.ones_like(X_binary.data)
    cells_expressing = M.dot(X_binary)
    cells_expressing_dense = np.asarray(cells_expressing.todense())
    
    n_cells_per_sample = np.asarray(M.sum(axis=1)).flatten()
    cell_expr_fraction_per_sample = cells_expressing_dense / n_cells_per_sample[:, np.newaxis]
    mean_cell_expr_fraction = cell_expr_fraction_per_sample.mean(axis=0)

    expressed_fraction = (pb_dense > 0).mean(axis=0)
    
    robust_mask = expressed_fraction >= min_sample_fraction
    pb_robust = pb_dense[:, robust_mask] + 1
    pseudo_reference = gmean(pb_robust, axis=0)

    ratios = pb_robust / pseudo_reference
    size_factors = np.median(ratios, axis=1, keepdims=True)
    size_factors[size_factors == 0] = 1.0
    pb_normalized = pb_dense / size_factors

    gene_mean = pb_normalized.mean(axis=0)
    gene_std = pb_normalized.std(axis=0)

    with np.errstate(divide='ignore', invalid='ignore'):
        gene_cv = gene_std / gene_mean
        
    stability_df = pd.DataFrame({
        'Expressed_Fraction': expressed_fraction,
        'Mean_Abs_Expression': gene_mean,
        'Std_Deviation': gene_std,
        'CV': gene_cv,
        'Mean_Cell_Expressed_Fraction': mean_cell_expr_fraction,
    }, index=adata.var_names)
    
    stability_df = stability_df.replace([np.inf, -np.inf], np.nan).dropna()
    stability_df = stability_df[stability_df['Expressed_Fraction'] >= min_sample_fraction]
    stability_df = stability_df[stability_df['Mean_Abs_Expression'] >= min_mean_expr]
    stability_df = stability_df.sort_values('CV', ascending=True)
    
    return stability_df
