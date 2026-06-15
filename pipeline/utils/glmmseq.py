import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
from typing import Union, Iterable
from anndata import AnnData
from rpy2.robjects import pandas2ri
from rpy2.robjects import conversion, default_converter
from pipeline.utils.env import find_env
from pipeline.utils.pseudobulk import pseudobulk
from pipeline.config.constants import CPU_CORE_COUNT

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

def pseudobulk_glmmseq_comp(
    adata: AnnData, 
    condition_col: str, group_test: str, group_control: str, 
    random_effect_col: str, # column name for random effect (e.g., patient ID)
    covariates: Union[str, Iterable[str]] = [],
    filter_dict: dict = {},
    min_cells: int = 30,
    test_method: str = "Wald"
) -> pd.DataFrame:
    if test_method not in {"Wald", "LRT"}:
        raise ValueError("test_method must be either 'Wald' or 'LRT'")
    
    group_keys = ["sample", condition_col]
    covariates = [covariates] if isinstance(covariates, str) else list(covariates)
    
    cond_mask = np.array(adata.obs[condition_col].isin([group_test, group_control]), dtype=bool)
    for col, vals in filter_dict.items():
        cond_mask &= adata.obs[col].isin(vals)
        
    counts_df, pb_meta = pseudobulk(adata[cond_mask], group_keys=group_keys, min_cells=min_cells)
    counts_df = pd.DataFrame(counts_df.sparse.to_dense())
    
    min_count_per_sample = 2
    min_samples = int(np.ceil(counts_df.shape[0] * 0.10))
    keep_genes = (counts_df >= min_count_per_sample).sum(axis=0) >= min_samples
    pb_counts = counts_df.loc[:, keep_genes]

    keep_genes_var = pb_counts.var(axis=0) > 0
    pb_counts = pb_counts.loc[:, keep_genes_var]
    
    assert isinstance(adata.obs, pd.DataFrame)
    obs_meta = adata.obs.drop_duplicates(subset=group_keys).copy()
    obs_meta["pseudobulk_id"] = obs_meta[group_keys].astype(str).agg("_".join, axis=1)
    meta = obs_meta.set_index("pseudobulk_id").loc[pb_counts.index]
    meta["n_cells"] = pb_meta["n_cells"]
    meta["library_size"] = pb_meta["library_size"]
    
    design_cols = list(set(covariates + [condition_col, random_effect_col]))
    ok = meta[design_cols].notna().all(axis=1)
    meta = meta.loc[ok]

    for col in design_cols:
        if meta[col].dtype == 'object' or meta[col].dtype.name == 'category':
            meta[col] = meta[col].astype(str).str.replace(r'[^a-zA-Z0-9_.]', '_', regex=True).astype('category')
    pb_counts = pb_counts.loc[meta.index]
 
    formula_str_fixed = "~ " + " + ".join(covariates + [condition_col]) if covariates else f"~ {condition_col}"
    formula_str_glmm = formula_str_fixed + f" + (1 | {random_effect_col})"
    formula_str_reduced = "~ " + " + ".join(covariates) if covariates else "~ 1"
    formula_str_reduced = formula_str_reduced + f" + (1 | {random_effect_col})"

    # In R, the count matrix should be in the format [genes x samples], so we apply Transpose (.T)
    count_data_r = pb_counts.T 

    with open(os.path.join(find_env("ROOT_DIR"), "R", "glmmseq.R"), "r") as f:
        r_script = f.read()

    with conversion.localconverter(default_converter + pandas2ri.converter):
        ro.globalenv['countdata'] = count_data_r
        ro.globalenv['meta_df'] = meta
        ro.globalenv['formula_str_fixed'] = formula_str_fixed
        ro.globalenv['formula_str_glmm'] = formula_str_glmm
        ro.globalenv['formula_str_reduced'] = formula_str_reduced
        ro.globalenv['condition_col'] = condition_col
        ro.globalenv['group_control'] = group_control
        ro.globalenv['group_test'] = group_test
        ro.globalenv['random_effect_col'] = random_effect_col
        ro.globalenv['n_cores'] = CPU_CORE_COUNT
        ro.globalenv['test_method'] = test_method
        
        ro.r(r_script)
        result_df = ro.globalenv['res_df']

    return result_df
