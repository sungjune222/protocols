import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Iterable, Tuple
from scipy.sparse import csr_matrix


def pseudobulk(
    adata: AnnData,
    group_keys: Iterable[str],
    min_cells: int = 15,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sep = "__"
    keys = list(group_keys)
    for k in keys:
        if k not in adata.obs.columns:
            raise ValueError(f"Key '{k}' not found in adata.obs")

    assert isinstance(adata.obs, pd.DataFrame)
    group_df = adata.obs[keys].astype(str)
    group_id = group_df.agg(sep.join, axis=1)

    sizes = group_id.value_counts()
    keep_groups = sizes[sizes >= min_cells].index
    keep_mask = group_id.isin(keep_groups).to_numpy()

    assert isinstance(adata.X, csr_matrix)
    X = adata.X[keep_mask]
    group_id = group_id[keep_mask]
    group_df = group_df.loc[keep_mask]

    cat = pd.Categorical(group_id)
    n_groups = len(cat.categories)

    rows = np.arange(X.shape[0])
    cols = cat.codes
    G = csr_matrix(
        (np.ones(X.shape[0], dtype=np.int8), (rows, cols)), shape=(X.shape[0], n_groups)
    )

    pb = (G.T @ X).astype(np.int64)

    counts_df = pd.DataFrame(
        pb.toarray(),
        index=pd.Index(cat.categories, name="pseudobulk_id"),
        columns=adata.var_names,
    )

    meta_df = pd.DataFrame(index=counts_df.index)
    split = meta_df.index.to_series().str.split(sep, expand=True)
    split.columns = keys
    for k in keys:
        meta_df[k] = split[k].values

    meta_df["n_cells"] = sizes.loc[meta_df.index].values
    meta_df["library_size"] = counts_df.sum(axis=1).values

    return counts_df, meta_df
