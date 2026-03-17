import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Iterable, Tuple
from scipy.sparse import csr_matrix

def pseudobulk(
    adata: AnnData,
    group_keys: Iterable[str],
    min_cells: int = 30,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    keys = list(group_keys)
    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        raise ValueError(f"Keys not found in adata.obs: {missing}")

    assert isinstance(adata.obs, pd.DataFrame)
    group_df = adata.obs.loc[:, keys]
    group_index = pd.MultiIndex.from_frame(group_df, names=keys)
    codes, uniques = group_index.factorize(sort=False)

    sizes = np.bincount(codes, minlength=len(uniques))
    keep_group = sizes >= min_cells
    keep_mask = keep_group[codes]

    assert isinstance(adata.X, csr_matrix)
    X = adata.X[keep_mask]
    codes = codes[keep_mask]

    kept_old_codes = np.flatnonzero(keep_group)
    remap = np.full(len(uniques), -1, dtype=np.int32)
    remap[kept_old_codes] = np.arange(len(kept_old_codes), dtype=np.int32)
    new_codes = remap[codes]
    kept_uniques = uniques[kept_old_codes]

    n_groups = len(kept_old_codes)
    n_cells = X.shape[0]

    G = csr_matrix(
        (
            np.ones(n_cells, dtype=np.int8),
            (new_codes, np.arange(n_cells, dtype=np.int64)),
        ),
        shape=(n_groups, n_cells),
    )
    pb = (G @ X).astype(np.int64, copy=False)

    assert isinstance(kept_uniques, pd.MultiIndex)
    meta_df = kept_uniques.to_frame(index=False)
    
    pseudobulk_id = pd.Index(
        meta_df.astype(str).agg("_".join, axis=1),
        name="pseudobulk_id",
    )

    meta_df.index = pseudobulk_id
    meta_df["n_cells"] = sizes[kept_old_codes]
    meta_df["library_size"] = np.asarray(pb.sum(axis=1)).ravel()

    counts_df = pd.DataFrame.sparse.from_spmatrix(
        pb,
        index=pseudobulk_id,
        columns=adata.var_names,
    )

    return counts_df, meta_df
