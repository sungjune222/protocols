import numpy as np
import pandas as pd
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from typing import Iterable, Union
from pipeline.utils.pseudobulk import pseudobulk
from pipeline.config.constants import CPU_CORE_COUNT


def pseudobulk_deseq2(
    adata: AnnData, group_keys: Union[str, Iterable[str]], min_cells: int = 15
) -> pd.DataFrame:
    sep = "__"
    group_keys = [group_keys] if isinstance(group_keys, str) else list(group_keys)

    for k in group_keys:
        if k not in adata.obs.columns:
            raise ValueError(f"Key '{k}' not found in adata.obs")

    sc_counts, sc_meta = pseudobulk(adata, ["sample", *group_keys], min_cells=1)
    sc_ncells = sc_meta["n_cells"]
    s_counts, s_meta = pseudobulk(adata, group_keys=["sample"], min_cells=1)
    s_ncells = s_meta["n_cells"]

    samples = adata.obs["sample"].unique().tolist()
    group_id_series = sc_meta[group_keys].astype(str).agg(sep.join, axis=1)
    group_ids = group_id_series.unique().tolist()

    results = []

    for group_id in group_ids:
        group_parts = group_id.split(sep)
        group_kv = dict(zip(group_keys, group_parts))

        rows, meta_rows, idx = [], [], []

        for sample in samples:
            key_s = f"{sample}"
            key_sc = f"{sample}{sep}{group_id}"

            if key_s not in s_counts.index:
                continue

            total_counts = s_counts.loc[key_s].to_numpy()
            total_ncells = int(s_ncells.loc[key_s])

            if key_sc in sc_counts.index:
                group_counts = sc_counts.loc[key_sc].to_numpy()
                group_ncells = int(sc_ncells.loc[key_sc])
            else:
                group_counts = np.zeros_like(total_counts)
                group_ncells = 0

            rest_counts = total_counts - group_counts
            rest_ncells = total_ncells - group_ncells

            if group_ncells < min_cells or rest_ncells < min_cells:
                continue

            rows.append(group_counts)
            meta_rows.append(
                {
                    "sample": sample,
                    "cluster": "group",
                    "batch": sample,
                    "cluster_id": group_id,
                    **group_kv,
                }
            )
            idx.append(f"{sample}_{group_id}")

            rows.append(rest_counts)
            meta_rows.append(
                {
                    "sample": sample,
                    "cluster": "rest",
                    "batch": sample,
                    "cluster_id": group_id,
                    **group_kv,
                }
            )
            idx.append(f"{sample}_{group_id}_rest")

        if len(rows) < 4:
            print("Not enough samples for group ", group_id)
            continue

        counts = pd.DataFrame(rows, index=idx, columns=adata.var_names)
        meta = pd.DataFrame(meta_rows, index=idx)
        meta["cluster"] = pd.Categorical(meta["cluster"], categories=["rest", "group"])
        meta["batch"] = pd.Categorical(meta["batch"])

        keep_genes = counts.sum(axis=0) > 0
        counts = counts.loc[:, keep_genes]

        dds = DeseqDataSet(
            counts=counts,
            metadata=meta,
            design="~batch + cluster",
            refit_cooks=True,
            n_cpus=CPU_CORE_COUNT,
        )
        dds.deseq2()

        ds = DeseqStats(
            dds, contrast=["cluster", "group", "rest"], n_cpus=CPU_CORE_COUNT
        )
        ds.summary()

        res = ds.results_df.copy()

        ds.lfc_shrink(coeff="cluster[T.group]", adapt=True)
        res_shrink = ds.results_df
        res["log2FoldChange_shrunk"] = res_shrink["log2FoldChange"]
        res["lfcSE_shrunk"] = res_shrink["lfcSE"]

        res["cluster_id"] = group_id
        res["contrast"] = "rest"
        res["group_keys"] = ",".join(group_keys)
        res = res.reset_index().rename(columns={"index": "gene"})
        results.append(res)

    return pd.concat(results, axis=0, ignore_index=True)
