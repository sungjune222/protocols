import numpy as np
import pandas as pd
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from typing import Iterable, Union
from pipeline.utils.pseudobulk import pseudobulk
from pipeline.config.constants import CPU_CORE_COUNT


def pseudobulk_deseq2(
    adata: AnnData,  group_keys: Union[str, Iterable[str]], covariates: Union[str, Iterable[str]], min_cells: int = 15
) -> pd.DataFrame:
    sep = "__"
    group_keys = [group_keys] if isinstance(group_keys, str) else list(group_keys)
    covariates = [covariates] if isinstance(covariates, str) else list(covariates)

    for k in group_keys:
        if k not in adata.obs.columns:
            raise ValueError(f"Key '{k}' not found in adata.obs")
        
    target_col = "comparison_group"
    formula_terms = covariates + [target_col]
    design_formula = "~" + " + ".join(formula_terms)

    sc_counts, sc_meta = pseudobulk(adata, group_keys=["sample", *group_keys], min_cells=1)
    sc_ncells = sc_meta["n_cells"]
    s_counts, s_meta = pseudobulk(adata, group_keys=["sample"], min_cells=1)
    s_ncells = s_meta["n_cells"]

    samples = adata.obs["sample"].unique().tolist()

    sample_cov_dict = dict()
    for sample in samples:
        sample_obs = adata.obs[adata.obs["sample"] == sample]
        sample_cov_dict[sample] = sample_obs[covariates].iloc[0].to_dict()

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

            if group_ncells >= min_cells:
                rows.append(group_counts)
                meta_rows.append(
                    {
                        target_col: "group",
                        "cluster": group_id,
                        **group_kv,
                        **sample_cov_dict[sample]
                    }
                )
                idx.append(f"{sample}_{group_id}")
            
            if rest_ncells >= min_cells:
                rows.append(rest_counts)
                meta_rows.append(
                    {
                        target_col: "rest",
                        "cluster": group_id,
                        **group_kv,
                        **sample_cov_dict[sample]
                    }
                )
                idx.append(f"{sample}_{group_id}_rest")

        if len(rows) < 4:
            print("Not enough samples for group ", group_id)
            continue

        counts = pd.DataFrame(rows, index=idx, columns=adata.var_names)
        meta = pd.DataFrame(meta_rows, index=idx)
        meta[target_col] = pd.Categorical(meta[target_col], categories=["rest", "group"])

        keep_genes = counts.sum(axis=0) > 0
        counts = counts.loc[:, keep_genes]

        dds = DeseqDataSet(
            counts=counts,
            metadata=meta,
            design=design_formula,
            refit_cooks=True,
            n_cpus=CPU_CORE_COUNT,
        )
        dds.deseq2()

        ds = DeseqStats(
            dds, contrast=[target_col, "group", "rest"], n_cpus=CPU_CORE_COUNT
        )
        ds.summary()

        res = ds.results_df.copy()
        ds.lfc_shrink(coeff=f"{target_col}[T.group]", adapt=True)
        res["log2FoldChange_shrunk"] = ds.results_df["log2FoldChange"]
        res["lfcSE_shrunk"] = ds.results_df["lfcSE"]

        res["cluster"] = group_id
        res["contrast"] = "rest"
        res["group_keys"] = ",".join(group_keys)
        res = res.reset_index().rename(columns={"index": "gene"})
        results.append(res)

    return pd.concat(results, axis=0, ignore_index=True)


def pseudobulk_deseq2_comp(
    adata: AnnData, condition_col: str, group_test: str, group_control: str, covariates: Union[str, Iterable[str]], min_cells: int = 15
) -> pd.DataFrame:
    covariates = [covariates] if isinstance(covariates, str) else list(covariates)
    formula_terms = covariates + [condition_col]
    design_formula = "~" + " + ".join(formula_terms)

    cond_mask = np.array(adata.obs[condition_col].isin([group_test, group_control]), dtype=bool)

    counts_df, pb_meta = pseudobulk(adata[cond_mask], group_keys=["sample"], min_cells=min_cells)
    pb_counts = counts_df.loc[:, counts_df.sum(axis=0) > 0]

    assert isinstance(adata.obs, pd.DataFrame)
    meta = adata.obs.drop_duplicates(subset=["sample"]).set_index("sample").loc[pb_counts.index]
    meta["n_cells"] = pb_meta["n_cells"]
    meta["library_size"] = pb_meta["library_size"]
    
    for col in meta.columns:
        if meta[col].dtype == 'object' or meta[col].dtype.name == 'category':
            meta[col] = meta[col].astype(str).astype('category')

    dds = DeseqDataSet(
        counts=pb_counts,
        metadata=meta,
        design=design_formula,
        refit_cooks=True,
        n_cpus=CPU_CORE_COUNT
    )
    dds.deseq2()
    
    ds = DeseqStats(
        dds, contrast=[condition_col, group_test, group_control], n_cpus=CPU_CORE_COUNT
    )
    ds.summary()
    
    res = ds.results_df.copy()
    ds.lfc_shrink(coeff=f"{condition_col}[T.{group_test}]", adapt=True)
    res["log2FoldChange_shrunk"] = ds.results_df["log2FoldChange"]
    res["lfcSE_shrunk"] = ds.results_df["lfcSE"]

    res['contrast'] = f"{group_test}_vs_{group_control}"
    return res.reset_index().rename(columns={'index': 'gene'})

