import gseapy as gp
import numpy as np
import os
import pandas as pd
from pipeline.utils.env import find_env_dir
from pipeline.config.constants import CPU_CORE_COUNT

def gsea(deg: pd.DataFrame, species:str, series: str):
    assert species in ["human", "mouse"], f"Species should be 'human' or 'mouse'"
    
    enrichment_dir = find_env_dir("ENRICHMENT")
    enrichment_dir = os.path.join(enrichment_dir, series)
    os.makedirs(enrichment_dir, exist_ok=True)

    if "gene" in deg.columns:
        deg = deg.set_index("gene")

    deg["stat"] = deg["log2FoldChange"] / deg["lfcSE"].replace(0, np.nan)
    clean_deg = deg.replace([np.inf, -np.inf], np.nan).dropna(subset=["stat"])

    contrast = clean_deg["contrast"].dropna().unique()
    assert len(contrast) == 1
    contrast = contrast[0]

    libraries = [
        "GO_Biological_Process_2026",
        "KEGG_2026",        
    ] # https://maayanlab.cloud/Enrichr/#libraries

    rank = clean_deg["stat"].sort_values(ascending=False).rename("score")

    results = []
    for lib in libraries:
        prerank = gp.prerank(
            rnk=rank,
            gene_sets=lib,
            organism=species,
            permutation_num=2000,
            min_size=10,
            max_size=500,
            seed=0,
            threads=CPU_CORE_COUNT,
        )

        assert prerank.res2d is not None
        res = prerank.res2d.sort_values("FDR q-val")
        res["Library"] = lib

        res.to_csv(
            os.path.join(
                enrichment_dir,
                f"{series}_gsea_{contrast}_{lib}.csv",
            )
        )
        results.append(res)
    
    return pd.concat(results, ignore_index=True)


def ora(deg: pd.DataFrame, species: str, series: str):
    assert species in ["human", "mouse"], f"Species should be 'human' or 'mouse'"
    
    enrichment_dir = find_env_dir("ENRICHMENT")
    enrichment_dir = os.path.join(enrichment_dir, series)
    os.makedirs(enrichment_dir, exist_ok=True)

    if "gene" in deg.columns:
        deg = deg.set_index("gene")

    clean_deg = deg.replace([np.inf, -np.inf], np.nan).dropna(subset=["padj", "log2FoldChange_shrunk"])

    contrast = clean_deg["contrast"].dropna().unique()
    assert len(contrast) == 1
    contrast = contrast[0]

    libraries = [
        "GO_Biological_Process_2026",
        "KEGG_2026",        
    ] # https://maayanlab.cloud/Enrichr/#libraries

    gene_lists = {
        "UP": clean_deg[
            (clean_deg["padj"] < 0.05)
            & (clean_deg["log2FoldChange_shrunk"] > 1)
        ].index.unique().tolist(),

        "DOWN": clean_deg[
            (clean_deg["padj"] < 0.05)
            & (clean_deg["log2FoldChange_shrunk"] < -1)
        ].index.unique().tolist(),
    }

    results = []
    
    for direction, gene_list in gene_lists.items():
        print(f"[ORA] {contrast} | {direction} genes: {len(gene_list)}")

        if len(gene_list) < 10:
            print("Not enough genes for enrichment. Skipping.")
            continue

        for lib in libraries:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism=species,
            )

            if enr.results is None or enr.results.empty: #type: ignore
                continue

            res = enr.results.sort_values("Adjusted P-value") #type: ignore
            res["Library"] = lib
            res["Direction"] = direction
            res["Contrast"] = contrast
            res["Series"] = series

            res.to_csv(
                os.path.join(
                    enrichment_dir,
                    f"{series}_ORA_{contrast}_{direction}_{lib}.csv",
                ),
                index=False,
            )

            results.append(res)

    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()