import matplotlib
matplotlib.use("Agg")

import argparse
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pathlib import Path
from sklearn.decomposition import PCA

class Args(argparse.Namespace):
    counts: str
    gtf: str
    output_dir: str
    project: str
    min_count: int
    n_hvg: int
    z_threshold: float

def parse_args() -> Args:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--min_count", type=int, default=5)
    parser.add_argument("--n_hvg", type=int, default=2000)
    parser.add_argument("--z_threshold", type=float, default=3.0)

    return parser.parse_args(namespace=Args())

def main():
    args = parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(args.counts, sep="\t", comment="#")

    meta = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    df = raw.drop(columns=meta[1:]).set_index("Geneid")
    df.columns = [Path(x).stem for x in df.columns]

    if df.shape[1] < 2:
        raise ValueError("At least two samples are required.")

    gene_map = {}
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            x = line.rstrip().split("\t")
            if len(x) < 9 or x[2] != "gene":
                continue

            gid = re.search(r'gene_id "([^"]+)"', x[8])
            name = re.search(r'(?:gene_name|gene) "([^"]+)"', x[8])
            if gid:
                gene_map[gid.group(1)] = (name.group(1) if name else gid.group(1))

    names = pd.Series(df.index, index=df.index).map(gene_map)
    names = names.fillna(pd.Series(df.index, index=df.index))

    mito = names.str.lower().str.startswith("mt-")
    ribo = names.str.lower().str.startswith(("rpl", "rps", "mrpl", "mrps"))

    # Basic QC
    library_size = df.sum()
    detected_genes = (df >= args.min_count).sum()
    mito_pct = df.loc[mito].sum() / library_size * 100
    ribo_pct = df.loc[ribo].sum() / library_size * 100

    df = df[(df >= args.min_count).sum(axis=1) >= 2]
    counts = df.T.astype(int)
    dds = DeseqDataSet(
        counts=counts,
        metadata=pd.DataFrame(index=counts.index),
        design="~1",
        quiet=True,
    )
    dds.vst(use_design=False)
    vst = dds.to_df(layer="vst_counts")

    hvg = vst.var().nlargest(min(args.n_hvg, vst.shape[1])).index
    X = vst[hvg]
    corr = X.T.corr()
    mean_corr = (corr.sum() - 1) / (len(corr) - 1)

    qc = pd.DataFrame({
        "library_size": library_size,
        "detected_genes": detected_genes,
        "mito_pct": mito_pct,
        "ribo_pct": ribo_pct,
        "mean_correlation": mean_corr,
    })
    qc.to_csv(out_dir / f"{args.project}_qc.tsv", sep="\t")

    # Correlation heatmap
    plt.figure(figsize=(8, 7))
    sns.heatmap(corr, cmap="vlag", square=True)
    plt.tight_layout()
    plt.savefig(out_dir / f"{args.project}_correlation.png", dpi=300)
    plt.close()

    # PCA
    pca = PCA(2)
    pcs = pca.fit_transform(X)

    pca_df = pd.DataFrame(
        pcs,
        index=X.index,
        columns=["PC1", "PC2"],
    )
    pca_df["condition"] = (
        pca_df.index.str.replace(r"\d+$", "", regex=True)
    )

    ev = pca.explained_variance_ratio_ * 100

    plt.figure(figsize=(8, 7))
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue="condition",
        s=100,
    )

    for sample, row in pca_df.iterrows():
        plt.text(row["PC1"], row["PC2"], sample, fontsize=8) #type: ignore

    plt.xlabel(f"PC1 ({ev[0]:.1f}%)")
    plt.ylabel(f"PC2 ({ev[1]:.1f}%)")
    plt.tight_layout()
    plt.savefig(out_dir / f"{args.project}_pca.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
