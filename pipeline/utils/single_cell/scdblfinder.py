import argparse
import os
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

class Args(argparse.Namespace):
    input: str
    singlets: str
    outdir: str
    sample_id: str


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description="Create clean h5ad from CellBender .h5 using singlet barcodes"
    )
    parser.add_argument(
        "--input", required=True, help="Path to CellBender output .h5 file"
    )
    parser.add_argument(
        "--singlets", required=True, help="CSV file containing singlet barcodes"
    )
    parser.add_argument(
        "--outdir", required=True, help="Directory to save output files"
    )
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    return parser.parse_args(namespace=Args())


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    print(f"Processing input: [{args.sample_id}]")

    adata = sc.read_10x_h5(args.input)
    adata.X = csr_matrix(adata.X)

    singlets = pd.read_csv(args.singlets)
    if "barcode" not in singlets.columns:
        raise ValueError(f"'barcode' column not found in {args.singlets}")

    singlet_barcodes = singlets["barcode"].astype(str)
    singlet_set = set(singlet_barcodes)

    keep_mask = adata.obs_names.astype(str).isin(singlet_set)
    adata_clean = adata[keep_mask].copy()

    if adata_clean.n_obs == 0:
        raise ValueError("No singlet barcodes matched the h5 barcodes.")

    output_path = os.path.join(args.outdir, f"{args.sample_id}_clean.h5ad")
    adata_clean.write_h5ad(output_path)
    

if __name__ == "__main__":
    main()