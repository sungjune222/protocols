import anndata
import argparse
import contextlib
import io
import os
import re
import sys
import scanpy as sc
import scipy.io
from scipy.sparse import csr_matrix
from sccomposite import RNA_modality

anndata.settings.allow_write_nullable_strings = True

class Args(argparse.Namespace):
    input: str
    outdir: str
    sample_id: str


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description="Run COMPOSITE for Doublet Detection (Optimized for Sparse Matrix)"
    )
    parser.add_argument(
        "--input", required=True, help="Path to CellBender output .h5 file"
    )
    parser.add_argument(
        "--outdir", required=True, help="Directory to save output files"
    )
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    return parser.parse_args(namespace=Args())


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.chdir(args.outdir)

    print(f"Processing input: [{args.sample_id}]")

    try:
        adata = sc.read_10x_h5(args.input)
        adata.var_names_make_unique()

        if not isinstance(adata.X, csr_matrix):
            print(f"[{args.sample_id}] Converting Dense matrix to Sparse CSR matrix...")
            adata.X = csr_matrix(adata.X)

        target_mtx = "RNA.mtx"
        print(
            f"[{args.sample_id}] Exporting {adata.n_obs} cells to Sparse Matrix Market format..."
        )
        scipy.io.mmwrite(target_mtx, adata.X.T)

        print(f"[{args.sample_id}] Running COMPOSITE model...")
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            doublet_class, consistency = RNA_modality.composite_rna(target_mtx)
        stdout_str = f.getvalue()
        print(stdout_str)

        # 4. GOF 점수 파싱
        gof_match = re.search(r"goodness-of-fit score is:\s*([0-9\.]+)", stdout_str)
        gof_score = float(gof_match.group(1)) if gof_match else 0.0

        with open(f"{args.sample_id}_gof_score.txt", "w") as file:
            file.write(str(gof_score))

        adata.obs["composite_class"] = doublet_class
        adata.obs["composite_consistency"] = consistency

        adata_clean = adata[adata.obs["composite_class"] == 0].copy()
        adata_clean.write_h5ad(f"{args.sample_id}_clean.h5ad", compression="gzip")

        print(
            f"[{args.sample_id}] Done. Clean: {adata_clean.n_obs}/{adata.n_obs}, GOF: {gof_score}"
        )

    except Exception as e:
        print(f"Error occurred during processing {args.sample_id}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
