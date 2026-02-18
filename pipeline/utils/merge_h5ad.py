import anndata
import argparse
import glob
import os
import sys
import scanpy as sc
from scipy.sparse import csr_matrix

anndata.settings.allow_write_nullable_strings = True

class Args(argparse.Namespace):
    input_dir: str
    project_id: str
    output_dir: str


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description="Merge clean h5ad files into a single project file"
    )
    parser.add_argument(
        "--input_dir", required=True, help="Directory containing sample subdirectories"
    )
    parser.add_argument(
        "--project_id", required=True, help="Project ID for naming the output file"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to save the merged file"
    )
    return parser.parse_args(namespace=Args())


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"=== Merging samples for Project: {args.project_id} ===")

    try:
        search_pattern = os.path.join(args.input_dir, "**", "*_clean.h5ad")
        files = sorted(glob.glob(search_pattern, recursive=True))

        if not files:
            print(f"Error: No '*_clean.h5ad' files found in {args.input_dir}")
            sys.exit(1)
        print(f"Found {len(files)} files to merge.")

        adatas = []
        for file_path in files:
            try:
                file_name = os.path.basename(file_path)
                sample_id = file_name.replace("_clean.h5ad", "")
                
                print(f"  -> Loading {sample_id} from {file_name}...")
                adata = sc.read_h5ad(file_path)
                
                adata.obs['sample_id'] = sample_id
                adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
                
                if not adata.obs_names[0].endswith(f"-{sample_id}"):
                    adata.obs_names = adata.obs_names + '-' + sample_id
                
                adatas.append(adata)
                
            except Exception as e:
                print(f"Warning: Failed to load {file_path}. Skipping. Error: {e}")

        if not adatas:
            print("Error: No valid AnnData objects loaded.")
            sys.exit(1)

        print("Concatenating objects...")
        combined_adata = anndata.concat(adatas, join='outer', merge="same")

        if not isinstance(combined_adata.X, csr_matrix):
            print(f"[{args.project_id}] Converting Dense matrix to Sparse CSR matrix...")
            combined_adata.X = csr_matrix(combined_adata.X)

        output_filename = f"{args.project_id}.h5ad"
        output_path = os.path.join(args.output_dir, output_filename)
        
        print(f"Saving merged file to {output_path}...")
        combined_adata.write_h5ad(output_path, compression="gzip")
        
        print(f"=== Merge Completed Successfully. Total Cells: {combined_adata.n_obs} ===")

    except Exception as e:
        print(f"Critical Error occurred during merging: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()