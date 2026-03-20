import anndata
import argparse
import glob
import os
import pandas as pd
import sys
import scanpy as sc
import xml.etree.ElementTree as ET
from scipy.sparse import csr_matrix

class Args(argparse.Namespace):
    input_dir: str
    project_id: str
    output_dir: str
    sra_xml: str


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
    parser.add_argument(
        "--sra_xml", required=True, help="Path to the SRA XML metadata file"
    )
    return parser.parse_args(namespace=Args())

def add_value(d: dict[str, str], key: str, value: str) -> None:
    value = (value or "").strip()
    if not key or not value:
        return

    if key not in d or d[key] == "":
        d[key] = value
    elif value not in d[key].split(" | "):
        d[key] += f" | {value}"

def parse_sra_xml(xml_path: str) -> dict[str, dict[str, str]]:
    tree = ET.parse(xml_path)
    root = tree.getroot()
    run_to_meta: dict[str, dict[str, str]] = {}

    for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
        base_meta: dict[str, str] = {}

        sample = pkg.find("./SAMPLE")
        if sample is not None:
            add_value(base_meta, "alias", sample.get("alias", ""))

            for attr in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
                tag = attr.findtext("TAG", default="")
                value = attr.findtext("VALUE", default="")
                add_value(base_meta, tag, value)

        for run in pkg.findall("./RUN_SET/RUN"):
            run_acc = (run.get("accession") or "").strip()
            if not run_acc:
                continue
            run_to_meta[run_acc] = base_meta

    return run_to_meta


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"=== Merging samples for Project: {args.project_id} ===")

    try:
        if not os.path.isfile(args.sra_xml):
            print(f"Error: SRA XML file not found: {args.sra_xml}")
            sys.exit(1)
        sra_meta = parse_sra_xml(args.sra_xml)

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
                run_id = file_name.replace("_clean.h5ad", "")
                
                print(f"  -> Loading {run_id} from {file_name}...")
                adata = sc.read_h5ad(file_path)
                
                assert isinstance(adata.obs, pd.DataFrame)
                adata.obs["run"] = run_id
                
                if not adata.obs_names[0].endswith(f"-{run_id}"):
                    adata.obs_names = (adata.obs_names.astype(str) + f"-{run_id}").tolist()
                
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
            
        assert isinstance(combined_adata.obs, pd.DataFrame)
        combined_adata.obs["run"] = combined_adata.obs["run"].astype(str)
        
        if sra_meta:
            meta_df = pd.DataFrame.from_dict(sra_meta, orient="index")
            meta_df.index = meta_df.index.astype(str)
            meta_df.index.name = "run"

            combined_adata.obs = combined_adata.obs.join(meta_df, on="run")

            loaded_runs = set(combined_adata.obs["run"].astype(str))
            xml_runs = set(meta_df.index.astype(str))
            matched = loaded_runs & xml_runs
            missing = sorted(loaded_runs - xml_runs)

            print(f"Matched XML metadata: {len(matched)} / {len(loaded_runs)} runs")
            if missing:
                print("Warning: metadata not found for these runs:")
                for x in missing[:20]:
                    print(f"  - {x}")
                if len(missing) > 20:
                    print(f"  ... and {len(missing) - 20} more")

        combined_adata.obs["run"] = combined_adata.obs["run"].astype("category")
        combined_adata.obs["sample"] = combined_adata.obs["alias"].astype("category")
        combined_adata.obs = combined_adata.obs.drop(columns=["alias"])

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