import anndata
import argparse
import glob
import os
import pandas as pd
import scanpy as sc
import xml.etree.ElementTree as ET
from scipy.sparse import csr_matrix

class Args(argparse.Namespace):
    input_dir: str
    project_id: str
    output_dir: str
    library_manifest: str
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
        "--library_manifest", default="", help="Path to the library manifest file"
    )
    parser.add_argument(
        "--sra_xml", default="", help="Path to the SRA metadata XML file"
    )
    return parser.parse_args(namespace=Args())

def add_value(d: dict[str, str], key: str, value: str) -> None:
    value = (value or "").strip()
    if not key or not value:
        return

    if key in d and d[key] != value:
        raise ValueError(
            f"Metadata conflict for key='{key}': existing='{d[key]}', new='{value}'"
        )

    d[key] = value

def parse_metadata_xml(xml_path: str) -> dict[str, dict[str, str]]:
    if not xml_path:
        return {}
    
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
            run_to_meta[run_acc] = dict(base_meta)

    return run_to_meta

def read_library_manifest(library_manifest_path: str):
    if not library_manifest_path:
        return {}

    df = pd.read_csv(library_manifest_path, sep="\t", dtype=str).fillna("")
    return {
        row["library_id"]: {
            "sample_id": row["sample_id"],
            "experiment": row["experiment"],
            "runs": row["runs"],
        }
        for _, row in df.iterrows()
    }

def merge_run_metadata(runs: str, run_to_meta: dict[str, dict[str, str]]) -> dict[str, str]:
    merged: dict[str, str] = {}

    for run_id in [x.strip() for x in runs.split(",") if x.strip()]:
        for key, value in run_to_meta.get(run_id, {}).items():
            add_value(merged, key, value)

    return merged

def check_sample_alias_mapping(
    sample_id: str,
    alias: str,
    library_id: str,
    sample_to_alias: dict[str, str],
    alias_to_sample: dict[str, str],
) -> None:
    if not alias:
        return

    if sample_id in sample_to_alias and sample_to_alias[sample_id] != alias:
        raise ValueError(
            f"Sample metadata conflict: same sample_id has different aliases. "
            f"library_id='{library_id}', sample_id='{sample_id}', "
            f"previous_alias='{sample_to_alias[sample_id]}', new_alias='{alias}'"
        )

    if alias in alias_to_sample and alias_to_sample[alias] != sample_id:
        raise ValueError(
            f"Sample metadata conflict: same alias maps to different sample_ids. "
            f"library_id='{library_id}', alias='{alias}', "
            f"previous_sample_id='{alias_to_sample[alias]}', new_sample_id='{sample_id}'"
        )

    sample_to_alias[sample_id] = alias
    alias_to_sample[alias] = sample_id

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"=== Merging samples for Project: {args.project_id} ===")

    library_manifest = read_library_manifest(args.library_manifest)
    sra_xml = parse_metadata_xml(args.sra_xml)

    search_pattern = os.path.join(args.input_dir, "**", "*_clean.h5ad")
    files = sorted(glob.glob(search_pattern, recursive=True))

    if not files:
        raise FileNotFoundError(f"No '*_clean.h5ad' files found in {args.input_dir}")
    print(f"Found {len(files)} files to merge.")

    adatas = []
    sample_to_alias: dict[str, str] = {}
    alias_to_sample: dict[str, str] = {}

    for file_path in files:
        file_name = os.path.basename(file_path)
        library_id = file_name.replace("_clean.h5ad", "")

        print(f"  -> Loading {library_id} from {file_name}...")
        adata = sc.read_h5ad(file_path)

        adata.obs["project_id"] = args.project_id
        adata.obs["library_id"] = library_id

        if library_id in library_manifest:
            sample_id = library_manifest[library_id]["sample_id"]
            experiment = library_manifest[library_id]["experiment"]
            runs = library_manifest[library_id]["runs"]

            adata.obs["sample"] = sample_id
            adata.obs["experiment"] = experiment

            meta = merge_run_metadata(runs, sra_xml)
            alias = meta.get("alias", "")

            check_sample_alias_mapping(
                sample_id=sample_id,
                alias=alias,
                library_id=library_id,
                sample_to_alias=sample_to_alias,
                alias_to_sample=alias_to_sample,
            )

            for key, value in meta.items():
                adata.obs[key] = value

        adata.obs_names = adata.obs_names.astype(str) + f"-{library_id}"
        adatas.append(adata)

    print("Concatenating h5ad objects...")
    combined_adata = anndata.concat(adatas, join="outer", merge="same")

    if not isinstance(combined_adata.X, csr_matrix):
        combined_adata.X = csr_matrix(combined_adata.X)

    combined_adata.obs["project_id"] = combined_adata.obs["project_id"].astype("category")
    combined_adata.obs["library_id"] = combined_adata.obs["library_id"].astype("category")

    if "sample" in combined_adata.obs.columns:
        combined_adata.obs["sample"] = combined_adata.obs["sample"].astype("category")

    if "experiment" in combined_adata.obs.columns:
        combined_adata.obs["experiment"] = combined_adata.obs["experiment"].astype("category")

    output_path = os.path.join(args.output_dir, f"{args.project_id}.h5ad")
    print(f"Saving merged file to {output_path}...")
    combined_adata.write_h5ad(output_path, compression="gzip")

    print(f"=== Merge completed successfully ===")
    print(f"Total cells: {combined_adata.n_obs}")
    print(f"Total genes: {combined_adata.n_vars}")
    print(f"Libraries: {combined_adata.obs['library_id'].nunique()}")

if __name__ == "__main__":
    main()
