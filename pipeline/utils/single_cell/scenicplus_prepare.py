import argparse
import pickle
import shutil
import anndata as ad
import pandas as pd
from pathlib import Path
from scenicplus.data_wrangling.adata_cistopic_wrangling import (  # type: ignore
    process_multiome_data,
)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--gex-in", type=Path, required=True)
    p.add_argument("--gex-out", type=Path, required=True)
    p.add_argument("--mudata-out", type=Path, required=True)
    p.add_argument("--cistopic", type=Path, required=True)
    p.add_argument("--celltype", type=str, required=True)
    p.add_argument("--regions-in-dir", type=Path, required=True)
    p.add_argument("--regions-out-dir", type=Path, required=True)
    p.add_argument("--cell-column", default="cell")
    p.add_argument("--cell-id-separator", default="___")
    return p.parse_args()

def read_bed(path: Path) -> set[str]:
    table = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2])
    return set(
        table[0].astype(str) + ":"
        + table[1].astype(str) + "-"
        + table[2].astype(str)
    )

def write_bed(regions: set[str], path: Path) -> None:
    with path.open("w") as f:
        for region in regions:
            chrom, coords = region.split(":")
            start, end = coords.split("-")
            f.write(f"{chrom}\t{start}\t{end}\n")

def prepare_regions(input_dir: Path, output_dir: Path) -> set[str]:
    top3k_dir = input_dir / "topics_top_3k"
    otsu_dir = input_dir / "topics_otsu"
    output_dir = output_dir / "topics_top3k_otsu_intersection"

    top3k_beds = sorted(top3k_dir.glob("*.bed"))
    if not top3k_beds:
        raise FileNotFoundError(f"No BED files found in {top3k_dir}")

    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True)

    all_regions: set[str] = set()
    for top3k_bed in top3k_beds:
        otsu_bed = otsu_dir / top3k_bed.name
        if not otsu_bed.is_file():
            raise FileNotFoundError(otsu_bed)

        overlap = read_bed(top3k_bed) & read_bed(otsu_bed)
        if overlap:
            write_bed(overlap, output_dir / top3k_bed.name)
            all_regions.update(overlap)

    if not all_regions:
        raise ValueError("No overlapping top-3k/Otsu regions found")
    return all_regions

def main() -> bool:
    args = parse_args()
    for path in (args.gex_in, args.cistopic, args.regions_in_dir):
        if not path.exists():
            raise FileNotFoundError(path)

    gex = ad.read_h5ad(args.gex_in, backed="r")
    gex.obs_names = gex.obs_names.astype(str)  # type: ignore
    if args.cell_column not in gex.obs:
        raise ValueError(
            f"{args.cell_column!r} is missing from GEX.obs"
            f"Available columns: {list(gex.obs.columns)}"
        )

    selected = gex.obs_names[
        gex.obs[args.cell_column].astype(str).eq(args.celltype)
    ]
    if len(selected) == 0:
        print(f"No cells found for cell type '{args.celltype}' in GEX.obs[{args.cell_column}]")
        return False

    with args.cistopic.open("rb") as handle:
        cistopic = pickle.load(handle)

    cistopic_cells = set(map(str, cistopic.cell_names))
    selected = selected[selected.isin(cistopic_cells)]

    region_set = prepare_regions(args.regions_in_dir, args.regions_out_dir)
    selected_regions = [r for r in cistopic.region_names if r in region_set]
    if len(selected_regions) != len(region_set):
        raise ValueError("Region sets and cisTopic peaks do not match")
    
    assert isinstance(gex.obs, pd.DataFrame)
    cistopic.add_cell_data(
        gex.obs.loc[selected].copy(),
        split_pattern=args.cell_id_separator
    )
    
    focused_gex = gex[selected].to_memory()
    gex.file.close()
    focused_gex.write_h5ad(args.gex_out, compression="gzip")

    mdata = process_multiome_data(
        GEX_anndata=focused_gex,
        cisTopic_obj=cistopic,
        use_raw_for_GEX_anndata=True,
        imputed_acc_kwargs={"selected_regions": selected_regions},
        bc_transform_func=lambda x: x,
    )
    mdata.write_h5mu(args.mudata_out)

    print(
        f"Prepared {args.celltype}: {len(selected)} cells, "
        f"{len(selected_regions)} regions"
    )
    return True

if __name__ == "__main__":
    valid_celltype = main()
    print(valid_celltype)
