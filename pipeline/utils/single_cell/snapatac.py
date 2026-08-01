import argparse
import anndata as ad
import numpy as np
import pandas as pd
import snapatac2 as snap
from pathlib import Path
from pipeline.config.constants import CPU_CORE_COUNT

def read_chrom_sizes(path: Path) -> dict[str, int]:
    table = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1])
    if table[0].duplicated().any():
        raise ValueError(f"Duplicated chromosomes in {path}")
    return dict(zip(table[0].astype(str), table[1].astype(int)))
        
def read_singlets(path: Path) -> list[str]:
    table = pd.read_csv(path, dtype=str)
    if "barcode" not in table:
        raise ValueError(f"`barcode` column is missing from {path}")
    barcodes = table["barcode"].str.strip()
    if barcodes.isna().any() or (barcodes == "").any() or barcodes.duplicated().any():
        raise ValueError(f"Invalid singlet barcodes in {path}")
    return barcodes.tolist()

def canonical_ids(barcodes: list[str], library_id: str, separator: str) -> list[str]:
    if separator in library_id or any(separator in barcode for barcode in barcodes):
        raise ValueError(
            f"Cell ID separator {separator!r} occurs in {library_id} or its barcodes"
        )
    return [f"{barcode}{separator}{library_id}" for barcode in barcodes]

def nonzero_axes(matrix) -> tuple[np.ndarray, np.ndarray]:
    cells = np.zeros(matrix.n_obs, dtype=bool)
    peaks = np.zeros(matrix.n_vars, dtype=bool)
    for chunk, start, end in matrix.chunked_X(1000):
        chunk = chunk.tocsr()
        cells[start:end] = np.diff(chunk.indptr) > 0
        peaks[np.unique(chunk.indices)] = True
    return cells, peaks

def write_bed(regions: list[str], path: Path) -> None:
    parsed: list[tuple[str, int, int]] = []
    for region in regions:
        chrom, coordinates = region.rsplit(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        start, end = int(start_text), int(end_text)
        parsed.append((chrom, start, end))

    with path.open("w") as handle:
        for chrom, start, end in parsed:
            handle.write(f"{chrom}\t{start}\t{end}\n")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Call de novo multiome peaks and export paired SCENIC+ inputs"
    )
    parser.add_argument("--cellranger-arc-dir", type=Path, required=True)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--library-list", type=Path, required=True)
    parser.add_argument("--scdblfinder-dir", type=Path, required=True)
    parser.add_argument("--rna-h5ad", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--gene-annotation", type=Path, required=True)
    parser.add_argument("--cell-id-separator", default="___")
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    required = [
        args.library_list,
        args.rna_h5ad,
        args.chrom_sizes,
        args.gene_annotation,
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    if not args.cell_id_separator:
        raise ValueError("Invalid cell ID separator")

    libraries = [line.strip() for line in args.library_list.read_text().splitlines() if line.strip()]
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    sample_dir = args.outdir / "samples"
    peak_call_dir = args.outdir / "peak_calls"
    scenic_dir = args.outdir / "scenic_input"
    temp_dir = args.outdir / "tmp"
    for path in (sample_dir, peak_call_dir, scenic_dir, temp_dir):
        path.mkdir(parents=True, exist_ok=True)

    sample_files: list[tuple[str, Path]] = []
    fragment_manifest: list[dict[str, str]] = []
    opened: list[tuple[str, object]] = []
    data = None
    peak_matrix = None

    for index, library_id in enumerate(libraries, start=1):
        print(f"[{index}/{len(libraries)}] Importing ATAC: {library_id}")
        fragment_file = args.cellranger_arc_dir.joinpath(library_id, "outs", "atac_fragments.tsv.gz")
        singlet_file = args.scdblfinder_dir.joinpath(library_id, f"{library_id}_singlet_barcodes.csv")
    
        if not fragment_file.is_file():
            raise FileNotFoundError(fragment_file)
        if not singlet_file.is_file():
            raise FileNotFoundError(singlet_file)

        barcodes = read_singlets(singlet_file)
        sample_file = sample_dir / f"{library_id}.h5ad"
        sample = snap.pp.import_fragments(
            fragment_file=fragment_file,
            chrom_sizes=chrom_sizes,
            file=sample_file,
            min_num_fragments=1,
            sorted_by_barcode=False,
            whitelist=barcodes,
            tempdir=temp_dir,
            n_jobs=CPU_CORE_COUNT,
        )
        
        snap.metrics.tsse(
            adata=sample,
            gene_anno=args.gene_annotation,
            n_jobs=CPU_CORE_COUNT,
        )
        before = sample.n_obs
        snap.pp.filter_cells(data=sample, n_jobs=CPU_CORE_COUNT)

        if sample.n_obs == 0:
            raise ValueError(f"No ATAC cells passed QC for {library_id}")

        original = np.asarray(sample.obs_names, dtype=str).tolist()
        sample.obs["orig_barcode"] = original
        sample.obs["library_id"] = [library_id] * sample.n_obs
        sample.obs_names = canonical_ids(original, library_id, args.cell_id_separator)

        snap.pp.add_tile_matrix(sample, n_jobs=CPU_CORE_COUNT)
        print(f"  ATAC QC: {before} -> {sample.n_obs} cells")

        sample_files.append((library_id, sample_file))
        fragment_manifest.append(
            {
                "library_id": library_id,
                "fragment_file": str(fragment_file.resolve()),
            }
        )

    opened = [(library_id, snap.read(path)) for library_id, path in sample_files]
    data = snap.AnnDataSet(
        adatas=opened,
        filename=args.outdir / "tile_matrix.h5ads",
        add_key="sample_id",
        use_absolute_path=True,
    )
    if not pd.Index(np.asarray(data.obs_names, dtype=str)).is_unique:
        raise ValueError("ATAC cell IDs are not unique")

    snap.pp.select_features(adata=data, n_jobs=CPU_CORE_COUNT)
    snap.tl.spectral(adata=data, num_threads=CPU_CORE_COUNT, random_state=0)

    mat = snap.pp.harmony(data, batch="sample_id", inplace=False, random_state=0)
    if mat.shape[0] != data.n_obs: #type: ignore
        mat = mat.T #type: ignore
    data.obsm["X_harmony"] = mat

    snap.pp.knn(adata=data, use_rep="X_harmony", random_state=0)
    snap.tl.leiden(adata=data, objective_function="CPM", random_state=0)

    peak_calls = snap.tl.macs3(
        adata=data,
        groupby="leiden",
        qvalue=0.05,
        tempdir=temp_dir,
        inplace=False,
        n_jobs=CPU_CORE_COUNT,
    )
    if not peak_calls:
        raise ValueError("MACS3 did not return any peaks")
    for group, table in peak_calls.items():
        table.write_parquet(peak_call_dir / f"{group}.parquet")

    consensus = snap.tl.merge_peaks(peaks=peak_calls, chrom_sizes=chrom_sizes)
    regions = consensus["Peaks"].to_list()
    if not regions:
        raise ValueError("The consensus peak set is empty")
    consensus.write_parquet(args.outdir / "consensus_peaks.parquet")

    internal_acc_path = args.outdir / "consensus_peak_matrix.snap.h5ad"
    peak_matrix = snap.pp.make_peak_matrix(
        adata=data,
        use_rep=regions,
        file=internal_acc_path,
        counting_strategy="fragment",
    )
    keep_cells, keep_peaks = nonzero_axes(peak_matrix)
    if not keep_cells.any() or not keep_peaks.any():
        raise ValueError("The consensus peak matrix is empty")
    if not keep_cells.all() or not keep_peaks.all():
        print(
            "  Removing "
            f"{int((~keep_cells).sum())} zero-peak cells and "
            f"{int((~keep_peaks).sum())} zero-count peaks"
        )
        peak_matrix.subset(
            np.flatnonzero(keep_cells),
            np.flatnonzero(keep_peaks),
        )

    atac_cells = pd.Index(np.asarray(peak_matrix.obs_names, dtype=str))
    if not atac_cells.is_unique:
        raise ValueError("Consensus ATAC matrix contains duplicated cell IDs")
    write_bed(
        np.asarray(peak_matrix.var_names, dtype=str).tolist(),
        scenic_dir / "consensus_peaks.bed",
    )

    rna = ad.read_h5ad(args.rna_h5ad)
    rna.obs_names = rna.obs_names.astype(str) #type: ignore
    if not {"orig_barcode", "library_id"}.issubset(rna.obs):
        raise ValueError("RNA is missing barcode or library metadata")

    missing_rna = atac_cells.difference(rna.obs_names)
    if len(missing_rna):
        raise ValueError(f"{len(missing_rna)} ATAC cells are missing from the RNA matrix")
    rna = rna[atac_cells].copy()
    if not rna.obs_names.equals(atac_cells):
        raise ValueError("RNA and ATAC cell order does not match")

    for column in rna.obs:
        peak_matrix.obs[column] = rna.obs[column].to_numpy()

    acc_path = scenic_dir / "ACC.h5ad"
    peak_matrix_memory = peak_matrix.to_memory()
    peak_matrix_memory.write_h5ad(acc_path, compression="gzip")
    rna.raw = rna.copy()
    gex_path = scenic_dir / "GEX.h5ad"
    rna.write_h5ad(gex_path, compression="gzip")

    final_libraries = set(rna.obs["library_id"].astype(str))
    final_fragments = [
        row for row in fragment_manifest
        if row["library_id"] in final_libraries
    ]
    pd.DataFrame(final_fragments).to_csv(
        scenic_dir / "fragments.tsv", sep="\t", index=False
    )

    print(
        f"Completed: {len(atac_cells)} paired cells, "
        f"{rna.n_vars} genes, {peak_matrix.n_vars} de novo peaks"
    )

if __name__ == "__main__":
    main()
