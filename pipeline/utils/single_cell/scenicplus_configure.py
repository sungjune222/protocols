import argparse
import yaml
import pandas as pd
import pyranges as pr
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    for name in (
        "config",
        "cistopic",
        "gex",
        "regions",
        "mudata",
        "ctx-db",
        "dem-db",
        "motif-annotation",
        "gtf",
        "chrom-sizes",
        "outdir",
    ):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--fdr", type=float, choices=(0.05, 1), required=True)
    parser.add_argument("--reference-genome", type=str, required=True)
    return parser.parse_args()

def write_genome_files(gtf: Path, fai: Path, outdir: Path) -> tuple[Path, Path]:
    annotation_out = outdir / "genome_annotation.tsv"
    chromsizes_out = outdir / "chromsizes.tsv"
    if annotation_out.is_file() and chromsizes_out.is_file():
        return annotation_out, chromsizes_out

    allowed_gene_types = {
        "protein_coding",
        "lncRNA",
        "miRNA",
        "snoRNA",
        "snRNA",
        "pseudogene"
    }
    annotation = pr.read_gtf(str(gtf), as_df=True)

    genes = annotation.loc[
        annotation["Feature"].eq("gene")
        & annotation["gene_type"].isin(allowed_gene_types)
    ].copy()
    transcripts = annotation.loc[annotation["Feature"].eq("transcript")].copy()
    transcripts["Transcription_Start_Site"] = transcripts["Start"].where(
        transcripts["Strand"].eq("+"), transcripts["End"] - 1
    )
    genes = genes.merge(transcripts[["gene_id", "Transcription_Start_Site"]], on="gene_id")

    gene_name = "gene_name" if "gene_name" in genes else "gene_id"
    genes["Gene"] = genes[gene_name]
    genes = genes[
        ["Chromosome", "Start", "End", "Strand", "Gene", "Transcription_Start_Site"]
    ].dropna().drop_duplicates(["Gene", "Transcription_Start_Site"])
    genes.to_csv(annotation_out, sep="\t", index=False)

    chromsizes = pd.read_csv(
        fai, sep="\t", header=None, usecols=[0, 1], names=["Chromosome", "End"]
    )
    chromsizes.insert(1, "Start", 0)
    chromsizes.to_csv(chromsizes_out, sep="\t", index=False)
    return annotation_out, chromsizes_out

def main() -> None:
    args = parse_args()
    for path in (
        args.config,
        args.cistopic,
        args.gex,
        args.regions,
        args.mudata,
        args.ctx_db,
        args.dem_db,
        args.motif_annotation,
        args.gtf,
        args.chrom_sizes,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    args.outdir.mkdir(parents=True, exist_ok=True)
    temp = args.outdir / "tmp"
    temp.mkdir(exist_ok=True)
    genome_annotation, chromsizes = write_genome_files(
        args.gtf, args.chrom_sizes, args.outdir
    )

    with args.config.open() as handle:
        config = yaml.safe_load(handle)

    config["input_data"].update(
        {
            "cisTopic_obj_fname": str(args.cistopic.resolve()),
            "GEX_anndata_fname": str(args.gex.resolve()),
            "region_set_folder": str(args.regions.resolve()),
            "ctx_db_fname": str(args.ctx_db.resolve()),
            "dem_db_fname": str(args.dem_db.resolve()),
            "path_to_motif_annotations": str(args.motif_annotation.resolve()),
        }
    )

    out = args.outdir
    result = out if args.fdr == 1 else out / "fdr_0.05"
    result.mkdir(exist_ok=True)
    config["output_data"].update(
        {
            "combined_GEX_ACC_mudata": str(args.mudata.resolve()),
            "dem_result_fname": str(out / "dem_results.hdf5"),
            "ctx_result_fname": str(out / "ctx_results.hdf5"),
            "output_fname_dem_html": str(out / "dem_results.html"),
            "output_fname_ctx_html": str(out / "ctx_results.html"),
            "cistromes_direct": str(out / "cistromes_direct.h5ad"),
            "cistromes_extended": str(out / "cistromes_extended.h5ad"),
            "tf_names": str(out / "tf_names.txt"),
            "genome_annotation": str(genome_annotation),
            "chromsizes": str(chromsizes),
            "search_space": str(out / "search_space.tsv"),
            "tf_to_gene_adjacencies": str(out / "tf_to_gene_adj.tsv"),
            "region_to_gene_adjacencies": str(out / "region_to_gene_adj.tsv"),
            "eRegulons_direct": str(result / "eRegulons_direct.tsv"),
            "eRegulons_extended": str(result / "eRegulons_extended.tsv"),
            "AUCell_direct": str(result / "AUCell_direct.h5mu"),
            "AUCell_extended": str(result / "AUCell_extended.h5mu"),
            "scplus_mdata": str(result / "scplus_mdata.h5mu"),
        }
    )

    config["params_general"].update(
        {"temp_dir": str(temp.resolve()), "n_cpu": args.threads, "seed": 0}
    )

    if args.reference_genome in ["GRCm39"]:
        data_prep_species = "mmusculus"
        motif_species = "mus_musculus"
    else:
        data_prep_species = "hsapiens"
        motif_species = "homo_sapiens"

    config["params_data_preparation"].update(
        {
            "bc_transform_func": '"lambda x: x"',
            "is_multiome": True,
            "species": data_prep_species,
        }
    )
    config["params_motif_enrichment"].update(
        {
            "species": motif_species,
            "annotation_version": "v10nr_clust",
            "dem_adj_pval_thr": 0.05,
        }
    )

    with args.config.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)

if __name__ == "__main__":
    main()
