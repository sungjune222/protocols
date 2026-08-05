import argparse
import gc
import logging
import pickle
import anndata as ad
import numpy as np
import pandas as pd
from gensim import utils #type: ignore
from pathlib import Path
from pycisTopic.cistopic_class import create_cistopic_object #type: ignore
from pycisTopic.lda_models import evaluate_models, run_cgs_models_mallet, LDAMallet #type: ignore
from pycisTopic.topic_binarization import binarize_topics #type: ignore
from pycisTopic.utils import region_names_to_coordinates #type: ignore
from scipy import sparse

# Alternative implementation of load_word_topics that reads the MALLET state file in chunks to reduce memory usage
def load_word_topics_streaming(self, chunksize=300_000):
    logger = logging.getLogger("LDAMalletWrapper")
    logger.info("loading assigned topics from %s in chunks", self.fstate())
    word_topics = np.zeros((self.num_topics, self.num_terms), dtype=np.float64)

    with utils.open(self.fstate(), "rb") as fin:
        next(fin)
        self.alpha = np.fromiter(next(fin).split()[2:], dtype=float)

    if len(self.alpha) != self.num_topics:
        raise ValueError("Mismatch between MALLET and requested topics")

    # MALLET state: 0=doc, 1=source, 2=pos, 3=typeindex, 4=type(region), 5=topic
    for chunk in pd.read_csv(
        self.fstate(),
        sep=r"\s+",
        header=None,
        skiprows=3,
        usecols=[4, 5],
        dtype={
            4: np.int32,
            5: np.int32,
        },
        chunksize=chunksize,
        compression="infer",
    ):
        regions = chunk.iloc[:, 0].to_numpy(copy=False)
        topics = chunk.iloc[:, 1].to_numpy(copy=False)
        np.add.at(word_topics, (topics, regions), 1)

    return word_topics

LDAMallet.load_word_topics = load_word_topics_streaming

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a selected cisTopic model and SCENIC+ region sets"
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--cell-id-separator", default="___")
    parser.add_argument("--additional-topics", action="store_true")
    return parser.parse_args()

def write_region_sets(region_sets: dict, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for path in outdir.glob("*.bed"):
        path.unlink()

    for name, regions in region_sets.items():
        coordinates = region_names_to_coordinates(regions.index).sort_values(
            ["Chromosome", "Start", "End"]
        )
        coordinates.to_csv(
            outdir / f"{name}.bed",
            sep="\t",
            header=False,
            index=False,
        )

def main() -> None:
    args = parse_args()
    gex_path = args.input_dir / "GEX.h5ad"
    acc_path = args.input_dir / "ACC.h5ad"
    fragments_path = args.input_dir / "fragments.tsv"
    for path in (gex_path, acc_path, fragments_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    separator = args.cell_id_separator
    gex = ad.read_h5ad(gex_path)
    acc = ad.read_h5ad(acc_path)

    gex.obs_names = gex.obs_names.astype(str) #type: ignore
    acc.obs_names = acc.obs_names.astype(str) #type: ignore
    if not gex.obs_names.equals(acc.obs_names):
        raise ValueError("GEX and ACC cell IDs or order do not match.")

    fragments = pd.read_csv(fragments_path, sep="\t", dtype=str)
    fragment_files = dict(zip(fragments["library_id"], fragments["fragment_file"]))
    if not all(Path(path).is_file() for path in fragment_files.values()):
        raise FileNotFoundError("A fragment file in fragments.tsv is missing")
    library_ids = set(gex.obs["library_id"].astype(str))
    if library_ids != set(fragment_files):
        raise ValueError("GEX libraries do not match the fragment manifest")

    if not args.additional_topics:
        args.outdir.mkdir(parents=True)
    matrix = sparse.csr_matrix(acc.X.T) #type: ignore
    cistopic = create_cistopic_object(
        fragment_matrix=matrix,
        cell_names=acc.obs_names.tolist(),
        region_names=acc.var_names.astype(str).tolist(), 
        path_to_fragments=fragment_files,
        project=args.project_id,
        tag_cells=False,
        split_pattern=separator,
    )
    cell_data = gex.obs.copy()
    cell_data["sample_id"] = cell_data["library_id"].astype(str) #type: ignore
    cistopic.add_cell_data(cell_data, split_pattern=separator)

    n_cells = int(gex.n_obs)
    n_genes = int(gex.n_vars)
    n_regions = int(acc.n_vars)

    del gex
    del matrix
    del acc
    del cell_data
    del fragments
    gc.collect()

    base_topics = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    extra_topics = [55, 60, 65, 70]
    topics = (
        base_topics + extra_topics
        if args.additional_topics
        else base_topics
    )

    model_dir = args.outdir / "models"
    if not args.additional_topics:
        model_dir.mkdir()

    models = []
    missing_topics = []
    for n_topic in topics:
        model_path = model_dir / f"Topic{n_topic}.pkl"
        if model_path.is_file():
            with model_path.open("rb") as handle:
                models.append(pickle.load(handle))
            print(f"Reusing {model_path.name}")
        else:
            missing_topics.append(n_topic)

    if missing_topics:
        print(f"Running missing topics: {missing_topics}")

        models.extend(
            run_cgs_models_mallet(
                cistopic,
                n_topics=missing_topics,
                n_cpu=args.threads,
                n_iter=150,
                random_state=0,
                save_path=str(model_dir),
                mallet_path="/opt/mallet/bin/mallet",
            )
        )
    models.sort(key=lambda model: model.n_topic)
        
    model = evaluate_models(
        models,
        select_model=None,
        return_model=True,
        plot=False,
        save=str(args.outdir / "model_selection.pdf")
    )
    cistopic.add_LDA_model(model)

    cistopic_path = args.outdir / "cistopic_obj.pkl"
    with cistopic_path.open("wb") as handle:
        pickle.dump(cistopic, handle)

    otsu = binarize_topics(cistopic, method="otsu", plot=False)
    top_n = min(3000, cistopic.fragment_matrix.shape[0])
    top = binarize_topics(cistopic, method="ntop", ntop=top_n, plot=False)
    region_root = args.outdir / "region_sets"
    write_region_sets(otsu, region_root / "topics_otsu")
    write_region_sets(top, region_root / "topics_top_3k")

    ready = {
        "cisTopic_obj_fname": "cistopic_obj.pkl",
        "GEX_anndata_fname": str(gex_path.resolve()),
        "region_set_folder": "region_sets",
        "selected_topics": int(model.n_topic),
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_regions": n_regions,
    }
    print(f"SCENIC+ inputs ready with {model.n_topic} selected topics")
    print(ready)

if __name__ == "__main__":
    main()
