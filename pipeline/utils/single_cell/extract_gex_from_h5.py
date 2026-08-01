import argparse
import h5py
import numpy as np
from pathlib import Path
from scipy import sparse

def extract_gex(raw_path: Path, out_path: Path) -> None:
    if out_path.exists():
        print(f"Reusing GEX-only H5: {out_path}")
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(raw_path, "r") as src:
        matrix_group = src["matrix"]
        matrix = sparse.csc_matrix(
            (
                matrix_group["data"][:], #type: ignore
                matrix_group["indices"][:], #type: ignore
                matrix_group["indptr"][:], #type: ignore
            ),
            shape=matrix_group["shape"][:], #type: ignore
        )
        feature_types = np.asarray(
            [x.decode() if isinstance(x, bytes) else str(x) for x in matrix_group["features/feature_type"][:]] #type: ignore
        )
        keep = np.flatnonzero(feature_types == "Gene Expression")
        if len(keep) == 0:
            raise ValueError(f"No Gene Expression features found in {raw_path}.")

        gex = matrix[keep, :].tocsc()
        with h5py.File(out_path, "x") as dst:
            dst.attrs.update(src.attrs)
            for name in src:
                if name != "matrix":
                    src.copy(name, dst)

            dst_matrix = dst.create_group("matrix")
            dst_matrix.attrs.update(matrix_group.attrs)
            compressed = {"compression": "gzip", "shuffle": True}
            for name, value in {
                "data": gex.data,
                "indices": gex.indices,
                "indptr": gex.indptr,
                "shape": gex.shape,
                "barcodes": matrix_group["barcodes"][:], #type: ignore
            }.items():
                kwargs = {} if name == "shape" else compressed
                dst_matrix.create_dataset(name, data=value, **kwargs)

            dst_features = dst_matrix.create_group("features")
            dst_features.attrs.update(matrix_group["features"].attrs) #type: ignore
            for name, dataset in matrix_group["features"].items(): #type: ignore
                value = dataset[:]
                if dataset.ndim and dataset.shape[0] == matrix.shape[0]: #type: ignore
                    value = value[keep]
                dst_features.create_dataset(name, data=value, **compressed)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    extract_gex(args.input, args.output)

if __name__ == "__main__":
    main()