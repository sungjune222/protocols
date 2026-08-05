import argparse
import pandas as pd
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct", type=Path, required=True)
    parser.add_argument("--extended", type=Path, required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    hits = []
    for source, path in (("direct", args.direct), ("extended", args.extended)):
        network = pd.read_csv(path, sep="\t")
        if not {"TF", "Gene"}.issubset(network):
            raise ValueError(f"TF/Gene columns are missing from {path}")

        upstream = network.loc[network["Gene"].eq(args.gene)].copy()
        upstream.insert(0, "direction", "TF_to_target")
        downstream = network.loc[
            network["TF"].eq(args.gene) & ~network["Gene"].eq(args.gene)
        ].copy()
        downstream.insert(0, "direction", "target_as_TF_to_gene")
        table = pd.concat([upstream, downstream], ignore_index=True)
        table.insert(0, "source", source)
        hits.append(table)

    result = pd.concat(hits, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(
        f"{args.gene}: "
        f"{result['direction'].eq('TF_to_target').sum()} upstream links, "
        f"{result['direction'].eq('target_as_TF_to_gene').sum()} target-as-TF links"
    )

if __name__ == "__main__":
    main()
