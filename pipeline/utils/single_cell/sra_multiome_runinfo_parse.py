import argparse
import csv
import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict

DELIM = re.compile(r"[,;\s]+")

def read_titles(path):
    titles = {}

    for item in ET.parse(path).iter("Item"):
        if item.get("Name") == "ExpXml":
            exp = ET.fromstring(f"<root>{item.text}</root>").find("Experiment")
            titles[exp.get("acc")] = exp.get("name", "") #type: ignore
    return titles

def sample_id(xs):
    return re.sub(r"[^A-Za-z0-9_]+", "_", "_".join(xs)).strip("_") or "sample"

def strategy(row):
    x = row.get("LibraryStrategy", "").strip().lower()
    return "gex" if x == "rna-seq" else "atac" if x == "atac-seq" else ""

def tokens(row, titles):
    x = titles.get(row.get("Experiment", ""), "") or (
         row.get("LibraryName", "").strip()
         or row.get("SampleName", "").strip()
         or row.get("Sample", "").strip()
     )
    if not x:
        return []

    x = re.sub(r"^\s*(?:GSM|SRX)\w+\s*:\s*", "", x, flags=re.I)
    x = x.split(";", 1)[0]
    x = re.sub(r",?\s*sc(?:RNA|ATAC)seq\s*$", "", x, flags=re.I) 
    return [t for t in DELIM.split(x) if t]

def read_experiments(path):
    experiments = defaultdict(list)

    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            run = row.get("Run", "").strip()
            exp = row.get("Experiment", "").strip()

            if run and exp:
                experiments[exp].append(row)

    return experiments

def build_pairs(experiments, titles):
    items = []

    for exp, rows in experiments.items():
        rows = sorted(rows, key=lambda r: r["Run"])
        row0 = rows[0]

        assay = strategy(row0)
        token = tokens(row0, titles)

        if assay not in {"gex", "atac"}:
            raise ValueError(f"Cannot detect assay from LibraryStrategy for experiment: {exp}")
        if len(token) < 2:
            raise ValueError(f"Unexpected sample naming structure for experiment: {exp}")

        items.append({
            "experiment": exp,
            "assay": assay,
            "tokens": token,
            "runs": rows,
        })

    prefix_map = defaultdict(lambda: {"gex": [], "atac": []})

    for i, item in enumerate(items):
        for n in range(1, len(item["tokens"]) + 1):
            prefix_map[tuple(item["tokens"][:n])][item["assay"]].append(i)

    candidates = [
        prefix
        for prefix, d in prefix_map.items()
        if len(d["gex"]) == 1 and len(d["atac"]) == 1
    ]

    assigned = set()
    pairs = []

    for prefix in sorted(candidates, key=lambda x: (-len(x), x)):
        gex = prefix_map[prefix]["gex"][0]
        atac = prefix_map[prefix]["atac"][0]

        if gex in assigned or atac in assigned:
            continue

        assigned.update([gex, atac])
        pairs.append({
            "sample_id": sample_id(prefix),
            "key": prefix,
            "gex": items[gex],
            "atac": items[atac],
        })

    if len(assigned) != len(items):
        bad = [x for i, x in enumerate(items) if i not in assigned]
        msg = "\n".join(f"  {x['experiment']}\t{x['assay']}" for x in bad)
        raise ValueError(
            "Unexpected sample naming structure. "
            "Could not make exact 1:1 GEX/ATAC pairs by delimiter-prefix matching:\n"
            + msg
        )

    return sorted(pairs, key=lambda p: p["sample_id"])

def write_outputs(pairs, args):
    seen = set()
    used = defaultdict(int)

    with open(args.library_manifest, "w", newline="") as f_lib, \
         open(args.run_manifest, "w", newline="") as f_run, \
         open(args.download_input_list, "w") as f_down:

        lib_writer = csv.writer(f_lib, delimiter="\t", lineterminator="\n")
        run_writer = csv.writer(f_run, delimiter="\t", lineterminator="\n")

        lib_writer.writerow([
            "sample_id",
            "gex_experiment",
            "gex_runs",
            "atac_experiment",
            "atac_runs",
        ])

        run_writer.writerow([
            "run",
            "sample_id",
            "assay",
            "experiment",
            "lane_id",
            "sra_file",
        ])

        for pair in pairs:
            base = pair["sample_id"]
            used[base] += 1
            sid = base if used[base] == 1 else f"{base}_{used[base]}"

            gex = pair["gex"]
            atac = pair["atac"]

            lib_writer.writerow([
                sid,
                gex["experiment"],
                ",".join(r["Run"] for r in gex["runs"]),
                atac["experiment"],
                ",".join(r["Run"] for r in atac["runs"]),
            ])

            for assay, item in [("gex", gex), ("atac", atac)]:
                for i, row in enumerate(item["runs"], start=1):
                    run = row["Run"]
                    sra_file = os.path.join(args.sra_dir, f"{run}.sra")

                    run_writer.writerow([
                        run,
                        sid,
                        assay,
                        item["experiment"],
                        f"L{i:03d}",
                        sra_file,
                    ])

                    url = row.get("download_path", "").strip()
                    if run not in seen and url and not os.path.exists(sra_file):
                        f_down.write(f"{url}\n")
                        f_down.write(f"  out={run}.sra\n")
                        f_down.write(f"  dir={args.sra_dir}\n")

                    seen.add(run)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runinfo_csv", required=True)
    p.add_argument("--sra_xml", required=True)
    p.add_argument("--sra_dir", required=True)
    p.add_argument("--library_manifest", required=True)
    p.add_argument("--run_manifest", required=True)
    p.add_argument("--download_input_list", required=True)
    args = p.parse_args()

    pairs = build_pairs(
        read_experiments(args.runinfo_csv),
        read_titles(args.sra_xml),
    )
    write_outputs(pairs, args)

    print(f"Detected {len(pairs)} paired multiome samples")


if __name__ == "__main__":
    main()