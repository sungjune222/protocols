import argparse
import csv
import re
from pathlib import Path

FASTQ_RE = re.compile(r"^(.+)_S\d+_L\d{3}_[RI]\d_001\.(?:fastq|fq)(?:\.\d+)?\.gz$")

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata_file", required=True)
    p.add_argument("--key", default="")
    p.add_argument("--library_list", required=True)
    return p.parse_args()

def cellranger_name(row):
    name = (row.get("Original File Name") or row.get("File Name") or "").strip()
    name = re.sub(r"\.(fastq|fq)\.\d+\.gz$", r".\1.gz", name)
    if not FASTQ_RE.match(name):
        raise ValueError(f"Not a 10x FASTQ name: {name}")
    return name

def main():
    args = parse_args()
    key = args.key.strip()
    libraries = set()

    with open(args.metadata_file, newline="", encoding="utf-8") as f:
        for meta in csv.DictReader(f, delimiter="\t"):
            if key and (meta.get("ARM Name") or "").strip() != key:
                continue

            match = FASTQ_RE.match(cellranger_name(meta))
            libraries.add(match.group(1)) #type: ignore

    if not libraries:
        raise SystemExit(f"No rows found for key={key or '<all>'}")

    Path(args.library_list).parent.mkdir(parents=True, exist_ok=True)
    with open(args.library_list, "w", encoding="utf-8") as f:
        for library_id in sorted(libraries):
            f.write(library_id + "\n")

    print(f"Wrote {len(libraries)} libraries for key={key or '<all>'}.")

if __name__ == "__main__":
    main()
