import argparse
import csv
import os
from collections import defaultdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runinfo_csv", required=True)
    p.add_argument("--sra_dir", required=True)
    p.add_argument("--library_manifest", required=True)
    p.add_argument("--run_manifest", required=True)
    p.add_argument("--download_input_list", required=True)
    args = p.parse_args()

    groups = defaultdict(list)

    with open(args.runinfo_csv, newline="") as f:
        for row in csv.DictReader(f):
            run = row["Run"].strip()
            exp = row["Experiment"].strip()
            sample = row["BioSample"].strip()
            url = row["download_path"].strip()
            library = f"{sample}__{exp}"

            groups[library].append({
                "run": run,
                "experiment": exp,
                "sample": sample,
                "url": url,
            })

    with open(args.library_manifest, "w", newline="") as f_lib, \
        open(args.run_manifest, "w", newline="") as f_run, \
        open(args.download_input_list, "w") as f_down:

        lib_writer = csv.writer(f_lib, delimiter="\t")
        run_writer = csv.writer(f_run, delimiter="\t")

        lib_writer.writerow(["library_id", "sample_id", "experiment", "runs"])
        run_writer.writerow(["run", "library_id", "sample_id", "experiment", "lane_id", "sra_file"])

        for library_id in sorted(groups):
            runs = sorted(groups[library_id], key=lambda x: x["run"])

            sample_id = runs[0]["sample"]
            exp = runs[0]["experiment"]

            lib_writer.writerow([
                library_id,
                sample_id,
                exp,
                ",".join(r["run"] for r in runs),
            ])

            for lane_idx, r in enumerate(runs, start=1):
                lane_id = f"L{lane_idx:03d}"
                sra_file = os.path.join(args.sra_dir, f"{r['run']}.sra")

                run_writer.writerow([
                    r["run"],
                    library_id,
                    sample_id,
                    exp,
                    lane_id,
                    sra_file,
                ])

                if r["url"] and not os.path.exists(sra_file):
                    f_down.write(r["url"] + "\n")
                    f_down.write(f"  out={r['run']}.sra\n")
                    f_down.write(f"  dir={args.sra_dir}\n")

if __name__ == "__main__":
    main()