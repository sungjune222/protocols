import argparse
import os
import shutil
import zipfile


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--zip_file", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    n = 0
    with zipfile.ZipFile(args.zip_file) as z:
        for info in z.infolist():
            base = os.path.basename(info.filename)

            if not (base.endswith(".fastq.gz") or base.endswith(".fq.gz") or
                    base.endswith(".fastq") or base.endswith(".fq")):
                continue

            dest = os.path.join(args.outdir, base)

            if os.path.exists(dest):
                raise RuntimeError(f"Duplicate FASTQ basename in zip: {base}")

            with z.open(info) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)

            n += 1

    if n == 0:
        raise RuntimeError(f"No FASTQ files found in zip: {args.zip_file}")

    print(f"Extracted {n} FASTQ files")


if __name__ == "__main__":
    main()