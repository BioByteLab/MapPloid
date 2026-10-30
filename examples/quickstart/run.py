#!/usr/bin/env python3
"""Run a small, reproducible 3x example without downloading data."""
from pathlib import Path
import argparse
import csv
import os
import sys
import time
import shlex

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")

import numpy as np
import pysam
from mapploid.cli import main


def make_inputs(root):
    """Deterministic random reference and coordinate-sorted synthetic 3x BAM."""
    root.mkdir(parents=True, exist_ok=True)
    ref, bam = root / "reference.fa", root / "synthetic_3x.bam"
    reference_rng = np.random.default_rng(20261030)
    sequences = [reference_rng.choice(list("ACGT"), 12000) for _ in range(3)]
    alternate = {"A": "C", "C": "G", "G": "T", "T": "A"}
    # Keep identities stable so the public pipeline can reuse completed results.
    if not ref.exists():
        ref.write_text("".join(f">chr{i+1}\n" + "".join(seq) + "\n"
                               for i, seq in enumerate(sequences)))
    if not bam.exists():
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": f"chr{i}", "LN": 12000} for i in range(1, 4)]}
        rng = np.random.default_rng(20261006)
        with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
            for chrom in range(3):
                for start in range(0, 12000, 1000):
                    bases = sequences[chrom][start:start + 1000]
                    alternatives = np.array([alternate[base] for base in bases])
                    for j in range(30):
                        read = pysam.AlignedSegment()
                        read.query_name = f"{chrom}_{start}_{j}"
                        read.query_sequence = "".join(np.where(rng.random(1000) < 1/3, alternatives, bases))
                        read.flag, read.reference_id, read.reference_start = 0, chrom, start
                        read.mapping_quality, read.cigar = 60, [(0, 1000)]
                        read.query_qualities = pysam.qualitystring_to_array("I" * 1000)
                        handle.write(read)
    if not Path(str(bam) + ".bai").exists():
        pysam.index(str(bam))
    return ref, bam


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "output")
    args = parser.parse_args()
    started = time.monotonic()
    root = args.out.resolve()
    print("[1/3] Create the small synthetic reference and indexed BAM (no downloads).", flush=True)
    ref, bam = make_inputs(root / "inputs")
    assets, result = root / "assets", root / "result"
    prepare_args = ["prepare", "--reference", str(ref), "--kmer", "50", "--threads", "1",
                    "--out", str(assets)]
    print("[2/3] Prepare reusable reference regions:", flush=True)
    print(shlex.join([sys.executable, str(ROOT / "mapploid"), *prepare_args]), flush=True)
    status = main(prepare_args)
    if status:
        return status
    call_args = ["call", "--bam", str(bam), "--assets", str(assets), "--out", str(result),
                 "--export-figures"]
    print("[3/3] Call ploidy and export PDF/PNG:", flush=True)
    print(shlex.join([sys.executable, str(ROOT / "mapploid"), *call_args]), flush=True)
    status = main(call_args)
    if status:
        return status
    with (result / "result.tsv").open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    if row["status"] != "CALLED" or row["predicted_ploidy"] != "3x":
        raise RuntimeError(f"Expected CALLED 3x, got {row}")
    print(f"Quickstart passed: CALLED 3x; {row['selected_sites']} fitting SNPs; "
          f"{time.monotonic() - started:.1f} seconds")
    print(f"Report: {result / 'report.html'}")
    print(f"Table:  {result / 'result.tsv'}")
    tutorial = Path(__file__).resolve().parent
    print(f"Standard workflow tutorial (English): {tutorial / 'README.md'}")
    print(f"Standard workflow tutorial (Simplified Chinese): {tutorial / 'README.zh-CN.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
