"""The two-command public interface."""
import argparse
import os
from pathlib import Path
import sys

from . import __version__


def build_parser():
    parser = argparse.ArgumentParser(prog="mapploid", description="Plant ploidy inference in high-mappability regions.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="prepare reusable reference regions once")
    prepare.add_argument("--reference", type=Path, required=True)
    prepare.add_argument("--read-length", type=int, default=150)
    prepare.add_argument("--out", type=Path, required=True)
    prepare.add_argument("--threads", type=int, default=8)
    prepare.add_argument("--exclude-bed", type=Path, action="append", default=[],
                         help="optional regions to exclude, such as organellar sequences")
    prepare.add_argument("--genmap-bedgraph", type=Path, help="reuse compatible GenMap scores")
    call = sub.add_parser("call", help="infer one sample's ploidy automatically")
    call.add_argument("--bam", type=Path, required=True)
    call.add_argument("--assets", type=Path, required=True)
    call.add_argument("--out", type=Path, required=True)
    call.add_argument("--sample", help="sample name; defaults to the BAM filename")
    filters = call.add_argument_group("site filters (optional)")
    filters.add_argument("--allele-min", type=float, default=0.15,
                         help="minimum frequency for each of the two leading alleles (default: 0.15)")
    filters.add_argument("--allele-max", type=float, default=0.85,
                         help="maximum frequency for each of the two leading alleles (default: 0.85)")
    filters.add_argument("--min-depth", type=int,
                         help="minimum quality-filtered A/C/G/T depth (default: automatic)")
    filters.add_argument("--max-depth", type=int,
                         help="maximum quality-filtered A/C/G/T depth (default: automatic)")
    filters.add_argument("--min-mapq", type=int, default=20,
                         help="minimum read mapping quality, 0-255 (default: 20)")
    filters.add_argument("--min-baseq", type=int, default=20,
                         help="minimum base quality, 0-255 (default: 20)")
    filters.add_argument("--min-allele-reads", type=int, default=2,
                         help="minimum supporting reads for each allele (default: 2)")
    filters.add_argument("--pair-purity", type=float, default=0.95,
                         help="minimum fraction explained by the two leading alleles, 0.5-1 (default: 0.95)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(key, "1")
    from .pipeline import PrepareConfig, CallConfig, prepare, call
    try:
        if args.command == "prepare":
            path = prepare(PrepareConfig(args.reference, args.out, args.read_length,
                                         args.threads, args.exclude_bed, args.genmap_bedgraph))
        else:
            path = call(CallConfig(
                args.bam, args.assets, args.out, args.sample,
                allele_min=args.allele_min, allele_max=args.allele_max,
                min_depth=args.min_depth, max_depth=args.max_depth,
                min_mapq=args.min_mapq, min_baseq=args.min_baseq,
                min_allele_reads=args.min_allele_reads, pair_purity=args.pair_purity,
            ))
        print(path)
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[MapPloid] ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("[MapPloid] interrupted; rerun the same command to resume", file=sys.stderr)
        return 130
