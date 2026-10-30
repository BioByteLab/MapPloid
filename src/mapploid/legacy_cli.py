"""Command-line interface for MapPloid data preparation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from . import __version__
from .bam_counts import BamCountConfig, build_bam_counts
from .ploidy_model import (
    DEFAULT_MAX_SITES,
    DEFAULT_MIN_SITES,
    PloidyModelConfig,
    run_ploidy_model,
    simulate_and_validate,
)
from .panel import PanelConfig, build_panels
from .reference_mask import ReferenceMaskConfig, build_reference_mask


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mapploid",
        description="Mapping-aware data preparation for sample ploidy inference.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    reference = subparsers.add_parser(
        "reference-mask",
        help="combine GenMap uniqueness, FASTA gaps and optional blacklists",
    )
    reference.add_argument("--reference", type=Path, required=True)
    reference.add_argument("--out-dir", type=Path, required=True)
    reference.add_argument("--kmer", dest="kmer_length", type=int, default=150)
    reference.add_argument("--min-region-length", type=int, default=1,
                           help="minimum continuous post-exclusion region length in bp (default: 1)")
    reference.add_argument("--mismatches", type=int, default=2)
    reference.add_argument("--threads", type=int, default=8)
    reference.add_argument("--genmap-executable", default="genmap")
    reference.add_argument(
        "--base-rule",
        choices=("strict", "start"),
        default="strict",
        help="strict requires every read start covering a retained base to be unique",
    )
    reference.add_argument(
        "--gap-padding",
        type=int,
        help="padding around ambiguous FASTA runs; default: k-mer length - 1",
    )
    reference.add_argument("--repeat-bed", action="append", type=Path, default=[])
    reference.add_argument("--repeat-gff", action="append", type=Path, default=[])
    reference.add_argument(
        "--low-complexity-bed", action="append", type=Path, default=[]
    )
    reference.add_argument(
        "--exclude-bed", action="append", type=Path, default=[],
        help="additional BED3 blacklist; may be repeated",
    )
    feature_group = reference.add_mutually_exclusive_group()
    feature_group.add_argument("--feature-bed", type=Path)
    feature_group.add_argument("--feature-gff", type=Path)
    reference.add_argument(
        "--feature-type", action="append",
        help="GFF feature to retain (default: gene); may be repeated",
    )
    reference.set_defaults(handler=_run_reference)

    panel = subparsers.add_parser(
        "panel-build",
        help="build nested, genome-balanced analysis panels from the master mask",
    )
    panel.add_argument("--mask", type=Path, required=True)
    panel.add_argument("--genome-sizes", type=Path, required=True)
    panel.add_argument("--out-dir", type=Path, required=True)
    panel.add_argument(
        "--target-mb", type=int, action="append",
        help="panel size in Mb; may be repeated (default: 20, 50 and 100)",
    )
    panel.add_argument("--min-run", type=int, default=300)
    panel.add_argument("--block-size", type=int, default=5_000_000)
    panel.add_argument("--chunk-size", type=int, default=10_000)
    panel.set_defaults(handler=_run_panel)

    counts = subparsers.add_parser(
        "bam-counts",
        help="extract per-site allele read counts and adaptive BAF filters",
    )
    counts.add_argument("--bam", type=Path, required=True)
    counts.add_argument("--reference", type=Path, required=True)
    counts.add_argument("--mask", type=Path, required=True)
    counts.add_argument("--out-dir", type=Path, required=True)
    counts.add_argument("--sample")
    counts.add_argument("--min-mapq", type=int, default=20)
    counts.add_argument("--min-baseq", type=int, default=20)
    counts.add_argument("--min-depth", type=int)
    counts.add_argument("--max-depth", type=int)
    counts.add_argument("--allele-min", type=float, default=0.15)
    counts.add_argument("--allele-max", type=float, default=0.85)
    counts.add_argument(
        "--pair-purity", type=float, default=0.95,
        help="minimum fraction of A/C/G/T depth explained by the top two alleles",
    )
    counts.add_argument(
        "--min-allele-reads", type=int, default=2,
        help="minimum read count for each retained allele (default: 2)",
    )
    counts.add_argument("--max-depth-cap", type=int, default=1_000_000)
    counts.add_argument(
        "--region", help="optional pilot region, for example chr1:1-1000000"
    )
    counts.add_argument("--samtools-executable", default="samtools")
    counts.add_argument(
        "--mode", choices=("panel", "full"), default="panel",
        help="panel is the scalable default; full scans the supplied master mask",
    )
    counts.set_defaults(handler=_run_bam_counts)

    model = subparsers.add_parser(
        "ploidy-model",
        help="fit aggregated beta-binomial 2x/3x/4x hypotheses",
    )
    model.add_argument("--counts", type=Path, required=True)
    model.add_argument("--out-dir", type=Path, required=True)
    model.add_argument("--sample")
    model.add_argument("--max-sites", type=int, default=DEFAULT_MAX_SITES)
    model.add_argument("--min-sites", type=int, default=DEFAULT_MIN_SITES)
    model.add_argument("--allele-min", type=float, default=0.15)
    model.add_argument("--allele-max", type=float, default=0.85)
    model.add_argument("--rho-min", type=float, default=1e-5)
    model.add_argument("--rho-max", type=float, default=0.25)
    model.add_argument(
        "--a-min", type=float, default=4 / 7,
        help="minimum M4P outer-peak fraction (default: 4/7)",
    )
    model.add_argument(
        "--sensitivity-a-min", type=float, default=0.5001,
        help="weaker M4P boundary used only for sensitivity analysis",
    )
    model.add_argument("--block-size", type=int, default=5_000_000)
    model.add_argument("--bootstrap", type=int, default=100)
    model.add_argument("--bootstrap-support", type=float, default=0.85, help="diagnostic warning threshold for bootstrap agreement; does not veto calls (default: 0.85)")
    model.add_argument("--delta-bic", type=float, default=10.0)
    model.add_argument("--min-blocks", type=int, default=5)
    model.add_argument("--chromosome-support", type=float, default=0.70, help="diagnostic warning threshold for chromosome agreement; does not veto calls (default: 0.70)")
    model.add_argument("--min-chromosome-sites", type=int, default=200)
    model.add_argument("--seed", type=int, default=20260803)
    model.set_defaults(handler=_run_ploidy_model)

    simulate = subparsers.add_parser(
        "simulate-model",
        help="simulate and validate the expected 2x/3x/4x allele-frequency peaks",
    )
    simulate.add_argument("--out-dir", type=Path, required=True)
    simulate.add_argument("--sites", type=int, default=20_000)
    simulate.add_argument("--mean-depth", type=float, default=30)
    simulate.add_argument("--rho", type=float, default=0.02)
    simulate.add_argument(
        "--contamination-fraction", type=float, default=0.0,
        help="optional uniform contamination used only for simulation stress tests",
    )
    simulate.add_argument(
        "--outer-fraction", type=float, default=0.65,
        help="4x fraction assigned jointly to the 0.25 and 0.75 peaks",
    )
    simulate.add_argument("--seed", type=int, default=20260803)
    simulate.set_defaults(handler=_run_simulation)
    return parser


def _run_reference(arguments: argparse.Namespace) -> int:
    outputs = build_reference_mask(
        ReferenceMaskConfig(
            reference=arguments.reference,
            output_dir=arguments.out_dir,
            kmer_length=arguments.kmer_length,
            min_region_length=arguments.min_region_length,
            mismatches=arguments.mismatches,
            threads=arguments.threads,
            genmap_executable=arguments.genmap_executable,
            base_rule=arguments.base_rule,
            gap_padding=arguments.gap_padding,
            repeat_beds=arguments.repeat_bed,
            repeat_gffs=arguments.repeat_gff,
            low_complexity_beds=arguments.low_complexity_bed,
            exclude_beds=arguments.exclude_bed,
            feature_bed=arguments.feature_bed,
            feature_gff=arguments.feature_gff,
            feature_types=tuple(arguments.feature_type or ["gene"]),
        )
    )
    for name, path in outputs.items():
        if path.exists():
            print(f"{name}\t{path}")
    return 0


def _run_panel(arguments: argparse.Namespace) -> int:
    outputs = build_panels(
        PanelConfig(
            mask=arguments.mask,
            genome_sizes=arguments.genome_sizes,
            output_dir=arguments.out_dir,
            target_mb=tuple(arguments.target_mb or (20, 50, 100)),
            min_run=arguments.min_run,
            block_size=arguments.block_size,
            chunk_size=arguments.chunk_size,
        )
    )
    for name, path in outputs.items():
        if path.exists():
            print(f"{name}\t{path}")
    return 0


def _run_bam_counts(arguments: argparse.Namespace) -> int:
    outputs = build_bam_counts(
        BamCountConfig(
            bam=arguments.bam,
            reference=arguments.reference,
            mask=arguments.mask,
            output_dir=arguments.out_dir,
            sample=arguments.sample,
            min_mapq=arguments.min_mapq,
            min_baseq=arguments.min_baseq,
            min_depth=arguments.min_depth,
            max_depth=arguments.max_depth,
            allele_min=arguments.allele_min,
            allele_max=arguments.allele_max,
            pair_purity=arguments.pair_purity,
            min_allele_reads=arguments.min_allele_reads,
            max_depth_cap=arguments.max_depth_cap,
            region=arguments.region,
            samtools_executable=arguments.samtools_executable,
            mode=arguments.mode,
        )
    )
    for name, path in outputs.items():
        if path.exists():
            print(f"{name}\t{path}")
    return 0


def _run_ploidy_model(arguments: argparse.Namespace) -> int:
    outputs = run_ploidy_model(
        PloidyModelConfig(
            counts=arguments.counts,
            output_dir=arguments.out_dir,
            sample=arguments.sample,
            max_sites=arguments.max_sites,
            min_sites=arguments.min_sites,
            allele_min=arguments.allele_min,
            allele_max=arguments.allele_max,
            rho_min=arguments.rho_min,
            rho_max=arguments.rho_max,
            tetraploid_a_min=arguments.a_min,
            sensitivity_a_min=arguments.sensitivity_a_min,
            block_size=arguments.block_size,
            bootstrap_replicates=arguments.bootstrap,
            bootstrap_support=arguments.bootstrap_support,
            delta_bic=arguments.delta_bic,
            min_blocks=arguments.min_blocks,
            chromosome_support=arguments.chromosome_support,
            min_chromosome_sites=arguments.min_chromosome_sites,
            seed=arguments.seed,
        )
    )
    for name, path in outputs.items():
        if path.exists():
            print(f"{name}\t{path}")
    return 0


def _run_simulation(arguments: argparse.Namespace) -> int:
    outputs = simulate_and_validate(
        output_dir=arguments.out_dir,
        sites=arguments.sites,
        mean_depth=arguments.mean_depth,
        rho=arguments.rho,
        contamination_fraction=arguments.contamination_fraction,
        outer_fraction=arguments.outer_fraction,
        seed=arguments.seed,
    )
    for name, path in outputs.items():
        if path.exists():
            print(f"{name}\t{path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        return arguments.handler(arguments)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
