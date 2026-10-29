"""Build a k-mer-aware reference-genome mappability mask.

The workflow has one responsibility:

1. run or reuse GenMap;
2. retain score-1 k-mer runs;
3. remove FASTA gaps and optional repeat/low-complexity regions;
4. write the final BED mask, excluded regions, statistics and a manifest.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
import json
import re
import shutil
import subprocess
import sys

from . import __version__
from .intervals import (
    annotate_bed,
    bed_stats,
    complement_bed,
    concatenate_files,
    copy_bed3,
    gff_to_bed,
    intersect_sorted_bed,
    iter_bed,
    merge_sorted_bed,
    pad_sorted_bed,
    sort_bed,
    subtract_sorted_bed,
)


@dataclass(slots=True)
class ReferenceMaskConfig:
    reference: Path
    output_dir: Path
    kmer_length: int = 150
    mismatches: int = 2
    threads: int = 8
    genmap_executable: str = "genmap"
    genmap_bedgraph: Path | None = None
    base_rule: str = "strict"
    gap_padding: int | None = None
    repeat_beds: list[Path] = field(default_factory=list)
    repeat_gffs: list[Path] = field(default_factory=list)
    low_complexity_beds: list[Path] = field(default_factory=list)
    exclude_beds: list[Path] = field(default_factory=list)
    feature_bed: Path | None = None
    feature_gff: Path | None = None
    feature_types: tuple[str, ...] = ("gene",)
    min_region_length: int = 1

    def validate(self) -> None:
        if self.kmer_length < 1:
            raise ValueError("k-mer length must be positive")
        if not 0 <= self.mismatches < self.kmer_length:
            raise ValueError("mismatches must be in [0, kmer_length)")
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if self.min_region_length < 1:
            raise ValueError("minimum continuous region length must be positive")
        if self.base_rule not in {"strict", "start"}:
            raise ValueError("base rule must be 'strict' or 'start'")
        if self.gap_padding is not None and self.gap_padding < 0:
            raise ValueError("gap padding cannot be negative")
        if self.feature_bed and self.feature_gff:
            raise ValueError("use only one of --feature-bed and --feature-gff")


def build_reference_mask(config: ReferenceMaskConfig) -> dict[str, Path]:
    config.validate()
    reference = config.reference.resolve()
    output_dir = config.output_dir.resolve()
    if not reference.is_file():
        raise FileNotFoundError(f"reference FASTA does not exist: {reference}")

    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "intermediate"
    work_dir.mkdir(exist_ok=True)
    outputs = _output_paths(output_dir)
    signature = _run_signature(config, reference)
    _check_or_write_signature(output_dir / "run_config.json", signature)

    required_outputs = [
        outputs["mappability_mask"],
        outputs["excluded_regions"],
        outputs["exclusion_sources"],
        outputs["mask_statistics"],
        outputs["manifest"],
    ]
    if config.feature_bed or config.feature_gff:
        required_outputs.append(outputs["feature_mask"])
    if all(path.is_file() for path in required_outputs):
        print(f"[MapPloid] completed result already exists; reusing {output_dir}")
        return outputs

    print("[MapPloid] scanning FASTA gaps and contig lengths", file=sys.stderr)
    raw_gaps = work_dir / "assembly_gaps.raw.bed"
    genome_lengths = scan_fasta(reference, raw_gaps)
    _write_genome_sizes(genome_lengths, output_dir / "genome.sizes")

    sorted_gaps = work_dir / "assembly_gaps.sorted.bed"
    merged_gaps = work_dir / "assembly_gaps.merged.bed"
    padded_gaps = work_dir / "assembly_gaps.padded.bed"
    sort_bed(raw_gaps, sorted_gaps)
    merge_sorted_bed(sorted_gaps, merged_gaps)
    padding = config.kmer_length - 1 if config.gap_padding is None else config.gap_padding
    pad_sorted_bed(merged_gaps, padded_gaps, genome_lengths, padding, padding)

    bedgraph, genmap_status = get_genmap_bedgraph(config, reference, work_dir)
    sequence_mask = work_dir / "genmap_unique_mask.bed"
    if not sequence_mask.is_file() or sequence_mask.stat().st_size == 0:
        print(f"[MapPloid] converting GenMap bedGraph: {bedgraph}", file=sys.stderr)
        raw_sequence_mask = work_dir / "genmap_unique_mask.raw.bed"
        unique_runs_to_mask(
            bedgraph,
            raw_sequence_mask,
            config.kmer_length,
            strict=config.base_rule == "strict",
        )
        sorted_sequence_mask = work_dir / "genmap_unique_mask.sorted.bed"
        sort_bed(raw_sequence_mask, sorted_sequence_mask)
        merge_sorted_bed(sorted_sequence_mask, sequence_mask)
    else:
        print(f"[MapPloid] reusing intermediate mask: {sequence_mask}", file=sys.stderr)
    validate_bed(sequence_mask, genome_lengths)

    blacklist_sources: list[tuple[Path, str]] = [(padded_gaps, "assembly_gap")]
    blacklist_sources += normalise_annotations(
        config.repeat_beds, config.repeat_gffs, work_dir, "repeat", genome_lengths
    )
    blacklist_sources += normalise_annotations(
        config.low_complexity_beds, [], work_dir, "low_complexity", genome_lengths
    )
    blacklist_sources += normalise_annotations(
        config.exclude_beds, [], work_dir, "user_excluded", genome_lengths
    )

    blacklist_concat = work_dir / "blacklist.concat.bed"
    blacklist_sorted = work_dir / "blacklist.sorted.bed"
    blacklist_merged = work_dir / "blacklist.merged.bed"
    concatenate_files([path for path, _ in blacklist_sources], blacklist_concat)
    sort_bed(blacklist_concat, blacklist_sorted)
    merge_sorted_bed(blacklist_sorted, blacklist_merged)
    after_exclusions = work_dir / "mask_after_exclusions.bed"
    subtract_sorted_bed(sequence_mask, blacklist_merged, after_exclusions)
    short_regions = work_dir / "short_continuous_regions.bed"
    filter_continuous_regions(
        after_exclusions, outputs["mappability_mask"], short_regions,
        config.min_region_length,
    )
    blacklist_sources.append((short_regions, "short_continuous_region"))
    complement_bed(
        outputs["mappability_mask"], genome_lengths, outputs["excluded_regions"]
    )

    reason_files: list[Path] = []
    low_map = work_dir / "low_mappability_or_edge.bed"
    complement_bed(sequence_mask, genome_lengths, low_map)
    for number, (path, reason) in enumerate(
        [(low_map, "low_mappability_or_edge"), *blacklist_sources], start=1
    ):
        annotated = work_dir / f"reason_{number:03d}_{reason}.bed"
        annotate_bed(path, annotated, reason)
        reason_files.append(annotated)
    reasons_concat = work_dir / "exclusion_sources.concat.bed"
    concatenate_files(reason_files, reasons_concat)
    sort_bed(reasons_concat, outputs["exclusion_sources"])

    feature_mask: Path | None = None
    if config.feature_bed or config.feature_gff:
        features = normalise_features(config, work_dir, genome_lengths)
        feature_mask = outputs["feature_mask"]
        intersect_sorted_bed(outputs["mappability_mask"], features, feature_mask)

    write_statistics(
        outputs["mask_statistics"], genome_lengths, sequence_mask,
        outputs["mappability_mask"], feature_mask, blacklist_sources
    )
    manifest = {
        "workflow": "reference-mask",
        "mapploid_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "run_config": signature,
        "genmap": {
            "bedgraph": str(bedgraph),
            "status": genmap_status,
        },
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    outputs["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"[MapPloid] reference mask completed: {outputs['mappability_mask']}",
        file=sys.stderr,
    )
    return outputs


def filter_continuous_regions(source: Path, kept: Path, removed: Path, minimum: int) -> None:
    """Filter complete post-exclusion runs, including exactly minimum bp."""
    if minimum < 1:
        raise ValueError("minimum continuous region length must be positive")
    with kept.open("wt", encoding="utf-8") as good, removed.open("wt", encoding="utf-8") as short:
        for interval in iter_bed(source):
            destination = good if interval.length >= minimum else short
            destination.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")


def get_genmap_bedgraph(
    config: ReferenceMaskConfig, reference: Path, work_dir: Path
) -> tuple[Path, str]:
    if config.genmap_bedgraph:
        bedgraph = config.genmap_bedgraph.resolve()
        if not bedgraph.is_file() or bedgraph.stat().st_size == 0:
            raise FileNotFoundError(f"GenMap bedGraph is missing or empty: {bedgraph}")
        return bedgraph, "provided"

    prefix = work_dir / f"genmap_k{config.kmer_length}_e{config.mismatches}"
    bedgraph = Path(f"{prefix}.bedgraph")
    if bedgraph.is_file() and bedgraph.stat().st_size > 0:
        print(f"[MapPloid] reusing completed GenMap bedGraph: {bedgraph}", file=sys.stderr)
        return bedgraph, "reused"

    executable = shutil.which(config.genmap_executable)
    if not executable:
        raise RuntimeError(f"GenMap executable not found: {config.genmap_executable}")
    index_dir = work_dir / "genmap_index"
    if not (index_dir / "index.sa.val").is_file():
        command = [
            executable, "index", "-F", str(reference), "-I", str(index_dir),
            "-A", "divsufsort",
        ]
        run_command(command, work_dir / "genmap_index.log")
    else:
        print(f"[MapPloid] reusing GenMap index: {index_dir}", file=sys.stderr)

    # GenMap -O is a file prefix, not an output directory.
    command = [
        executable, "map", "-K", str(config.kmer_length), "-E", str(config.mismatches),
        "-I", str(index_dir), "-O", str(prefix), "-bg", "-T",
        str(config.threads),
    ]
    run_command(command, work_dir / "genmap_map.log")
    if not bedgraph.is_file() or bedgraph.stat().st_size == 0:
        found = sorted(str(path) for path in work_dir.glob("*.bedgraph"))
        raise RuntimeError(
            f"GenMap finished but expected output is absent: {bedgraph}. "
            f"bedGraph files found: {found or 'none'}"
        )
    return bedgraph, "generated"


def run_command(command: list[str], log_path: Path) -> None:
    print(f"[MapPloid] running: {' '.join(command)}", file=sys.stderr)
    with log_path.open("wt", encoding="utf-8") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise RuntimeError(
            f"command failed with exit code {result.returncode}; see {log_path}"
        )


def unique_runs_to_mask(
    bedgraph: Path, output: Path, kmer_length: int, strict: bool = True
) -> None:
    """Stream a GenMap bedGraph and retain continuous score-1 runs.

    GenMap output can contain tens of millions of rows. Bytes and an exact
    string comparison to ``1`` avoid expensive float parsing.
    """

    trim = kmer_length - 1 if strict else 0
    chrom: bytes | None = None
    run_start = 0
    run_end = 0
    line_count = 0

    def flush(handle) -> None:
        nonlocal chrom
        if chrom is not None and run_end > run_start + trim:
            handle.write(
                chrom + b"\t" + str(run_start + trim).encode() + b"\t"
                + str(run_end).encode() + b"\n"
            )
        chrom = None

    with bedgraph.open("rb", buffering=8 * 1024 * 1024) as source, output.open(
        "wb", buffering=8 * 1024 * 1024
    ) as destination:
        for line in source:
            line_count += 1
            if line_count % 10_000_000 == 0:
                print(
                    f"[MapPloid] parsed {line_count:,} bedGraph rows",
                    file=sys.stderr,
                )
            last_tab = line.rfind(b"\t")
            if last_tab < 0 or line[last_tab + 1 :].strip() != b"1":
                flush(destination)
                continue
            fields = line[:last_tab].split(b"\t")
            if len(fields) != 3:
                raise ValueError(f"invalid GenMap bedGraph line {line_count}")
            current_chrom, start_text, end_text = fields
            start = int(start_text)
            end = int(end_text)
            if chrom == current_chrom and start <= run_end:
                run_end = max(run_end, end)
            else:
                flush(destination)
                chrom = current_chrom
                run_start = start
                run_end = end
        flush(destination)
    # A reference can legitimately have no uniquely mappable regions.


def scan_fasta(reference: Path, gaps_bed: Path) -> dict[str, int]:
    """Write every run of non-ACGT symbols and return contig lengths."""

    invalid_bases = re.compile(r"[^ACGTacgt]+")
    lengths: dict[str, int] = {}
    name: str | None = None
    position = 0
    gap_start: int | None = None

    def finish_gap(handle, end: int) -> None:
        nonlocal gap_start
        if gap_start is not None:
            handle.write(f"{name}\t{gap_start}\t{end}\n")
            gap_start = None

    with reference.open("rt", encoding="utf-8") as source, gaps_bed.open(
        "wt", encoding="utf-8"
    ) as gaps:
        for raw in source:
            if raw.startswith(">"):
                if name is not None:
                    finish_gap(gaps, position)
                    lengths[name] = position
                name = raw[1:].strip().split()[0]
                if not name or name in lengths:
                    raise ValueError(f"invalid or duplicate FASTA name: {name!r}")
                position = 0
                gap_start = None
                continue
            if name is None:
                if raw.strip():
                    raise ValueError("FASTA sequence appears before its header")
                continue
            sequence = raw.strip()
            search_from = 0

            # Continue or close an ambiguous run from the previous FASTA line.
            if gap_start is not None:
                prefix = invalid_bases.match(sequence)
                if prefix and prefix.end() == len(sequence):
                    position += len(sequence)
                    continue
                if prefix:
                    finish_gap(gaps, position + prefix.end())
                    search_from = prefix.end()
                else:
                    finish_gap(gaps, position)

            for match in invalid_bases.finditer(sequence, search_from):
                start = position + match.start()
                end = position + match.end()
                if match.end() == len(sequence):
                    gap_start = start
                else:
                    gaps.write(f"{name}\t{start}\t{end}\n")
            position += len(sequence)
        if name is None:
            raise ValueError("reference FASTA contains no sequences")
        finish_gap(gaps, position)
        lengths[name] = position
    return lengths


def normalise_annotations(
    beds: list[Path], gffs: list[Path], work_dir: Path, reason: str,
    genome_lengths: dict[str, int]
) -> list[tuple[Path, str]]:
    results: list[tuple[Path, str]] = []
    annotations = [(path, "bed") for path in beds]
    annotations.extend((path, "gff") for path in gffs)
    for number, (path, kind) in enumerate(annotations, start=1):
        if not path.is_file():
            raise FileNotFoundError(f"annotation does not exist: {path}")
        raw = work_dir / f"{reason}_{number:03d}.raw.bed"
        copy_bed3(path, raw) if kind == "bed" else gff_to_bed(path, raw)
        sorted_path = work_dir / f"{reason}_{number:03d}.sorted.bed"
        merged = work_dir / f"{reason}_{number:03d}.merged.bed"
        sort_bed(raw, sorted_path)
        merge_sorted_bed(sorted_path, merged)
        validate_bed(merged, genome_lengths)
        results.append((merged, reason))
    return results


def normalise_features(
    config: ReferenceMaskConfig, work_dir: Path, genome_lengths: dict[str, int]
) -> Path:
    raw = work_dir / "features.raw.bed"
    if config.feature_bed:
        copy_bed3(config.feature_bed, raw)
    else:
        assert config.feature_gff is not None
        gff_to_bed(config.feature_gff, raw, set(config.feature_types))
    sorted_path = work_dir / "features.sorted.bed"
    merged = work_dir / "features.merged.bed"
    sort_bed(raw, sorted_path)
    merge_sorted_bed(sorted_path, merged)
    validate_bed(merged, genome_lengths)
    return merged


def validate_bed(path: Path, genome_lengths: dict[str, int]) -> None:
    for interval in iter_bed(path):
        if interval.chrom not in genome_lengths:
            raise ValueError(f"{path}: unknown contig {interval.chrom!r}")
        if interval.end > genome_lengths[interval.chrom]:
            raise ValueError(f"{path}: interval exceeds {interval.chrom} length")


def write_statistics(
    output: Path, genome_lengths: dict[str, int], sequence_mask: Path, final_mask: Path,
    feature_mask: Path | None, blacklist_sources: list[tuple[Path, str]]
) -> None:
    genome_bases = sum(genome_lengths.values())
    rows = [
        ("reference_genome", len(genome_lengths), genome_bases),
        ("genmap_unique_mask", *bed_stats(sequence_mask)),
        ("mappability_mask", *bed_stats(final_mask)),
    ]
    if feature_mask:
        rows.append(("mappability_feature_mask", *bed_stats(feature_mask)))
    for path, reason in blacklist_sources:
        rows.append((f"excluded_{reason}", *bed_stats(path)))
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("stage\tintervals\tbases\tgenome_fraction\n")
        for stage, intervals, bases in rows:
            handle.write(
                f"{stage}\t{intervals}\t{bases}\t{bases / genome_bases:.8f}\n"
            )


def _run_signature(config: ReferenceMaskConfig, reference: Path) -> dict[str, object]:
    stat = reference.stat()
    return {
        "reference": str(reference),
        "reference_size": stat.st_size,
        "reference_mtime_ns": stat.st_mtime_ns,
        "kmer_length": config.kmer_length,
        "min_region_length": config.min_region_length,
        "mismatches": config.mismatches,
        "base_rule": config.base_rule,
        "gap_padding": config.gap_padding,
        "genmap_bedgraph": (
            str(config.genmap_bedgraph.resolve()) if config.genmap_bedgraph else None
        ),
        "repeat_beds": [str(path.resolve()) for path in config.repeat_beds],
        "repeat_gffs": [str(path.resolve()) for path in config.repeat_gffs],
        "low_complexity_beds": [str(path.resolve()) for path in config.low_complexity_beds],
        "exclude_beds": [str(path.resolve()) for path in config.exclude_beds],
        "feature_bed": str(config.feature_bed.resolve()) if config.feature_bed else None,
        "feature_gff": str(config.feature_gff.resolve()) if config.feature_gff else None,
        "feature_types": list(config.feature_types),
    }


def _check_or_write_signature(path: Path, signature: dict[str, object]) -> None:
    if path.is_file():
        previous = json.loads(path.read_text(encoding="utf-8"))
        if previous != signature:
            raise ValueError(
                f"existing output directory was created with different inputs or parameters: {path.parent}"
            )
        return
    path.write_text(json.dumps(signature, indent=2) + "\n", encoding="utf-8")


def _write_genome_sizes(lengths: dict[str, int], output: Path) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        for name, length in lengths.items():
            handle.write(f"{name}\t{length}\n")


def _output_paths(output_dir: Path) -> dict[str, Path]:
    return {
        "mappability_mask": output_dir / "mappability_mask.bed",
        "excluded_regions": output_dir / "excluded_regions.bed",
        "exclusion_sources": output_dir / "exclusion_sources.bed",
        "mask_statistics": output_dir / "mask_statistics.tsv",
        "feature_mask": output_dir / "mappability_feature_mask.bed",
        "manifest": output_dir / "run_manifest.json",
    }
