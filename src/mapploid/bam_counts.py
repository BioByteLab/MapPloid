"""Stream BAM pileups into compact candidate allele-count tables.

The default workflow never materialises counts for every reference position.
It records the depth distribution while retaining only sites with exactly two
alleles inside the requested frequency truncation.  Coordinates in TSV output
are one-based.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import gzip
import json
import math
import re
import shutil
import subprocess
import sys

from . import __version__
from .intervals import iter_bed


BASES = ("A", "C", "G", "T")
BASE_INDEX = {base: index for index, base in enumerate(BASES)}
BASE_COUNT_HEADER = (
    "chrom", "position", "ref", "mpileup_depth", "depth", "A", "C", "G", "T"
)
SITE_HEADER = BASE_COUNT_HEADER + (
    "pair_depth", "pair_purity", "residual_count",
    "allele1", "allele1_count", "allele1_frequency",
    "allele2", "allele2_count", "allele2_frequency", "minor_baf",
)


@dataclass(slots=True)
class BamCountConfig:
    bam: Path
    reference: Path
    mask: Path
    output_dir: Path
    sample: str | None = None
    min_mapq: int = 20
    min_baseq: int = 20
    min_depth: int | None = None
    max_depth: int | None = None
    allele_min: float = 0.15
    allele_max: float = 0.85
    pair_purity: float = 0.95
    min_allele_reads: int = 2
    max_depth_cap: int = 1_000_000
    region: str | None = None
    samtools_executable: str = "samtools"
    mode: str = "panel"

    def validate(self) -> None:
        if not 0 <= self.min_mapq <= 255:
            raise ValueError("minimum MAPQ must be in [0, 255]")
        if not 0 <= self.min_baseq <= 255:
            raise ValueError("minimum base quality must be in [0, 255]")
        if self.min_depth is not None and self.min_depth < 1:
            raise ValueError("minimum depth must be positive")
        if self.max_depth is not None and self.max_depth < 1:
            raise ValueError("maximum depth must be positive")
        if self.min_depth and self.max_depth and self.min_depth > self.max_depth:
            raise ValueError("minimum depth cannot exceed maximum depth")
        if not 0 <= self.allele_min < self.allele_max <= 1:
            raise ValueError("allele truncation must satisfy 0 <= min < max <= 1")
        if self.allele_min >= 0.5 or self.allele_max <= 0.5:
            raise ValueError("allele truncation must contain 0.5")
        if not 0.5 <= self.pair_purity <= 1:
            raise ValueError("pair purity must be in [0.5, 1]")
        if self.min_allele_reads < 1:
            raise ValueError("minimum allele reads must be positive")
        if self.max_depth_cap < 1:
            raise ValueError("maximum pileup depth must be positive")
        if self.mode not in {"panel", "full"}:
            raise ValueError("BAM count mode must be 'panel' or 'full'")


def build_bam_counts(config: BamCountConfig) -> dict[str, Path]:
    """Extract counts, estimate depth limits, and retain biallelic sites."""

    config.validate()
    bam = config.bam.resolve()
    reference = config.reference.absolute()
    mask = config.mask.resolve()
    output_dir = config.output_dir.resolve()
    for path, label in ((bam, "BAM"), (reference, "reference"), (mask, "mask")):
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")

    sample = config.sample or _sample_from_bam(bam)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
        raise ValueError("sample name may contain only letters, digits, '.', '_' and '-'")
    executable = shutil.which(config.samtools_executable)
    if executable is None:
        raise RuntimeError(f"samtools executable not found: {config.samtools_executable}")

    mask_bases, warnings = validate_inputs(bam, reference, mask, config.region)
    if config.mode == "full":
        warnings.append("full-mask mode used; intended for validation or uncertain samples")
    output_dir.mkdir(parents=True, exist_ok=True)
    intermediate = output_dir / "intermediate"
    intermediate.mkdir(exist_ok=True)
    outputs = _output_paths(output_dir, sample)

    signature = _signature(config, sample, bam, reference, mask)
    _check_or_write_signature(output_dir / "run_config.json", signature)
    required = ["depth_distribution", "thresholds", "biallelic_sites",
                "filter_summary", "manifest"]
    if all(outputs[name].is_file() for name in required):
        print(f"[MapPloid] completed BAM counts already exist: {output_dir}")
        return outputs

    state_path = intermediate / "extraction_summary.json"
    candidates = intermediate / f"{sample}.candidate_sites.tsv.gz"
    reusable = candidates.is_file() and state_path.is_file()
    if reusable:
        extraction = json.loads(state_path.read_text(encoding="utf-8"))
        print("[MapPloid] reusing candidate-count extraction", file=sys.stderr)
    else:
        extraction = extract_base_counts(
            config, executable, bam, reference, mask, candidates, intermediate,
        )
        state_path.write_text(json.dumps(extraction, indent=2) + "\n", encoding="utf-8")

    observed_histogram = {int(depth): int(count) for depth, count in
                          dict(extraction["depth_histogram"]).items()}
    histogram = dict(observed_histogram)
    if config.region is None:
        histogram[0] = histogram.get(0, 0) + max(0, mask_bases - int(extraction["pileup_sites"]))
    thresholds = estimate_depth_thresholds(histogram, config.min_depth, config.max_depth)
    thresholds.update({
        "allele_min": config.allele_min,
        "allele_max": config.allele_max,
        "pair_purity": config.pair_purity,
        "min_allele_reads": config.min_allele_reads,
        "mask_bases": mask_bases,
        "zero_depth_inferred": config.region is None,
    })
    write_depth_distribution(outputs["depth_distribution"], histogram)
    write_thresholds(outputs["thresholds"], thresholds)
    filtering = filter_sites(config, candidates, outputs["biallelic_sites"], thresholds)
    minimum_depth = int(thresholds["min_depth"])
    maximum_depth = int(thresholds["max_depth"])
    filtering.update({
        "low_depth": sum(count for depth, count in observed_histogram.items()
                         if depth < minimum_depth),
        "high_depth": sum(count for depth, count in observed_histogram.items()
                          if depth > maximum_depth),
    })
    category_histograms = dict(extraction["category_depth_histograms"])
    for category in ("not_biallelic", "multiallelic", "low_pair_purity"):
        category_histogram = {
            int(depth): int(count)
            for depth, count in dict(category_histograms[category]).items()
        }
        filtering[category] = sum(
            count for depth, count in category_histogram.items()
            if minimum_depth <= depth <= maximum_depth
        )
    write_filter_summary(outputs["filter_summary"], extraction, filtering)

    manifest = {
        "workflow": "bam-counts-stream",
        "mapploid_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "run_config": signature,
        "warnings": warnings,
        "extraction": {
            key: value for key, value in extraction.items()
            if key not in {"depth_histogram", "category_depth_histograms"}
        },
        "thresholds": thresholds,
        "filtering": filtering,
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n",
                                   encoding="utf-8")
    print(
        f"[MapPloid] {sample}: retained {filtering['retained_sites']:,} "
        f"biallelic sites ({thresholds['min_depth']} <= depth <= "
        f"{thresholds['max_depth']})",
        file=sys.stderr,
    )
    return outputs


def validate_inputs(
    bam: Path, reference: Path, mask: Path, region: str | None, *, allow_empty: bool = False
) -> tuple[int, list[str]]:
    """Check BAM/reference compatibility and return total mask length."""

    import pysam

    try:
        pysam.quickcheck("-q", str(bam))
    except pysam.SamtoolsError as exc:
        raise ValueError(f"BAM failed quickcheck: {bam}") from exc

    warnings: list[str] = []
    with pysam.AlignmentFile(bam, "rb") as alignment:
        header = alignment.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError("BAM header does not declare coordinate sorting")
        if not alignment.has_index():
            raise ValueError(f"BAM index is missing: {bam}")
        bam_lengths = dict(zip(alignment.references, alignment.lengths))
    with pysam.FastaFile(str(reference)) as fasta:
        reference_lengths = dict(zip(fasta.references, fasta.lengths))
    if bam_lengths != reference_lengths:
        raise ValueError("BAM and reference sequence dictionaries differ")

    previous: tuple[str, int, int] | None = None
    mask_bases = 0
    for interval in iter_bed(mask):
        if interval.chrom not in bam_lengths:
            raise ValueError(f"mask contig is absent from BAM: {interval.chrom}")
        if interval.end > bam_lengths[interval.chrom]:
            raise ValueError(f"mask interval exceeds {interval.chrom} length")
        # The reference-mask workflow sorts BED with LC_ALL=C, so validate
        # against that same lexical contig order rather than BAM @SQ order.
        key = (interval.chrom, interval.start, interval.end)
        if previous and (key[0] < previous[0] or
                         (key[0] == previous[0] and key[1] < previous[1])):
            raise ValueError("mask BED is not coordinate sorted")
        if previous and key[0] == previous[0] and key[1] < previous[2]:
            raise ValueError("mask BED contains overlapping intervals")
        previous = key
        mask_bases += interval.length

    if mask_bases == 0 and not allow_empty:
        raise ValueError("mask BED contains no bases")
    index = _find_bam_index(bam)
    if index and index.stat().st_mtime < bam.stat().st_mtime:
        warnings.append("BAM index timestamp is older than BAM")
    if not any("M5" in row for row in header.get("SQ", [])):
        warnings.append("BAM @SQ has no M5 checksum; names and lengths were compared")
    if region is not None:
        warnings.append("pilot region used; zero-depth mask bases were not inferred")
    return mask_bases, warnings


def extract_base_counts(
    config: BamCountConfig,
    executable: str,
    bam: Path,
    reference: Path,
    mask: Path,
    candidates: Path,
    intermediate: Path,
) -> dict[str, object]:
    """Stream pileup once, retaining only frequency-qualified candidates."""

    command = [
        executable, "mpileup",
        "-B",
        "-f", str(reference),
        "-l", str(mask),
        "-q", str(config.min_mapq),
        "-Q", str(config.min_baseq),
        "-d", str(config.max_depth_cap),
        "--ff", "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
        "--no-output-ins", "--no-output-ins",
        "--no-output-del", "--no-output-del",
        "--no-output-ends",
    ]
    if config.region:
        command.extend(("-r", config.region))
    command.append(str(bam))

    candidate_temporary = Path(f"{candidates}.tmp")
    log_path = intermediate / "samtools_mpileup.log"
    histogram: Counter[int] = Counter()
    category_histograms: dict[str, Counter[int]] = {
        "candidate": Counter(),
        "not_biallelic": Counter(),
        "multiallelic": Counter(),
        "low_pair_purity": Counter(),
    }
    pileup_sites = 0
    acgt_observations = 0
    candidate_sites = 0
    print(f"[MapPloid] running samtools mpileup for {bam.name}", file=sys.stderr)
    try:
        with log_path.open("wt", encoding="utf-8") as log, subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=log,
            text=True,
            encoding="utf-8",
            bufsize=1024 * 1024,
        ) as process, gzip.open(
            candidate_temporary, "wt", encoding="utf-8", compresslevel=1
        ) as candidate_handle:
            assert process.stdout is not None
            candidate_handle.write("\t".join(BASE_COUNT_HEADER) + "\n")
            try:
                for line_number, line in enumerate(process.stdout, start=1):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        raise RuntimeError(f"malformed mpileup row {line_number}")
                    try:
                        mpileup_depth = int(fields[3])
                    except ValueError as exc:
                        raise RuntimeError(
                            f"invalid mpileup depth at row {line_number}"
                        ) from exc
                    ref = fields[2].upper()
                    counts = count_simple_pileup_bases(fields[4], ref)
                    depth = sum(counts)
                    pileup_sites += 1
                    acgt_observations += depth
                    histogram[depth] += 1
                    category = classify_site(
                        counts,
                        depth,
                        config.allele_min,
                        config.allele_max,
                        config.pair_purity,
                    )
                    category_histograms[category][depth] += 1
                    if category == "candidate":
                        row = (
                            f"{fields[0]}\t{fields[1]}\t{ref}\t{mpileup_depth}\t"
                            f"{depth}\t" + "\t".join(str(value) for value in counts)
                            + "\n"
                        )
                    if category == "candidate":
                        candidate_handle.write(row)
                        candidate_sites += 1
                    if line_number % 10_000_000 == 0:
                        print(
                            f"[MapPloid] counted {line_number:,} panel positions; "
                            f"{candidate_sites:,} candidates",
                            file=sys.stderr,
                        )
                return_code = process.wait()
            except BaseException:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                raise
    except BaseException:
        candidate_temporary.unlink(missing_ok=True)
        raise
    if return_code:
        candidate_temporary.unlink(missing_ok=True)
        raise RuntimeError(f"samtools mpileup failed; see {log_path}")
    candidate_temporary.replace(candidates)
    return {
        "pileup_sites": pileup_sites,
        "acgt_observations": acgt_observations,
        "candidate_sites": candidate_sites,
        "depth_histogram": dict(sorted(histogram.items())),
        "category_depth_histograms": {
            category: dict(sorted(values.items()))
            for category, values in category_histograms.items()
        },
        "command": command,
    }


def classify_site(
    counts: tuple[int, int, int, int],
    depth: int,
    allele_min: float,
    allele_max: float,
    pair_purity: float = 0.95,
) -> str:
    if depth <= 0:
        return "not_biallelic"
    ordered = sorted(counts, reverse=True)
    if ordered[1] <= 0:
        return "not_biallelic"
    significant = sum(count / depth >= allele_min for count in counts)
    if significant > 2:
        return "multiallelic"
    pair_depth = ordered[0] + ordered[1]
    if pair_depth / depth < pair_purity:
        return "low_pair_purity"
    second_fraction = ordered[1] / pair_depth
    if not allele_min <= second_fraction <= allele_max:
        return "not_biallelic"
    return "candidate"


def count_pileup_bases(read_bases: str, reference_base: str) -> tuple[int, int, int, int]:
    """Count A/C/G/T in a SAMtools pileup read-bases string."""

    counts = [0, 0, 0, 0]
    index = 0
    while index < len(read_bases):
        symbol = read_bases[index]
        if symbol == "^":                 # start marker plus encoded MAPQ
            index += 2
            continue
        if symbol == "$":                 # end-of-read marker
            index += 1
            continue
        if symbol in "+-":                # indel length followed by sequence
            index += 1
            end = index
            while end < len(read_bases) and read_bases[end].isdigit():
                end += 1
            if end == index:
                raise ValueError("invalid indel in mpileup read-bases column")
            index = end + int(read_bases[index:end])
            continue
        base = reference_base if symbol in ".," else symbol.upper()
        if base in BASE_INDEX:
            counts[BASE_INDEX[base]] += 1
        index += 1
    return tuple(counts)


def count_simple_pileup_bases(
    read_bases: str, reference_base: str
) -> tuple[int, int, int, int]:
    """Count simplified mpileup symbols using C-level ``str.count`` loops."""

    counts = [
        read_bases.count("A") + read_bases.count("a"),
        read_bases.count("C") + read_bases.count("c"),
        read_bases.count("G") + read_bases.count("g"),
        read_bases.count("T") + read_bases.count("t"),
    ]
    reference_count = read_bases.count(".") + read_bases.count(",")
    reference_index = BASE_INDEX.get(reference_base)
    if reference_index is not None:
        counts[reference_index] += reference_count
    return tuple(counts)


def estimate_depth_thresholds(
    histogram: dict[int, int], minimum: int | None, maximum: int | None
) -> dict[str, int | float | str]:
    """Use positive-depth median/MAD and quantiles for robust cutoffs."""

    positive = {depth: count for depth, count in histogram.items()
                if depth > 0 and count > 0}
    if not positive:
        return {"depth_median": 0, "depth_mad": 0, "depth_q05": 0, "depth_q995": 0,
                "min_depth": minimum or 1, "max_depth": maximum or minimum or 1,
                "min_depth_source": "NO_COVERAGE", "max_depth_source": "NO_COVERAGE"}
    median = weighted_quantile(positive, 0.5)
    q05 = weighted_quantile(positive, 0.05)
    q995 = weighted_quantile(positive, 0.995)
    deviations: Counter[int] = Counter()
    for depth, count in positive.items():
        deviations[abs(depth - median)] += count
    mad = weighted_quantile(deviations, 0.5)
    robust_upper = median + 6 * 1.4826 * mad if mad else q995
    robust_lower = median - 1.4826 * mad if mad else q05
    automatic_min = max(1, q05, int(math.ceil(robust_lower)))
    automatic_max = max(automatic_min, int(math.ceil(min(q995, robust_upper))))
    chosen_min = minimum if minimum is not None else automatic_min
    chosen_max = maximum if maximum is not None else automatic_max
    min_source = "user" if minimum is not None else "max(Q05, median-1.4826*MAD)"
    max_source = "user" if maximum is not None else "min(Q99.5, median+6*1.4826*MAD)"
    # A user-specified bound takes precedence over the automatic opposite bound.
    if chosen_min > chosen_max:
        if minimum is None:
            chosen_min = chosen_max
            min_source += ";clamped_to_user_max"
        elif maximum is None:
            chosen_max = chosen_min
            max_source += ";clamped_to_user_min"
        else:
            raise ValueError("minimum depth cannot exceed maximum depth")
    return {
        "depth_median": median,
        "depth_mad": mad,
        "depth_q05": q05,
        "depth_q995": q995,
        "min_depth": chosen_min,
        "max_depth": chosen_max,
        "min_depth_source": min_source,
        "max_depth_source": max_source,
    }


def filter_sites(
    config: BamCountConfig,
    source: Path,
    output: Path,
    thresholds: dict[str, object],
) -> dict[str, int]:
    """Keep depth-PASS sites whose two leading alleles explain the pileup."""

    minimum_depth = int(thresholds["min_depth"])
    maximum_depth = int(thresholds["max_depth"])
    counters: Counter[str] = Counter()
    temporary = Path(f"{output}.tmp")
    with gzip.open(source, "rt", encoding="utf-8") as incoming, gzip.open(
        temporary, "wt", encoding="utf-8", compresslevel=1
    ) as retained:
        header = incoming.readline().rstrip("\n").split("\t")
        indexes = {name: index for index, name in enumerate(header)}
        retained.write("\t".join(SITE_HEADER) + "\n")
        for line in incoming:
            fields = line.rstrip("\n").split("\t")
            depth = int(fields[indexes["depth"]])
            if depth < minimum_depth:
                counters["LOW_DEPTH"] += 1
                continue
            if depth > maximum_depth:
                counters["HIGH_DEPTH"] += 1
                continue
            counts = [int(fields[indexes[base]]) for base in BASES]
            ordered = sorted(range(len(BASES)), key=lambda index: (-counts[index], index))
            first, second = ordered[:2]
            if counts[second] <= 0:
                counters["NOT_BIALLELIC"] += 1
                continue
            significant = sum(count / depth >= config.allele_min for count in counts)
            if significant > 2:
                counters["MULTIALLELIC"] += 1
                continue
            pair_depth = counts[first] + counts[second]
            pair_purity = pair_depth / depth
            if pair_purity < config.pair_purity:
                counters["LOW_PAIR_PURITY"] += 1
                continue
            first_frequency = counts[first] / pair_depth
            second_frequency = counts[second] / pair_depth
            if not (
                config.allele_min <= first_frequency <= config.allele_max
                and config.allele_min <= second_frequency <= config.allele_max
            ):
                counters["NOT_BIALLELIC"] += 1
                continue
            if counts[second] < config.min_allele_reads:
                counters["LOW_ALLELE_SUPPORT"] += 1
                continue
            row = fields + [
                str(pair_depth), f"{pair_purity:.6f}", str(depth - pair_depth),
                BASES[first], str(counts[first]), f"{first_frequency:.6f}",
                BASES[second], str(counts[second]), f"{second_frequency:.6f}",
                f"{min(first_frequency, second_frequency):.6f}",
            ]
            retained.write("\t".join(row) + "\n")
            counters["PASS"] += 1
    temporary.replace(output)
    return {
        "retained_sites": counters["PASS"],
        "low_depth": counters["LOW_DEPTH"],
        "high_depth": counters["HIGH_DEPTH"],
        "not_biallelic": counters["NOT_BIALLELIC"],
        "multiallelic": counters["MULTIALLELIC"],
        "low_pair_purity": counters["LOW_PAIR_PURITY"],
        "low_allele_support": counters["LOW_ALLELE_SUPPORT"],
    }


def weighted_quantile(histogram: dict[int, int] | Counter[int], probability: float) -> int:
    if not histogram or not 0 <= probability <= 1:
        raise ValueError("invalid histogram or quantile probability")
    target = max(1, math.ceil(sum(histogram.values()) * probability))
    cumulative = 0
    for value in sorted(histogram):
        cumulative += histogram[value]
        if cumulative >= target:
            return int(value)
    return int(max(histogram))


def write_depth_distribution(output: Path, histogram: dict[int, int]) -> None:
    total = sum(histogram.values())
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("depth\tsite_count\tfraction\n")
        for depth in sorted(histogram):
            count = histogram[depth]
            handle.write(f"{depth}\t{count}\t{count / total if total else 0:.8f}\n")


def write_thresholds(output: Path, thresholds: dict[str, object]) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in thresholds.items():
            handle.write(f"{key}\t{value}\n")


def write_filter_summary(
    output: Path, extraction: dict[str, object], filtering: dict[str, int]
) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tcount\n")
        handle.write(f"pileup_sites\t{extraction['pileup_sites']}\n")
        for key, value in filtering.items():
            handle.write(f"{key}\t{value}\n")


def _sample_from_bam(bam: Path) -> str:
    return bam.name[:-4] if bam.name.lower().endswith(".bam") else bam.name


def _find_bam_index(bam: Path) -> Path | None:
    candidates = (Path(f"{bam}.bai"), bam.with_suffix(".bai"), Path(f"{bam}.csi"))
    return next((path for path in candidates if path.is_file()), None)


def _file_stamp(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _signature(
    config: BamCountConfig, sample: str, bam: Path, reference: Path, mask: Path
) -> dict[str, object]:
    return {
        "sample": sample,
        "bam": _file_stamp(bam),
        "reference": _file_stamp(reference),
        "mask": _file_stamp(mask),
        "min_mapq": config.min_mapq,
        "min_baseq": config.min_baseq,
        "min_depth": config.min_depth,
        "max_depth": config.max_depth,
        "allele_min": config.allele_min,
        "allele_max": config.allele_max,
        "pair_purity": config.pair_purity,
        "min_allele_reads": config.min_allele_reads,
        "max_depth_cap": config.max_depth_cap,
        "region": config.region,
        "backend": config.samtools_executable,
        "mode": config.mode,
        "mapploid_version": __version__,
        "filter_algorithm": "strict-pair-depth-mac-v4-configurable",
    }


def _check_or_write_signature(path: Path, signature: dict[str, object]) -> None:
    if path.is_file():
        previous = json.loads(path.read_text(encoding="utf-8"))
        if previous != signature:
            raise ValueError(f"output directory has different inputs or parameters: {path.parent}")
        return
    path.write_text(json.dumps(signature, indent=2) + "\n", encoding="utf-8")


def _output_paths(output_dir: Path, sample: str) -> dict[str, Path]:
    return {
        "depth_distribution": output_dir / f"{sample}.depth_distribution.tsv",
        "thresholds": output_dir / f"{sample}.thresholds.tsv",
        "biallelic_sites": output_dir / f"{sample}.biallelic_sites.tsv.gz",
        "filter_summary": output_dir / f"{sample}.filter_summary.tsv",
        "manifest": output_dir / "run_manifest.json",
    }
