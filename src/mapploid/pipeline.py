"""Reference preparation and a single, bounded region-to-ploidy workflow."""
from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import fcntl
import hashlib
import json
import os
import re
import shutil
import sys
import time
import traceback
import uuid

from . import __version__
from .reporting import OUTPUT_SCHEMA, RESULT_FIELDS, clean_json, finite, write_sample_report
from .bam_counts import BamCountConfig, build_bam_counts, validate_inputs, _find_bam_index
from .panel import PanelConfig, build_panels
from .ploidy_model import DEFAULT_MAX_SITES, DEFAULT_MIN_SITES, PloidyModelConfig, run_ploidy_model
from .reference_mask import ReferenceMaskConfig, build_reference_mask


@dataclass
class PrepareConfig:
    reference: Path
    output_dir: Path
    kmer_length: int = 150
    threads: int = 8
    exclude_beds: list[Path] = field(default_factory=list)
    genmap_bedgraph: Path | None = None
    genmap_executable: str = "genmap"
    min_region_length: int = 10_000


@dataclass
class CallConfig:
    bam: Path
    assets: Path
    output_dir: Path
    sample: str | None = None
    samtools_executable: str = "samtools"
    allele_min: float = 0.15
    allele_max: float = 0.85
    min_depth: int | None = None
    max_depth: int | None = None
    min_mapq: int = 20
    min_baseq: int = 20
    min_allele_reads: int = 2
    pair_purity: float = 0.95

    keep_details: bool = False
    export_figures: bool = False

    def count_filters(self):
        return {name: getattr(self, name) for name in (
            "allele_min", "allele_max", "min_depth", "max_depth",
            "min_mapq", "min_baseq", "min_allele_reads", "pair_purity",
        )}

    def model_filters(self):
        return {name: getattr(self, name) for name in (
            "allele_min", "allele_max", "min_allele_reads", "pair_purity",
        )}

    def validate(self):
        BamCountConfig(self.bam, self.assets / "reference.fa",
                       self.assets / "mappable_regions.bed", self.output_dir,
                       **self.count_filters()).validate()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def identity(path):
    path = Path(path).resolve(strict=True)
    info = path.stat()
    return dict(path=str(path), size=info.st_size, mtime_ns=info.st_mtime_ns)


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temp.replace(path)


def write_table(path, rows):
    path = Path(path)
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


@contextmanager
def output_lock(root):
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".mapploid.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"another MapPloid process is using {root}") from exc
        yield


def check_signature(root, signature):
    path = root / "details" / "request.json"
    if path.exists():
        if read_json(path) != signature:
            raise ValueError(f"inputs or settings changed; use a new output directory: {root}")
    else:
        if any(p.name != ".mapploid.lock" for p in root.iterdir()):
            raise ValueError(f"output directory is not an empty MapPloid run: {root}")
        atomic_json(path, signature)


def output_checksums(root, files):
    return {str(Path(path).relative_to(root)): sha256(path) for path in files}


def checked_path(root, relative):
    # Assets and checkpoints never direct writes/reads outside their own directory.
    path = root / relative
    if Path(relative).is_absolute() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"invalid asset/checkpoint path: {relative}")
    return path


def verify_checksums(root, hashes):
    for relative, expected in hashes.items():
        path = checked_path(root, relative)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"output is missing or modified: {path}")


def prepare(config: PrepareConfig):
    """Build reusable k-mer-aware regions, including small references."""
    if config.kmer_length < 3 or config.threads < 1:
        raise ValueError("k-mer length must be at least 3; threads must be positive")
    if config.min_region_length < 1:
        raise ValueError("minimum continuous region length must be positive")
    reference = config.reference.resolve(strict=True)
    root = config.output_dir.resolve()
    signature = dict(workflow="prepare", version=__version__, reference=identity(reference),
                     reference_sha256=sha256(reference), kmer_length=config.kmer_length,
                     min_region_length=config.min_region_length,
                     mismatches=2, genmap_executable=config.genmap_executable,
                     exclude_beds=[dict(identity(p), sha256=sha256(p)) for p in config.exclude_beds],
                     genmap_bedgraph=(dict(identity(config.genmap_bedgraph), sha256=sha256(config.genmap_bedgraph))
                                      if config.genmap_bedgraph else None))
    with output_lock(root):
        check_signature(root, signature)
        if (root / "assets.json").exists():
            load_assets(root)
            return root / "assets.json"
        print("[MapPloid] preparing high-mappability reference regions", file=sys.stderr)
        with (root / "details" / "prepare.log").open("a") as log, redirect_stdout(log), redirect_stderr(log):
            mask_outputs = build_reference_mask(ReferenceMaskConfig(
                reference=reference, output_dir=root / "details" / "mask",
                kmer_length=config.kmer_length, mismatches=2, threads=config.threads,
                exclude_beds=config.exclude_beds, genmap_bedgraph=config.genmap_bedgraph,
                genmap_executable=config.genmap_executable,
                min_region_length=config.min_region_length,
            ))
            panel_outputs = build_panels(PanelConfig(
                mask=mask_outputs["mappability_mask"],
                genome_sizes=root / "details" / "mask" / "genome.sizes",
                output_dir=root / "details" / "panels",
            ))
            reference_link = root / "reference.fa"
            if not reference_link.is_symlink():
                if reference_link.exists():
                    raise ValueError(f"unexpected reference file: {reference_link}")
                reference_link.symlink_to(reference)
            elif reference_link.resolve() != reference:
                raise ValueError("reference link changed")
            import pysam
            # Index the asset-local link; never create/replace an index beside the source FASTA.
            pysam.faidx(str(reference_link))
        panel_manifest = read_json(panel_outputs["manifest"])
        stages = []
        for name, info in panel_manifest["panels"].items():
            bed = Path(info["path"])
            bases = sum(int(line.split()[2]) - int(line.split()[1])
                        for line in bed.read_text().splitlines() if line.strip())
            if not stages or bases != stages[-1]["bases"]:
                stages.append(dict(name=name, bed=str(bed.relative_to(root)), bases=bases))
        good = mask_outputs["mappability_mask"]
        shutil.copyfile(good, root / "mappable_regions.bed")
        hashes = output_checksums(root, [root / "mappable_regions.bed", root / "reference.fa.fai",
                                        *[root / row["bed"] for row in stages]])
        manifest = dict(schema=2, version=__version__, reference=signature["reference"],
                        reference_sha256=signature["reference_sha256"], kmer_length=config.kmer_length,
                        min_region_length=config.min_region_length,
                        mismatches=2, high_mappability_bases=panel_manifest["source"]["source_bases"],
                        eligible_bases=panel_manifest["source"]["eligible_bases"],
                        stages=stages, sha256=hashes)
        atomic_json(root / "assets.json", manifest)
    print(f"[MapPloid] Reference directory ready: {root}; use --assets {root} for call", file=sys.stderr)
    return root / "assets.json"


def load_assets(root):
    if not (root / "assets.json").is_file():
        raise ValueError(
            f"--assets must point to the reference directory created by prepare --out: {root}; "
            "first run mapploid prepare --reference genome.fa --out reference_assets"
        )
    manifest = read_json(root / "assets.json")
    if manifest.get("schema") not in (1, 2) or manifest.get("version") not in ("0.10.1", "0.10.2", "0.10.3", __version__):
        raise ValueError("reference assets were prepared by a different MapPloid version")
    if manifest["schema"] == 1:
        # Existing assets remain readable, without reinterpreting them as 10kb-filtered.
        manifest.setdefault("kmer_length", manifest.get("read_length"))
        manifest.setdefault("min_region_length", 300)
    if not isinstance(manifest.get("kmer_length"), int) or manifest["kmer_length"] < 3:
        raise ValueError("reference assets contain an invalid k-mer length")
    if not isinstance(manifest.get("min_region_length"), int) or manifest["min_region_length"] < 1:
        raise ValueError("reference assets contain an invalid minimum continuous region length")
    if identity(root / "reference.fa") != manifest["reference"]:
        raise ValueError("reference FASTA changed; prepare a new asset directory")
    verify_checksums(root, manifest["sha256"])
    previous = -1
    for stage in manifest["stages"]:
        bed = checked_path(root, stage["bed"])
        if stage["bed"] not in manifest["sha256"] or stage["bases"] <= previous:
            raise ValueError("invalid or repeated analysis region stage")
        previous = stage["bases"]
    if not manifest["stages"]:
        raise ValueError("reference assets contain no analysis stage record")
    return manifest


def check_kmer_length(bam, kmer_length):
    import pysam
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        seen = 0
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            length = read.query_length
            if length and length < kmer_length:
                raise ValueError(
                    f"BAM contains reads of length {length}, shorter than asset k-mer length {kmer_length}; "
                    "prepare assets using a k-mer length compatible with the trimmed reads"
                )
            seen += 1
            if seen >= 10_000:
                break


def empty_evidence(reason):
    return dict(predicted_ploidy="UNCERTAIN", selected_sites=0, available_sites=0,
                bootstrap_support="NA", informative_chromosomes=0, uncertainty_reasons=reason,
                warnings="NONE")


def run_request(root, request):
    """Check ownership without overwriting a completed or incompatible run."""
    record_path = root / "details" / "run.json"
    if record_path.exists():
        previous = read_json(record_path)
        if previous.get("output_schema") != OUTPUT_SCHEMA or previous.get("request") != request:
            raise ValueError(f"inputs or settings changed; use a new output directory: {root}")
        return previous
    if any(p.name != ".mapploid.lock" for p in root.iterdir()):
        raise ValueError(f"output directory is not an empty MapPloid run: {root}")
    return {}


def cleanup_stages(root, keep_details):
    if keep_details:
        return
    for p in (root / "details").glob("stage_*"):
        if p.is_dir() and re.fullmatch(r"stage_[0-9]+", p.name) and not p.is_symlink():
            try:
                shutil.rmtree(p)
            except OSError as exc:
                print(f"[MapPloid] Result is complete; could not clean {p}: {exc}", file=sys.stderr)


def compact_row(sample, evidence=None, stage=None, error=None):
    evidence = evidence or {}
    if error is not None:
        return dict(sample=sample, predicted_ploidy="NA", status="FAILED",
                    selected_sites="NA", delta_bic="NA", analyzed_region_mb="NA",
                    reason=str(error).replace("\n", " ").replace("\t", " "), warnings="NONE")
    predicted = evidence["predicted_ploidy"]
    delta = evidence.get("delta_bic")
    return dict(sample=sample, predicted_ploidy=predicted,
                status="UNCERTAIN" if predicted == "UNCERTAIN" else "CALLED",
                selected_sites=evidence["selected_sites"],
                delta_bic=delta if finite(delta) else "NA",
                analyzed_region_mb=stage["bases"] / 1e6,
                reason=evidence["uncertainty_reasons"], warnings=evidence.get("warnings", "NONE"))


def call(config: CallConfig):
    """Call with compact durable artifacts, resumable stage checkpoints and offline reports."""
    config.validate()
    root = config.output_dir.resolve()
    bam, assets = config.bam.resolve(), config.assets.resolve()
    sample = config.sample or bam.stem
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample) or sample in (".", ".."):
        raise ValueError("sample name must contain only letters, digits, '.', '_' or '-'")
    request = dict(workflow="call", version=__version__, output_schema=OUTPUT_SCHEMA,
                   sample=sample, bam=str(bam), assets=str(assets),
                   samtools=config.samtools_executable, filters=config.count_filters(),
                   model_site_limits=dict(min_sites=DEFAULT_MIN_SITES, max_sites=DEFAULT_MAX_SITES),
                   output_options=dict(keep_details=config.keep_details, export_figures=config.export_figures))
    if config.export_figures:
        request["output_options"]["figure_formats"] = ["pdf", "png"]
    with output_lock(root):
        previous = run_request(root, request)
        complete = root / ".mapploid.complete.json"
        # Verify final artifacts before touching them, including the run record and log.
        if complete.exists():
            verify_checksums(root, read_json(complete)["sha256"])
        record = dict(previous)
        history = []
        start = time.monotonic()
        model_config = PloidyModelConfig(Path("unused"), Path("unused"), sample=sample,
                                        min_sites=DEFAULT_MIN_SITES, max_sites=DEFAULT_MAX_SITES,
                                        **config.model_filters())
        model_settings = {k: v for k, v in asdict(model_config).items()
                          if k not in ("counts", "output_dir", "sample")}
        # A mismatch is a reuse error, not a new FAILED result.
        if previous.get("input_identity"):
            old = previous["input_identity"]
            try:
                index = _find_bam_index(bam)
                current = dict(bam=identity(bam), bam_index=dict(identity(index), sha256=sha256(index)),
                               assets=dict(identity(assets / "assets.json"), sha256=sha256(assets / "assets.json")))
            except (OSError, TypeError) as exc:
                raise ValueError("inputs or settings changed; use a new output directory") from exc
            if current != old:
                raise ValueError("inputs or settings changed; use a new output directory")
        if complete.exists():
            load_assets(assets)
            cleanup_stages(root, config.keep_details)
            return root / "result.tsv"
        record.update(version=__version__, output_schema=OUTPUT_SCHEMA, request=request,
                      state="RUNNING", model_settings=model_settings,
                      cache_policy="retain" if config.keep_details else "clean after verified completion",
                      stages=[], started_at=previous.get("started_at", datetime.now(timezone.utc).isoformat()))
        record.pop("result", None)
        record.pop("error", None)
        atomic_json(root / "details/run.json", clean_json(record))
        log_path = root / "details/run.log"
        try:
            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"MapPloid {__version__}: starting/resuming {sample}\n")
                log.flush()
            manifest = load_assets(assets)
            reference = assets / "reference.fa"
            index = _find_bam_index(bam)
            if index is None:
                raise ValueError(f"BAM index is missing: {bam}")
            record["input_identity"] = dict(
                bam=identity(bam), bam_index=dict(identity(index), sha256=sha256(index)),
                assets=dict(identity(assets / "assets.json"), sha256=sha256(assets / "assets.json")))
            first_bed = assets / manifest["stages"][0]["bed"]
            _, input_warnings = validate_inputs(bam, reference, first_bed, None, allow_empty=True)
            check_kmer_length(bam, manifest["kmer_length"])
            record["input_warnings"] = input_warnings
            record["reference"] = {k: v for k, v in manifest.items() if k not in ("stages", "sha256")}
            atomic_json(root / "details/run.json", clean_json(record))
            for number, stage in enumerate(manifest["stages"], 1):
                directory = root / "details" / f"stage_{number}"
                directory.mkdir(exist_ok=True)
                bed = assets / stage["bed"]
                checkpoint = directory / "complete.json"
                started = time.monotonic()
                print(f"[MapPloid] {sample}: analyzing {stage['bases'] / 1e6:.3f} Mb", file=sys.stderr)
                if checkpoint.exists():
                    saved = read_json(checkpoint)
                    verify_checksums(root, saved["sha256"])
                    summary = saved["summary"]
                    evidence = summary["evidence"]
                else:
                    files, counts_info, model_info, chromosomes = [], {}, {}, []
                    with log_path.open("a", encoding="utf-8") as log, redirect_stdout(log), redirect_stderr(log):
                        if stage["bases"] == 0:
                            evidence = empty_evidence("NO_ELIGIBLE_REGIONS")
                        else:
                            counts = build_bam_counts(BamCountConfig(
                                bam=bam, reference=reference, mask=bed, output_dir=directory / "counts",
                                sample=sample, samtools_executable=config.samtools_executable,
                                **config.count_filters()))
                            files.extend(counts.values())
                            counts_info = read_json(counts["manifest"])
                            if not counts_info["filtering"]["retained_sites"]:
                                evidence = empty_evidence("NO_COVERAGE" if not counts_info["extraction"]["acgt_observations"]
                                                          else "NO_VALID_SITES")
                            else:
                                model = run_ploidy_model(PloidyModelConfig(
                                    counts=counts["biallelic_sites"], output_dir=directory / "model", sample=sample,
                                    min_sites=DEFAULT_MIN_SITES, max_sites=DEFAULT_MAX_SITES,
                                    **config.model_filters()))
                                model_info = read_json(model["manifest"])
                                evidence = model_info["result"]
                                files.extend(model.values())
                                if "chromosomes" in model:
                                    with Path(model["chromosomes"]).open() as handle:
                                        chromosomes = list(csv.DictReader(handle, delimiter="\t"))
                    summary = clean_json(dict(
                        stage=number, region_bases=stage["bases"], evidence=evidence,
                        walltime_s=time.monotonic() - started,
                        counts={k: v for k, v in counts_info.items() if k != "outputs"},
                        model={k: v for k, v in model_info.items() if k != "outputs"},
                        chromosomes=chromosomes))
                    atomic_json(checkpoint, dict(summary=summary, sha256=output_checksums(root, files)))
                history.append(summary)
                record["stages"] = history
                atomic_json(root / "details/run.json", clean_json(record))
                if evidence["predicted_ploidy"] != "UNCERTAIN":
                    break
            row = compact_row(sample, evidence, stage)
            record.update(state="COMPLETE", result=row,
                          completed_at=datetime.now(timezone.utc).isoformat(),
                          walltime_s=sum(s["walltime_s"] for s in history),
                          current_invocation_walltime_s=time.monotonic() - start)
            write_table(root / "result.tsv", [{k: row[k] for k in RESULT_FIELDS}])
            temporary = root / "regions.bed.tmp"
            shutil.copyfile(bed, temporary)
            temporary.replace(root / "regions.bed")
            artifacts = write_sample_report(root, record, config.export_figures)
            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"Complete: {row['predicted_ploidy']}; {row['status']}\n")
            atomic_json(root / "details/run.json", clean_json(record))
            artifacts += [root / "result.tsv", root / "regions.bed", root / "details/run.json", log_path]
            hashes = output_checksums(root, artifacts)
            verify_checksums(root, hashes)
            atomic_json(complete, dict(version=__version__, output_schema=OUTPUT_SCHEMA, sha256=hashes))
        except BaseException as exc:
            # Keep checkpoints and logs. Only a verified COMPLETE marker authorizes cleanup.
            record.update(state="INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "FAILED",
                          stages=history, error=str(exc),
                          result=compact_row(sample, error="INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else exc))
            with log_path.open("a", encoding="utf-8") as log:
                log.write(traceback.format_exc())
            atomic_json(root / "details/run.json", clean_json(record))
            write_table(root / "result.tsv", [{k: record["result"][k] for k in RESULT_FIELDS}])
            try:
                write_sample_report(root, record)
            except Exception as report_error:
                with log_path.open("a", encoding="utf-8") as log:
                    log.write(f"Could not render failure report: {report_error}\n")
            raise
        cleanup_stages(root, config.keep_details)
    print(f"[MapPloid] {sample}: {row['predicted_ploidy']}; report: {root / 'report.html'}", file=sys.stderr)
    return root / "result.tsv"
