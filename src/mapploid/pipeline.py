"""Reference preparation and a single, bounded region-to-ploidy workflow."""
from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from dataclasses import dataclass, field
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
from .bam_counts import BamCountConfig, build_bam_counts, validate_inputs, _find_bam_index
from .panel import PanelConfig, build_panels
from .ploidy_model import PloidyModelConfig, run_ploidy_model
from .reference_mask import ReferenceMaskConfig, build_reference_mask


@dataclass
class PrepareConfig:
    reference: Path
    output_dir: Path
    read_length: int = 150
    threads: int = 8
    exclude_beds: list[Path] = field(default_factory=list)
    genmap_bedgraph: Path | None = None
    genmap_executable: str = "genmap"


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
    """Build reusable read-length-aware regions, including small references."""
    if config.read_length < 3 or config.threads < 1:
        raise ValueError("read length must be at least 3; threads must be positive")
    reference = config.reference.resolve(strict=True)
    root = config.output_dir.resolve()
    signature = dict(workflow="prepare", version=__version__, reference=identity(reference),
                     reference_sha256=sha256(reference), read_length=config.read_length,
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
                read_length=config.read_length, mismatches=2, threads=config.threads,
                exclude_beds=config.exclude_beds, genmap_bedgraph=config.genmap_bedgraph,
                genmap_executable=config.genmap_executable,
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
        manifest = dict(schema=1, version=__version__, reference=signature["reference"],
                        reference_sha256=signature["reference_sha256"], read_length=config.read_length,
                        mismatches=2, high_mappability_bases=panel_manifest["source"]["source_bases"],
                        eligible_bases=panel_manifest["source"]["eligible_bases"],
                        stages=stages, sha256=hashes)
        atomic_json(root / "assets.json", manifest)
    print(f"[MapPloid] reference assets ready: {root}", file=sys.stderr)
    return root / "assets.json"


def load_assets(root):
    manifest = read_json(root / "assets.json")
    if manifest.get("schema") != 1 or manifest.get("version") != __version__:
        raise ValueError("reference assets were prepared by a different MapPloid version")
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


def check_read_length(bam, read_length):
    import pysam
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        seen = 0
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            length = read.query_length
            if length and length < read_length:
                raise ValueError(
                    f"BAM contains reads of length {length}, shorter than asset read length {read_length}; "
                    "prepare assets using a read length compatible with the trimmed reads"
                )
            seen += 1
            if seen >= 10_000:
                break


def empty_evidence(reason):
    return dict(predicted_ploidy="UNCERTAIN", selected_sites=0, available_sites=0,
                bootstrap_support="NA", informative_chromosomes=0, uncertainty_reasons=reason,
                warnings="NONE")


def call(config: CallConfig):
    """Count and fit 20/50/100 Mb of qualifying regions, stopping when supported."""
    config.validate()
    root = config.output_dir.resolve()
    assets = config.assets.resolve(strict=True)
    bam = config.bam.resolve(strict=True)
    reference = assets / "reference.fa"
    manifest = load_assets(assets)
    sample = config.sample or bam.stem
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
        raise ValueError("sample name must contain only letters, digits, '.', '_' or '-'")
    index = _find_bam_index(bam)
    if index is None:
        raise ValueError(f"BAM index is missing: {bam}")
    first_bed = assets / manifest["stages"][0]["bed"]
    _, input_warnings = validate_inputs(bam, reference, first_bed, None, allow_empty=True)
    check_read_length(bam, manifest["read_length"])
    signature = dict(workflow="call", version=__version__, sample=sample, bam=identity(bam),
                     bam_index=dict(identity(index), sha256=sha256(index)), assets=str(assets),
                     assets_sha256=sha256(assets / "assets.json"), samtools=config.samtools_executable,
                     filter_schema=1, filters=config.count_filters())
    with output_lock(root):
        check_signature(root, signature)
        complete = root / "details" / "complete.json"
        if complete.exists():
            verify_checksums(root, read_json(complete)["sha256"])
            return root / "result.tsv"
        history = []
        start = time.monotonic()
        atomic_json(root / "details" / "status.json", dict(state="RUNNING", sample=sample))
        try:
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
                    evidence = saved["evidence"]
                    elapsed = saved["walltime_s"]
                else:
                    files = []
                    with (directory / "run.log").open("a") as log, redirect_stdout(log), redirect_stderr(log):
                        if stage["bases"] == 0:
                            evidence = empty_evidence("NO_ELIGIBLE_REGIONS")
                        else:
                            counts = build_bam_counts(BamCountConfig(
                                bam=bam, reference=reference, mask=bed, output_dir=directory / "counts",
                                sample=sample, samtools_executable=config.samtools_executable,
                                **config.count_filters(),
                            ))
                            files.extend(counts.values())
                            info = read_json(counts["manifest"])
                            if not info["filtering"]["retained_sites"]:
                                evidence = empty_evidence("NO_COVERAGE" if not info["extraction"]["acgt_observations"]
                                                          else "NO_VALID_SITES")
                            else:
                                model = run_ploidy_model(PloidyModelConfig(
                                    counts=counts["biallelic_sites"], output_dir=directory / "model", sample=sample,
                                    **config.model_filters(),
                                ))
                                evidence = read_json(model["manifest"])["result"]
                                files.extend(model.values())
                    elapsed = time.monotonic() - started
                    atomic_json(checkpoint, dict(evidence=evidence, walltime_s=elapsed,
                                                sha256=output_checksums(root, files)))
                history.append(dict(stage=number, region_bases=stage["bases"],
                                    predicted_ploidy=evidence["predicted_ploidy"],
                                    selected_sites=evidence["selected_sites"],
                                    uncertainty_reasons=evidence["uncertainty_reasons"], walltime_s=round(elapsed, 3)))
                write_table(root / "details" / "stages.tsv", history)
                if evidence["predicted_ploidy"] != "UNCERTAIN":
                    break
            called = evidence["predicted_ploidy"] != "UNCERTAIN"
            row = dict(sample=sample, predicted_ploidy=evidence["predicted_ploidy"],
                       status="CALLED" if called else "UNCERTAIN",
                       mappable_region_mb=manifest["high_mappability_bases"] / 1e6,
                       eligible_region_mb=manifest["eligible_bases"] / 1e6,
                       analyzed_region_mb=stage["bases"] / 1e6, selected_sites=evidence["selected_sites"],
                       bootstrap_support=evidence["bootstrap_support"],
                       informative_chromosomes=evidence["informative_chromosomes"],
                       reason=evidence["uncertainty_reasons"], warnings=evidence.get("warnings", "NONE"))
            write_table(root / "result.tsv", [row])
            temporary = root / "regions.bed.tmp"
            shutil.copyfile(bed, temporary)
            temporary.replace(root / "regions.bed")
            atomic_json(complete, dict(version=__version__, sample=sample,
                                      completed_at=datetime.now(timezone.utc).isoformat(),
                                      stop_reason="EVIDENCE_PASSED" if called else "REGION_BUDGET_EXHAUSTED",
                                      input_warnings=input_warnings, walltime_s=time.monotonic() - start,
                                      sha256=output_checksums(root, [root / "result.tsv", root / "regions.bed"])))
            atomic_json(root / "details" / "status.json", dict(state="COMPLETE", result=row))
        except BaseException as exc:
            atomic_json(root / "details" / "status.json", dict(state="FAILED", error=str(exc)))
            (root / "details" / "error.log").write_text(traceback.format_exc())
            raise
    print(f"[MapPloid] {sample}: {row['predicted_ploidy']}; result: {root / 'result.tsv'}", file=sys.stderr)
    return root / "result.tsv"
