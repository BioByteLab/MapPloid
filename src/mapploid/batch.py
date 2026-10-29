"""Sequential BAM-list calls with independent sample outputs and live summary."""
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime, timezone
import csv
import json
import re
import sys

from . import __version__
from .pipeline import (CallConfig, atomic_json, call, check_signature, identity,
                       load_assets, output_lock, sha256, write_table)
from .ploidy_model import DEFAULT_MAX_SITES, DEFAULT_MIN_SITES, PloidyModelConfig

from .reporting import RESULT_FIELDS
FIELDS = (*RESULT_FIELDS, 'result_file')


@dataclass(frozen=True)
class BatchResult:
    summary: Path
    failed: int


def read_bam_list(path):
    """Whole lines are paths; preserve embedded spaces and resolve relative to list."""
    path=Path(path).resolve(strict=True)
    entries=[];seen_paths={};seen_names={}
    with path.open(encoding='utf-8-sig') as handle:
        for number,line in enumerate(handle,1):
            value=line.strip()
            if not value or value.startswith('#'):
                continue
            bam=Path(value).expanduser()
            if not bam.is_absolute():bam=path.parent/bam
            bam=bam.resolve()
            if bam.suffix.lower()!='.bam':
                raise ValueError(f'{path}, line {number}: expected one BAM file path: {value}')
            sample=bam.stem
            if not re.fullmatch(r'[A-Za-z0-9_.-]+',sample) or sample in ('.', '..'):
                raise ValueError(f'{path}, line {number}: invalid sample name from BAM filename: {sample}; use letters, digits, dots, underscores or hyphens')
            if bam in seen_paths:
                raise ValueError(f'{path}, line {number}: duplicate BAM path, first listed on line {seen_paths[bam]}: {bam}')
            if sample in seen_names:
                raise ValueError(f'{path}, line {number}: duplicate sample name: {sample}; use distinct BAM filenames (already used on line {seen_names[sample]})')
            seen_paths[bam]=number;seen_names[sample]=number
            entries.append(dict(sample=sample,bam=bam))
    if not entries:raise ValueError(f'BAM list contains no paths: {path}; provide one BAM file path per line')
    return entries


def call_batch(bam_list, assets, output_dir, **filters):
    """Continue after sample errors, with live compact tables and offline reports."""
    from .pipeline import run_request, compact_row, verify_checksums, output_checksums
    from .reporting import OUTPUT_SCHEMA, clean_json, write_batch_report
    import traceback
    template = CallConfig(Path("unused.bam"), Path(assets), Path(output_dir), **filters)
    template.validate()
    bam_list = Path(bam_list).resolve(strict=True)
    entries = read_bam_list(bam_list)
    assets = Path(assets).resolve(strict=True)
    load_assets(assets)
    root = Path(output_dir).resolve()
    request = dict(workflow="batch-call", version=__version__, output_schema=OUTPUT_SCHEMA,
                   bam_list=dict(identity(bam_list), sha256=sha256(bam_list)),
                   assets=str(assets), assets_sha256=sha256(assets / "assets.json"),
                   filters=template.count_filters(), samtools=template.samtools_executable,
                   model_site_limits=dict(min_sites=DEFAULT_MIN_SITES, max_sites=DEFAULT_MAX_SITES),
                   output_options=dict(keep_details=template.keep_details, export_figures=template.export_figures),
                   samples=[dict(sample=x["sample"], bam=str(x["bam"])) for x in entries])
    if template.export_figures:
        request["output_options"]["figure_formats"] = ["pdf", "png"]
    summary, rows, failed = root / "result.tsv", [], 0
    with output_lock(root):
        run_request(root, request)
        complete = root / ".mapploid.complete.json"
        if complete.exists():
            verify_checksums(root, json.loads(complete.read_text())["sha256"])
            complete.unlink()
        record = dict(version=__version__, output_schema=OUTPUT_SCHEMA, request=request,
                      state="RUNNING", total=len(entries), completed=0, failed=0)
        atomic_json(root / "details/run.json", record)
        log_path = root / "details/run.log"

        def update(state):
            record.update(state=state, completed=len(rows), failed=failed, results=rows)
            atomic_json(root / "details/run.json", clean_json(record))
            if rows:
                write_table(summary, [{k: row[k] for k in FIELDS} for row in rows])
            write_batch_report(root, rows, len(entries), state, __version__)

        update("RUNNING")
        try:
            for number, entry in enumerate(entries, 1):
                sample, bam = entry["sample"], entry["bam"]
                message = f"[{number}/{len(entries)}] Analyzing {sample}"
                print(f"[MapPloid] {message}", file=sys.stderr)
                with log_path.open("a") as log:
                    log.write(message + "\n")
                try:
                    result = call(CallConfig(bam, assets, root / "samples" / sample, sample, **filters))
                    with result.open() as handle:
                        sample_row = next(csv.DictReader(handle, delimiter="\t"))
                    if sample_row.get("sample") != sample or sample_row.get("status") not in ("CALLED", "UNCERTAIN"):
                        raise ValueError(f"Invalid sample result format: {result}")
                    row = {field: sample_row[field] for field in RESULT_FIELDS}
                    row["result_file"] = str(result.relative_to(root))
                except (OSError, ValueError, RuntimeError) as exc:
                    failed += 1
                    row = compact_row(sample, error=exc)
                    result = root / "samples" / sample / "result.tsv"
                    row["result_file"] = str(result.relative_to(root)) if result.exists() else "NA"
                    with log_path.open("a") as log:
                        log.write(traceback.format_exc())
                    print(f"[MapPloid] {sample} failed: {exc}", file=sys.stderr)
                rows.append(row)
                update("RUNNING")
            state = "COMPLETE_WITH_ERRORS" if failed else "COMPLETE"
            record["completed_at"] = datetime.now(timezone.utc).isoformat()
            update(state)
            # Only all-success batches get a final reuse checksum marker; failed samples are retryable.
            if not failed:
                artifacts = [summary, root / "report.html", root / "details/run.json", log_path]
                hashes = output_checksums(root, artifacts)
                verify_checksums(root, hashes)
                atomic_json(complete, dict(version=__version__, output_schema=OUTPUT_SCHEMA, sha256=hashes))
        except BaseException as exc:
            record["error"] = str(exc)
            update("INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "FAILED")
            with log_path.open("a") as log:
                log.write(traceback.format_exc())
            raise
    print(f"[MapPloid] Batch complete: {len(entries)} samples, {failed} failed; report: {root / 'report.html'}", file=sys.stderr)
    return BatchResult(summary, failed)
