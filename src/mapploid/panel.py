"""Build deterministic, nested analysis panels from a master mappability mask."""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import json

from . import __version__
from .intervals import Interval, iter_bed


@dataclass(slots=True)
class PanelConfig:
    mask: Path
    genome_sizes: Path
    output_dir: Path
    target_mb: tuple[int, ...] = (20, 50, 100)
    min_run: int = 300
    block_size: int = 5_000_000
    chunk_size: int = 10_000

    def validate(self) -> None:
        if not self.target_mb or any(value < 1 for value in self.target_mb):
            raise ValueError("panel targets must be positive megabase values")
        if len(set(self.target_mb)) != len(self.target_mb):
            raise ValueError("panel targets must be unique")
        if self.min_run < 1:
            raise ValueError("minimum panel run must be positive")
        if self.block_size < 1:
            raise ValueError("panel block size must be positive")
        if self.chunk_size < self.min_run:
            raise ValueError("panel chunk size must be at least min-run")


def build_panels(config: PanelConfig) -> dict[str, Path]:
    """Create spatially balanced panels whose selected bases are nested."""

    config.validate()
    mask = config.mask.resolve()
    genome_sizes_path = config.genome_sizes.resolve()
    output_dir = config.output_dir.resolve()
    for path, label in ((mask, "mask"), (genome_sizes_path, "genome sizes")):
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")

    genome_sizes = read_genome_sizes(genome_sizes_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = tuple(sorted(config.target_mb))
    signature = _signature(config, mask, genome_sizes_path, targets)
    config_path = output_dir / "run_config.json"
    _check_or_write_signature(config_path, signature)
    outputs = {
        f"panel_{target}mb": output_dir / f"analysis_panel_{target}mb.bed"
        for target in targets
    }
    outputs.update({
        "statistics": output_dir / "panel_statistics.tsv",
        "blocks": output_dir / "panel_blocks.tsv",
        "manifest": output_dir / "run_manifest.json",
    })
    chunks_by_block, source_stats = prepare_chunks(
        mask, genome_sizes, config.min_run, config.block_size, config.chunk_size
    )
    ordered_chunks = balanced_chunk_order(chunks_by_block)
    available_bases = sum(interval.length for interval in ordered_chunks)
    # Keep only distinct attainable stages; short/empty masks are valid assets.
    distinct_targets = []
    last_bases = None
    for target in targets:
        actual = min(target * 1_000_000, available_bases)
        if actual != last_bases:
            distinct_targets.append(target)
            last_bases = actual
    targets = tuple(distinct_targets)
    outputs = {key: path for key, path in outputs.items()
               if not key.startswith("panel_") or key in {f"panel_{t}mb" for t in targets}}

    selected: list[Interval] = []
    cumulative = 0
    target_index = 0
    target_stats: list[dict[str, int | float]] = []
    while target_index < len(targets):
        target_bases = min(targets[target_index] * 1_000_000, available_bases)
        while cumulative < target_bases:
            interval = ordered_chunks[len(selected)]
            selected.append(interval)
            cumulative += interval.length
        merged = merge_intervals(sorted(selected, key=lambda row: (row.chrom, row.start, row.end)))
        panel_path = outputs[f"panel_{targets[target_index]}mb"]
        write_panel(panel_path, merged)
        target_stats.append(panel_statistics(
            targets[target_index], target_bases, merged, config.block_size
        ))
        target_index += 1

    write_statistics(outputs["statistics"], source_stats, target_stats)
    write_block_statistics(outputs["blocks"], targets, outputs, config.block_size)
    manifest = {
        "workflow": "panel-build",
        "mapploid_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "run_config": signature,
        "source": source_stats,
        "panels": {
            f"{target}mb": {
                "path": str(outputs[f"panel_{target}mb"]),
                **target_stats[index],
            }
            for index, target in enumerate(targets)
        },
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    outputs["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return outputs


def read_genome_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected contig and length")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid contig length") from exc
            if not fields[0] or length < 1 or fields[0] in sizes:
                raise ValueError(f"{path}:{line_number}: invalid or duplicate contig")
            sizes[fields[0]] = length
    if not sizes:
        raise ValueError(f"genome sizes file is empty: {path}")
    return sizes


def prepare_chunks(
    mask: Path,
    genome_sizes: dict[str, int],
    min_run: int,
    block_size: int,
    chunk_size: int,
) -> tuple[dict[tuple[str, int], list[Interval]], dict[str, int | float]]:
    """Filter short source runs and split eligible bases into balanced chunks."""

    chunks_by_block: dict[tuple[str, int], list[tuple[int, Interval]]] = defaultdict(list)
    source_intervals = source_bases = eligible_intervals = eligible_bases = 0
    previous: Interval | None = None
    for interval in iter_bed(mask):
        if interval.chrom not in genome_sizes:
            raise ValueError(f"mask contig is absent from genome sizes: {interval.chrom}")
        if interval.end > genome_sizes[interval.chrom]:
            raise ValueError(f"mask interval exceeds {interval.chrom} length")
        if previous and (interval.chrom < previous.chrom or (
            interval.chrom == previous.chrom and interval.start < previous.end
        )):
            raise ValueError("mask BED must be lexically sorted and non-overlapping")
        previous = interval
        source_intervals += 1
        source_bases += interval.length
        if interval.length < min_run:
            continue
        eligible_intervals += 1
        eligible_bases += interval.length
        cursor = interval.start
        while cursor < interval.end:
            block_number = cursor // block_size
            block_end = min(interval.end, (block_number + 1) * block_size)
            piece_length = block_end - cursor
            if piece_length >= min_run:
                for chunk in split_piece(
                    Interval(interval.chrom, cursor, block_end), min_run, chunk_size
                ):
                    chunks_by_block[(interval.chrom, block_number)].append(
                        (interval.length, chunk)
                    )
            cursor = block_end

    for block, rows in chunks_by_block.items():
        rows.sort(key=lambda row: (-row[0], row[1].start, row[1].end))
        chunks_by_block[block] = [row[1] for row in rows]
    return chunks_by_block, {
        "source_intervals": source_intervals,
        "source_bases": source_bases,
        "eligible_intervals": eligible_intervals,
        "eligible_bases": eligible_bases,
        "removed_short_intervals": source_intervals - eligible_intervals,
        "removed_short_bases": source_bases - eligible_bases,
        "min_run": min_run,
    }


def split_piece(piece: Interval, min_run: int, chunk_size: int) -> list[Interval]:
    lengths: list[int] = []
    remaining = piece.length
    while remaining > chunk_size:
        lengths.append(chunk_size)
        remaining -= chunk_size
    if remaining:
        if remaining < min_run and lengths:
            lengths[-1] += remaining
        elif remaining >= min_run:
            lengths.append(remaining)
    result: list[Interval] = []
    cursor = piece.start
    for length in lengths:
        result.append(Interval(piece.chrom, cursor, cursor + length))
        cursor += length
    return result


def balanced_chunk_order(
    chunks_by_block: dict[tuple[str, int], list[Interval]],
) -> list[Interval]:
    queues = {
        block: deque(chunks)
        for block, chunks in sorted(chunks_by_block.items())
        if chunks
    }
    ordered: list[Interval] = []
    while queues:
        exhausted: list[tuple[str, int]] = []
        for block in list(queues):
            queue = queues[block]
            ordered.append(queue.popleft())
            if not queue:
                exhausted.append(block)
        for block in exhausted:
            del queues[block]
    return ordered


def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    merged: list[Interval] = []
    for interval in intervals:
        if merged and interval.chrom == merged[-1].chrom and interval.start <= merged[-1].end:
            previous = merged[-1]
            merged[-1] = Interval(previous.chrom, previous.start, max(previous.end, interval.end))
        else:
            merged.append(interval)
    return merged


def write_panel(path: Path, intervals: list[Interval]) -> None:
    temporary = Path(f"{path}.tmp")
    with temporary.open("wt", encoding="utf-8") as handle:
        for interval in intervals:
            handle.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")
    temporary.replace(path)


def panel_statistics(
    target_mb: int, target_bases: int, intervals: list[Interval], block_size: int
) -> dict[str, int | float]:
    bases = sum(interval.length for interval in intervals)
    blocks = {(row.chrom, row.start // block_size) for row in intervals}
    chromosomes = {row.chrom for row in intervals}
    return {
        "target_mb": target_mb,
        "target_bases": target_bases,
        "actual_bases": bases,
        "intervals": len(intervals),
        "chromosomes": len(chromosomes),
        "blocks": len(blocks),
    }


def write_statistics(
    path: Path,
    source: dict[str, int | float],
    panels: list[dict[str, int | float]],
) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        handle.write("panel\ttarget_bases\tactual_bases\tintervals\tchromosomes\tblocks\n")
        for row in panels:
            handle.write(
                f"{row['target_mb']}mb\t{row['target_bases']}\t{row['actual_bases']}\t"
                f"{row['intervals']}\t{row['chromosomes']}\t{row['blocks']}\n"
            )
        handle.write(
            f"source\t{source['source_bases']}\t{source['eligible_bases']}\t"
            f"{source['eligible_intervals']}\t.\t.\n"
        )


def write_block_statistics(
    path: Path,
    targets: tuple[int, ...],
    outputs: dict[str, Path],
    block_size: int,
) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        handle.write("panel\tchrom\tblock_start\tblock_end\tselected_bases\n")
        for target in targets:
            totals: dict[tuple[str, int], int] = defaultdict(int)
            for interval in iter_bed(outputs[f"panel_{target}mb"]):
                cursor = interval.start
                while cursor < interval.end:
                    block_number = cursor // block_size
                    end = min(interval.end, (block_number + 1) * block_size)
                    totals[(interval.chrom, block_number)] += end - cursor
                    cursor = end
            for (chrom, block_number), bases in sorted(totals.items()):
                start = block_number * block_size
                handle.write(
                    f"{target}mb\t{chrom}\t{start}\t{start + block_size}\t{bases}\n"
                )


def _file_stamp(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _signature(
    config: PanelConfig,
    mask: Path,
    genome_sizes: Path,
    targets: tuple[int, ...],
) -> dict[str, object]:
    return {
        "mask": _file_stamp(mask),
        "genome_sizes": _file_stamp(genome_sizes),
        "target_mb": list(targets),
        "min_run": config.min_run,
        "block_size": config.block_size,
        "chunk_size": config.chunk_size,
        "selection": "round-robin-longest-run-v1",
    }


def _check_or_write_signature(path: Path, signature: dict[str, object]) -> None:
    if path.is_file():
        previous = json.loads(path.read_text(encoding="utf-8"))
        if previous != signature:
            raise ValueError(f"output directory has different inputs or parameters: {path.parent}")
        return
    path.write_text(json.dumps(signature, indent=2) + "\n", encoding="utf-8")
