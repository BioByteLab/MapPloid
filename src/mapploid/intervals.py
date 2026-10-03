"""Streaming BED interval utilities used by the mask workflow.

All coordinates are zero-based, half-open.  Large files are sorted with the
system ``sort`` command and are otherwise processed one chromosome at a time,
so a whole plant genome does not have to be stored as Python objects.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import gzip
import os
import shutil
import subprocess
from typing import Iterable, Iterator, Mapping, Sequence


@dataclass(frozen=True, slots=True)
class Interval:
    chrom: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if not self.chrom:
            raise ValueError("BED chromosome name cannot be empty")
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"Invalid BED interval: {self}")

    @property
    def length(self) -> int:
        return self.end - self.start


def iter_bed(path: str | Path) -> Iterator[Interval]:
    """Yield valid BED3 intervals, rejecting malformed coordinates."""

    bed_path = Path(path)
    with _open_text(bed_path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected at least 3 BED columns"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{bed_path}:{line_number}: BED start/end must be integers"
                ) from exc
            interval = Interval(fields[0], start, end)
            if interval.length:
                yield interval


def write_bed(intervals: Iterable[Interval], destination: str | Path) -> None:
    with Path(destination).open("wt", encoding="utf-8", newline="") as handle:
        for interval in intervals:
            handle.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")


def sort_bed(source: str | Path, destination: str | Path) -> None:
    """Sort a BED-like file lexically by contig and numerically by coordinates."""

    sort_executable = shutil.which("sort")
    if sort_executable is None:
        # Small-platform fallback. The production Linux workflow uses external
        # sort so very large repeat tracks can spill to disk.
        with Path(source).open("rt", encoding="utf-8") as handle:
            lines = [line for line in handle if line.strip() and not line.startswith("#")]
        lines.sort(key=_bed_line_key)
        with Path(destination).open("wt", encoding="utf-8") as handle:
            handle.writelines(lines)
        return

    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    with Path(destination).open("wb") as output_handle:
        completed = subprocess.run(
            [sort_executable, "-k1,1", "-k2,2n", "-k3,3n", str(source)],
            stdout=output_handle,
            stderr=subprocess.PIPE,
            env=environment,
            check=False,
        )
    if completed.returncode:
        message = completed.stderr.decode("utf-8", errors="replace")
        raise RuntimeError(f"BED sort failed: {message.strip()}")


def _bed_line_key(line: str) -> tuple[str, int, int]:
    fields = line.split()
    return fields[0], int(fields[1]), int(fields[2])


def merge_sorted_bed(source: str | Path, destination: str | Path) -> None:
    """Merge overlapping or directly adjacent intervals from a sorted BED."""

    iterator = iter_bed(source)
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        try:
            current = next(iterator)
        except StopIteration:
            return

        for interval in iterator:
            if interval.chrom < current.chrom or (
                interval.chrom == current.chrom and interval.start < current.start
            ):
                raise ValueError(f"BED is not coordinate sorted: {source}")
            if interval.chrom == current.chrom and interval.start <= current.end:
                current = Interval(
                    current.chrom, current.start, max(current.end, interval.end)
                )
                continue
            output.write(f"{current.chrom}\t{current.start}\t{current.end}\n")
            current = interval
        output.write(f"{current.chrom}\t{current.start}\t{current.end}\n")


def concatenate_files(sources: Sequence[str | Path], destination: str | Path) -> None:
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        for source in sources:
            with Path(source).open("rt", encoding="utf-8") as input_handle:
                for line in input_handle:
                    if line.strip():
                        output.write(line)


def grouped_bed(path: str | Path) -> Iterator[tuple[str, list[Interval]]]:
    current_chrom: str | None = None
    current_intervals: list[Interval] = []
    previous_start = -1
    for interval in iter_bed(path):
        if current_chrom is None:
            current_chrom = interval.chrom
        if interval.chrom < current_chrom:
            raise ValueError(f"BED is not lexically sorted by contig: {path}")
        if interval.chrom != current_chrom:
            yield current_chrom, current_intervals
            current_chrom = interval.chrom
            current_intervals = []
            previous_start = -1
        if interval.start < previous_start:
            raise ValueError(f"BED is not numerically sorted within {current_chrom}: {path}")
        current_intervals.append(interval)
        previous_start = interval.start
    if current_chrom is not None:
        yield current_chrom, current_intervals


class _GroupCursor:
    def __init__(self, groups: Iterator[tuple[str, list[Interval]]]):
        self._groups = groups
        self.current: tuple[str, list[Interval]] | None = next(groups, None)

    def advance(self) -> None:
        self.current = next(self._groups, None)


def subtract_sorted_bed(
    minuend: str | Path, subtrahend: str | Path, destination: str | Path
) -> None:
    """Write A minus B for two sorted, internally merged BED files."""

    b_cursor = _GroupCursor(grouped_bed(subtrahend))
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        for chrom, a_intervals in grouped_bed(minuend):
            while b_cursor.current is not None and b_cursor.current[0] < chrom:
                b_cursor.advance()
            b_intervals: list[Interval] = []
            if b_cursor.current is not None and b_cursor.current[0] == chrom:
                b_intervals = b_cursor.current[1]
                b_cursor.advance()
            for interval in subtract_lists(a_intervals, b_intervals):
                output.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")


def subtract_lists(
    minuend: Sequence[Interval], subtrahend: Sequence[Interval]
) -> Iterator[Interval]:
    b_index = 0
    for a_interval in minuend:
        cursor = a_interval.start
        while b_index < len(subtrahend) and subtrahend[b_index].end <= cursor:
            b_index += 1
        scan_index = b_index
        while scan_index < len(subtrahend):
            b_interval = subtrahend[scan_index]
            if b_interval.start >= a_interval.end:
                break
            if b_interval.start > cursor:
                yield Interval(
                    a_interval.chrom, cursor, min(b_interval.start, a_interval.end)
                )
            cursor = max(cursor, b_interval.end)
            if cursor >= a_interval.end:
                break
            scan_index += 1
        if cursor < a_interval.end:
            yield Interval(a_interval.chrom, cursor, a_interval.end)


def intersect_sorted_bed(
    first: str | Path, second: str | Path, destination: str | Path
) -> None:
    b_cursor = _GroupCursor(grouped_bed(second))
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        for chrom, a_intervals in grouped_bed(first):
            while b_cursor.current is not None and b_cursor.current[0] < chrom:
                b_cursor.advance()
            b_intervals: list[Interval] = []
            if b_cursor.current is not None and b_cursor.current[0] == chrom:
                b_intervals = b_cursor.current[1]
                b_cursor.advance()
            for interval in intersect_lists(a_intervals, b_intervals):
                output.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")


def intersect_lists(
    first: Sequence[Interval], second: Sequence[Interval]
) -> Iterator[Interval]:
    first_index = 0
    second_index = 0
    while first_index < len(first) and second_index < len(second):
        left = first[first_index]
        right = second[second_index]
        start = max(left.start, right.start)
        end = min(left.end, right.end)
        if start < end:
            yield Interval(left.chrom, start, end)
        if left.end <= right.end:
            first_index += 1
        else:
            second_index += 1


def complement_bed(
    mask: str | Path,
    genome_lengths: Mapping[str, int],
    destination: str | Path,
) -> None:
    """Write the complement of a merged mask across all reference contigs."""

    mask_groups = iter(grouped_bed(mask))
    current_group = next(mask_groups, None)
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        for chrom in sorted(genome_lengths):
            cursor = 0
            chrom_length = genome_lengths[chrom]
            while current_group is not None and current_group[0] < chrom:
                current_group = next(mask_groups, None)
            intervals: Sequence[Interval] = []
            if current_group is not None and current_group[0] == chrom:
                intervals = current_group[1]
                current_group = next(mask_groups, None)
            for interval in intervals:
                if interval.start > cursor:
                    output.write(f"{chrom}\t{cursor}\t{interval.start}\n")
                cursor = max(cursor, interval.end)
            if cursor < chrom_length:
                output.write(f"{chrom}\t{cursor}\t{chrom_length}\n")


def pad_sorted_bed(
    source: str | Path,
    destination: str | Path,
    genome_lengths: Mapping[str, int],
    left: int,
    right: int,
) -> None:
    """Pad intervals, clamp them to contig limits, and merge the result."""

    previous: Interval | None = None
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        for interval in iter_bed(source):
            if interval.chrom not in genome_lengths:
                raise ValueError(
                    f"BED contig {interval.chrom!r} is absent from the reference FASTA"
                )
            padded = Interval(
                interval.chrom,
                max(0, interval.start - left),
                min(genome_lengths[interval.chrom], interval.end + right),
            )
            if previous is None:
                previous = padded
            elif padded.chrom == previous.chrom and padded.start <= previous.end:
                previous = Interval(
                    previous.chrom, previous.start, max(previous.end, padded.end)
                )
            else:
                output.write(
                    f"{previous.chrom}\t{previous.start}\t{previous.end}\n"
                )
                previous = padded
        if previous is not None:
            output.write(f"{previous.chrom}\t{previous.start}\t{previous.end}\n")


def bed_stats(path: str | Path) -> tuple[int, int]:
    interval_count = 0
    base_count = 0
    for interval in iter_bed(path):
        interval_count += 1
        base_count += interval.length
    return interval_count, base_count


def annotate_bed(
    source: str | Path, destination: str | Path, reason: str
) -> None:
    safe_reason = reason.replace("\t", "_").replace("\n", "_")
    with Path(destination).open("wt", encoding="utf-8", newline="") as output:
        for interval in iter_bed(source):
            output.write(
                f"{interval.chrom}\t{interval.start}\t{interval.end}\t{safe_reason}\n"
            )


def copy_bed3(source: str | Path, destination: str | Path) -> None:
    write_bed(iter_bed(source), destination)


def gff_to_bed(
    source: str | Path,
    destination: str | Path,
    feature_types: set[str] | None = None,
) -> None:
    """Convert GFF/GTF 1-based closed coordinates to BED coordinates."""

    source_path = Path(source)
    with _open_text(source_path) as input_handle, Path(
        destination
    ).open("wt", encoding="utf-8", newline="") as output:
        for line_number, raw in enumerate(input_handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"{source_path}:{line_number}: expected at least 5 GFF columns"
                )
            if feature_types and fields[2] not in feature_types:
                continue
            try:
                start = int(fields[3]) - 1
                end = int(fields[4])
            except ValueError as exc:
                raise ValueError(
                    f"{source_path}:{line_number}: invalid GFF coordinates"
                ) from exc
            interval = Interval(fields[0], start, end)
            if interval.length:
                output.write(
                    f"{interval.chrom}\t{interval.start}\t{interval.end}\n"
                )


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")
