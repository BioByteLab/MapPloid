from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from mapploid.intervals import (
    intersect_sorted_bed,
    merge_sorted_bed,
    sort_bed,
    subtract_sorted_bed,
)


class IntervalTests(unittest.TestCase):
    def test_merge_subtract_and_intersect(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw_a = root / "a.raw.bed"
            raw_a.write_text(
                "chr2\t0\t10\nchr1\t20\t40\nchr1\t0\t10\nchr1\t8\t25\n",
                encoding="utf-8",
            )
            sorted_a = root / "a.sorted.bed"
            merged_a = root / "a.bed"
            sort_bed(raw_a, sorted_a)
            merge_sorted_bed(sorted_a, merged_a)
            self.assertEqual(
                merged_a.read_text(encoding="utf-8"),
                "chr1\t0\t40\nchr2\t0\t10\n",
            )

            raw_b = root / "b.raw.bed"
            raw_b.write_text("chr1\t5\t8\nchr1\t15\t30\n", encoding="utf-8")
            sorted_b = root / "b.sorted.bed"
            merged_b = root / "b.bed"
            sort_bed(raw_b, sorted_b)
            merge_sorted_bed(sorted_b, merged_b)

            subtracted = root / "subtracted.bed"
            subtract_sorted_bed(merged_a, merged_b, subtracted)
            self.assertEqual(
                subtracted.read_text(encoding="utf-8"),
                "chr1\t0\t5\nchr1\t8\t15\nchr1\t30\t40\nchr2\t0\t10\n",
            )

            intersected = root / "intersected.bed"
            intersect_sorted_bed(merged_a, merged_b, intersected)
            self.assertEqual(
                intersected.read_text(encoding="utf-8"),
                "chr1\t5\t8\nchr1\t15\t30\n",
            )


if __name__ == "__main__":
    unittest.main()

