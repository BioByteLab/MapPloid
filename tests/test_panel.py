from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

from mapploid.intervals import iter_bed
from mapploid.panel import PanelConfig, build_panels


class PanelTests(unittest.TestCase):
    def test_nested_balanced_panels_and_minimum_run(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            sizes = root / "genome.sizes"
            sizes.write_text("chr1\t2000000\nchr2\t2000000\n", encoding="utf-8")
            mask = root / "mask.bed"
            mask.write_text(
                "chr1\t0\t900000\n"
                "chr1\t900100\t900200\n"
                "chr1\t900300\t2000000\n"
                "chr2\t0\t900000\n"
                "chr2\t900100\t900200\n"
                "chr2\t900300\t2000000\n",
                encoding="utf-8",
            )
            outputs = build_panels(
                PanelConfig(
                    mask=mask,
                    genome_sizes=sizes,
                    output_dir=root / "panels",
                    target_mb=(1, 2),
                    min_run=300,
                    block_size=500_000,
                    chunk_size=10_000,
                )
            )
            panel1 = list(iter_bed(outputs["panel_1mb"]))
            panel2 = list(iter_bed(outputs["panel_2mb"]))
            self.assertGreaterEqual(sum(row.length for row in panel1), 1_000_000)
            self.assertGreaterEqual(sum(row.length for row in panel2), 2_000_000)
            self.assertTrue(all(row.length >= 300 for row in panel1 + panel2))
            self.assertEqual({row.chrom for row in panel1}, {"chr1", "chr2"})
            for small in panel1:
                self.assertTrue(any(
                    small.chrom == large.chrom
                    and large.start <= small.start
                    and large.end >= small.end
                    for large in panel2
                ))
            manifest = json.loads(outputs["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["source"]["removed_short_intervals"], 2)
            self.assertEqual(manifest["run_config"]["min_run"], 300)


if __name__ == "__main__":
    unittest.main()
