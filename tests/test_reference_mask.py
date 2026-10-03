from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

from mapploid.reference_mask import ReferenceMaskConfig, build_reference_mask


class ReferenceMaskTests(unittest.TestCase):
    def test_reference_workflow_with_existing_bedgraph(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.fa"
            # The N run crosses a FASTA line boundary: [40, 45).
            reference.write_text(
                ">chr1\n" + "A" * 40 + "NN\n" + "NNN" + "C" * 55 + "\n",
                encoding="utf-8",
            )
            bedgraph = root / "reference.genmap.bedgraph"
            bedgraph.write_text(
                "chr1\t0\t30\t1\n"
                "chr1\t30\t35\t0.5\n"
                "chr1\t35\t91\t1\n",
                encoding="utf-8",
            )
            repeats = root / "repeats.bed"
            repeats.write_text("chr1\t60\t70\n", encoding="utf-8")
            features = root / "features.gff3"
            features.write_text(
                "##gff-version 3\n"
                "chr1\ttest\tgene\t1\t80\t.\t+\t.\tID=gene1\n",
                encoding="utf-8",
            )
            results = build_reference_mask(
                ReferenceMaskConfig(
                    reference=reference,
                    output_dir=root / "output",
                    read_length=10,
                    mismatches=1,
                    genmap_bedgraph=bedgraph,
                    gap_padding=0,
                    repeat_beds=[repeats],
                    feature_gff=features,
                )
            )

            self.assertEqual(
                results["mappability_mask"].read_text(encoding="utf-8"),
                "chr1\t9\t30\nchr1\t45\t60\nchr1\t70\t91\n",
            )
            self.assertEqual(
                results["feature_mask"].read_text(encoding="utf-8"),
                "chr1\t9\t30\nchr1\t45\t60\nchr1\t70\t80\n",
            )
            self.assertEqual(
                results["excluded_regions"].read_text(encoding="utf-8"),
                "chr1\t0\t9\nchr1\t30\t45\nchr1\t60\t70\nchr1\t91\t100\n",
            )
            reasons = results["exclusion_sources"].read_text(encoding="utf-8")
            self.assertIn("assembly_gap", reasons)
            self.assertIn("repeat", reasons)
            self.assertIn("low_mappability_or_edge", reasons)
            manifest = json.loads(results["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["workflow"], "reference-mask")
            self.assertEqual(manifest["run_config"]["read_length"], 10)

    def test_reuses_genmap_output_prefix_from_interrupted_run(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
            output = root / "output"
            intermediate = output / "intermediate"
            intermediate.mkdir(parents=True)
            # GenMap -O is this file prefix, not a directory.
            reused = intermediate / "genmap_k10_e1.bedgraph"
            reused.write_text("chr1\t0\t91\t1\n", encoding="utf-8")

            results = build_reference_mask(
                ReferenceMaskConfig(
                    reference=reference,
                    output_dir=output,
                    read_length=10,
                    mismatches=1,
                    genmap_executable="definitely-not-installed",
                )
            )
            self.assertEqual(
                results["mappability_mask"].read_text(encoding="utf-8"),
                "chr1\t9\t91\n",
            )
            manifest = json.loads(results["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["genmap"]["status"], "reused")
            self.assertEqual(Path(manifest["genmap"]["bedgraph"]), reused)

            # A completed result is safe to invoke again with identical parameters.
            second = build_reference_mask(
                ReferenceMaskConfig(
                    reference=reference,
                    output_dir=output,
                    read_length=10,
                    mismatches=1,
                    genmap_executable="definitely-not-installed",
                )
            )
            self.assertEqual(second["mappability_mask"], results["mappability_mask"])


if __name__ == "__main__":
    unittest.main()
