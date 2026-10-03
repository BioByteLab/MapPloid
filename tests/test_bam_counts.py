from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import shutil
import unittest

import pysam

from mapploid.bam_counts import (
    BamCountConfig,
    build_bam_counts,
    classify_site,
    count_pileup_bases,
    count_simple_pileup_bases,
    estimate_depth_thresholds,
    filter_sites,
    weighted_quantile,
)


class BamCountUnitTests(unittest.TestCase):
    def test_parse_mpileup_bases(self):
        # Six reference A observations, two G observations and ignored indels.
        read_bases = "^].,.,.,Gg+2TT-1a$*<>"
        self.assertEqual(count_pileup_bases(read_bases, "A"), (6, 0, 2, 0))
        self.assertEqual(count_simple_pileup_bases(".,.,Gg*<>", "A"), (4, 0, 2, 0))

    def test_weighted_quantile_and_depth_limits(self):
        histogram = {10: 2, 20: 5, 30: 3}
        self.assertEqual(weighted_quantile(histogram, 0.5), 20)
        thresholds = estimate_depth_thresholds(histogram, None, None)
        self.assertEqual(thresholds["depth_median"], 20)
        self.assertEqual(thresholds["min_depth"], 10)
        self.assertGreaterEqual(thresholds["max_depth"], 20)

    def test_candidate_classification(self):
        self.assertEqual(
            classify_site((10, 0, 10, 0), 20, 0.15, 0.85, 0.95),
            "candidate",
        )
        self.assertEqual(
            classify_site((7, 7, 6, 0), 20, 0.15, 0.85, 0.95),
            "multiallelic",
        )
        self.assertEqual(
            classify_site((19, 0, 1, 0), 20, 0.15, 0.85, 0.95),
            "not_biallelic",
        )
        self.assertEqual(
            classify_site((10, 0, 8, 2), 20, 0.15, 0.85, 0.95),
            "low_pair_purity",
        )
        self.assertEqual(
            classify_site((10, 0, 9, 1), 20, 0.15, 0.85, 0.95),
            "candidate",
        )

    def test_minimum_allele_reads_is_applied_after_frequency_filter(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "candidate.tsv.gz"
            output = root / "retained.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("chrom\tposition\tref\tmpileup_depth\tdepth\tA\tC\tG\tT\n")
                handle.write("chr1\t1\tA\t6\t6\t5\t0\t1\t0\n")
                handle.write("chr1\t2\tA\t6\t6\t4\t0\t2\t0\n")
            config = BamCountConfig(
                bam=root / "unused.bam",
                reference=root / "unused.fa",
                mask=root / "unused.bed",
                output_dir=root / "unused-output",
            )
            summary = filter_sites(
                config, source, output, {"min_depth": 1, "max_depth": 100}
            )
            self.assertEqual(summary["retained_sites"], 1)
            self.assertEqual(summary["low_allele_support"], 1)
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = handle.read().splitlines()
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[1].split("\t")[1], "2")

    def test_minimum_allele_reads_must_be_positive(self):
        config = BamCountConfig(
            bam=Path("sample.bam"),
            reference=Path("reference.fa"),
            mask=Path("mask.bed"),
            output_dir=Path("output"),
            min_allele_reads=0,
        )
        with self.assertRaisesRegex(ValueError, "minimum allele reads"):
            config.validate()


@unittest.skipUnless(shutil.which("samtools"), "samtools is not installed")
class BamCountIntegrationTests(unittest.TestCase):
    def test_observed_zero_acgt_and_uncovered_positions_are_both_counted(self):
        with TemporaryDirectory() as temporary:
            root=Path(temporary)
            reference=root/'ref.fa'
            reference.write_text('>chr1\n'+'A'*100+'\n')
            pysam.faidx(str(reference))
            mask=root/'mask.bed';mask.write_text('chr1\t0\t100\n')
            bam=root/'sample.bam'
            header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100}]}
            with pysam.AlignmentFile(str(bam),'wb',header=header) as handle:
                read=self._read('ambiguous','A',60,35,0)
                read.query_sequence='N'*20
                read.query_qualities=pysam.qualitystring_to_array('I'*20)
                handle.write(read)
            pysam.index(str(bam))
            result=build_bam_counts(BamCountConfig(bam,reference,mask,root/'out'))
            self.assertIn('0\t100\t1.00000000',result['depth_distribution'].read_text())

    def test_small_bam_end_to_end(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
            pysam.faidx(str(reference))
            mask = root / "mask.bed"
            mask.write_text("chr1\t0\t100\n", encoding="utf-8")
            bam = root / "sample.bam"
            self._write_bam(bam)
            pysam.index(str(bam))

            outputs = build_bam_counts(
                BamCountConfig(
                    bam=bam,
                    reference=reference,
                    mask=mask,
                    output_dir=root / "output",
                    sample="synthetic",
                    min_mapq=20,
                    min_baseq=20,
                    min_depth=4,
                    max_depth=100,
                )
            )
            with gzip.open(outputs["biallelic_sites"], "rt", encoding="utf-8") as handle:
                rows = [line.rstrip("\n").split("\t") for line in handle]
            self.assertEqual(len(rows), 2)
            values = dict(zip(rows[0], rows[1]))
            self.assertEqual(values["position"], "21")
            self.assertEqual(values["A"], "6")
            self.assertEqual(values["G"], "6")
            self.assertEqual(values["C"], "0")
            self.assertEqual(values["T"], "0")
            self.assertEqual(values["pair_depth"], "12")
            self.assertEqual(values["pair_purity"], "1.000000")
            self.assertEqual(values["residual_count"], "0")
            self.assertEqual(values["allele1_frequency"], "0.500000")
            self.assertEqual(values["minor_baf"], "0.500000")
            self.assertNotIn("base_counts", outputs)
            self.assertTrue(
                (root / "output" / "intermediate" / "synthetic.candidate_sites.tsv.gz").is_file()
            )

            depth_rows = outputs["depth_distribution"].read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertIn("0\t80\t0.80000000", depth_rows)
            thresholds = outputs["thresholds"].read_text(encoding="utf-8")
            self.assertIn("min_allele_reads\t2", thresholds)
            filtering = outputs["filter_summary"].read_text(encoding="utf-8")
            self.assertIn("low_allele_support\t0", filtering)

    @staticmethod
    def _write_bam(path: Path) -> None:
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 100}],
        }
        records = []
        for allele, count in (("A", 6), ("G", 6)):
            for index in range(count):
                records.append(
                    BamCountIntegrationTests._read(
                        f"{allele}{index}", allele, 60, 35, 16 if index % 2 else 0
                    )
                )
        records.extend(
            [
                BamCountIntegrationTests._read("low_mq", "C", 10, 35, 0),
                BamCountIntegrationTests._read("low_bq", "T", 60, 10, 0),
                BamCountIntegrationTests._read("duplicate", "C", 60, 35, 1024),
                BamCountIntegrationTests._read("secondary", "T", 60, 35, 256),
                BamCountIntegrationTests._read("supplementary", "C", 60, 35, 2048),
            ]
        )
        with pysam.AlignmentFile(path, "wb", header=header) as output:
            for record in records:
                output.write(record)

    @staticmethod
    def _read(name: str, allele: str, mapq: int, baseq: int, flag: int):
        read = pysam.AlignedSegment()
        read.query_name = name
        sequence = list("A" * 20)
        sequence[10] = allele
        read.query_sequence = "".join(sequence)
        read.flag = flag
        read.reference_id = 0
        read.reference_start = 10
        read.mapping_quality = mapq
        read.cigar = ((0, 20),)
        read.query_qualities = pysam.qualitystring_to_array(chr(baseq + 33) * 20)
        return read


if __name__ == "__main__":
    unittest.main()
