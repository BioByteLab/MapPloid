from contextlib import redirect_stderr, redirect_stdout
from dataclasses import replace
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

from mapploid.cli import build_parser, main
from mapploid.pipeline import PrepareConfig, prepare, load_assets
from mapploid.reference_mask import ReferenceMaskConfig, build_reference_mask


class PrepareOptionsTests(unittest.TestCase):
    def fixture(self, root):
        reference = root / 'ref.fa'
        reference.write_text('>chr1\n' + 'A' * 40000 + '\n')
        scores = root / 'scores.bedgraph'
        # With k=5, strict coverage removes four bp at each score-1 run start.
        # The resulting continuous intervals have lengths 9999, 10000, 10001.
        scores.write_text('chr1\t0\t10003\t1\n'
                          'chr1\t10003\t12000\t0.5\n'
                          'chr1\t12000\t22004\t1\n'
                          'chr1\t22004\t24000\t0.5\n'
                          'chr1\t24000\t34005\t1\n')
        return reference, scores

    def test_cli_defaults_and_complete_prepare_flow(self):
        args = build_parser().parse_args(['prepare', '--reference', 'r.fa', '--out', 'a'])
        self.assertEqual(args.kmer, 150)
        self.assertEqual(args.min_region_length, 10000)
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, scores = self.fixture(root)
            assets = root / 'assets'
            command = ['prepare', '--reference', str(reference), '--out', str(assets),
                       '--kmer', '5', '--genmap-bedgraph', str(scores)]
            with redirect_stdout(StringIO()) as output, redirect_stderr(StringIO()):
                self.assertEqual(main(command), 0)
            self.assertEqual(output.getvalue().strip(), str(assets))
            manifest = load_assets(assets)
            self.assertEqual(manifest['schema'], 2)
            self.assertEqual(manifest['kmer_length'], 5)
            self.assertEqual(manifest['min_region_length'], 10000)
            self.assertNotIn('read_length', manifest)
            self.assertEqual((assets / 'mappable_regions.bed').read_text(),
                             'chr1\t12004\t22004\nchr1\t24004\t34005\n')
            stage = assets / manifest['stages'][0]['bed']
            self.assertEqual(stage.read_text(), (assets / 'mappable_regions.bed').read_text())
            request = json.loads((assets / 'details/request.json').read_text())
            self.assertEqual(request['min_region_length'], 10000)
            self.assertEqual(request['kmer_length'], 5)
            modified = (assets / 'assets.json').stat().st_mtime_ns
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                self.assertEqual(main(command), 0)
                self.assertEqual((assets / 'assets.json').stat().st_mtime_ns, modified)
                error = StringIO()
                with redirect_stderr(error):
                    self.assertEqual(main(command + ['--min-region-length', '9999']), 1)
                self.assertIn('settings changed', error.getvalue())
                with redirect_stderr(StringIO()):
                    self.assertEqual(main(command + ['--kmer', '6']), 1)

    def test_threshold_applies_after_exclusions_split_a_long_run(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, scores = self.fixture(root)
            excluded = root / 'exclude.bed'
            excluded.write_text('chr1\t25000\t25001\n')
            out = build_reference_mask(ReferenceMaskConfig(
                reference, root / 'mask', kmer_length=5, min_region_length=10000,
                genmap_bedgraph=scores, exclude_beds=[excluded]))
            self.assertEqual(out['mappability_mask'].read_text(), 'chr1\t12004\t22004\n')
            reasons = out['exclusion_sources'].read_text()
            self.assertIn('chr1\t24004\t25000\tshort_continuous_region', reasons)
            self.assertIn('chr1\t25001\t34005\tshort_continuous_region', reasons)
            self.assertIn('user_excluded', reasons)

    def test_invalid_parameters_fail_before_opening_reference(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            base = PrepareConfig(root / 'missing.fa', root / 'assets')
            for config, message in [(replace(base, min_region_length=0), 'region length'),
                                    (replace(base, min_region_length=-1), 'region length'),
                                    (replace(base, kmer_length=2), 'k-mer length')]:
                with self.subTest(config=config), self.assertRaisesRegex(ValueError, message):
                    prepare(config)
                self.assertFalse(config.output_dir.exists())

    def test_legacy_assets_keep_their_original_threshold(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, scores = self.fixture(root)
            assets = root / 'assets'
            prepare(PrepareConfig(reference, assets, kmer_length=5,
                                  genmap_bedgraph=scores, min_region_length=300))
            manifest = json.loads((assets / 'assets.json').read_text())
            manifest['schema'] = 1
            manifest['read_length'] = manifest.pop('kmer_length')
            manifest.pop('min_region_length')
            (assets / 'assets.json').write_text(json.dumps(manifest))
            loaded = load_assets(assets)
            self.assertEqual(loaded['kmer_length'], 5)
            self.assertEqual(loaded['min_region_length'], 300)

    def test_threshold_with_no_qualifying_runs_produces_empty_valid_assets(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference, scores = self.fixture(root)
            assets = root / 'assets'
            prepare(PrepareConfig(reference, assets, kmer_length=5,
                                  genmap_bedgraph=scores, min_region_length=20000))
            manifest = load_assets(assets)
            self.assertEqual(manifest['eligible_bases'], 0)
            self.assertEqual((assets / 'mappable_regions.bed').read_text(), '')
            self.assertEqual(manifest['stages'][0]['bases'], 0)


if __name__ == '__main__':
    unittest.main()
