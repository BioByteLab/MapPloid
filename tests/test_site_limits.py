from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

from mapploid.legacy_cli import build_parser
from mapploid.ploidy_model import PloidyModelConfig, run_ploidy_model


class SiteLimitTests(unittest.TestCase):
    def test_lower_boundary_and_upper_sampling_cap(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            for n in (24_999, 25_000, 100_001):
                counts = root/f'counts_{n}.tsv'
                with counts.open('w') as handle:
                    handle.write('chrom\tposition\tdepth\tpair_depth\tpair_purity\tallele1_count\tallele2_count\n')
                    for i in range(n):
                        handle.write(f'chr1\t{i+1}\t30\t30\t1.0\t15\t15\n')
                outputs = run_ploidy_model(PloidyModelConfig(
                    counts, root/f'model_{n}', sample='test', bootstrap_replicates=0))
                manifest = json.loads(outputs['manifest'].read_text())
                result = manifest['result']
                self.assertEqual(result['available_sites'], n)
                self.assertEqual(result['selected_sites'], min(n,100_000))
                self.assertEqual('LOW_SITE_COUNT' in result['uncertainty_reasons'], n<25_000)
                self.assertEqual(result['predicted_ploidy'], 'UNCERTAIN' if n < 25_000 else '2x')
                self.assertIn('INSUFFICIENT_CHROMOSOME_EVIDENCE', result['warnings'])

    def test_legacy_cli_default_and_explicit_limits(self):
        parser = build_parser()
        args = ['ploidy-model','--counts','counts.tsv','--out-dir','model']
        default = parser.parse_args(args)
        self.assertEqual((default.min_sites,default.max_sites),(25_000,100_000))
        custom = parser.parse_args(args+['--min-sites','1000','--max-sites','2000'])
        self.assertEqual((custom.min_sites,custom.max_sites),(1000,2000))
