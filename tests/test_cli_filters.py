from genmap_fixture import mock_genmap
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import replace
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import gzip
import json
import math
import shutil
import unittest

import numpy as np
import pysam
from scipy.stats import betabinom

from mapploid.cli import main
from mapploid.bam_counts import estimate_depth_thresholds
import mapploid.pipeline as pipeline
import mapploid.ploidy_model as model


class FilterModelTests(unittest.TestCase):
    def test_custom_mac_likelihood_and_shape_match_independent_probabilities(self):
        # Explicit enumeration catches a filter reaching counts but not likelihood/QC.
        for mac in (1, 3, 4):
            for n in (2 * mac, 2 * mac + 3, 20):
                for lo, hi in ((.10, .90), (.10, .70)):
                    retained = [x for x in range(n + 1) if min(x, n-x) >= mac
                                and lo <= x/n <= hi and lo <= (n-x)/n <= hi]
                    for ploidy in (2, 3, 4):
                        expected = np.zeros(len(retained))
                        for mean, weight in zip(model.EXPECTED_MEANS[ploidy],
                                                model.component_weights(ploidy, .7)):
                            a, b = model.beta_parameters(mean, .03)
                            raw = betabinom.pmf(retained, n, a, b)
                            expected += weight * raw / raw.sum()
                        observed = [math.exp(model.truncated_mixture_log_likelihood(
                            np.array([n]), np.array([x]), np.array([1]),
                            ploidy, .03, .7, lo, hi, min_allele_reads=mac))
                            for x in retained]
                        np.testing.assert_allclose(observed, expected, rtol=1e-10)
                        histogram = Counter({(n, x): probability * 1e6
                                             for x, probability in zip(retained, expected)})
                        fit = model.ModelFit(ploidy, 0, 0, .03, .7, True, 1)
                        observed_profile, fitted_profile, _ = model.fitted_frequency_profile(
                            histogram, fit, lo, hi, min_allele_reads=mac)
                        np.testing.assert_allclose(observed_profile, fitted_profile, atol=1e-12)

    def test_custom_mac_reaches_global_bootstrap_chromosome_and_sensitivity_fits(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / 'counts.tsv.gz'
            model.simulate_allele_counts(counts, 3, sites=60, mean_depth=1,
                                        chromosomes=3, min_allele_reads=1)
            # Every site has depth 2 or 3; an accidental MAC=2 anywhere is invalid.
            outputs = model.run_ploidy_model(model.PloidyModelConfig(
                counts, root/'model', min_sites=1, min_chromosome_sites=1,
                bootstrap_replicates=2, min_allele_reads=1))
            data = json.loads(outputs['manifest'].read_text())
            self.assertEqual(data['result']['selected_sites'], 60)
            self.assertEqual(data['result']['bootstrap_valid'], 2)
            self.assertTrue(all(math.isfinite(fit['bic']) for fit in data['fits'].values()))
            self.assertTrue(math.isfinite(data['sensitivity_4x_fit']['bic']))
            with outputs['chromosomes'].open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(math.isfinite(float(row['bic_2x'])) for row in rows))

    def test_purity_uses_exact_counts_and_configured_threshold(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root/'counts.tsv'
            counts.write_text('chrom\tposition\tdepth\tpair_depth\tpair_purity\t'
                              'allele1_count\tallele2_count\n'
                              'chr1\t1\t18\t17\t0.944444\t12\t5\n')
            args = (counts, 's', 100, .15, .85, 5000, 1)
            self.assertEqual(model.load_and_aggregate_sites(*args).selected_sites, 0)
            self.assertEqual(model.load_and_aggregate_sites(
                *args, pair_purity=.9444444).selected_sites, 1)
            self.assertEqual(model.load_and_aggregate_sites(
                *args, pair_purity=.9444445).selected_sites, 0)

    def test_single_depth_bound_keeps_user_value_and_ordered_limits(self):
        for low, high, expected in ((None, None, (30, 30)), (5, 50, (5, 50)),
                                    (40, None, (40, 40)), (None, 10, (10, 10))):
            values = estimate_depth_thresholds({0: 50, 30: 100}, low, high)
            self.assertEqual((values['min_depth'], values['max_depth']), expected)

    def test_invalid_cli_filters_fail_before_opening_inputs_or_creating_output(self):
        with TemporaryDirectory() as temporary:
            out = Path(temporary)/'out'
            base = ['call', '--bam', '/missing.bam', '--assets', '/missing-assets', '--out', str(out)]
            cases = [(['--allele-min', '.5'], 'allele'),
                     (['--allele-max', '1.1'], 'allele'),
                     (['--allele-min', 'nan'], 'allele'),
                     (['--min-depth', '0'], 'depth'),
                     (['--min-depth', '30', '--max-depth', '10'], 'depth'),
                     (['--max-depth', '-1'], 'depth'),
                     (['--min-mapq', '256'], 'MAPQ'),
                     (['--min-baseq', '-1'], 'base quality'),
                     (['--min-allele-reads', '0'], 'allele reads'),
                     (['--pair-purity', '.4'], 'purity')]
            for flags, message in cases:
                with self.subTest(flags=flags), redirect_stderr(StringIO()) as error:
                    self.assertEqual(main(base + flags), 1)
                    self.assertIn(message, error.getvalue())
                    self.assertNotIn('missing', error.getvalue())
                    self.assertFalse(out.exists())


@unittest.skipUnless(shutil.which('samtools'), 'samtools required')
class FilterCliIntegrationTests(unittest.TestCase):
    def test_all_eight_cli_options_affect_real_bam_and_model_and_protect_resume(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root/'ref.fa'
            reference.write_text('>chr1\n'+'A'*2000+'\n')
            scores = root/'scores.bedgraph'
            scores.write_text('chr1\t0\t1951\t1\n')
            self.enterContext(mock_genmap(scores))
            assets = root/'assets'
            pipeline.prepare(pipeline.PrepareConfig(reference, assets, 50, 1, min_region_length=300))
            bam = root/'sample.bam'
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                      'SQ': [{'SN': 'chr1', 'LN': 2000}]}
            groups = [({'A':9, 'G':1}, 60, 35), ({'A':24, 'G':4, 'C':2}, 60, 35),
                      ({'A':8, 'G':2}, 15, 35), ({'A':8, 'G':2}, 60, 15),
                      ({'A':18, 'G':2}, 60, 35), ({'A':1, 'G':1}, 60, 35),
                      ({'A':30, 'G':30}, 60, 35)]
            with pysam.AlignmentFile(str(bam), 'wb', header=header) as handle:
                for group, (alleles, mapq, baseq) in enumerate(groups):
                    for allele, count in alleles.items():
                        for i in range(count):
                            read = pysam.AlignedSegment()
                            read.query_name = f'{group}_{allele}_{i}'
                            read.query_sequence = 'A'*25 + allele + 'A'*24
                            read.flag = 0
                            read.reference_id = 0
                            read.reference_start = 100 + 200*group
                            read.mapping_quality = mapq
                            read.cigar = [(0, 50)]
                            read.query_qualities = pysam.qualitystring_to_array(chr(baseq+33)*50)
                            handle.write(read)
            pysam.index(str(bam))
            filters = dict(allele_min=.10, allele_max=.90, min_depth=4, max_depth=40,
                           min_mapq=10, min_baseq=10, min_allele_reads=1, pair_purity=.90)
            output = root/'out'
            args = ['call', '--bam', str(bam), '--assets', str(assets), '--out', str(output), '--keep-details']
            for key, value in filters.items():
                args.extend(['--'+key.replace('_', '-'), str(value)])
            with redirect_stdout(StringIO()):
                self.assertEqual(main(args), 0)
            stage = output/'details/stage_1'
            with gzip.open(stage/'counts/sample.biallelic_sites.tsv.gz', 'rt') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual([int(row['position']) for row in rows], [126, 326, 526, 726, 926])
            model_config = json.loads((stage/'model/run_config.json').read_text())
            self.assertEqual(model_config['min_sites'], 25_000)
            self.assertEqual(model_config['max_sites'], 100_000)
            request_path = output/'details/run.json'
            run_record = json.loads(request_path.read_text())
            request = run_record['request']
            self.assertEqual(request['model_site_limits'], dict(min_sites=25_000, max_sites=100_000))
            count_config = json.loads((stage/'counts/run_config.json').read_text())
            for name, value in filters.items():
                self.assertEqual(count_config[name], value)
                if name in ('allele_min', 'allele_max', 'min_allele_reads', 'pair_purity'):
                    self.assertEqual(model_config[name], value)
            evidence = json.loads((stage/'model/run_manifest.json').read_text())
            self.assertEqual(evidence['result']['selected_sites'], 5)
            self.assertTrue(all(math.isfinite(f['bic']) for f in evidence['fits'].values()))
            result = output/'result.tsv'
            stamp = result.stat().st_mtime_ns
            with redirect_stdout(StringIO()):
                self.assertEqual(main(args), 0)
            self.assertEqual(result.stat().st_mtime_ns, stamp)
            config = pipeline.CallConfig(bam, assets, output, keep_details=True, **filters)
            changes = dict(allele_min=.12, allele_max=.88, min_depth=5, max_depth=50,
                           min_mapq=11, min_baseq=11, min_allele_reads=2, pair_purity=.91)
            for key, value in changes.items():
                with self.subTest(changed=key), self.assertRaisesRegex(ValueError, 'settings changed'):
                    pipeline.call(replace(config, **{key:value}))
            # A completed run made with earlier limits must not be silently reused.
            request['model_site_limits'] = dict(min_sites=5_000, max_sites=50_000)
            request_path.write_text(json.dumps(run_record))
            with self.assertRaisesRegex(ValueError, 'settings changed'):
                pipeline.call(config)
            request.pop('model_site_limits')
            request_path.write_text(json.dumps(run_record))
            with self.assertRaisesRegex(ValueError, 'settings changed'):
                pipeline.call(config)


if __name__ == '__main__':
    unittest.main()
