from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import json
import math
import unittest

import numpy as np
from scipy.stats import betabinom
import mapploid.ploidy_model as model


class ModelCorrections(unittest.TestCase):
    def test_low_depth_mac_likelihood_matches_independent_enumeration(self):
        for n in (4, 5, 6, 7, 8, 20):
            for lo, hi in ((.15, .85), (.10, .70)):
                retained = [x for x in range(n + 1) if min(x, n-x) >= 2
                            and lo <= x/n <= hi and lo <= (n-x)/n <= hi]
                for ploidy in (2, 3, 4):
                    weights = model.component_weights(ploidy, .7)
                    expected = np.zeros(len(retained))
                    for mean, weight in zip(model.EXPECTED_MEANS[ploidy], weights):
                        alpha, beta = model.beta_parameters(mean, .03)
                        raw = betabinom.pmf(retained, n, alpha, beta)
                        expected += weight * raw / raw.sum()
                    observed = []
                    for x in retained:
                        ll = model.truncated_mixture_log_likelihood(
                            np.array([n]), np.array([x]), np.array([1]), ploidy, .03, .7, lo, hi)
                        observed.append(math.exp(ll))
                    np.testing.assert_allclose(observed, expected, rtol=1e-10)
                    self.assertAlmostEqual(sum(observed), 1., places=10)
        for p in (2, 3, 4):
            ll = model.truncated_mixture_log_likelihood(
                np.array([4]), np.array([2]), np.array([100]), p, .03, .7, .15, .85)
            self.assertAlmostEqual(ll, 0., places=10)

    def test_chromosome_diagnostics_record_shortage_and_conflict(self):
        for n in (0, 1, 2):
            self.assertIn('INSUFFICIENT_CHROMOSOME_EVIDENCE', model.chromosome_evidence_warnings(n, 1.))
        self.assertEqual(model.chromosome_evidence_warnings(10, .70), [])
        self.assertIn('CHROMOSOME_HETEROGENEITY', model.chromosome_evidence_warnings(1, 0.))
        self.assertIn('CHROMOSOME_HETEROGENEITY', model.chromosome_evidence_warnings(10, .69))

    def test_failed_models_do_not_vote_in_bootstrap_or_chromosomes(self):
        fits = {p: model.ModelFit(p, -100., 200.+p, .03, .7 if p==4 else None, p!=4, 1)
                for p in (2, 3, 4)}
        blocks = {'chr1:1-5000': Counter({(30, 15): 1000})}
        with patch.object(model, 'fit_all_ploidies', return_value=fits):
            rows, votes = model.block_bootstrap(blocks, fits, 3, .15, .85, (1e-5,.25), 4/7, 71)
            self.assertEqual(sum(votes.values()), 0)
            self.assertTrue(all(r['status']=='UNRESOLVED' and r['best_ploidy'] is None for r in rows))
            rows, votes = model.chromosome_model_support(blocks, .15, .85, (1e-5,.25), 4/7, 200, 0)
            self.assertEqual(sum(votes.values()), 0)
            self.assertEqual(rows[0]['status'], 'UNRESOLVED')

    def test_duplicate_counts_and_filter_contract_mismatch_are_rejected(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root/'counts.tsv'
            header = 'chrom\tposition\tpair_depth\tpair_purity\tallele1_count\tallele2_count\n'
            counts.write_text(header+'chr1\t1\t10\t1\t5\t5\n'*2)
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                model.load_and_aggregate_sites(counts, 's', 100, .15,.85,5000,1)
            (root/'run_config.json').write_text(json.dumps(dict(allele_min=.15,allele_max=.85,
                                                             min_allele_reads=1,pair_purity=.95)))
            with self.assertRaisesRegex(ValueError, 'min_allele_reads'):
                model.validate_counts_contract(counts, model.PloidyModelConfig(counts,root/'out'))

    def test_low_depth_simulator_uses_the_same_support(self):
        import gzip,csv
        with TemporaryDirectory() as temporary:
            path=Path(temporary)/'counts.tsv.gz'
            model.simulate_allele_counts(path,3,100,1,.02,seed=8)
            with gzip.open(path,'rt') as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertTrue(all(int(r['allele1_count'])>=2 and int(r['allele2_count'])>=2 for r in rows))
            self.assertTrue(any(int(r['pair_depth'])==4 for r in rows))


if __name__ == '__main__':
    unittest.main()
