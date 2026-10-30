"""Exercise the full result writer with diagnostic failures and remaining gates."""
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import csv
import json
import unittest

import mapploid.ploidy_model as model


class CallPolicyTests(unittest.TestCase):
    def run_case(self, root, *, support=1., chromosomes=None,
                 selected=25_000, delta=100., converged=True):
        counts = root / 'counts.tsv'
        counts.write_text('test input; site aggregation is supplied by this fixture\n')
        histogram = Counter({(60, 30): selected})
        data = model.AggregatedSites(histogram, {'chr1:1-5000000': histogram},
                                     Counter(chr1=selected), selected, selected)
        fits = {p: model.ModelFit(p, -100., 100. + (p-2)*delta, .03,
                                .7 if p == 4 else None, converged or p != 4, 1)
                for p in (2, 3, 4)}
        if chromosomes is None:
            chromosomes = Counter({2: 3})
        with ExitStack() as stack:
            stack.enter_context(patch.object(model, 'load_and_aggregate_sites', return_value=data))
            stack.enter_context(patch.object(model, 'fit_all_ploidies', return_value=fits))
            stack.enter_context(patch.object(model, 'fit_ploidy', return_value=fits[4]))
            stack.enter_context(patch.object(model, 'block_bootstrap', return_value=([], Counter({2: round(100*support)}))))
            stack.enter_context(patch.object(model, 'chromosome_model_support', return_value=([], chromosomes)))
            outputs = model.run_ploidy_model(model.PloidyModelConfig(counts, root / 'out', sample='policy'))
        manifest = json.loads(outputs['manifest'].read_text())
        with outputs['result'].open() as handle:
            serialized = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(serialized['predicted_ploidy'], manifest['result']['predicted_ploidy'])
        self.assertEqual(serialized['uncertainty_reasons'], manifest['result']['uncertainty_reasons'])
        self.assertEqual(serialized['warnings'], manifest['result']['warnings'])
        self.assertEqual(manifest['mapploid_version'], '0.10.5')
        return manifest, outputs

    def test_supporting_diagnostic_failures_are_warning_only_individually_and_together(self):
        cases = [
            ({'support': .10}, {'LOW_BOOTSTRAP_SUPPORT'}),
            ({'chromosomes': Counter()}, {'INSUFFICIENT_CHROMOSOME_EVIDENCE'}),
            ({'chromosomes': Counter({2: 1, 3: 2})}, {'CHROMOSOME_HETEROGENEITY'}),
            ({'support': .10, 'chromosomes': Counter({3: 1})},
             {'LOW_BOOTSTRAP_SUPPORT', 'INSUFFICIENT_CHROMOSOME_EVIDENCE',
              'CHROMOSOME_HETEROGENEITY'}),
        ]
        for kwargs, expected in cases:
            with self.subTest(kwargs=kwargs), TemporaryDirectory() as temporary:
                manifest, _ = self.run_case(Path(temporary), **kwargs)
                result = manifest['result']
                self.assertEqual(result['predicted_ploidy'], '2x')
                self.assertEqual(result['uncertainty_reasons'], 'NONE')
                self.assertEqual(result['qc_status'], 'WARNING')
                self.assertTrue(expected <= set(result['warnings'].split(';')))

    def test_remaining_three_gates_still_make_calls_uncertain(self):
        cases = [({'selected': 24_999}, {'LOW_SITE_COUNT'}),
                 ({'delta': 9.}, {'LOW_DELTA_BIC'}),
                 ({'converged': False}, {'MODEL_NOT_CONVERGED'}),
                 ({'selected': 24_999, 'delta': 9., 'converged': False},
                  {'LOW_SITE_COUNT', 'LOW_DELTA_BIC', 'MODEL_NOT_CONVERGED'})]
        for kwargs, expected in cases:
            with self.subTest(kwargs=kwargs), TemporaryDirectory() as temporary:
                manifest, _ = self.run_case(Path(temporary), support=.10,
                                           chromosomes=Counter({3: 1}), **kwargs)
                result = manifest['result']
                self.assertEqual(result['predicted_ploidy'], 'UNCERTAIN')
                self.assertEqual(result['qc_status'], 'FAIL')
                self.assertEqual(set(result['uncertainty_reasons'].split(';')), expected)
                self.assertIn('LOW_BOOTSTRAP_SUPPORT', result['warnings'])

    def test_planned_denominator_and_two_display_only_diagnostics(self):
        from mapploid.reporting import diagnostics
        with TemporaryDirectory() as temporary:
            manifest, _ = self.run_case(Path(temporary), support=.20)
            result = manifest['result']
            self.assertEqual(result['bootstrap_planned'], 100)
            self.assertEqual(result['bootstrap_valid'], 20)
            self.assertEqual(result['bootstrap_unresolved'], 80)
            self.assertEqual(result['bootstrap_support'], .20)
            self.assertNotIn('bootstrap_support_valid', result)
            self.assertNotIn('fit_js_divergence', result)
            gates, rows = diagnostics(dict(evidence=result, model=manifest), {})
            self.assertEqual([row[-1] for row in gates], ['PASS'] * 3)
            self.assertEqual([row[0] for row in rows], ['Bootstrap stability', 'Chromosome consistency'])
            self.assertEqual([row[-1] for row in rows], ['WARNING', 'CONSISTENT'])
            self.assertEqual(result['predicted_ploidy'], '2x')

    def test_missing_diagnostics_are_not_shown_as_consistent(self):
        from mapploid.reporting import diagnostics
        _, rows = diagnostics(dict(evidence=dict(bootstrap_planned=100)), {})
        self.assertEqual([row[-1] for row in rows], ['NOT_CALCULATED', 'NOT_ASSESSABLE'])

    def test_boundary_values_pass_remaining_gates(self):
        with TemporaryDirectory() as temporary:
            manifest, _ = self.run_case(Path(temporary), selected=25_000, delta=10.)
            self.assertEqual(manifest['result']['predicted_ploidy'], '2x')

    def test_v0102_model_outputs_cannot_be_reused(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, outputs = self.run_case(root)
            signature_path = root / 'out/run_config.json'
            signature = json.loads(signature_path.read_text())
            signature['model_version'] = 'beta-binomial-v7-configurable-filters'
            signature_path.write_text(json.dumps(signature))
            before = outputs['result'].read_bytes()
            with self.assertRaisesRegex(ValueError, 'different inputs or parameters'):
                model.run_ploidy_model(model.PloidyModelConfig(root / 'counts.tsv', root / 'out', sample='policy'))
            self.assertEqual(outputs['result'].read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
