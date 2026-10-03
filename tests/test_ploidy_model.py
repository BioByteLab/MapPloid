from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import gzip
import unittest

import numpy as np

from mapploid.ploidy_model import (
    DEFAULT_TETRAPLOID_A_MIN,
    PloidyModelConfig,
    beta_binomial_logpmf,
    beta_parameters,
    fit_all_ploidies,
    fit_ploidy,
    proportional_quotas,
    run_ploidy_model,
    simulate_allele_counts,
    truncation_log_normalizers,
)


class PloidyModelUnitTests(unittest.TestCase):
    def test_truncated_component_is_normalized(self):
        depth = np.asarray([20])
        low = np.asarray([3])
        high = np.asarray([17])
        alpha, beta = beta_parameters(1 / 3, 0.02)
        normalizer = truncation_log_normalizers(depth, low, high, alpha, beta)[0]
        counts = np.arange(3, 18)
        depths = np.full(counts.shape, 20)
        probabilities = np.exp(
            beta_binomial_logpmf(counts, depths, alpha, beta) - normalizer
        )
        self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=10)

    def test_proportional_quotas_have_exact_total(self):
        quotas = proportional_quotas(Counter({"chr1": 100, "chr2": 50, "chr3": 25}), 80)
        self.assertEqual(sum(quotas.values()), 80)
        self.assertGreater(quotas["chr1"], quotas["chr2"])
        self.assertGreater(quotas["chr2"], quotas["chr3"])

    def test_exact_theoretical_peaks_select_correct_ploidy(self):
        histograms = {
            2: Counter({(60, 30): 3_000}),
            3: Counter({(60, 20): 1_500, (60, 40): 1_500}),
            4: Counter({(60, 15): 1_050, (60, 30): 900, (60, 45): 1_050}),
        }
        for expected, histogram in histograms.items():
            fits = fit_all_ploidies(histogram, 0.15, 0.85, (1e-5, 0.25))
            observed = min(fits.values(), key=lambda fit: fit.bic).ploidy
            self.assertEqual(observed, expected)

    def test_constrained_m4p_cannot_collapse_onto_pure_half_peak(self):
        histogram = Counter({(60, 30): 3_000})
        fits = fit_all_ploidies(
            histogram, 0.15, 0.85, (1e-5, 0.25)
        )
        self.assertEqual(min(fits.values(), key=lambda fit: fit.bic).ploidy, 2)
        self.assertGreater(fits[4].bic - fits[2].bic, 100)
        self.assertGreaterEqual(
            fits[4].outer_fraction, DEFAULT_TETRAPLOID_A_MIN
        )

    def test_m4p_accepts_multiple_biological_outer_peak_weights(self):
        total = 7_000
        for outer_fraction in (4 / 7, 0.70, 0.90):
            side = round(total * outer_fraction / 2)
            middle = total - 2 * side
            histogram = Counter({
                (80, 20): side,
                (80, 40): middle,
                (80, 60): side,
            })
            fits = fit_all_ploidies(
                histogram, 0.15, 0.85, (1e-5, 0.25)
            )
            self.assertEqual(
                min(fits.values(), key=lambda fit: fit.bic).ploidy, 4
            )


class PloidyModelIntegrationTests(unittest.TestCase):
    def test_simulator_and_uncertain_bootstrap_workflow(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "simulated_3x.tsv.gz"
            simulate_allele_counts(
                counts,
                ploidy=3,
                sites=1_000,
                mean_depth=35,
                rho=0.01,
                contamination_fraction=0.0,
                chromosomes=10,
                seed=71,
            )
            with gzip.open(counts, "rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1_000)
            frequencies = [
                int(row["allele1_count"]) / int(row["pair_depth"])
                for row in rows
            ]
            near_expected = sum(
                min(abs(value - 1 / 3), abs(value - 2 / 3)) < 0.15
                for value in frequencies
            )
            self.assertGreater(near_expected / len(frequencies), 0.90)

            outputs = run_ploidy_model(
                PloidyModelConfig(
                    counts=counts,
                    output_dir=root / "model",
                    sample="sim3",
                    max_sites=1_000,
                    min_sites=2_000,
                    block_size=500_000,
                    bootstrap_replicates=3,
                    bootstrap_support=0.5,
                    delta_bic=0,
                    min_blocks=2,
                    seed=71,
                )
            )
            result = outputs["result"].read_text(encoding="utf-8").splitlines()
            values = dict(zip(result[0].split("\t"), result[1].split("\t")))
            self.assertEqual(values["best_hypothesis"], "3x")
            self.assertEqual(values["predicted_ploidy"], "UNCERTAIN")
            self.assertEqual(values["qc_status"], "FAIL")
            self.assertIn("LOW_SITE_COUNT", values["uncertainty_reasons"])
            self.assertEqual(values["best_model"], "M3")
            bootstrap = outputs["bootstrap"].read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(bootstrap), 4)

    def test_2x_scientific_limit_is_warning_not_uncertain(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "simulated_2x.tsv.gz"
            simulate_allele_counts(
                counts,
                ploidy=2,
                sites=2_000,
                mean_depth=35,
                rho=0.01,
                contamination_fraction=0.0,
                chromosomes=10,
                seed=73,
            )
            outputs = run_ploidy_model(
                PloidyModelConfig(
                    counts=counts,
                    output_dir=root / "model",
                    sample="sim2",
                    max_sites=2_000,
                    min_sites=1_000,
                    block_size=500_000,
                    bootstrap_replicates=3,
                    bootstrap_support=0.5,
                    delta_bic=0,
                    min_blocks=1,
                    min_chromosome_sites=1,
                    max_fit_js=0.10,
                    seed=73,
                )
            )
            lines = outputs["result"].read_text(encoding="utf-8").splitlines()
            result = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            self.assertEqual(result["predicted_ploidy"], "2x")
            self.assertEqual(result["qc_status"], "WARNING")
            self.assertEqual(result["uncertainty_reasons"], "NONE")
            self.assertIn("DISOMIC_4X_NOT_EXCLUDED", result["warnings"])

    def test_parameter_boundary_is_warning_only(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "simulated_2x.tsv.gz"
            simulate_allele_counts(
                counts,
                ploidy=2,
                sites=2_000,
                mean_depth=40,
                rho=0.08,
                contamination_fraction=0.0,
                chromosomes=10,
                seed=79,
            )
            outputs = run_ploidy_model(
                PloidyModelConfig(
                    counts=counts,
                    output_dir=root / "model",
                    sample="boundary2",
                    max_sites=2_000,
                    min_sites=1_000,
                    rho_max=0.08,
                    block_size=500_000,
                    bootstrap_replicates=0,
                    delta_bic=0,
                    min_blocks=1,
                    min_chromosome_sites=1,
                    max_fit_js=0.10,
                    seed=79,
                )
            )
            lines = outputs["result"].read_text(encoding="utf-8").splitlines()
            result = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            self.assertEqual(result["predicted_ploidy"], "2x")
            self.assertEqual(result["uncertainty_reasons"], "NONE")
            self.assertIn("RHO_AT_BOUNDARY", result["warnings"])

    def test_poor_theoretical_fit_is_uncertain(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "simulated_3x.tsv.gz"
            simulate_allele_counts(
                counts,
                ploidy=3,
                sites=2_000,
                mean_depth=35,
                rho=0.01,
                contamination_fraction=0.0,
                chromosomes=10,
                seed=83,
            )
            outputs = run_ploidy_model(
                PloidyModelConfig(
                    counts=counts,
                    output_dir=root / "model",
                    sample="fitfail3",
                    max_sites=2_000,
                    min_sites=1_000,
                    block_size=500_000,
                    bootstrap_replicates=0,
                    delta_bic=0,
                    min_blocks=1,
                    min_chromosome_sites=1,
                    max_fit_js=0.0,
                    seed=83,
                )
            )
            lines = outputs["result"].read_text(encoding="utf-8").splitlines()
            result = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            self.assertEqual(result["predicted_ploidy"], "UNCERTAIN")
            self.assertEqual(result["qc_status"], "FAIL")
            self.assertIn("POOR_THEORETICAL_FIT", result["uncertainty_reasons"])


if __name__ == "__main__":
    unittest.main()
