"""Fast beta-binomial mixture models for 2x/3x/polysomic-like 4x inference.

The expected allele-frequency means are fixed by allele dosage.  Only a
sample-level beta-binomial overdispersion and (for 4x) the proportion of
0.25/0.75 outer-peak sites are estimated.  The primary 4x model constrains
that outer proportion so it cannot collapse into the 2x model.  Likelihoods
are conditioned on the allele-frequency truncation.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import heapq
import json
import math
import re

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, gammaln, logsumexp

from . import __version__
MODEL_VERSION = "beta-binomial-v9-three-call-gates-output4"
MIN_ALLELE_READS = 2
DEFAULT_MIN_SITES = 25_000
DEFAULT_MAX_SITES = 100_000
DEFAULT_TETRAPLOID_A_MIN = 4 / 7
DEFAULT_SENSITIVITY_A_MIN = 0.5001
DEFAULT_FIT_BIN_WIDTH = 0.005
PLOIDIES = (2, 3, 4)
EXPECTED_MEANS = {
    2: (0.5,),
    3: (1 / 3, 2 / 3),
    4: (0.25, 0.5, 0.75),
}
MODEL_LABELS = {2: "M2", 3: "M3", 4: "M4P"}
PREDICTION_LABELS = {2: "2x", 3: "3x", 4: "4x"}


@dataclass(frozen=True, slots=True)
class ModelSite:
    chrom: str
    position: int
    depth: int
    allele_count: int


@dataclass(slots=True)
class AggregatedSites:
    histogram: Counter[tuple[int, int]]
    block_histograms: dict[str, Counter[tuple[int, int]]]
    chromosome_counts: Counter[str]
    available_sites: int
    selected_sites: int


@dataclass(slots=True)
class ModelFit:
    ploidy: int
    log_likelihood: float
    bic: float
    rho: float
    outer_fraction: float | None
    converged: bool
    iterations: int


@dataclass(slots=True)
class PloidyModelConfig:
    counts: Path
    output_dir: Path
    sample: str | None = None
    max_sites: int = DEFAULT_MAX_SITES
    min_sites: int = DEFAULT_MIN_SITES
    allele_min: float = 0.15
    allele_max: float = 0.85
    rho_min: float = 1e-5
    rho_max: float = 0.25
    tetraploid_a_min: float = DEFAULT_TETRAPLOID_A_MIN
    sensitivity_a_min: float = DEFAULT_SENSITIVITY_A_MIN
    block_size: int = 5_000_000
    bootstrap_replicates: int = 100
    bootstrap_support: float = 0.85
    delta_bic: float = 10.0
    min_blocks: int = 5
    chromosome_support: float = 0.70
    min_chromosome_sites: int = 200
    chromosome_delta_bic: float = 10.0
    min_informative_chromosomes: int = 3
    seed: int = 20260803
    min_allele_reads: int = MIN_ALLELE_READS
    pair_purity: float = 0.95

    def validate(self) -> None:
        if self.max_sites < 1 or self.min_sites < 1:
            raise ValueError("maximum and minimum site counts must be positive")
        if not 0 <= self.allele_min < 0.5 < self.allele_max <= 1:
            raise ValueError("allele truncation must contain 0.5")
        if self.min_allele_reads < 1:
            raise ValueError("minimum allele reads must be positive")
        if not 0.5 <= self.pair_purity <= 1:
            raise ValueError("pair purity must be in [0.5, 1]")
        if not 0 < self.rho_min < self.rho_max < 1:
            raise ValueError("rho bounds must satisfy 0 < min < max < 1")
        if not 0.5 < self.tetraploid_a_min < 1:
            raise ValueError("tetraploid a-min must be in (0.5, 1)")
        if not 0.5 < self.sensitivity_a_min <= self.tetraploid_a_min:
            raise ValueError("sensitivity a-min must be in (0.5, tetraploid a-min]")
        if self.block_size < 1 or self.bootstrap_replicates < 0:
            raise ValueError("block size must be positive and bootstrap count nonnegative")
        if not 0 < self.bootstrap_support <= 1:
            raise ValueError("bootstrap support threshold must be in (0, 1]")
        if self.delta_bic < 0 or self.min_blocks < 1:
            raise ValueError("delta BIC must be nonnegative and minimum blocks positive")
        if not 0 < self.chromosome_support <= 1:
            raise ValueError("chromosome support threshold must be in (0, 1]")
        if self.chromosome_delta_bic < 0 or self.min_informative_chromosomes < 1:
            raise ValueError("invalid chromosome evidence thresholds")
        if self.min_chromosome_sites < 1:
            raise ValueError("minimum chromosome site count must be positive")


def run_ploidy_model(config: PloidyModelConfig) -> dict[str, Path]:
    """Fit M2/M3/M4P and calculate block and chromosome support."""

    config.validate()
    counts = config.counts.resolve()
    if not counts.is_file():
        raise FileNotFoundError(f"allele-count table does not exist: {counts}")
    sample = config.sample or _sample_from_counts(counts)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
        raise ValueError("sample name may contain only letters, digits, '.', '_' and '-'")
    output_dir = config.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = _output_paths(output_dir, sample)
    signature = _signature(config, counts, sample)
    _check_or_write_signature(output_dir / "run_config.json", signature)
    required = (
        "comparison", "sensitivity", "bootstrap", "blocks", "chromosomes",
        "result", "manifest",
    )
    if all(outputs[name].is_file() for name in required):
        return outputs

    validate_counts_contract(counts, config)
    data = load_and_aggregate_sites(
        counts=counts,
        sample=sample,
        max_sites=config.max_sites,
        allele_min=config.allele_min,
        allele_max=config.allele_max,
        block_size=config.block_size,
        seed=config.seed,
        min_allele_reads=config.min_allele_reads,
        pair_purity=config.pair_purity,
    )
    if data.selected_sites == 0:
        raise RuntimeError("no sites passed model input validation and truncation")

    fits = fit_all_ploidies(
        data.histogram,
        allele_min=config.allele_min,
        allele_max=config.allele_max,
        rho_bounds=(config.rho_min, config.rho_max),
        tetraploid_a_min=config.tetraploid_a_min,
        min_allele_reads=config.min_allele_reads,
    )
    ordered = sorted(fits.values(), key=lambda fit: fit.bic)
    best = ordered[0]
    delta_bic = ordered[1].bic - best.bic

    sensitivity_4x = fit_ploidy(
        data.histogram,
        4,
        config.allele_min,
        config.allele_max,
        (config.rho_min, config.rho_max),
        tetraploid_a_min=config.sensitivity_a_min,
        min_allele_reads=config.min_allele_reads,
    )
    sensitivity_fits = {2: fits[2], 3: fits[3], 4: sensitivity_4x}
    sensitivity_best = min(sensitivity_fits.values(), key=lambda fit: fit.bic)

    bootstrap_rows, votes = block_bootstrap(
        data.block_histograms,
        fits,
        replicates=config.bootstrap_replicates,
        allele_min=config.allele_min,
        allele_max=config.allele_max,
        rho_bounds=(config.rho_min, config.rho_max),
        tetraploid_a_min=config.tetraploid_a_min,
        seed=config.seed,
        min_allele_reads=config.min_allele_reads,
    )
    support = (
        votes.get(best.ploidy, 0) / config.bootstrap_replicates
        if config.bootstrap_replicates else 1.0
    )
    chromosome_rows, chromosome_votes = chromosome_model_support(
        data.block_histograms,
        allele_min=config.allele_min,
        allele_max=config.allele_max,
        rho_bounds=(config.rho_min, config.rho_max),
        tetraploid_a_min=config.tetraploid_a_min,
        min_sites=config.min_chromosome_sites,
        min_delta_bic=config.chromosome_delta_bic,
        min_allele_reads=config.min_allele_reads,
    )
    informative_chromosomes = sum(chromosome_votes.values())
    chromosome_agreement = (
        chromosome_votes.get(best.ploidy, 0) / informative_chromosomes
        if informative_chromosomes else 0.0
    )

    reasons: list[str] = []
    warnings: list[str] = []
    if data.selected_sites < config.min_sites:
        reasons.append("LOW_SITE_COUNT")
    if len(data.block_histograms) < config.min_blocks:
        warnings.append("TOO_FEW_BLOCKS")
    if delta_bic < config.delta_bic:
        reasons.append("LOW_DELTA_BIC")
    if config.bootstrap_replicates and support < config.bootstrap_support:
        warnings.append("LOW_BOOTSTRAP_SUPPORT")
    rho_tolerance = 0.01 * (config.rho_max - config.rho_min)
    if (
        best.rho <= config.rho_min + rho_tolerance
        or best.rho >= config.rho_max - rho_tolerance
    ):
        warnings.append("RHO_AT_BOUNDARY")
    if best.ploidy == 4 and best.outer_fraction is not None:
        a_tolerance = max(1e-4, 0.01 * (1 - config.tetraploid_a_min))
        if best.outer_fraction <= config.tetraploid_a_min + a_tolerance:
            warnings.append("A_AT_BOUNDARY")
        if best.outer_fraction >= 1 - 1e-3:
            warnings.append("A_AT_UPPER_BOUNDARY")
    if sensitivity_best.ploidy != best.ploidy:
        warnings.append("MODEL_ASSUMPTION_SENSITIVE")
    warnings.extend(chromosome_evidence_warnings(
        informative_chromosomes, chromosome_agreement,
        config.min_informative_chromosomes, config.chromosome_support,
    ))
    if not all(fit.converged for fit in fits.values()):
        reasons.append("MODEL_NOT_CONVERGED")

    predicted = "UNCERTAIN" if reasons else PREDICTION_LABELS[best.ploidy]
    qc_status = "FAIL" if reasons else ("WARNING" if warnings else "PASS")

    write_model_comparison(outputs["comparison"], fits, best.bic)
    write_sensitivity(
        outputs["sensitivity"], config, best, sensitivity_4x, sensitivity_best
    )
    write_bootstrap(outputs["bootstrap"], bootstrap_rows)
    write_block_summary(outputs["blocks"], data)
    write_chromosome_support(outputs["chromosomes"], chromosome_rows)
    result = {
        "sample": sample,
        "predicted_ploidy": predicted,
        "best_hypothesis": f"{best.ploidy}x",
        "best_model": MODEL_LABELS[best.ploidy],
        "bootstrap_support": support,
        "bootstrap_planned": config.bootstrap_replicates,
        "bootstrap_valid": sum(votes.values()),
        "bootstrap_unresolved": config.bootstrap_replicates - sum(votes.values()),
        "delta_bic": delta_bic,
        "selected_sites": data.selected_sites,
        "available_sites": data.available_sites,
        "blocks": len(data.block_histograms),
        "informative_chromosomes": informative_chromosomes,
        "chromosome_agreement": chromosome_agreement,
        "sensitivity_best_hypothesis": f"{sensitivity_best.ploidy}x",
        "qc_status": qc_status,
        "uncertainty_reasons": "NONE" if not reasons else ";".join(reasons),
        "warnings": "NONE" if not warnings else ";".join(warnings),
    }
    write_result(outputs["result"], result, votes, config.bootstrap_replicates)
    profiles = {p: fitted_frequency_profile(
        data.histogram, fits[p], config.allele_min, config.allele_max,
        config.min_allele_reads) for p in PLOIDIES}
    observed, _, edges = profiles[best.ploidy]
    manifest = {
        "frequency_profile": {"edges": edges, "observed": observed,
                              "fitted_models": {str(p): profiles[p][1] for p in PLOIDIES},
                              "sites": data.selected_sites,
                              "definition": "minor count / sum of the two leading allele counts",
                              "model": f"{best.ploidy}x"},
        "workflow": "ploidy-model",
        "model_version": MODEL_VERSION,
        "mapploid_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "run_config": signature,
        "result": result,
        "fits": {str(ploidy): _fit_dict(fit) for ploidy, fit in fits.items()},
        "sensitivity_4x_fit": _fit_dict(sensitivity_4x),
        "bootstrap_votes": {f"{p}x": votes.get(p, 0) for p in PLOIDIES},
        "chromosome_votes": {
            f"{p}x": chromosome_votes.get(p, 0) for p in PLOIDIES
        },
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
    outputs["manifest"].write_text(json.dumps(manifest, indent=2) + "\n",
                                   encoding="utf-8")
    return outputs


def retained_count_bounds(depth, allele_min: float, allele_max: float,
                          min_allele_reads: int = MIN_ALLELE_READS):
    """Common support for both allele filters, MAC, likelihood and simulation."""
    lower = max(allele_min, 1 - allele_max)
    upper = min(allele_max, 1 - allele_min)
    if isinstance(depth, (int, np.integer)):
        return (max(min_allele_reads, math.ceil(lower * depth - 1e-12)),
                min(depth - min_allele_reads, math.floor(upper * depth + 1e-12)))
    n = np.asarray(depth)
    return (np.maximum(min_allele_reads, np.ceil(lower * n - 1e-12).astype(int)),
            np.minimum(n - min_allele_reads, np.floor(upper * n + 1e-12).astype(int)))


def chromosome_evidence_warnings(informative, agreement, minimum=3, support=0.70):
    """Report chromosome diagnostics without vetoing the global BIC call."""
    warnings = []
    if informative < minimum:
        warnings.append("INSUFFICIENT_CHROMOSOME_EVIDENCE")
    if informative and agreement < support:
        warnings.append("CHROMOSOME_HETEROGENEITY")
    return warnings


def validate_counts_contract(counts: Path, config: PloidyModelConfig):
    """Check the producer contract when importing a pipeline count table."""
    manifest = counts.parent / "run_config.json"
    if manifest.exists():
        saved = json.loads(manifest.read_text())
        expected = {"allele_min": config.allele_min, "allele_max": config.allele_max,
                    "min_allele_reads": config.min_allele_reads, "pair_purity": config.pair_purity}
        for key, value in expected.items():
            if key not in saved or saved[key] != value:
                raise ValueError(f"count/model filter mismatch: {key}")


def load_and_aggregate_sites(
    counts: Path,
    sample: str,
    max_sites: int,
    allele_min: float,
    allele_max: float,
    block_size: int,
    seed: int,
    min_allele_reads: int = MIN_ALLELE_READS,
    pair_purity: float = 0.95,
) -> AggregatedSites:
    """Deterministically sample sites across chromosomes and aggregate counts."""

    chromosome_available: Counter[str] = Counter()
    for site in _iter_model_sites(counts, sample, allele_min, allele_max, seed,
                                  min_allele_reads, pair_purity):
        chromosome_available[site.chrom] += 1
    total = sum(chromosome_available.values())
    if total == 0:
        return AggregatedSites(Counter(), {}, Counter(), 0, 0)
    target = min(max_sites, total)
    quotas = proportional_quotas(chromosome_available, target)

    heaps: dict[str, list[tuple[int, int, int, int]]] = {
        chrom: [] for chrom, quota in quotas.items() if quota
    }
    for site in _iter_model_sites(counts, sample, allele_min, allele_max, seed,
                                  min_allele_reads, pair_purity):
        quota = quotas.get(site.chrom, 0)
        if quota == 0:
            continue
        token = _site_hash(sample, site.chrom, site.position, seed, "site-sample")
        item = (-token, site.position, site.depth, site.allele_count)
        heap = heaps[site.chrom]
        if len(heap) < quota:
            heapq.heappush(heap, item)
        elif token < -heap[0][0]:
            heapq.heapreplace(heap, item)

    histogram: Counter[tuple[int, int]] = Counter()
    block_histograms: dict[str, Counter[tuple[int, int]]] = {}
    chromosome_counts: Counter[str] = Counter()
    for chrom, heap in heaps.items():
        for _, position, depth, allele_count in heap:
            key = (depth, allele_count)
            histogram[key] += 1
            chromosome_counts[chrom] += 1
            block_start = ((position - 1) // block_size) * block_size + 1
            block_end = block_start + block_size - 1
            block = f"{chrom}:{block_start}-{block_end}"
            block_histograms.setdefault(block, Counter())[key] += 1
    return AggregatedSites(
        histogram=histogram,
        block_histograms=block_histograms,
        chromosome_counts=chromosome_counts,
        available_sites=total,
        selected_sites=sum(histogram.values()),
    )


def proportional_quotas(counts: Counter[str], target: int) -> dict[str, int]:
    """Allocate an exact deterministic sample size proportional to chromosomes."""

    total = sum(counts.values())
    if target >= total:
        return dict(counts)
    raw = {chrom: target * count / total for chrom, count in counts.items()}
    quotas = {chrom: min(counts[chrom], int(math.floor(value)))
              for chrom, value in raw.items()}
    remaining = target - sum(quotas.values())
    order = sorted(counts, key=lambda chrom: (-(raw[chrom] - quotas[chrom]), chrom))
    for chrom in order:
        if remaining == 0:
            break
        if quotas[chrom] < counts[chrom]:
            quotas[chrom] += 1
            remaining -= 1
    return quotas


def fit_all_ploidies(
    histogram: Counter[tuple[int, int]],
    allele_min: float,
    allele_max: float,
    rho_bounds: tuple[float, float],
    tetraploid_a_min: float = DEFAULT_TETRAPLOID_A_MIN,
    initial_fits: dict[int, ModelFit] | None = None,
    min_allele_reads: int = MIN_ALLELE_READS,
) -> dict[int, ModelFit]:
    return {
        ploidy: fit_ploidy(
            histogram,
            ploidy,
            allele_min,
            allele_max,
            rho_bounds,
            tetraploid_a_min,
            initial_fits.get(ploidy) if initial_fits else None,
            min_allele_reads=min_allele_reads,
        )
        for ploidy in PLOIDIES
    }


def fit_ploidy(
    histogram: Counter[tuple[int, int]],
    ploidy: int,
    allele_min: float,
    allele_max: float,
    rho_bounds: tuple[float, float],
    tetraploid_a_min: float = DEFAULT_TETRAPLOID_A_MIN,
    initial_fit: ModelFit | None = None,
    min_allele_reads: int = MIN_ALLELE_READS,
) -> ModelFit:
    """Maximum-likelihood fit for one fixed-ploidy hypothesis."""

    if ploidy not in PLOIDIES:
        raise ValueError(f"unsupported ploidy: {ploidy}")
    depth, allele_count, frequency = histogram_arrays(histogram)
    observations = int(frequency.sum())
    bounds = [rho_bounds]
    if ploidy == 4:
        bounds.append((tetraploid_a_min, 1 - 1e-4))

    if initial_fit is not None:
        initial = [initial_fit.rho]
        if ploidy == 4:
            initial.append(max(tetraploid_a_min, initial_fit.outer_fraction or 0.65))
        starts = [initial]
    else:
        mid_rho = math.sqrt(rho_bounds[0] * rho_bounds[1])
        starts = [
            [mid_rho] + ([tetraploid_a_min] if ploidy == 4 else []),
            [min(0.02, rho_bounds[1])]
            + ([max(tetraploid_a_min, 0.65)] if ploidy == 4 else []),
            [min(0.08, rho_bounds[1])]
            + ([max(tetraploid_a_min, 0.90)] if ploidy == 4 else []),
        ]

    def objective(parameters: np.ndarray) -> float:
        return -truncated_mixture_log_likelihood(
            depth,
            allele_count,
            frequency,
            ploidy,
            rho=float(parameters[0]),
            outer_fraction=float(parameters[1]) if ploidy == 4 else None,
            allele_min=allele_min,
            allele_max=allele_max,
            min_allele_reads=min_allele_reads,
        )

    results = [
        minimize(objective, np.asarray(start, dtype=float), method="L-BFGS-B",
                 bounds=bounds, options={"maxiter": 250, "ftol": 1e-10})
        for start in starts
    ]
    # A failed warm start receives the same finite, predetermined retry set.
    if initial_fit is not None and not any(item.success and math.isfinite(float(item.fun)) for item in results):
        return fit_ploidy(histogram, ploidy, allele_min, allele_max, rho_bounds,
                          tetraploid_a_min, initial_fit=None, min_allele_reads=min_allele_reads)
    successful = [item for item in results if item.success and math.isfinite(float(item.fun))]
    result = min(successful or results,
                 key=lambda item: float(item.fun) if math.isfinite(float(item.fun)) else math.inf)
    log_likelihood = -float(result.fun)
    parameters = result.x
    parameter_count = 2 if ploidy == 4 else 1
    bic = -2 * log_likelihood + parameter_count * math.log(observations)
    return ModelFit(
        ploidy=ploidy,
        log_likelihood=log_likelihood,
        bic=bic,
        rho=float(parameters[0]),
        outer_fraction=float(parameters[1]) if ploidy == 4 else None,
        converged=bool(result.success and math.isfinite(log_likelihood)),
        iterations=int(result.nit),
    )


def truncated_mixture_log_likelihood(
    depth: np.ndarray,
    allele_count: np.ndarray,
    frequency: np.ndarray,
    ploidy: int,
    rho: float,
    outer_fraction: float | None,
    allele_min: float,
    allele_max: float,
    min_allele_reads: int = MIN_ALLELE_READS,
) -> float:
    """Weighted log likelihood, conditional on allele-frequency truncation."""

    means = EXPECTED_MEANS[ploidy]
    if ploidy == 2:
        weights = np.asarray((1.0,))
    elif ploidy == 3:
        weights = np.asarray((0.5, 0.5))
    else:
        outer = DEFAULT_TETRAPLOID_A_MIN if outer_fraction is None else outer_fraction
        weights = np.asarray((outer / 2, 1 - outer, outer / 2))

    component_logs = []
    unique_depth, inverse = np.unique(depth, return_inverse=True)
    low, high = retained_count_bounds(unique_depth, allele_min, allele_max, min_allele_reads)
    row_low, row_high = retained_count_bounds(depth, allele_min, allele_max, min_allele_reads)
    if np.any(allele_count < row_low) or np.any(allele_count > row_high):
        return -math.inf
    for mean in means:
        alpha, beta = beta_parameters(mean, rho)
        raw = beta_binomial_logpmf(allele_count, depth, alpha, beta)
        normalizer = truncation_log_normalizers(
            unique_depth, low, high, alpha, beta
        )
        component_logs.append(raw - normalizer[inverse])
    stacked = np.vstack(component_logs)
    mixture = logsumexp(stacked + np.log(weights)[:, None], axis=0)
    return float(np.dot(frequency, mixture))


def fitted_frequency_profile(histogram, fit, allele_min, allele_max,
                             min_allele_reads=MIN_ALLELE_READS, bin_width=DEFAULT_FIT_BIN_WIDTH):
    """Folded minor-allele frequency probabilities for the exact fitted site set."""
    lower = max(0.0, min(allele_min, 1 - allele_max))
    bins = max(1, int(math.ceil((0.5 - lower) / bin_width)))
    observed = np.zeros(bins, dtype=float)
    predicted = np.zeros(bins, dtype=float)
    depth_counts: Counter[int] = Counter()

    def index(value: float) -> int:
        return min(bins - 1, max(0, int((value - lower) / bin_width)))

    for (depth, count), frequency in histogram.items():
        observed[index(min(count, depth - count) / depth)] += frequency
        depth_counts[depth] += frequency

    means = EXPECTED_MEANS[fit.ploidy]
    weights = component_weights(fit.ploidy, fit.outer_fraction)
    for depth, depth_frequency in depth_counts.items():
        low, high = (int(x) for x in retained_count_bounds(depth, allele_min, allele_max, min_allele_reads))
        counts = np.arange(low, high + 1, dtype=int)
        depths = np.full(counts.shape, depth, dtype=int)
        component_logs = []
        for mean in means:
            alpha, beta = beta_parameters(mean, fit.rho)
            log_probability = beta_binomial_logpmf(
                counts, depths, alpha, beta
            )
            component_logs.append(log_probability - logsumexp(log_probability))
        probabilities = np.exp(
            logsumexp(
                np.vstack(component_logs) + np.log(weights)[:, None], axis=0
            )
        )
        for count, probability in zip(counts, probabilities):
            minor_baf = min(int(count), depth - int(count)) / depth
            predicted[index(minor_baf)] += depth_frequency * probability

    if observed.sum() == 0 or predicted.sum() == 0:
        return observed.tolist(), predicted.tolist(), []
    observed /= observed.sum()
    predicted /= predicted.sum()
    edges = [lower + i * bin_width for i in range(bins)] + [0.5]
    return observed.tolist(), predicted.tolist(), edges


def component_weights(
    ploidy: int, outer_fraction: float | None = None
) -> np.ndarray:
    if ploidy == 2:
        return np.asarray((1.0,))
    if ploidy == 3:
        return np.asarray((0.5, 0.5))
    outer = DEFAULT_TETRAPLOID_A_MIN if outer_fraction is None else outer_fraction
    return np.asarray((outer / 2, 1 - outer, outer / 2))


def truncation_log_normalizers(
    unique_depth: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
    alpha: float,
    beta: float,
) -> np.ndarray:
    """Exact beta-binomial mass inside the retained count interval.

    Directly summing short PMF vectors is much faster than repeatedly calling
    the generic SciPy beta-binomial CDF during numerical optimization.
    """

    maximum = int(unique_depth.max())
    if maximum <= 5_000:
        counts = np.arange(maximum + 1, dtype=int)[None, :]
        depths = unique_depth[:, None]
        valid = (
            (counts >= low[:, None])
            & (counts <= high[:, None])
            & (counts <= depths)
        )
        safe_counts = np.minimum(counts, depths)
        log_probabilities = beta_binomial_logpmf(
            safe_counts, depths, alpha, beta
        )
        return logsumexp(np.where(valid, log_probabilities, -np.inf), axis=1)

    # Avoid a potentially huge matrix for pathological ultra-high-depth data.
    values = np.empty(len(unique_depth), dtype=float)
    for index, depth in enumerate(unique_depth):
        retained = np.arange(low[index], high[index] + 1, dtype=int)
        depths = np.full(retained.shape, int(depth), dtype=int)
        values[index] = logsumexp(
            beta_binomial_logpmf(retained, depths, alpha, beta)
        )
    return values


def beta_parameters(mean: float, rho: float) -> tuple[float, float]:
    concentration = (1 - rho) / rho
    return mean * concentration, (1 - mean) * concentration


def beta_binomial_logpmf(
    count: np.ndarray, depth: np.ndarray, alpha: float, beta: float
) -> np.ndarray:
    return (
        gammaln(depth + 1)
        - gammaln(count + 1)
        - gammaln(depth - count + 1)
        + betaln(count + alpha, depth - count + beta)
        - betaln(alpha, beta)
    )


def histogram_arrays(
    histogram: Counter[tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = sorted(histogram.items())
    if not rows:
        raise ValueError("cannot fit an empty histogram")
    depth = np.asarray([key[0] for key, _ in rows], dtype=int)
    count = np.asarray([key[1] for key, _ in rows], dtype=int)
    frequency = np.asarray([value for _, value in rows], dtype=float)
    return depth, count, frequency


def block_bootstrap(
    block_histograms: dict[str, Counter[tuple[int, int]]],
    global_fits: dict[int, ModelFit],
    replicates: int,
    allele_min: float,
    allele_max: float,
    rho_bounds: tuple[float, float],
    tetraploid_a_min: float,
    seed: int,
    min_allele_reads: int = MIN_ALLELE_READS,
) -> tuple[list[dict[str, object]], Counter[int]]:
    """Resample physical chromosome blocks and refit all three hypotheses."""

    if replicates == 0:
        return [], Counter()
    names = sorted(block_histograms)
    if not names:
        raise ValueError("no blocks are available for bootstrap")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    votes: Counter[int] = Counter()
    for replicate in range(1, replicates + 1):
        selected = rng.integers(0, len(names), size=len(names))
        histogram: Counter[tuple[int, int]] = Counter()
        for index in selected:
            histogram.update(block_histograms[names[int(index)]])
        fits = fit_all_ploidies(
            histogram,
            allele_min,
            allele_max,
            rho_bounds,
            tetraploid_a_min,
            initial_fits=global_fits,
            min_allele_reads=min_allele_reads,
        )
        valid = all(fit.converged and math.isfinite(fit.bic) for fit in fits.values())
        winner = min(fits.values(), key=lambda fit: fit.bic).ploidy if valid else None
        if valid:
            votes[winner] += 1
        rows.append({
            "replicate": replicate,
            "status": "PASS" if valid else "UNRESOLVED",
            "sites": sum(histogram.values()),
            "best_ploidy": winner,
            **{f"bic_{ploidy}x": fits[ploidy].bic for ploidy in PLOIDIES},
        })
    return rows, votes


def chromosome_model_support(
    block_histograms: dict[str, Counter[tuple[int, int]]],
    allele_min: float,
    allele_max: float,
    rho_bounds: tuple[float, float],
    tetraploid_a_min: float,
    min_sites: int,
    min_delta_bic: float,
    min_allele_reads: int = MIN_ALLELE_READS,
) -> tuple[list[dict[str, object]], Counter[int]]:
    """Fit the primary models per chromosome to expose regional conflict."""

    chromosome_histograms: dict[str, Counter[tuple[int, int]]] = {}
    for block, histogram in block_histograms.items():
        chrom = block.split(":", 1)[0]
        chromosome_histograms.setdefault(chrom, Counter()).update(histogram)
    rows: list[dict[str, object]] = []
    votes: Counter[int] = Counter()
    for chrom in sorted(chromosome_histograms):
        histogram = chromosome_histograms[chrom]
        sites = sum(histogram.values())
        if sites < min_sites:
            rows.append({
                "chrom": chrom,
                "sites": sites,
                "best_ploidy": None,
                "delta_bic": math.nan,
                "bic_2x": math.nan,
                "bic_3x": math.nan,
                "bic_4x": math.nan,
                "status": "TOO_FEW_SITES",
            })
            continue
        fits = fit_all_ploidies(
            histogram,
            allele_min,
            allele_max,
            rho_bounds,
            tetraploid_a_min,
            min_allele_reads=min_allele_reads,
        )
        ordered = sorted(fits.values(), key=lambda fit: fit.bic)
        winner = ordered[0].ploidy
        delta_bic = ordered[1].bic - ordered[0].bic
        valid = all(fit.converged and math.isfinite(fit.bic) for fit in fits.values())
        status = ("UNRESOLVED" if not valid else
                  "PASS" if delta_bic >= min_delta_bic else "LOW_DELTA_BIC")
        if status == "PASS":
            votes[winner] += 1
        rows.append({
            "chrom": chrom,
            "sites": sites,
            "best_ploidy": winner,
            "delta_bic": delta_bic,
            **{f"bic_{ploidy}x": fits[ploidy].bic for ploidy in PLOIDIES},
            "status": status,
        })
    return rows, votes


def simulate_allele_counts(
    output: Path,
    ploidy: int,
    sites: int = 20_000,
    mean_depth: float = 30,
    rho: float = 0.02,
    contamination_fraction: float = 0.0,
    outer_fraction: float = 0.65,
    allele_min: float = 0.15,
    allele_max: float = 0.85,
    chromosomes: int = 10,
    seed: int = 20260803,
    min_allele_reads: int = MIN_ALLELE_READS,
) -> Path:
    """Simulate a filtered biallelic count table from the inference model."""

    if ploidy not in PLOIDIES:
        raise ValueError("simulator supports only 2x, 3x and 4x")
    if not 0 <= allele_min < 0.5 < allele_max <= 1:
        raise ValueError("allele truncation must contain 0.5")
    if sites < 1 or mean_depth <= 0 or chromosomes < 1:
        raise ValueError("sites, mean depth and chromosome count must be positive")
    if min_allele_reads < 1:
        raise ValueError("minimum allele reads must be positive")
    if not 0 <= contamination_fraction < 1:
        raise ValueError("simulation contamination fraction must be in [0, 1)")
    if ploidy == 4 and not DEFAULT_TETRAPLOID_A_MIN <= outer_fraction < 1:
        raise ValueError("4x outer fraction must be in [4/7, 1)")
    rng = np.random.default_rng(seed + ploidy)
    means = np.asarray(EXPECTED_MEANS[ploidy])
    if ploidy == 2:
        weights = np.asarray((1.0,))
    elif ploidy == 3:
        weights = np.asarray((0.5, 0.5))
    else:
        weights = np.asarray((outer_fraction / 2, 1 - outer_fraction,
                              outer_fraction / 2))

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{output}.tmp")
    with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=1) as handle:
        handle.write(
            "chrom\tposition\tdepth\tpair_depth\tpair_purity\t"
            "allele1_count\tallele2_count\n"
        )
        for index in range(sites):
            depth = max(2 * min_allele_reads, int(rng.poisson(mean_depth)))
            low, high = retained_count_bounds(depth, allele_min, allele_max, min_allele_reads)
            while low > high:
                depth += 1
                low, high = retained_count_bounds(depth, allele_min, allele_max, min_allele_reads)
            if rng.random() < contamination_fraction:
                count = int(rng.integers(low, high + 1))
            else:
                mean = float(rng.choice(means, p=weights))
                alpha, beta = beta_parameters(mean, rho)
                while True:
                    probability = float(rng.beta(alpha, beta))
                    count = int(rng.binomial(depth, probability))
                    if low <= count <= high:
                        break
            chrom_index = index % chromosomes + 1
            position = (index // chromosomes) * 10_000 + 1
            handle.write(
                f"sim_chr_{chrom_index}\t{position}\t{depth}\t{depth}\t1.0\t"
                f"{count}\t{depth-count}\n"
            )
    temporary.replace(output)
    return output


def simulate_and_validate(
    output_dir: Path,
    sites: int = 20_000,
    mean_depth: float = 30,
    rho: float = 0.02,
    contamination_fraction: float = 0.0,
    outer_fraction: float = 0.65,
    seed: int = 20260803,
) -> dict[str, Path]:
    """Generate 2x/3x/4x data and verify the best BIC hypothesis."""

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = output_dir / "simulation_validation.tsv"
    expected = output_dir / "expected_peaks.tsv"
    peak_validation = output_dir / "peak_validation.tsv"
    with expected.open("wt", encoding="utf-8") as handle:
        handle.write("ploidy\texpected_allele_frequencies\n")
        for ploidy in PLOIDIES:
            handle.write(f"{ploidy}x\t{','.join(f'{value:.6f}' for value in EXPECTED_MEANS[ploidy])}\n")

    paths: dict[str, Path] = {
        "summary": summary,
        "expected_peaks": expected,
        "peak_validation": peak_validation,
    }
    rows = []
    peak_rows: list[tuple[int, float, float, int, float]] = []
    for ploidy in PLOIDIES:
        count_path = output_dir / f"simulated_{ploidy}x.tsv.gz"
        simulate_allele_counts(
            count_path, ploidy, sites, mean_depth, rho, contamination_fraction,
            outer_fraction, seed=seed,
        )
        data = load_and_aggregate_sites(
            count_path, f"simulated_{ploidy}x", sites, 0.15, 0.85,
            5_000_000, seed,
        )
        fits = fit_all_ploidies(
            data.histogram,
            0.15,
            0.85,
            (1e-5, 0.25),
            DEFAULT_TETRAPLOID_A_MIN,
        )
        ordered = sorted(fits.values(), key=lambda fit: fit.bic)
        rows.append((ploidy, ordered[0].ploidy, ordered[1].bic - ordered[0].bic,
                     data.selected_sites))
        peak_rows.extend(summarize_simulated_peaks(count_path, ploidy))
        paths[f"simulated_{ploidy}x"] = count_path
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("true_ploidy\tpredicted_ploidy\tdelta_bic\tsites\tstatus\n")
        for true, predicted, delta, selected in rows:
            status = "PASS" if true == predicted else "FAIL"
            handle.write(f"{true}x\t{predicted}x\t{delta:.6f}\t{selected}\t{status}\n")
    with peak_validation.open("wt", encoding="utf-8") as handle:
        handle.write(
            "ploidy\texpected_peak\tobserved_mean\tassigned_sites\t"
            "absolute_error\tstatus\n"
        )
        for ploidy, expected_mean, observed, assigned, error in peak_rows:
            status = "PASS" if error <= 0.03 else "FAIL"
            handle.write(
                f"{ploidy}x\t{expected_mean:.6f}\t{observed:.6f}\t{assigned}\t"
                f"{error:.6f}\t{status}\n"
            )
    return paths


def summarize_simulated_peaks(
    path: Path, ploidy: int
) -> list[tuple[int, float, float, int, float]]:
    """Assign simulated frequencies to the nearest theoretical dosage peak."""

    means = EXPECTED_MEANS[ploidy]
    sums = [0.0] * len(means)
    counts = [0] * len(means)
    with _open_text(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        indexes = {name: index for index, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            depth = int(fields[indexes.get("pair_depth", indexes["depth"])])
            value = int(fields[indexes["allele1_count"]]) / depth
            nearest = min(range(len(means)), key=lambda index: abs(value - means[index]))
            sums[nearest] += value
            counts[nearest] += 1
    rows = []
    for index, expected in enumerate(means):
        observed = sums[index] / counts[index] if counts[index] else math.nan
        error = abs(observed - expected) if counts[index] else math.inf
        rows.append((ploidy, expected, observed, counts[index], error))
    return rows


def _iter_model_sites(
    path: Path,
    sample: str,
    allele_min: float,
    allele_max: float,
    seed: int,
    min_allele_reads: int = MIN_ALLELE_READS,
    pair_purity: float = 0.95,
):
    with _open_text(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        indexes = {name: index for index, name in enumerate(header)}
        required = (
            "chrom", "position", "pair_depth", "pair_purity",
            "allele1_count", "allele2_count",
        )
        missing = [name for name in required if name not in indexes]
        if missing:
            raise ValueError(f"count table is missing columns: {', '.join(missing)}")
        last_position = {}
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            try:
                position = int(fields[indexes["position"]])
                depth = int(fields[indexes["pair_depth"]])
                purity = float(fields[indexes["pair_purity"]])
                first = int(fields[indexes["allele1_count"]])
                second = int(fields[indexes["allele2_count"]])
                total_depth = int(fields[indexes["depth"]]) if "depth" in indexes else None
            except (ValueError, IndexError) as exc:
                raise ValueError(f"invalid count row at {path}:{line_number}") from exc
            if (
                depth <= 0
                or position <= 0
                or min(first, second) < 0
                or first + second != depth
                or not 0.5 <= purity <= 1
            ):
                raise ValueError(f"invalid model counts at {path}:{line_number}")
            chrom = fields[indexes["chrom"]]
            if not chrom or position <= last_position.get(chrom, 0):
                raise ValueError(f"duplicate or unordered coordinate at {path}:{line_number}")
            last_position[chrom] = position
            # Count tables round purity for display; use exact counts if available.
            if total_depth is not None:
                if total_depth < depth:
                    raise ValueError(f"invalid total depth at {path}:{line_number}")
                purity = depth / total_depth
            if min(first, second) < min_allele_reads or purity < pair_purity:
                continue
            token = _site_hash(sample, chrom, position, seed,
                               "allele-choice")
            count = (first, second)[token & 1]
            if count > depth:
                continue
            fraction = count / depth
            low, high = retained_count_bounds(depth, allele_min, allele_max, min_allele_reads)
            if low <= count <= high:
                yield ModelSite(fields[indexes["chrom"]], position, depth, count)


def _site_hash(sample: str, chrom: str, position: int, seed: int, purpose: str) -> int:
    token = f"{purpose}|{seed}|{sample}|{chrom}|{position}".encode()
    return int.from_bytes(hashlib.blake2b(token, digest_size=8).digest(), "big")


def write_model_comparison(
    output: Path, fits: dict[int, ModelFit], best_bic: float
) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write(
            "model\tploidy\tlog_likelihood\tbic\tdelta_bic\trho\t"
            "outer_fraction\tconverged\titerations\n"
        )
        for ploidy in PLOIDIES:
            fit = fits[ploidy]
            outer = "NA" if fit.outer_fraction is None else f"{fit.outer_fraction:.8f}"
            handle.write(
                f"{MODEL_LABELS[ploidy]}\t{ploidy}x\t{fit.log_likelihood:.8f}\t"
                f"{fit.bic:.8f}\t"
                f"{fit.bic-best_bic:.8f}\t{fit.rho:.8f}\t"
                f"{outer}\t"
                f"{str(fit.converged).lower()}\t{fit.iterations}\n"
            )


def write_sensitivity(
    output: Path,
    config: PloidyModelConfig,
    primary_best: ModelFit,
    sensitivity_4x: ModelFit,
    sensitivity_best: ModelFit,
) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"primary_a_min\t{config.tetraploid_a_min:.8f}\n")
        handle.write(f"sensitivity_a_min\t{config.sensitivity_a_min:.8f}\n")
        handle.write(f"primary_best_model\t{MODEL_LABELS[primary_best.ploidy]}\n")
        handle.write(f"sensitivity_best_model\t{MODEL_LABELS[sensitivity_best.ploidy]}\n")
        handle.write(f"sensitivity_4x_bic\t{sensitivity_4x.bic:.8f}\n")
        handle.write(
            "sensitivity_4x_outer_fraction\t"
            f"{sensitivity_4x.outer_fraction:.8f}\n"
        )


def write_bootstrap(output: Path, rows: list[dict[str, object]]) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("replicate\tsites\tbest_ploidy\tbic_2x\tbic_3x\tbic_4x\tstatus\n")
        for row in rows:
            handle.write(
                f"{row['replicate']}\t{row['sites']}\t"
                + (f"{row['best_ploidy']}x" if row['best_ploidy'] else "NA") + "\t"
                f"{float(row['bic_2x']):.8f}\t{float(row['bic_3x']):.8f}\t"
                f"{float(row['bic_4x']):.8f}\t{row['status']}\n"
            )


def write_block_summary(output: Path, data: AggregatedSites) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write("block\tchrom\tselected_sites\n")
        for block in sorted(data.block_histograms):
            chrom = block.split(":", 1)[0]
            handle.write(f"{block}\t{chrom}\t{sum(data.block_histograms[block].values())}\n")


def write_chromosome_support(
    output: Path, rows: list[dict[str, object]]
) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write(
            "chrom\tsites\tbest_model\tbest_hypothesis\tdelta_bic\t"
            "bic_2x\tbic_3x\tbic_4x\tstatus\n"
        )
        for row in rows:
            ploidy = row["best_ploidy"]
            model = "NA" if ploidy is None else MODEL_LABELS[int(ploidy)]
            hypothesis = "NA" if ploidy is None else f"{int(ploidy)}x"
            values = []
            for key in ("delta_bic", "bic_2x", "bic_3x", "bic_4x"):
                value = float(row[key])
                values.append("NA" if math.isnan(value) else f"{value:.8f}")
            handle.write(
                f"{row['chrom']}\t{row['sites']}\t{model}\t{hypothesis}\t"
                + "\t".join(values)
                + f"\t{row['status']}\n"
            )


def write_result(
    output: Path,
    result: dict[str, object],
    votes: Counter[int],
    replicates: int,
) -> None:
    with output.open("wt", encoding="utf-8") as handle:
        handle.write(
            "sample\tpredicted_ploidy\tbest_model\tbest_hypothesis\t"
            "bootstrap_support\tdelta_bic\tselected_sites\tavailable_sites\t"
            "blocks\tinformative_chromosomes\t"
            "chromosome_agreement\tsensitivity_best_hypothesis\t"
            "bootstrap_2x\tbootstrap_3x\t"
            "bootstrap_4x\tqc_status\tuncertainty_reasons\twarnings\n"
        )
        support_values = [votes.get(ploidy, 0) / replicates if replicates else 0
                          for ploidy in PLOIDIES]
        handle.write(
            f"{result['sample']}\t{result['predicted_ploidy']}\t"
            f"{result['best_model']}\t{result['best_hypothesis']}\t"
            f"{float(result['bootstrap_support']):.6f}\t"
            f"{float(result['delta_bic']):.6f}\t{result['selected_sites']}\t"
            f"{result['available_sites']}\t"
            f"{result['blocks']}\t"
            f"{result['informative_chromosomes']}\t"
            f"{float(result['chromosome_agreement']):.6f}\t"
            f"{result['sensitivity_best_hypothesis']}\t"
            + "\t".join(f"{value:.6f}" for value in support_values)
            + f"\t{result['qc_status']}\t{result['uncertainty_reasons']}\t"
            + f"{result['warnings']}\n"
        )


def _fit_dict(fit: ModelFit) -> dict[str, object]:
    return {
        "ploidy": fit.ploidy,
        "log_likelihood": fit.log_likelihood,
        "bic": fit.bic,
        "rho": fit.rho,
        "outer_fraction": fit.outer_fraction,
        "converged": fit.converged,
        "iterations": fit.iterations,
    }


def _sample_from_counts(path: Path) -> str:
    name = path.name
    for suffix in (".biallelic_sites.tsv.gz", ".tsv.gz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("rt", encoding="utf-8")


def _file_stamp(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _signature(
    config: PloidyModelConfig,
    counts: Path,
    sample: str,
) -> dict[str, object]:
    return {
        "model_version": MODEL_VERSION,
        "sample": sample,
        "counts": _file_stamp(counts),
        "max_sites": config.max_sites,
        "min_sites": config.min_sites,
        "allele_min": config.allele_min,
        "allele_max": config.allele_max,
        "rho_min": config.rho_min,
        "rho_max": config.rho_max,
        "tetraploid_a_min": config.tetraploid_a_min,
        "sensitivity_a_min": config.sensitivity_a_min,
        "block_size": config.block_size,
        "bootstrap_replicates": config.bootstrap_replicates,
        "bootstrap_support": config.bootstrap_support,
        "delta_bic": config.delta_bic,
        "min_blocks": config.min_blocks,
        "chromosome_support": config.chromosome_support,
        "min_chromosome_sites": config.min_chromosome_sites,
        "chromosome_delta_bic": config.chromosome_delta_bic,
        "min_informative_chromosomes": config.min_informative_chromosomes,
        "min_allele_reads": config.min_allele_reads,
        "pair_purity": config.pair_purity,
        "seed": config.seed,
    }


def _check_or_write_signature(path: Path, signature: dict[str, object]) -> None:
    if path.is_file():
        previous = json.loads(path.read_text(encoding="utf-8"))
        if previous != signature:
            raise ValueError(f"output directory has different inputs or parameters: {path.parent}")
        return
    path.write_text(json.dumps(signature, indent=2) + "\n", encoding="utf-8")


def _output_paths(output_dir: Path, sample: str) -> dict[str, Path]:
    return {
        "comparison": output_dir / f"{sample}.model_comparison.tsv",
        "sensitivity": output_dir / f"{sample}.model_sensitivity.tsv",
        "bootstrap": output_dir / f"{sample}.bootstrap.tsv",
        "blocks": output_dir / f"{sample}.block_summary.tsv",
        "chromosomes": output_dir / f"{sample}.chromosome_support.tsv",
        "result": output_dir / f"{sample}.ploidy_result.tsv",
        "manifest": output_dir / "run_manifest.json",
    }
