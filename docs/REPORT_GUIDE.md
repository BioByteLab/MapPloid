# Report and diagnostic guide

**English** | [简体中文](REPORT_GUIDE.zh-CN.md)

[Main documentation](../README.md) · [Command-line reference](CLI.md)

Open `report.html` in a browser; it works offline. The report, plot labels, logs and machine-readable identifiers use English. This guide provides the same explanations in English and Chinese. Read `result.tsv` for downstream analysis and `details/run.json` for complete diagnostic values.

## Start with the result

| Label | Meaning |
|---|---|
| CALLED | The three decision requirements passed; the reported ploidy is 2x, 3x or 4x |
| UNCERTAIN | No definitive call; any displayed best candidate remains a candidate |
| FAILED | Runtime failure; inspect `details/run.log` and the reason in `details/run.json` |
| Best candidate | Model with the best available BIC; not necessarily a definitive call |
| Analyzed region (Mb) | Final region size, not genome size |
| Available qualifying SNPs | Available sites that passed site filtering |
| Fitting SNPs | Sites actually used by the main fit, after any maximum-site sampling |
| Not calculated / NOT_CALCULATED | The metric was not calculated; do not interpret it as zero |
| Not assessable | Insufficient evidence to evaluate the metric, such as no informative chromosomes |

The main TSV contains `sample`, `predicted_ploidy`, `status`, `selected_sites`, `delta_bic` and `analyzed_region_mb`. The batch TSV adds `result_file`. Missing values are `NA`; non-finite JSON numbers are `null`. Reasons and warnings are under `result` in the sample's `details/run.json`.

## Decision requirements

| Report label | Requirement |
|---|---|
| Fitting SNPs | At least 25,000 SNPs used in fitting |
| BIC separation | Second-best BIC minus best BIC at least 10 |
| Main-model convergence | All three main models (2x, 3x, 4x) converged |
| PASS / FAIL | The requirement passed / failed |

The analysis expands from about 20 Mb to 50 Mb and then 100 Mb when the requirements fail, subject to available regions. Runtime failures differ from inconclusive analyses. None of the diagnostic warnings below adds another decision gate.

## Diagnostic evidence

| Report label | Interpretation |
|---|---|
| Bootstrap stability | Agreement with the global candidate across the planned block-bootstrap replicates; unresolved replicates remain in the denominator |
| Valid replicates | Resolved bootstrap replicates relative to the planned count |
| Chromosome consistency | Number of informative reference sequences and their agreement with the global candidate |
| Fit shape (JS divergence) | Difference between observed and fitted distribution shapes; lower means closer agreement |
| WARNING | A diagnostic threshold was not met; it does not veto the call |

Bootstrap support is not the probability that a call is correct. A reference sequence is informative with at least 200 fitting SNPs, converged models with finite BIC and delta BIC at least 10. No informative sequences means agreement cannot be assessed, not 0% agreement.

## Three-row evidence figure

Each row shows one candidate (2x, 3x, 4x) using the same fitting SNPs, bins and axis scales.

| Plot label | Meaning |
|---|---|
| Observed and fitted distributions | Comparison of observed site frequencies with model predictions |
| Observed | Gray observed distribution |
| Fitted model | Colored fitted distribution |
| Minor allele fraction | Minor allele count divided by the sum of the two leading allele counts |
| Sites (%) | Percentage of fitting sites in each bin |
| ΔBIC | This model's BIC minus the best model's BIC; the best model has zero |
| BEST | Best candidate among fitted models; not an independent declaration of a successful call |
| not converged | The corresponding model did not converge |

Figure ΔBIC values differ from the single TSV `delta_bic`, which compares the best two models. Model curves use each model's fitted parameters, observed depth distribution and the likelihood's discrete truncation; no smoothing is applied. `--export-figures` exports PDF and 300-dpi PNG. No standalone SVG is exported; offline HTML embeds vector markup internally. With no usable model evidence, no empty figure is generated.

## Reasons and warnings

### Reasons for uncertainty

| Identifier | Explanation and next check |
|---|---|
| LOW_SITE_COUNT | Too few fitting SNPs; inspect coverage, retained-site counts and usable reference regions |
| LOW_DELTA_BIC | The best two models are insufficiently separated; inspect the frequency distribution and regional consistency |
| MODEL_NOT_CONVERGED | At least one main model did not converge; inspect model diagnostics before interpreting a candidate |
| NO_ELIGIBLE_REGIONS | No usable analysis regions; check reference preparation and exclusions |
| NO_COVERAGE | No valid A/C/G/T observations in the analysis regions; check the BAM, reference and coverage |
| NO_VALID_SITES | Coverage exists, but no sites passed the filters; inspect actual depth and allele filters |

### Priority diagnostic warnings

| Identifier | Explanation |
|---|---|
| LOW_BOOTSTRAP_SUPPORT | Agreement below 0.85 across 100 planned block-bootstrap replicates |
| INSUFFICIENT_CHROMOSOME_EVIDENCE | Fewer than 3 informative reference sequences |
| CHROMOSOME_HETEROGENEITY | Fewer than 70% of informative reference sequences support the global candidate |
| POOR_THEORETICAL_FIT | JS divergence exceeds 0.05 or cannot be calculated as a finite value |

### Additional diagnostic warnings

| Identifier | Explanation |
|---|---|
| TOO_FEW_BLOCKS | Fewer than five occupied 5 Mb blocks for resampling |
| RHO_AT_BOUNDARY | Fitted dispersion is near an allowed boundary |
| A_AT_BOUNDARY | Tetraploid outer-peak weight is near its lower bound |
| A_AT_UPPER_BOUNDARY | Tetraploid central-peak weight is near zero |
| MODEL_ASSUMPTION_SENSITIVE | Relaxing the tetraploid outer-peak constraint changes the best candidate |

## Detailed sections and batch overview

| Report label | Contents |
|---|---|
| Detailed diagnostics | Reasons and warnings, plus per-reference-sequence fits when available |
| Reference-sequence fits | Sequence, SNP count, candidate, delta BIC and status |
| Inputs and analysis settings | BAM, reference assets, k-mer length, minimum region length, eligible region size, seed and requested filters |
| Actual site filters | Filters actually applied in the final analysis stage |
| Input notices | Input-specific messages, when present |
| Result table / Run record / Run log / Analysis regions | Links to TSV, detailed JSON, execution log and BED regions |
| Sample overview | Batch progress, status counts and called ploidy counts |
| Search sample / All statuses | Search sample names and filter by status |

Batch counts describe predictions, not accuracy. Click a sample name to open its report. An individual failure does not stop remaining samples. See [resuming and exit codes](CLI.md#exit-codes-and-resuming) for rerun behavior.
