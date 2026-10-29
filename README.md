# MapPloid 0.10.4

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/brand/mapploid-wordmark-white.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/brand/mapploid-wordmark-transparent.svg">
    <img src="assets/brand/mapploid-wordmark.svg" width="580" alt="MapPloid">
  </picture>
</p>

**English** | [简体中文](README.zh-CN.md)

Infer plant ploidy from allele frequencies in high-mappability reference regions.

[Logo assets and usage](assets/brand/README.md)

## Installation

```bash
conda env create -f environment.yml
conda activate mapploid
python -m pip install --no-deps .
mapploid --version
```

The environment includes SAMtools, GenMap and the Python dependencies.
License: [MIT](LICENSE). [Availability](docs/availability.md).

## Quick start

```bash
conda activate mapploid
./mapploid prepare --reference genome.fa --kmer 125 --out reference_assets
./mapploid call --bam sample.bam --assets reference_assets --out sample_result
./mapploid call --bam-list bams.txt --assets reference_assets --out batch_results
```

Open `report.html` in a browser. Reports work offline and contain their own figures.
Use `result.tsv` for downstream analysis.

Python 3.10+, NumPy, SciPy, pysam and matplotlib are required, plus SAMtools and
GenMap on PATH. Install this directory with `python -m pip install .`, or use its
`mapploid` launcher with an existing environment. `--version` reports 0.10.4.

## One-command example (no downloads)

From this version directory, with the MapPloid environment active:

```bash
python examples/quickstart/run.py
```

The example generates a 36 kb toy reference and a sorted/indexed synthetic 3x
BAM, prepares assets using supplied toy mappability scores, and runs the public
`call` command with the default filters, 25,000-SNP minimum and 100 bootstrap
replicates. It checks for `CALLED 3x` and exits successfully. No GenMap run or
network access is needed for this example; Python dependencies and SAMtools
must already be installed. The supplied scores are a fixture, not measured
mappability of the toy all-A reference. This checks the workflow, not biological
accuracy or GenMap.

Results: `examples/quickstart/output/result/report.html` and `result.tsv`;
PDF/PNG: `examples/quickstart/output/result/figures/`. Repeat the same command
to reuse the completed result, or add `--out /path/to/new-example` for a fresh run.
See the [quickstart tutorial](examples/quickstart/README.md) for the quick test,
step-by-step reference preparation, single-sample calls and batch analysis.

## English and Chinese documentation

All documentation is available in English and Simplified Chinese; use the links
at the top of each page to switch versions. Commands, code and code comments
use English. CLI help, runtime logs, report labels, fields and diagnostic
identifiers also use English for consistent lookup and programmatic use.
Both languages explain them in the [command-line reference](docs/CLI.md)
and [report guide](docs/REPORT_GUIDE.md).

## Results

```text
sample_result/
  result.tsv
  report.html
  regions.bed
  details/
    run.json
    run.log
```

A small hidden completion marker and lock support result verification and resuming.

| result.tsv column | Meaning |
|---|---|
| sample | Sample name |
| predicted_ploidy | 2x, 3x, 4x, UNCERTAIN; NA for a failed run |
| status | CALLED, UNCERTAIN or FAILED |
| selected_sites | SNPs actually used in the main model fit |
| delta_bic | Second-best BIC minus best BIC |
| analyzed_region_mb | Final analysis region size, in Mb |

Missing/uncomputed table values are NA; JSON uses null for non-finite values.
Machine-readable numeric values retain their precision. Output schema is **3**,
separate from the software version. Old result directories must not be reused.

The report shows:
- The final call and three decision requirements: at least 25,000 fitting SNPs,
  delta BIC at least 10, and convergence of all three main models.
- Bootstrap stability, informative chromosome count and chromosome agreement,
  and observed-versus-fitted shape (JS divergence).
- One figure with three rows: observed versus fitted minor-allele frequency
  distributions for 2x, 3x and 4x, with relative BIC and the best candidate marked.
  All rows use the actual fitting site set, common bins and common axis scales.
- Expandable detailed warnings, per-reference-sequence fits and actual filters.

Diagnostic warnings do not veto a call. Bootstrap support is a resampling
agreement rate, not a probability that the call is correct. An UNCERTAIN result
may have a best candidate; the report labels it as a candidate only. When no
sites or no coverage are available, the report explains why without empty plots.
No informative chromosomes means agreement is not assessable, rather than 0%.

`details/run.json` records input identities, requested and actual settings,
random seed, stage summaries, model fits, diagnostic metrics and plot data.
`reason` and `warnings` are retained in `details/run.json` under `result`,
and in the report’s collapsed detailed diagnostics, rather than the main TSV.
`details/run.log` contains execution messages and exception details.
BAM identity uses its path, size and modification time; indexes and the asset
manifest also have SHA256 checksums. This does not claim a full BAM checksum.

### Optional output

```bash
# Preserve stage tables, bootstrap replicates and intermediate site data.
./mapploid call --bam sample.bam --assets reference_assets --out detailed --keep-details

# Also export the three-row evidence figure in PDF and 300-dpi PNG.
./mapploid call --bam sample.bam --assets reference_assets --out exported --export-figures
```

Default compact runs retain stage checkpoints while running. Only after final
outputs are written and verified are stage directories removed. Interrupted and
failed runs retain their checkpoints. Successfully completed UNCERTAIN analyses
also use compact cleanup. `--keep-details` retains these directories.

Rerun the same command to resume. Completed results are checked and reused;
incomplete stages continue from retained checkpoints. Changed inputs, settings,
output options or software version require a new output directory. Only one
process may use an output directory at a time. Reference assets are never
removed by sample-result cleanup.

## Batch analysis

`bams.txt` contains one BAM path per line. Blank lines and lines beginning with
`#` are ignored. Relative paths resolve from the list's directory. Samples must
share a reference and have distinct BAM filenames. Sample names can contain
letters, digits, dots, underscores and hyphens.

```text
batch_results/
  result.tsv       Same six columns, plus relative result_file
  report.html      Search/filter samples and open individual reports
  samples/
    sample_A/
    sample_B/
  details/
    run.json
    run.log
```

The overview updates after each sample. A runtime failure is recorded as FAILED
and remaining samples continue; exit code is 1 if any sample failed. UNCERTAIN
is a completed analysis, not a runtime failure. An interrupted command exits
with code 130. Missing BAMs can be supplied and failed samples retried by rerunning
the same batch command. The report shows prediction counts, not accuracy.

## Inputs and filters

BAMs must be coordinate-sorted and indexed with a BAI or CSI. Reference sequence
names and lengths must match. `--assets` is the directory produced by
`prepare --out`, not a FASTA or BED file.

Reference preparation defaults: `--kmer 150`, two GenMap mismatches,
`--min-region-length 10000`, `--threads 8`. Optional `--exclude-bed` may be
repeated. `--genmap-bedgraph` can supply existing matching GenMap scores.
Strict score-1 regions are retained after gaps and exclusions; short regions
are not joined across gaps. Use a new asset directory for changed preparation
settings. The source FASTA must remain accessible.

This release preserves the manuscript v20 source's read-length check: up to
10,000 primary alignments are inspected, and reads shorter than the asset k-mer
length cause an input error.

| call option | Default |
|---|---|
| --allele-min / --allele-max | 0.15 / 0.85 |
| --min-depth / --max-depth | Automatic, estimated per analysis stage |
| --min-mapq / --min-baseq | 20 / 20 |
| --min-allele-reads | 2 |
| --pair-purity | 0.95 |

Allele frequencies use the sum of the two leading allele counts. Both must
pass the requested frequency interval. Site depth uses total A/C/G/T counts
after quality filtering and is not whole-sample sequencing depth. Depth bounds
are inclusive. Actual bounds are in the report and run record.

## Decision and diagnostic rules

The analysis begins with about 20 Mb, expanding to about 50 and 100 Mb if the
three decision requirements fail. Smaller references use available regions.
The main fit uses at least 25,000 and at most 100,000 SNPs. Above the maximum,
deterministic sampling allocates quotas proportional to chromosome counts.

Three decision failures produce uncertainty: LOW_SITE_COUNT, LOW_DELTA_BIC,
MODEL_NOT_CONVERGED. NO_ELIGIBLE_REGIONS, NO_COVERAGE and NO_VALID_SITES also
produce UNCERTAIN without a usable model fit. Runtime errors produce FAILED.

Priority diagnostics:
- LOW_BOOTSTRAP_SUPPORT: agreement below 0.85 across 100 planned block
  bootstrap replicates; unresolved replicates remain in the denominator.
- INSUFFICIENT_CHROMOSOME_EVIDENCE: fewer than 3 informative reference sequences.
- CHROMOSOME_HETEROGENEITY: agreement below 70% among informative sequences.
  An informative sequence has at least 200 fitting SNPs, converged models
  with finite BIC, and delta BIC at least 10.
- POOR_THEORETICAL_FIT: JS divergence above 0.05 or not finite.

Detailed diagnostics retain TOO_FEW_BLOCKS (fewer than five occupied 5 Mb
blocks), RHO_AT_BOUNDARY, A_AT_BOUNDARY, A_AT_UPPER_BOUNDARY and
MODEL_ASSUMPTION_SENSITIVE. These are warnings only.

## Compatibility and development

Reference assets from 0.10.1, 0.10.2 and 0.10.3 remain readable. Sample outputs
require new directories. Version 0.10.4 changes presentation and removes the
automatic 2x interpretation warning; the likelihoods, site selection and three
decision gates are unchanged. See [release notes](RELEASE_NOTES_0.10.4.md)
and [validation record](validation/README.md).

```bash
PYTHONPATH=src python -B -m unittest discover -s tests -v
```
