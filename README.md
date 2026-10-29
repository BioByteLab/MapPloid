# MapPloid 0.10.4

<p align="center">
  <img src="assets/brand/mapploid-wordmark.svg" width="580" alt="MapPloid">
</p>

**English** | [简体中文](README.zh-CN.md)

Infer plant ploidy from allele frequencies in high-mappability reference regions.

## Installation

Linux (tested on Ubuntu), Python 3.10+. The environment below includes NumPy 1.26+,
SciPy 1.11+, pysam 0.22+, matplotlib 3.8+, SAMtools and GenMap.

```bash
conda env create -f environment.yml
conda activate mapploid
python -m pip install --no-deps .
mapploid --version
```

With an existing environment, install with `python -m pip install .` and ensure
SAMtools and GenMap are on PATH. License: [MIT](LICENSE), including non-academic use.

## Try the example

From the repository root, with the environment active:

```bash
python examples/quickstart/run.py
```

This generates synthetic inputs and checks for `CALLED 3x`, without downloads or
a GenMap run. It uses supplied toy mappability scores and checks the software
workflow, not biological accuracy. Open `examples/quickstart/output/result/report.html`;
the same directory contains `result.tsv` and PDF/PNG figures under `figures/`.

## Analyze your data

Use coordinate-sorted, indexed BAMs with matching reference sequence names and lengths.
Choose a k-mer length no greater than the read lengths; 125 is an example.

```bash
mapploid prepare --reference genome.fa --kmer 125 --out reference_assets
mapploid call --bam sample.bam --assets reference_assets --out sample_result
mapploid call --bam-list bams.txt --assets reference_assets --out batch_results
```

`bams.txt` contains one BAM path per line. Relative paths resolve from the list's
directory; samples must share a reference and have distinct BAM filenames.
Reuse reference assets for the same reference and preparation settings.

## Read the results

| File | Use |
|---|---|
| `report.html` | Offline report; batch reports link to individual samples |
| `result.tsv` | Ploidy, status, fitting SNP count, delta BIC and analyzed region size |
| `regions.bed` | Analysis regions for a sample |
| `details/run.json` | Settings, model fits, reasons and warnings |
| `details/run.log` | Execution messages and errors |

`CALLED` reports 2x, 3x or 4x; `UNCERTAIN` means no definitive call; `FAILED` means
a runtime error. Add `--export-figures` for PDF/PNG, or `--keep-details` to retain
intermediate data. Repeat the same command to resume or reuse verified results.
Changed inputs, settings, output options or software version require a new output directory.

## Documentation

- [Quickstart and standard workflow](examples/quickstart/README.md)
- [Commands, input requirements and filters](docs/CLI.md)
- [Results, decision rules and diagnostics](docs/REPORT_GUIDE.md)
- [Changes and compatibility](CHANGELOG.md)

Documentation is available in English and Simplified Chinese. Commands, code,
logs and report identifiers use English. Output schema is 3; use new sample/batch
output directories for older results. Reference assets from 0.10.1–0.10.3 remain supported.

## Verification

The 2026-10-29 release checks passed 58 tests with zero skips, a fresh synthetic
`CALLED 3x` quickstart and wheel packaging. Reproduce the checks with:

```bash
PYTHONPATH=src python -B -m unittest discover -s tests -v
python examples/quickstart/run.py
python scripts/validate_output.py
```

The output check exercises successful, uncertain and failed samples, batch
continuation and result reuse. Generated data and validation records stay local.
