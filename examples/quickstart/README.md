# Quickstart: quick test and standard workflow

**English** | [简体中文](README.zh-CN.md)

Run all commands below from the repository root. First activate an environment with the dependencies installed:

```bash
conda activate mapploid
```

Requirements: Python 3.10+, NumPy, SciPy, pysam, matplotlib, SAMtools and GenMap. For a first installation, run `python -m pip install .` from this version directory.

## 1. One-command test (no downloads)

```bash
python examples/quickstart/run.py
```

The script creates small inputs, prepares reference regions, calls ploidy and prints the command for each step.
Expect `Quickstart passed: CALLED 3x`.
Runtime depends on the machine; repeat the command to reuse completed results.

After running the example, open these local files under `examples/quickstart/`:

- Example report (`output/result/report.html`): ploidy result and observed and model-fitted distributions.
- Main result table (`output/result/result.tsv`): six columns for downstream analysis.
- PDF (`output/result/figures/evidence.pdf`) / PNG (`output/result/figures/evidence.png`): exported figure; the PNG is 300 dpi.
- Detailed record (`output/result/details/run.json`): settings, diagnostic reasons and warnings.

Repeating the same command reuses completed results. Choose a new directory to run the full calculation again:

```bash
python examples/quickstart/run.py --out /tmp/mapploid_demo_new
```

The example uses fixed seeds to generate a 36 kb random reference and a matching 30x synthetic BAM, and runs GenMap to calculate mappability. The script constructs the BAM directly; it does not perform real sequencing or alignment. This tests the workflow, not biological accuracy. It uses default filters, the 25,000-site minimum and 100 bootstrap replicates. The small reference may produce diagnostic warnings such as `TOO_FEW_BLOCKS`.

## 2. Standard workflow: prepare the reference once → analyze samples → read results

### Step 1: prepare inputs

Use the same `genome.fa` used for BAM alignment, and a coordinate-sorted BAM with a BAI/CSI index. Reference sequence names and lengths must match. To sort and index an unprepared BAM:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

The examples below use `sample.sorted.bam`; substitute the path to an existing valid BAM if available.

### Step 2: prepare reference regions (once per reference and settings)

This example assumes 150 bp reads. `--kmer` must not exceed the read lengths actually used; account for shorter reads after trimming.

```bash
python mapploid prepare --reference genome.fa --kmer 150 --threads 8 --out reference_assets
```

This step runs GenMap to calculate mappability and create analysis regions. Runtime depends on reference size. `reference_assets/` becomes `--assets` in the next step and can be reused for samples sharing the reference. Keep the source FASTA available.

### Step 3: analyze one sample

```bash
python mapploid call --bam sample.sorted.bam --assets reference_assets --out sample_result --export-figures
```

Open `sample_result/report.html` or read `sample_result/result.tsv`.
`--export-figures` also creates `figures/evidence.pdf` and `figures/evidence.png`.

| Key field | Interpretation |
|---|---|
| predicted_ploidy | 2x / 3x / 4x; UNCERTAIN for an unresolved call; NA for a failed run |
| status | CALLED = a call was made; UNCERTAIN = insufficient evidence; FAILED = runtime failure |
| selected_sites | SNPs actually used in fitting |
| delta_bic | Second-best BIC minus best BIC; larger values mean stronger separation, not accuracy |

Reasons and warnings are in `details/run.json`; execution errors are in `details/run.log`. Use a new output directory when inputs or settings change. Rerun the same command to resume after interruption. Default outputs are compact; add `--keep-details` to retain intermediate tables.

### Multiple samples: reuse the reference and provide a BAM list

Create `bams.txt` with one BAM path per line, for example:

```text
/data/sample_A.bam
/data/sample_B.bam
```

```bash
python mapploid call --bam-list bams.txt --assets reference_assets --out batch_results --export-figures
```

Samples must share the reference and have distinct BAM filenames. Relative paths resolve from the list file's directory.
The summary table and overview are in `batch_results/`; individual results are in `batch_results/samples/SAMPLE/`.

### Practice manually with the generated example

After the one-command test, run these commands separately to see how `prepare --out` connects to `call --assets`. They reuse existing outputs:

```bash
python mapploid prepare --reference examples/quickstart/output/inputs/reference.fa --kmer 50 --threads 1 --out examples/quickstart/output/assets
python mapploid call --bam examples/quickstart/output/inputs/synthetic_3x.bam --assets examples/quickstart/output/assets --out examples/quickstart/output/result --export-figures
```


Further reading: [command-line reference](../../docs/CLI.md) · [report guide](../../docs/REPORT_GUIDE.md).
