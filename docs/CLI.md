# Command-line reference

**English** | [简体中文](CLI.zh-CN.md)

[Main documentation](../README.md) · [Quickstart](../examples/quickstart/README.md)

Run these commands from the version directory, or use the installed `mapploid` command without the `python` prefix. 

```bash
python mapploid --help
python mapploid prepare --help
python mapploid call --help
python mapploid --version
python examples/quickstart/run.py --help
```

## General options

| Option | Meaning |
|---|---|
| -h, --help | Show help for the current command and exit |
| --version | Show the software version and exit; use before the subcommand |

## Prepare reference regions

```bash
python mapploid prepare --reference genome.fa --kmer 150 --out reference_assets
```

| Option | Required / default | Meaning |
|---|---|---|
| --reference FASTA | Required | Uncompressed reference FASTA, matching the BAM alignment reference |
| --out DIR | Required | Reference asset directory; pass it to `call --assets` |
| --kmer BP | 150 | GenMap k-mer length in bp; must not exceed the read lengths used |
| --min-region-length BP | 10000 | Minimum continuous high-mappability region length after exclusions |
| --threads N | 8 | Threads for reference preparation |
| --exclude-bed BED | None | Exclude regions in a BED file; may be repeated |

Preparation runs GenMap to calculate mappability scores; resuming the same task automatically reuses generated scores. Keep the original FASTA accessible. Reuse assets only for the same reference and preparation settings; use a new directory when these change.

## Call one sample or a batch

```bash
python mapploid call --bam sample.bam --assets reference_assets --out sample_result
python mapploid call --bam-list bams.txt --assets reference_assets --out batch_results
```

| Option | Required / default | Meaning |
|---|---|---|
| --bam FILE | Exactly one of `--bam` / `--bam-list` | One coordinate-sorted BAM with a BAI/CSI index |
| --bam-list FILE | Exactly one of `--bam` / `--bam-list` | One BAM path per line; each BAM must be sorted and indexed |
| --assets DIR | Required | Directory created by `prepare --out` |
| --out DIR | Required | Single-sample output directory or batch summary directory |
| --sample NAME | BAM filename without `.bam` | Single-sample name; not allowed with `--bam-list` |
| --keep-details | Off by default | Retain stage tables and intermediate site data after completion |
| --export-figures | Off by default | Export the report figure as PDF and 300-dpi PNG when model evidence is available |

Blank lines and lines starting with `#` in BAM lists are ignored. Relative paths resolve from the list's directory. BAM filenames must be distinct and samples must share a reference. Batch results are under `samples/SAMPLE/`; the batch TSV adds a `result_file` column.

### Optional site filters

| Option | Default | Meaning |
|---|---|---|
| --allele-min FLOAT | 0.15 | Minimum frequency of each of the two leading alleles |
| --allele-max FLOAT | 0.85 | Maximum frequency of each of the two leading alleles |
| --min-depth N | Automatic | Minimum site depth after quality filtering; boundary included |
| --max-depth N | Automatic | Maximum site depth after quality filtering; boundary included |
| --min-mapq N | 20 | Minimum read mapping quality, in the range 0–255 |
| --min-baseq N | 20 | Minimum base quality, in the range 0–255 |
| --min-allele-reads N | 2 | Minimum supporting reads for each retained allele |
| --pair-purity FLOAT | 0.95 | Minimum fraction of valid base counts explained by the two leading alleles |

Frequencies use the sum of the two leading allele counts. Site depth uses all quality-filtered A/C/G/T counts; it is not whole-sample sequencing depth. Automatic depth bounds are estimated per analysis stage and recorded in the report. The public `call` command uses the built-in model decision settings; these are described in the [report guide](REPORT_GUIDE.md).

## Quickstart script

```bash
python examples/quickstart/run.py
python examples/quickstart/run.py --out /tmp/mapploid_demo_new
```

`--out DIR` selects the example working directory (default: `examples/quickstart/output/`). The script creates inputs, runs GenMap to prepare assets and runs a synthetic 3x call. Repeat the same command to reuse its result; use a new directory for a fresh calculation. `-h` / `--help` displays script help.

## Exit codes and resuming

| Exit code | Meaning |
|---|---|
| 0 | Completed successfully, including a completed UNCERTAIN analysis |
| 1 | Runtime failure; for a batch, at least one sample failed |
| 2 | Invalid command-line arguments rejected by the argument parser |
| 130 | Interrupted during execution |

Rerun the same command to resume or reuse verified results. Changed inputs, settings, output options or software version require a new output directory. Failed/interrupted samples retain checkpoints. For a missing BAM, supply the file and rerun the same batch command. Execution details are in `details/run.log`; structured reasons and warnings are in `details/run.json`. See the [report guide](REPORT_GUIDE.md) for their interpretation.

## Input checks

BAMs must be coordinate-sorted and indexed with a BAI or CSI. Reference sequence
names and lengths must match. `--assets` is the directory produced by
`prepare --out`, not a FASTA or BED file.

Reference preparation defaults: `--kmer 150`, two GenMap mismatches,
`--min-region-length 10000`, `--threads 8`. Optional `--exclude-bed` may be
repeated to exclude specified regions.
Strict score-1 regions are retained after gaps and exclusions; short regions
are not joined across gaps. Use a new asset directory for changed preparation
settings. The source FASTA must remain accessible.

Read-length check: up to 10,000 primary alignments are inspected, and reads shorter than the asset k-mer
length cause an input error.

Both leading alleles must pass the frequency interval. Sample names may contain letters, digits, dots, underscores and hyphens.

Successful compact runs remove stage directories only after final outputs are verified; completed UNCERTAIN runs also use compact cleanup. `--keep-details` retains stage data. Only one process may use an output directory at a time. Sample cleanup never removes reference assets.
