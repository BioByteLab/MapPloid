# MapPloid

Mapping-aware allele counts for plant ploidy screening from whole-genome resequencing alignments.

[中文说明](README.zh-CN.md)

MapPloid builds reusable high-mappability reference regions, extracts allele counts from BAM files, and compares diploid, triploid and polysomic-like tetraploid models. It expands the analysed region when evidence is insufficient and reports `UNCERTAIN` when the configured evidence requirements are not met.

## Availability and requirements

| Item | Description |
| --- | --- |
| Project | MapPloid |
| Version | 0.10.1 |
| Intended project URL | https://github.com/BioByteLab/MapPloid |
| Operating system | Linux; tested on Ubuntu |
| Language | Python |
| Python requirements | Python >=3.10, NumPy >=1.26, SciPy >=1.11, pysam >=0.22 |
| External tools | GenMap and SAMtools, available on `PATH` |
| License | [MIT](LICENSE) |
| Restrictions on non-academic use | None beyond the MIT license terms |
| Archived version / DOI | Not yet assigned |

The repository URL is the intended publication destination until the repository is created. A Git tag or GitHub URL is not a DOI. See [availability details](docs/availability.md).

## Installation

Use a separate environment. The environment file includes the external tools:

```bash
conda env create -f environment.yml
conda activate mapploid
python -m pip install --no-deps .
mapploid --version
```

Alternatively, install the Python package into an existing compatible environment with `python -m pip install .` and separately provide GenMap and SAMtools on `PATH`. SAMtools must support the `mpileup --no-output-ins`, `--no-output-del` and `--no-output-ends` options used by the counting backend.

## Quick start

```bash
mapploid prepare --reference genome.fa --read-length 125 --threads 8 --out reference_assets
mapploid call --bam sample.bam --assets reference_assets --out sample_result
```

The reference must be an uncompressed FASTA. The BAM must be coordinate-sorted, indexed and use the same reference sequence dictionary. Choose the read length for your actual reads; do not reuse PE150 assets for PE125 reads. Exclude organellar or other unwanted regions with `prepare --exclude-bed` when appropriate.

Reference preparation is performed once per compatible reference/read-length configuration. Sample analysis proceeds through approximately 20, 50 and 100 Mb of spatially distributed high-mappability sequence, stopping when evidence is sufficient. This public interface does not fall back to the complete mask.

Outputs:

- `result.tsv`: predicted ploidy, status and supporting evidence.
- `regions.bed`: the final analysed regions.
- `details/`: configuration, stage results, diagnostics and logs.

Use a new output directory after changing inputs, software version or filtering parameters. Completed stages are reused only after consistency checks. Do not use the same output directory concurrently.

## Validation and minimal example

Run the existing test suite from the repository root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Generate a small deterministic reference and artificial BAM examples, then run real GenMap, SAMtools and the default MapPloid model:

```bash
python scripts/validate_synthetic.py
```

Expected final message: `SYNTHETIC_END_TO_END_OK: 2x, 3x, 4x`.
The generated inputs and output tables are placed in `validation/synthetic_e2e/`.
These deliberately simplified examples test software integration, not biological accuracy or realistic sequencing depth limits. See [validation notes](validation/README.md).

## Interpretation and limitations

Labels `2x`, `3x` and `4x` indicate compatibility with the implemented allele-count models; they do not directly count chromosomes. A `2x` result cannot exclude disomic tetraploidy, and `4x` does not establish autopolyploid origin. Higher ploidies and unmodelled dosage patterns may receive misleading confident labels, rather than necessarily returning `UNCERTAIN`.

The default evidence rules require at least three informative reference chromosomes. Fragmented assemblies and limited genomic regions have not been validated as equivalent substitutes. Confidence measures describe model stability and are not posterior probabilities of biological ploidy. Independent cytometric or cytological evidence remains important.

Software tests do not establish that this public release reproduces every historical manuscript result. Earlier experiments used frozen versions and must retain their actual version identifiers. Broader simulator and species-matched validation is ongoing; no universal depth threshold is claimed here.

## Citation and support

The manuscript is in preparation; no published article DOI is assigned. Report the software version and commit used in analyses. Once the repository is public, use its Issues page for reproducible bug reports without uploading private or large sequencing datasets.
