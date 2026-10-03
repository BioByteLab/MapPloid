# Validation

The unit-test suite covers counting, intervals, reference masks, model behaviour, CLI filters and pipeline checkpoints. Tests requiring SAMtools may be skipped if it is unavailable; report skipped tests explicitly.

The synthetic example creates six 1,500 bp artificial chromosomes and short pre-aligned BAMs at approximately 30 reads per site. Alleles are generated from prescribed dosage distributions. This is an integration fixture, not simulated real whole-genome resequencing and not a depth-boundary experiment.

Run from a clean output directory:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/validate_synthetic.py
```

Expected labels are `2x`, `3x`, `4x` for the respective synthetic inputs. The script asserts these labels and writes `validation/synthetic_e2e/verification.json` on success. Generated files are excluded from Git.

Release check on 2026-10-03: all 39 unit tests passed with zero skips. The synthetic end-to-end example passed for 2x, 3x and 4x using real GenMap 1.3.0 and SAMtools/HTSlib 1.24. Historical v0.10.0 experiments are not relabelled as v0.10.1 validation.
