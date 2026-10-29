# MapPloid 0.10.4 validation

**English** | [简体中文](README.zh-CN.md)

The public source was checked from a clean staging directory on 2026-10-29.
Results are recorded in [unit tests](github_tests.log) and
[fresh quickstart](github_quickstart.log). These are software checks, not
validation of biological accuracy or universal sequencing-depth requirements.

From the repository root, with the documented environment active:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python examples/quickstart/run.py
python scripts/validate_output.py
```

The suite covers model behavior, reference preparation, BAM counts, pipeline
reuse, sample calls and batch handling. The quickstart uses a synthetic 3x BAM
and supplied toy mappability scores; it checks `CALLED 3x` and exports PDF/PNG
figures and an offline HTML report. It does not run GenMap.

The optional output validator includes an empty BAM and an intentionally missing
BAM to check `UNCERTAIN`, `FAILED`, batch continuation and output reuse.
Generated inputs and reports are created locally and excluded from Git.
The old `scripts/validate_synthetic.py` entrypoint is superseded by the current
quickstart and output validator; historical v0.10.1 checks are not new v0.10.4 results.

Earlier local validation logs and audit records are retained in this directory
for provenance; absolute paths in them refer to the original development runs.
The GitHub preparation results above are the portable-source checks.
See [logo review](../assets/brand/review/REVIEW.md) for the completed visual review.

58 tests passed with zero skips in 61.340 seconds. A fresh quickstart returned
`CALLED 3x` with 35,255 fitting SNPs in 12.4 seconds. Wheel packaging passed;
see [packaging log](github_wheel.log).
