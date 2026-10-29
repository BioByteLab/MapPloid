# Changelog

**English** | [简体中文](CHANGELOG.zh-CN.md)

## Repository cleanup — 2026-10-29

- Keep one color SVG wordmark and focused bilingual user documentation.
- Remove design drafts, visual reviews, generated outputs, logs and internal
  manifests from the current repository; preserve them in the local release archive.
- Consolidate release history and installation information. Retain source,
  tests, the quickstart and output verification script.
- Allow the output verification script to create its missing parent directory.
- Software version remains 0.10.4. The published v0.10.4 tag retains its original
  commit; this repository cleanup follows that tag. Model behavior is unchanged.

## 0.10.4 — 2026-10-29

Based on the manuscript v20 MapPloid 0.10.3 source.

- Six-column sample result table with global delta BIC; reasons and warnings
  remain in `details/run.json` and collapsed report diagnostics.
- Self-contained offline HTML reports and searchable/filterable batch overview.
- Three decision gates displayed separately from bootstrap, chromosome and shape diagnostics.
- One three-row observed/fitted figure for 2x, 3x and 4x from the fitting site set,
  with optional PDF/PNG export; no separate BIC bar chart.
- One-command synthetic 3x example, bilingual documentation and a peak-shaped M wordmark.
- Remove the automatic 2x interpretation warning and its generating branch.
- Compact successful outputs retain a run record and log; `--keep-details`
  preserves stage data. Failed/interrupted runs retain recovery checkpoints.
- Output schema 3 requires new sample/batch output directories. Reference assets
  from 0.10.1–0.10.3 remain supported.
- Model likelihoods, filtering, sampling and decision gates are unchanged.
  Historical software and experiment outputs retain their original versions.
- Release checks: 58 tests passed, zero skips; fresh quickstart returned
  `CALLED 3x` with 35,255 fitting SNPs; wheel packaging passed.

## 0.10.1

Initial public distribution of the v19 manuscript companion software.

- Rename output labels: `2x-like` → `2x`, `3x` → `3x`,
  `4x-polysomic-like` → `4x`; `UNCERTAIN` is unchanged.
  This label change does not alter models, filters, thresholds or region selection.
- Include installation instructions, MIT license, environment definition and a
  deterministic synthetic example generator; fix its parent-directory creation.
- Historical checks: 39 tests passed with zero skips; the synthetic example
  passed with real GenMap and SAMtools for all three expected labels; wheel packaging passed.

A 2x label cannot exclude a disomically inherited tetraploid, and a 4x label
does not establish autopolyploid origin. These are interpretation limits.
Old analysis files keep their original labels and software versions; use new
output directories after version changes. Software checks do not establish
universal ploidy accuracy or sequencing-depth requirements.
