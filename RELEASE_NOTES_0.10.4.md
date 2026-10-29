# MapPloid 0.10.4

**English** | [简体中文](RELEASE_NOTES_0.10.4.zh-CN.md)

Updated 2026-10-29. Based on the manuscript v20 MapPloid 0.10.3 source.

- Six-column sample result table, including global delta BIC; reasons and warnings
  are retained in details/run.json and collapsed report diagnostics.
- Self-contained offline HTML reports and a searchable/filterable batch overview.
- Three decision gates shown separately from bootstrap, chromosome and shape diagnostics.
- One three-row observed/fitted figure for 2x, 3x and 4x from the fitted site set;
  optional PDF/PNG export. No separate BIC bar chart.
- One-command local synthetic 3x example: python examples/quickstart/run.py.
- Removed the automatic 2x interpretation warning and its generating branch.
- Compact successful outputs retain a run record and log; optional --keep-details
  preserves full stage data. Failed/interrupted runs retain recovery checkpoints.
- Output schema 3; use new sample/batch output directories. Older reference assets
  remain supported.
- No change to model likelihoods, filtering, sampling or decision gates.
- Existing v0.10.3 releases and experiment outputs are untouched.
- Main documentation, quickstart, command-line reference, report guide, release
  notes and validation record are available in English and Simplified Chinese.
  Code, code comments and example commands use English.
- Peak-M project icon and outlined wordmark, with transparent and monochrome SVG/PNG assets and theme-aware README headers.
