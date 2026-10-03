# MapPloid v0.10.1

Initial public distribution of the v19 manuscript companion software.

- Core Python source is identical to the frozen v19 companion v0.10.1.
- Output labels are 2x, 3x, 4x and UNCERTAIN; see LABEL_CHANGE.md for biological interpretation.
- Includes installation instructions, MIT license, environment definition and a deterministic synthetic example generator.
- Fixes parent-directory creation in the standalone synthetic example script; no model, threshold or filtering changes.
- All 39 existing tests pass, with zero skips. The synthetic example passes with real GenMap and SAMtools for all three expected labels. Wheel packaging succeeds.

This is a software release, not a claim of universal ploidy accuracy. Historical experiments retain their frozen versions. Depth and species-matched validation are ongoing. No publication DOI or software-archive DOI has been assigned.
