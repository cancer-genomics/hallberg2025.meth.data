# OsMethExpData

Companion data package for Hallberg et al. 2025 (*Cancer Research Communications*),
methylation arm. Provides the methylation `SummarizedExperiment` consumed by the
manuscript pipeline in `cancer-genomics/hallberg2025`. All sample identifiers are CG lab
IDs (e.g. `CGOV167T`).

This package has **zero coupling** to `OsSeqExpData` (the sequencing/FACETS/trellis
companion package) in either direction.

## Status

- Installable as of card `MET-01`: package skeleton, `methylation_se()` accessor, and
  tests exist.
- `inst/extdata/methylation_se.rds` does **not** exist yet. Calling `methylation_se()`
  now raises an informative "no file found" error rather than returning `NULL`. Serving
  the real data is card `MET-02a`.
- The published `methylation_se` is, for now, still served from
  `ovarian.subtypes/data/methylation_se.rda` in the base package. Migrating callers to
  this package is a three-card sequence (`MET-02a` -> `MAN-03` -> `BASE-03`); both copies
  ship side by side in the middle so `Rscript tests/verify_snapshot.R` stays green
  throughout.

## Generating pipeline

`data-raw/` contains the targets-based pipeline that produces `methylation_se` from JHU
IDATs. See `data-raw/CLAUDE.md` for the DAG, frozen leaves, and tolerance regime, and
`OsMethExpData/CLAUDE.md` for the arm-level orientation (known issues, card sequence,
zero-coupling boundary with `OsSeqExpData`).
