# OsMethExpData — methylation companion data package

**This arm has zero coupling to `OsSeqExpData` in either direction.** Nothing here needs
trellis, FACETS, CNV, or structural-variant context — do not go read it. The only shared
dependency is the base package `ovarian.subtypes`.

## What this package is (and is not, yet)

It is a **package shell**: `DESCRIPTION` and `data-raw/`, and nothing else. There is no
`R/`, `NAMESPACE`, `man/`, `inst/extdata/`, `tests/`, or README, so it is **not
`R CMD INSTALL`-able** and CI cannot check it. Card `MET-01` fixes that.

Consequently the published `methylation_se` is still served from
`ovarian.subtypes/data/methylation_se.rda`, not from here. Migrating it is a three-card
sequence — `MET-02a` (add serving here) → `MAN-03` (repoint the manuscript pipeline) →
`BASE-03` (remove the old copy) — with a green `verify_snapshot.R` at each boundary. Ship
both copies in the middle; do not try to do it in one step.

Like `OsSeqExpData`, this directory now has its own nested `.git` (`MET-09`, private repo
`cancer-genomics/OsMethExpData`) in addition to being tracked directly in the outer repo.

## Provenance record

`DESCRIPTION:12-13` points at `provenance/methylation_provenance_verification_plan.md` —
the primary provenance record for this arm, tracked and fresh-clone-reachable since
card `DOC-04` moved it there from the gitignored `code/archive/`. Do not delete the
sentence or the file: nothing else records the full provenance DAG and step sequence.

## The generating pipeline

`data-raw/_targets.R` is a single unified DAG from JHU IDATs down to the published
`methylation_se`. See `data-raw/CLAUDE.md` for the DAG, the frozen leaves, the tolerance
regime, and how to run it.
