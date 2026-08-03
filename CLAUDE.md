# hallberg2025.meth.data — methylation companion data package

**This arm has zero coupling to `hallberg2025.seq.data` in either direction.** Nothing here needs
trellis, FACETS, CNV, or structural-variant context — do not go read it. The only shared
dependency is the base package `hallberg2025.base`.

## What this package is (and is not, yet)

It is a **package shell**: `DESCRIPTION` and `data-raw/`, and nothing else. There is no
`R/`, `NAMESPACE`, `man/`, `inst/extdata/`, `tests/`, or README, so it is **not
`R CMD INSTALL`-able** and CI cannot check it. Card `MET-01` fixes that.

Consequently the published `methylation_se` is still served from
`hallberg2025.base/data/methylation_se.rda`, not from here. Migrating it is a three-card
sequence — `MET-02a` (add serving here) → `MAN-03` (repoint the manuscript pipeline) →
`BASE-03` (remove the old copy) — with a green `verify_snapshot.R` at each boundary. Ship
both copies in the middle; do not try to do it in one step.

Like `hallberg2025.seq.data`, this directory now has its own nested `.git` (`MET-09`, private repo
`cancer-genomics/hallberg2025.meth.data`) in addition to being tracked directly in the outer repo.

## Provenance record

`DESCRIPTION:12-13` points at `provenance/methylation_provenance_verification_plan.md` —
the primary provenance record for this arm, tracked and fresh-clone-reachable since
card `DOC-04` moved it there from the gitignored `code/archive/`. Do not delete the
sentence or the file: nothing else records the full provenance DAG and step sequence.

## The generating pipeline

`data-raw/_targets.R` is a single unified DAG from JHU IDATs down to the published
`methylation_se`. See `data-raw/CLAUDE.md` for the DAG, the frozen leaves, the tolerance
regime, and how to run it.
