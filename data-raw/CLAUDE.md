# hallberg2025.meth.data/data-raw — the methylation provenance pipeline

One `targets` DAG, one store, expressing the whole chain from JHU IDATs down to the
published `methylation_se`. **No trellis / FACETS / CNV context is needed here.**

The authoritative narrative version of this lives in the file `_targets.R:1-58` header —
read that before changing the graph. What follows is the operating summary.

## The DAG

```
[batch1 IDATs] --build_batch1_matrices--> bVals/mVals/targets      --check_*_b1
[batch2 IDATs] --build_batch2_matrices--> bVals081820/mVals081820  --check_*_b2
                     |                          |
                     +--------- build_orse -----+--> orse          --check_orse
                                                       |
[bValsselect.rds]  (branch C, frozen leaf) --+         |
[combmetadata.rds] (branch C, frozen leaf) --+-- build_se_lab_tcga --> se_lab_tcga
                                             |        --check_se_lab_tcga
read_methylation_data … drop_signet ---------+--> methylation_se   --verification
```

**Frozen leaves** — `format = "file"` nodes with no code producer. This is the honest
reproducibility floor, and the graph documents it by giving these nodes no upstream: the
two IDAT directories, `methylationmanifest.rds`, and the branch-C TCGA files
`bValsselect.rds` / `combmetadata.rds` (Step 3 closed as *not regenerable*). Frozen inputs
are staged in `<repo root>/extdata/`.

## Tolerance regime

Deepens toward the roots. The **terminal gate governs acceptance** regardless of any
assay-level check below it.

| Layer | Tolerance | Why |
|---|---|---|
| matrices | `compare_matrix` — set overlap + finite max diff | minfi funnorm is not bit-reproducible across versions |
| `orse` | `compare_orse(tol = 2)` | deterministic arithmetic, but integer-encoded (beta ×1000 / M ×100), so sub-permille matrix drift flips an integer by ±1 |
| `se_lab_tcga` | `tol = 2e-3`; colData labels expected **exact** | JHU beta cols inherit ≤2/1000 drift from `orse`; TCGA cols are exact (frozen); labels are name-based — and labels are the Change-002 surface |
| leaf | `verify_against_baseline` vs the committed `methylation_se.rds` (`inst/extdata/`), then `Rscript tests/verify_snapshot.R` at the repo root | the terminal gate |

**All `*_check` / `verification` targets are diagnostics.** They return structures; they do
not `stop()`. A tolerance miss is a reported value, not a build error — that is deliberate,
so every check still reports when a deeper layer drifts. Do not "fix" this by making them
throw.

## Running it

Cold runs are **cluster-only** (IDATs, minfi, memory):

```bash
sbatch run_targets.sh                     # on JHPCE
```

**`ENV-18` (2026-08-05): the R/package layer now runs inside a pinned Apptainer
image** (`hallberg2025-meth-data-raw.sif`, built from `docker/Dockerfile.meth-data-raw`,
`ENV-17`) instead of `module load conda_R/4.5.x` — closes the gap this project's own
`REPRODUCIBILITY_PLAN.qmd` flagged as "a reproducibility claim resting on a log line, not
a committed artifact." SLURM submission itself is unchanged. JHPCE ships the tool as
`singularity`, not `apptainer` (same `docker://`-derived image/CLI surface). Confirmed
end-to-end on the real cluster tree, real IDATs, both batches: the containerized
recompute of `batch1_matrices`/`batch2_matrices` from raw IDATs reproduces every
diagnostic number (`finite_maxdiff`, correlations, `check_orse`/`check_se_lab_tcga`
diffs, `verification`) **exactly** against the last real module-based run
(`slurm-34256940.out`, 2026-07-19) — containerization introduced zero new drift. The
non-identical-but-highly-correlated tolerance results themselves (`within_tol: FALSE`
throughout) are a pre-existing, already-documented condition (see "Tolerance regime"
above) unrelated to this cutover, not something this card resolved or needs to.

Two Apptainer/JHPCE-specific gotchas, both real and non-obvious, worth knowing before
touching `run_targets.sh` again:
- The `singularity` module (version 3.11.4) calls `os.getenv("HOSTNAME")`, unset in a
  non-login batch shell by default — `run_targets.sh` exports it first.
- Singularity doesn't auto-bind any cluster storage mount (`/dcs07`, `/dcs11`,
  `/dcl01`, all needed here) — every one must be passed via explicit `--bind`.
- `data-raw/.Rprofile` re-sources the project root's `.Rprofile` (see its own header),
  which activates `uvr`'s host-managed library — invisible to a normal `Rscript`
  call, but since R sources `.Rprofile` from cwd regardless of container boundaries,
  this silently overrode the container's own baked-in packages until
  `run_targets.sh` set `SINGULARITYENV_R_PROFILE_USER=/dev/null`.

If Apptainer/registry access ever becomes unreliable on JHPCE, `module load
conda_R/4.5.x` (the invocation this replaced) remains a documented fallback — see git
history on `run_targets.sh` for the pre-`ENV-18` version.

Downstream-only re-runs work anywhere once a cluster run has populated `_targets/`:

```r
targets::tar_make(callr_function = NULL, reporter = "verbose")
```

`callr_function = NULL` makes the current R session the worker, avoiding a
callr-subprocess segfault seen on macOS in this project. Keep it.

A local cold run additionally needs the frozen inputs staged in `<root>/extdata/` and the
IDATs present, or the `format = "file"` root targets error out.

## Hardcoded paths — known, carded, do not "fix" opportunistically

`run_targets.sh`'s absolute-path hardcoding and `_targets.R`'s lack of an environment
override were both **fixed by `MET-05`** (`Sys.getenv("BATCH1_IDAT_DIR"/"BATCH2_IDAT_DIR",
<default>)`). `MET-06` then replaced the *default* each falls back to with a repo-relative
indirection layer instead of a one-user, one-mount absolute path baked into this file —
see `data-raw/idat_links/`. The `Sys.getenv()` escape hatch from `MET-05` is unchanged.

**Mechanism updated by `MET-09`:** `MET-06` originally implemented the indirection as two
tracked OS symlinks (paths since renamed away — see below). Publishing this
package as its own repo (`MET-09`) found that `R CMD build`'s initial staging copy
(`cp -pLR`, dereferencing every symlink) aborts on dangling symlinks *before*
`.Rbuildignore` is ever read — and these are dangling by design off-cluster. The symlinks
are now tracked plain-text pointer files (`idat_links/batch1_idats.path` /
`batch2_idats.path`, one absolute path per line), read via `read_idat_link()` in
`_targets.R`. Same one-source-of-truth indirection, no filesystem symlink for `R CMD build`
to walk.

**Batch2 history** (preserved, extended — do not delete): for a long time the pipeline
deliberately did **not** use the data-warehouse copy of batch2, because that copy carried
a GenomeStudio `JHU_EST_1941_GSFile/` subdir owned by skoul with mode `drwxr-----` (group
lacked the traverse bit), which made `format = "file"` fail with "no read permission"
(job 34217423) — and its numbered-sentrix layout was *suspected*, but never actually
checked, not to match this sample sheet's `Basename` paths. The pipeline used the `dcl01`
scratch copy instead (Step 4 validated end-to-end).

**2026-07-27 update (`MET-04`/`MET-06`):** Rob asked Shashi to `chmod -R g+rX
JHU_EST_1941_GSFile/` on the warehouse copy, and Shashi did so — confirmed via `namei -l`
showing a group traverse (`x`) bit the full length of the path. `MET-06` then re-verified
the layout claim directly: `minfi::read.metharray.sheet()` against the warehouse copy
(`/dcs07/scharpf/data/data-warehouse/methyl-epic/ovarian-subtypes-2020-08-17`) resolves
all 48/48 Grn and all 48/48 Red idat files through its
`<Sentrix_ID>/<Sentrix_ID>_<Sentrix_Position>` subdirectory layout — the original layout
concern does not hold (or no longer does). `md5sum` of all 96 files against the `dcl01`
copy found them byte-identical, with matching `Sample_Name` assignment for all 48 samples.
**Batch2 is now repointed to the durable warehouse copy** (`idat_links/batch2_idats.path`);
the `dcl01` copy remains on disk (untouched) but is no longer the default.

## Change 002

`change002_delta_targets.R` + `R/historical_buggy_variant.R` run the historical buggy
variant in parallel as an A/B diagnostic. This is **not** part of the forward analysis and
is separate on purpose. Its store is `_targets_change002_delta/`, which is gitignored —
and note that `scripts/build_public.sh` uses rsync, which does not consult `.gitignore`,
so that store ships unless explicitly excluded (card `REL-02`).

The underlying defect: 15 "additional JHU" Normal samples were mislabeled in the shipped
`methylation_se.rda`. Root cause fixed and artifacts regenerated (CHANGES.md Change 002),
but the latent risk statement stands — any future analysis using `methylation_se` for
those 15 patients would silently substitute the wrong normal reference.
