## Methylation provenance-verification pipeline -- UNIFIED single DAG.
## (code/methylation_provenance_verification_plan.md)
##
## ONE pipeline, one store (_targets/), expressing the whole chain from the JHU
## IDATs down to the published `methylation_se`. targets makes every dependency
## explicit (tar_visnetwork()) and caches every intermediate, so re-verifying a
## single layer after a code change re-runs only that layer's descendants -- the
## cached upstream objects (matrices, orse) are reused without recomputation.
## This supersedes the former per-step files (_targets_step{4,5}.R): those froze
## the frontier one generation at a time as a bring-up scaffold; caching now gives
## the same "don't recompute upstream" behavior for free. The Change-002 A/B delta
## diagnostic stays separate on purpose (change002_delta_targets.R) -- it runs the
## historical buggy variant in parallel and is not part of the forward analysis.
##
## THE DAG (roots at top; each computed layer has a direct intermediate check):
##
##   [batch1 IDATs] --build_batch1_matrices--> bVals/mVals/targets  --check_*_b1
##   [batch2 IDATs] --build_batch2_matrices--> bVals081820/mVals081820 --check_*_b2
##                        |                         |
##                        +--------- build_orse ----+--> orse  --check_orse
##                                                          |
##   [bValsselect.rds] (branch C, frozen leaf) --+          |
##   [combmetadata.rds] (branch C, frozen leaf) -+-- build_se_lab_tcga --> se_lab_tcga
##                                               |          --check_se_lab_tcga
##   read_methylation_data...drop_signet --------+--> methylation_se --verification
##
## Frozen leaves (format="file", no code producer -- the honest reproducibility
## floor): the two IDAT dirs (roots of the JHU arm), methylationmanifest.rds, and
## the branch-C TCGA files bValsselect.rds / combmetadata.rds (Step 3 closed as
## not-regenerable; see the provenance plan). The graph itself now documents that
## floor: these nodes simply have no upstream.
##
## TOLERANCE REGIME (deepens toward the roots; the TERMINAL gate governs
## acceptance regardless of any assay-level check below):
##   matrices : minfi funnorm is NOT bit-reproducible across versions -> tolerance
##              (compare_matrix). Reports set overlap + finite max diff.
##   orse     : deterministic arithmetic on the matrices, BUT integer-encoded
##              (beta x1000 / M x100), so sub-permille matrix drift can flip an
##              integer by +/-1 -> compare_orse(tol = 2).
##   se_lab_tcga : JHU beta cols = orse beta / 1000, so they inherit <= 2/1000
##              drift; TCGA cols are exact (branch-C frozen) -> tol = 2e-3. colData
##              labels (the Change-002 surface) are name-based -> expected exact.
##   leaf     : verify_against_baseline vs the committed methylation_se.rds, then
##              the TERMINAL gate `Rscript tests/verify_snapshot.R` at the repo root.
##
## All *_check / verification targets are DIAGNOSTICS: they return structures, they
## do not stop() the pipeline. A tolerance miss is a reported value, not a build
## error, so every check reports even when a deeper layer drifts.
##
## RUN: cluster-only for a cold run (IDATs + minfi + memory). Set the two IDAT
## dirs below, then `sbatch run_targets.sh` on JHPCE. A local run needs the frozen
## inputs staged in <root>/extdata/ (see the former locate_step4_inputs.sh) and
## the IDATs present; without them the format="file" root targets error. Once a
## cluster run has populated _targets/, downstream-only re-runs work anywhere.
##
##   Rscript -e "targets::tar_make(callr_function = NULL, reporter = 'verbose')"
## (callr_function = NULL: the current R session is the worker, avoiding a
## callr-subprocess segfault seen on macOS in this project -- see CLAUDE.md.)

library(targets)
library(here)

## hallberg2025.base is not installed -- load from source (callr_function = NULL
## means the current R session is the worker, so load_all() suffices).
pkgload::load_all(here("..", "..", "hallberg2025.base"), quiet = TRUE)

tar_option_set(packages = c("SummarizedExperiment", "S4Vectors",
                            "tibble", "dplyr", "readr", "magrittr",
                            "minfi",
                            "IlluminaHumanMethylationEPICanno.ilm10b2.hg19"))

tar_source("R/generate.R")
tar_source("R/env_check.R")

epath <- function(...) here("..", "..", "extdata", ...)

## Reads the absolute path out of a tracked idat_links/*.path pointer file
## (MET-09: real OS symlinks here broke `R CMD build`'s cp -pLR staging copy --
## it dereferences every symlink before .Rbuildignore is ever consulted, and
## aborts on the dangling ones these deliberately are off-cluster. A tracked
## plain-text pointer carries the same one-source-of-truth path without being
## a filesystem symlink R's build step has to walk).
## Takes a resolved path, not a name -- idat_dir() below has already established
## that the file exists.
read_idat_link <- function(f) {
  ## Skip blank lines and #-comments rather than blindly taking line 1: the
  ## .path.example template leads with instructions, so a file copied from it and
  ## edited anywhere below would otherwise read a comment as the path.
  ln <- trimws(readLines(f, warn = FALSE))
  ln <- ln[nzchar(ln) & !startsWith(ln, "#")]
  if (!length(ln)) {
    stop("No path found in ", f, " -- every line is blank or a comment.\n",
         "Put the absolute IDAT directory on a line of its own.", call. = FALSE)
  }
  ln[1L]
}

## Absolute path to one IDAT directory: the environment variable if set, else the
## idat_links/ pointer file, else a self-describing sentinel.
##
## Two things this has to get right, both learned the hard way (MET-14):
##
## 1. `Sys.getenv(var, read_idat_link(name))` does NOT short-circuit. Sys.getenv()
##    calls as.character() on `unset` before it consults the environment, which
##    forces the promise -- so the fallback ran, and stopped, even when the
##    variable was set. The override was documented in three places and worked in
##    none of them. Read the variable first, explicitly.
##
## 2. This runs at FILE SCOPE, so it runs every time _targets.R is sourced --
##    including under tar_validate(), a static DAG check that reads no IDAT and
##    touches no cluster. Stopping here left a fresh clone unable to validate at
##    all, which failed CI at step 12 of 23 and skipped the eleven steps after it.
##    An unconfigured directory must fail at tar_make() time, where the
##    format="file" target actually needs it, not at definition time.
idat_dir <- function(var, name) {
  from_env <- Sys.getenv(var, unset = "")
  if (nzchar(from_env)) return(from_env)

  f <- here("idat_links", name)
  if (file.exists(f)) return(read_idat_link(f))

  message("No IDAT directory configured for ", var, ".\n",
          "MET-13: the pointer files are gitignored -- they held absolute cluster ",
          "paths and this repo is public. Only *.path.example ships.\n",
          "The pipeline can still be validated; it cannot be run. Fix either way:\n",
          "  cp ", f, ".example ", f, "   # then edit in the local IDAT directory\n",
          "  export ", var, "=...                      # bypasses that file entirely")
  sprintf("<unconfigured IDAT directory: set %s, or fill in %s>", var, f)
}

## >>> IDAT-DIR PATHS. <<<
## Both default to a repo-relative pointer file in idat_links/ rather than a
## hardcoded absolute cluster path -- MET-06, mechanism updated by MET-09. The
## Sys.getenv() override from MET-05 is kept: any deployment off this cluster (or
## onto a different mount of the same data) can still override without touching
## this file. MET-14 is what made that override actually work.
##
## MET-13: the pointer files themselves are NO LONGER tracked -- this repo is
## public and their entire content was an absolute lab cluster path. Only
## *.path.example ships. A fresh clone therefore needs one setup step before
## data-raw can RUN: copy an .example and fill it in, or set the two env vars.
## It needs no setup at all to be VALIDATED -- idat_dir() above explains itself
## and returns a sentinel rather than stopping, so tar_validate() works on a bare
## clone (MET-14; stopping here had been failing CI at step 12 of 23).
##
## batch1: durable Azure-backed data-warehouse copy (numbered sentrix dirs only;
##   no GenomeStudio folder, so format="file" hashes it cleanly). An alternate
##   copy on dcl01 scratch, under dhallber's tree, was also verified (48 EPIC).
BATCH1_IDAT_DIR <- idat_dir("BATCH1_IDAT_DIR", "batch1_idats.path")
## batch2: history -- NOT the data-warehouse copy for a long time, because that
##   copy carried a GenomeStudio JHU_EST_1941_GSFile/ subdir owned by skoul, mode
##   drwxr----- (group had no traverse bit), which made format="file" fail with
##   "no read permission" (cf. job 34217423), and its numbered-sentrix layout was
##   suspected (never actually checked) not to match this sample sheet's Basename
##   paths. The pipeline instead used the dcl01 copy Step 4 validated end-to-end
##   (JHU_EST_1941_idats/ layout, group-readable) -- dcl01 is scratch.
##
##   2026-07-27 (MET-04/MET-06): Rob asked Shashi to `chmod -R g+rX
##   JHU_EST_1941_GSFile/` on the warehouse copy, and Shashi did so -- confirmed
##   via `namei -l` showing a group traverse (x) bit the length of the path.
##   MET-06 then re-checked the layout claim directly: minfi::read.metharray.sheet()
##   against the batch2 warehouse copy (the path idat_links/batch2_idats.path
##   points at) resolves all 48/48 Grn and all 48/48 Red idat files via its
##   <Sentrix_ID>/<Sentrix_ID>_<Sentrix_Position> subdirectory layout -- the
##   layout claim was wrong, or true of an earlier state that no longer holds.
##   md5sum comparison of all 96 files against the dcl01 copy found them
##   byte-identical, and Sample_Name assignment matched for all 48 samples.
##   Batch2 is therefore now repointed to the durable warehouse copy; dcl01
##   remains on disk but is no longer the default.
BATCH2_IDAT_DIR <- idat_dir("BATCH2_IDAT_DIR", "batch2_idats.path")

list(

  ## == Frozen roots (format="file", no code producer) =======================

  ## JHU IDAT directories -- the roots of the array arm.
  tar_target(batch1_idat_dir, BATCH1_IDAT_DIR, format = "file"),
  tar_target(batch2_idat_dir, BATCH2_IDAT_DIR, format = "file"),

  ## Cross-reactive probe list (450K, Chen et al. 2013). Vendored from
  ## methylationArrayAnalysis inst/extdata (md5 3bec72842c8c24fae1f9c80516c71397,
  ## 29233 probes) so reproduction needs no heavyweight Bioconductor install.
  tar_target(xreactive_file,
             here("..", "..", "hallberg2025.base", "inst", "extdata",
                  "48639-non-specific-probes-Illumina450k.csv"),
             format = "file"),

  ## Irreducible: no code producer found in any tree (provenance plan).
  tar_target(methmanifest_file, epath("methylationmanifest.rds"), format = "file"),
  ## Batch-2 phenotype CSV used by build_orse().
  tar_target(methdat_file,      epath("methdat082620.csv"),        format = "file"),

  ## Branch-C TCGA frozen leaves (Step 3 closed as not-regenerable; the honest
  ## reproducibility floor). They enter build_se_lab_tcga() with no upstream.
  tar_target(bValsselect_file,  epath("bValsselect.rds"),   format = "file"),
  tar_target(combmetadata_file, epath("combmetadata.rds"),  format = "file"),

  ## == Frozen twins (compared-against only; NOT build inputs) ================
  ## The published matrices / se.rds / se_lab_tcga.rds, used by the per-layer
  ## direct intermediate checks (safeguard B/D).
  tar_target(bVals_frozen_file,        epath("bVals.rds"),        format = "file"),
  tar_target(mVals_frozen_file,        epath("mVals.rds"),        format = "file"),
  tar_target(targets1_frozen_file,     epath("targets.rds"),      format = "file"),
  tar_target(bVals081820_frozen_file,  epath("bVals081820.rds"),  format = "file"),
  tar_target(mVals081820_frozen_file,  epath("mVals081820.rds"),  format = "file"),
  tar_target(se_jhu_frozen_file,
             here("..", "..", "output", "methylation.Rmd", "se.rds"),
             format = "file"),
  tar_target(se_lab_tcga_frozen_file,
             epath("se_lab_tcga.rds"),
             format = "file"),

  ## == Step-1 static inputs =================================================
  tar_target(match_table_file,
             here("..", "..", "hallberg2025.base", "inst", "extdata", "match_table.csv"),
             format = "file"),
  tar_target(signet_file,
             here("..", "..", "hallberg2025.base", "inst", "extdata", "stomach_muc_signet.csv"),
             format = "file"),
  tar_target(baseline_file,
             here("..", "inst", "extdata", "methylation_se.rds"),
             format = "file"),

  ## build_se_lab_tcga() matches the "additional JHU" columns of orse against a
  ## manifest, including 12 WGS-only Normal samples (e.g. CGOV463N) dropped from
  ## the public manifest during the facets join. Use the raw, pre-facets-join
  ## manifest so those samples resolve.
  tar_target(raw_manifest_file,
             here("..", "..", "hallberg2025.base", "inst", "extdata", "manifest.rds"),
             format = "file"),
  tar_target(raw_manifest, readRDS(raw_manifest_file)),
  tar_target(pkg_manifest,   { data(manifest,   package = "hallberg2025.base"); manifest   }),
  tar_target(pkg_discordant, { data(discordant, package = "hallberg2025.base"); discordant }),

  ## == Layer: IDATs -> bVals/mVals matrices (minfi) =========================
  tar_target(batch1_matrices,
             hallberg2025.base:::build_batch1_matrices(batch1_idat_dir, xreactive_file)),
  tar_target(batch2_matrices,
             hallberg2025.base:::build_batch2_matrices(batch2_idat_dir, xreactive_file)),

  ## Direct intermediate checks vs frozen matrices (tolerance -- minfi drift).
  tar_target(check_bVals_b1, compare_matrix(batch1_matrices$bVals, bVals_frozen_file)),
  tar_target(check_mVals_b1, compare_matrix(batch1_matrices$mVals, mVals_frozen_file)),
  tar_target(check_bVals_b2, compare_matrix(batch2_matrices$bVals, bVals081820_frozen_file)),
  tar_target(check_mVals_b2, compare_matrix(batch2_matrices$mVals, mVals081820_frozen_file)),

  ## Per-sample beta correlation vs frozen (assurance the version-drift is small):
  ## expect all r > 0.99 (documented: batch 1 >= 0.992, batch 2 >= 0.9997).
  tar_target(cor_bVals_b1, beta_cor(batch1_matrices$bVals, bVals_frozen_file)),
  tar_target(cor_bVals_b2, beta_cor(batch2_matrices$bVals, bVals081820_frozen_file)),

  ## == Layer: matrices -> orse (JHU SummarizedExperiment) ===================
  ## build_orse accepts matrices/targets1 as in-memory objects (rp() shim);
  ## baseDir for batch-2 targets is the batch-2 IDAT dir (holds the sample sheet).
  tar_target(orse_regenerated,
             hallberg2025.base:::build_orse(batch2_matrices$bVals, batch2_matrices$mVals,
                                           batch1_matrices$bVals, batch1_matrices$mVals,
                                           batch1_matrices$targets,
                                           methmanifest_file, batch2_idat_dir, methdat_file)),
  ## orse is integer-encoded; tol = 2 absorbs a +/-1 integer flip from matrix drift.
  tar_target(check_orse, compare_orse(orse_regenerated, se_jhu_frozen_file, tol = 2)),

  ## == Layer: (orse + branch-C TCGA) -> se_lab_tcga =========================
  tar_target(se_lab_tcga_regenerated,
             hallberg2025.base:::build_se_lab_tcga(bValsselect_file, combmetadata_file,
                                                  orse_regenerated, raw_manifest)),
  ## JHU beta cols inherit <= 2/1000 drift from orse; TCGA cols exact; labels exact.
  tar_target(check_se_lab_tcga,
             compare_se_lab_tcga(se_lab_tcga_regenerated, se_lab_tcga_frozen_file, tol = 2e-3)),

  ## == Layer: se_lab_tcga -> methylation_se (leaf) ==========================
  tar_target(methylation_se_regenerated,
             regenerate_methylation_se_from_obj(se_lab_tcga_regenerated,
                                                combmetadata_file,
                                                pkg_manifest, pkg_discordant,
                                                match_table_file, signet_file)),

  ## Anchor: verify against the committed baseline. The TERMINAL gate
  ## (tests/verify_snapshot.R at the repo root) is the acceptance authority.
  tar_target(verification,
             verify_against_baseline(methylation_se_regenerated, baseline_file)),

  ## ENV-03: report drift against the captured cluster environment
  ## (tests/snapshots/env_osmethexpdata.rds) -- minfi, the EPIC annotation and
  ## manifest packages, R, and the Bioconductor release are versioned by
  ## nothing else for this arm. A diagnostic like every other *_check target:
  ## it reports, it does not stop() the pipeline.
  tar_target(env_check, check_environment_drift())

)
