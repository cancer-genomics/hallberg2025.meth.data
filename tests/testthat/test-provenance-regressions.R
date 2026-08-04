## Card MET-03 -- per-layer regressions for the methylation provenance chain.
##
## The chain (hallberg2025.meth.data/data-raw/_targets.R, data-raw/CLAUDE.md "The DAG")
## was cleared one layer at a time -- matrices, orse, se_lab_tcga, leaf -- but
## each clearance exists only as a diagnostic `*_check`/`verification` target
## inside the pipeline (deliberately non-enforcing: they return structures,
## they never stop() -- see _targets.R:46-48) and a note in the plan. Nothing
## re-checks them, so a future change could silently un-clear a layer and
## nothing would fail. These tests are that regression net: they call the
## SAME compare_*()/verify_against_baseline() functions the pipeline itself
## uses (sourced from data-raw/R/generate.R, never reimplemented here) against
## real values cached in the data-raw/_targets/ store, and assert the
## DOCUMENTED TOLERANCE for each layer (data-raw/CLAUDE.md "Tolerance
## regime") -- not bit-exactness.
##
## They must never trigger a rebuild: a test that runs a cluster pipeline is
## not a test. Each layer skips cleanly, with an informative reason, when the
## cached target it needs isn't available locally -- covering two distinct
## cases: (1) no _targets/ store at all (fresh laptop clone, CI), and (2) a
## store that exists but predates the current DAG and so lacks specific
## target names (see "Evidence" below -- true of this repo's checked-out
## store for the matrices/orse layers as of 2026-07-28).

## generate.R defines compare_matrix()/compare_orse()/compare_se_lab_tcga()/
## beta_cor()/verify_against_baseline() and is tar_source()'d by _targets.R --
## it is pipeline code, not a package, so there is nothing to library(); pull
## the definitions in directly so these tests exercise the pipeline's own
## comparison logic rather than a reimplementation of it. This only defines
## functions (safe even when the store is absent/targets isn't installed).
source(testthat::test_path("..", "..", "data-raw", "R", "generate.R"))

STORE_PATH <- testthat::test_path("..", "..", "data-raw", "_targets")

skip_unless_store_available <- function() {
  testthat::skip_if_not(
    requireNamespace("targets", quietly = TRUE),
    "targets package not installed -- provenance regressions need it to read the cached store"
  )
  testthat::skip_if_not(
    dir.exists(STORE_PATH),
    paste0(
      "hallberg2025.meth.data/data-raw/_targets/ store not found at ", STORE_PATH, " -- ",
      "provenance regressions compare against a cluster-populated cache (see ",
      "data-raw/CLAUDE.md 'Running it'); this is expected on a laptop clone ",
      "or in CI and is not a failure. Skipping."
    )
  )
  ## Warm-up: dim()/nrow()/ncol() on a SummarizedExperiment that was just
  ## deserialized by tar_read_raw() in a session where SummarizedExperiment's
  ## namespace has not yet been touched can return NULL on the FIRST S4
  ## dispatch and the correct value on every call after (reproduced 3/3 tries
  ## against this store, 2026-07-28) -- a cold-dispatch quirk of the object,
  ## not of the data. Forcing the namespace to load before any compare_*()
  ## call avoids a spurious "same_dim = FALSE" that has nothing to do with an
  ## actual regression.
  requireNamespace("SummarizedExperiment", quietly = TRUE)
}

## Reads a target from the explicit store path; returns NULL (never errors)
## when the name isn't cached, so callers can skip with a specific reason
## instead of failing the whole file on the first missing target.
read_cached_target <- function(name) {
  tryCatch(
    targets::tar_read_raw(name, store = STORE_PATH),
    error = function(e) NULL
  )
}

skip_if_missing <- function(values, names) {
  missing <- names[vapply(values, is.null, logical(1))]
  testthat::skip_if(
    length(missing) > 0,
    paste0(
      "target(s) ", paste(missing, collapse = ", "), " not cached in ",
      STORE_PATH, " -- this store predates the current unified DAG for this ",
      "layer (data-raw/_targets.R's header documents the old per-step -> ",
      "unified-DAG rename) and needs a fresh cluster-side tar_make() to ",
      "populate; not a regression. Skipping this layer here."
    )
  )
}

## ── Layer 1: matrices (minfi funnorm output) ───────────────────────────────
## Tolerance: compare_matrix() -- set overlap on probes/samples plus a finite
## max abs diff, NOT bit-exactness (minfi funnorm is not bit-reproducible
## across package versions; data-raw/CLAUDE.md). Do not tighten this to
## diff == 0 -- that fails on the next minfi release for a reason that means
## nothing.
test_that("matrices layer: batch1 bVals/mVals stay within compare_matrix()'s documented tolerance", {
  skip_unless_store_available()
  b1 <- read_cached_target("batch1_matrices")
  bVals_frozen <- read_cached_target("bVals_frozen_file")
  mVals_frozen <- read_cached_target("mVals_frozen_file")
  skip_if_missing(list(b1, bVals_frozen, mVals_frozen),
                  c("batch1_matrices", "bVals_frozen_file", "mVals_frozen_file"))

  bVals_res <- compare_matrix(b1$bVals, bVals_frozen)
  mVals_res <- compare_matrix(b1$mVals, mVals_frozen)

  expect_gt(bVals_res$n_common_probes, 0)
  expect_gt(bVals_res$n_common_samples, 0)
  expect_true(bVals_res$within_tol,
              info = paste0("batch1 bVals finite_maxdiff=", bVals_res$finite_maxdiff))
  expect_true(mVals_res$within_tol,
              info = paste0("batch1 mVals finite_maxdiff=", mVals_res$finite_maxdiff,
                             " n_nonfinite_mismatch=", mVals_res$n_nonfinite_mismatch))
})

test_that("matrices layer: batch2 bVals/mVals stay within compare_matrix()'s documented tolerance", {
  skip_unless_store_available()
  b2 <- read_cached_target("batch2_matrices")
  bVals_frozen <- read_cached_target("bVals081820_frozen_file")
  mVals_frozen <- read_cached_target("mVals081820_frozen_file")
  skip_if_missing(list(b2, bVals_frozen, mVals_frozen),
                  c("batch2_matrices", "bVals081820_frozen_file", "mVals081820_frozen_file"))

  bVals_res <- compare_matrix(b2$bVals, bVals_frozen)
  mVals_res <- compare_matrix(b2$mVals, mVals_frozen)

  expect_gt(bVals_res$n_common_probes, 0)
  expect_gt(bVals_res$n_common_samples, 0)
  expect_true(bVals_res$within_tol,
              info = paste0("batch2 bVals finite_maxdiff=", bVals_res$finite_maxdiff))
  expect_true(mVals_res$within_tol,
              info = paste0("batch2 mVals finite_maxdiff=", mVals_res$finite_maxdiff,
                             " n_nonfinite_mismatch=", mVals_res$n_nonfinite_mismatch))
})

## ── Layer 2: orse (JHU SummarizedExperiment from the matrices) ────────────
## Tolerance: compare_orse(tol = 2) -- orse is integer-encoded (beta x1000,
## M x100), so sub-permille drift inherited from the matrices layer can flip
## an integer by +/-1; tol = 2 absorbs that without masking a real regression.
test_that("orse layer: compare_orse(tol = 2) stays within the documented integer tolerance", {
  skip_unless_store_available()
  orse <- read_cached_target("orse_regenerated")
  frozen <- read_cached_target("se_jhu_frozen_file")
  skip_if_missing(list(orse, frozen), c("orse_regenerated", "se_jhu_frozen_file"))

  res <- compare_orse(orse, frozen, tol = 2)

  expect_true(res$same_dim, info = "orse dims differ from the frozen se.rds")
  expect_true(res$same_rownames_set)
  expect_true(res$same_colnames_set)
  expect_true(res$assays_equal,
              info = paste0("orse beta/M max abs diff exceeds tol=2: ",
                             paste(names(res$assay_max_absdiff), res$assay_max_absdiff,
                                   sep = "=", collapse = ", ")))
})

## ── Layer 3: se_lab_tcga (orse + frozen branch-C TCGA) ─────────────────────
## Tolerance: compare_se_lab_tcga(tol = 2e-3) on assay VALUES (JHU beta cols
## inherit <= 2/1000 drift from orse; TCGA cols are exact/frozen) -- but
## colData LABELS are asserted EXACT. This is deliberately the strict half of
## this test: labels are the Change-002 surface (CHANGES.md) -- a positional
## (not name-based) colData assignment silently swapped diagnosis/tumor/study
## labels across "additional JHU" samples while every value stayed internally
## consistent-looking. A value tolerance would never have caught that; only
## an exact label check does.
test_that("se_lab_tcga layer: values within tol = 2e-3, colData labels exact (Change-002 surface)", {
  skip_unless_store_available()
  regen <- read_cached_target("se_lab_tcga_regenerated")
  frozen <- read_cached_target("se_lab_tcga_frozen_file")
  skip_if_missing(list(regen, frozen), c("se_lab_tcga_regenerated", "se_lab_tcga_frozen_file"))

  res <- compare_se_lab_tcga(regen, frozen, tol = 2e-3)

  expect_true(res$same_dim)
  expect_true(res$same_rownames_set)
  expect_true(res$same_colnames_set)
  expect_true(res$assays_equal,
              info = paste0("se_lab_tcga beta max abs diff exceeds tol=2e-3: ",
                             res$assay_max_absdiff))
  expect_true(res$coldata_equal,
              info = paste(
                "se_lab_tcga colData labels (lab_id/diagnosis/tumor/study) no longer",
                "match the frozen twin EXACTLY -- this is exactly the Change-002 failure",
                "mode (a positional, not name-based, label assignment); see CHANGES.md."
              ))
})

## ── Layer 4: leaf (methylation_se, the published object) ──────────────────
## Tolerance: verify_against_baseline() against the committed
## hallberg2025.base/data/methylation_se.rda -- the terminal gate. Accept
## either bit-identical (the expected/normal case) or, failing that,
## all.equal()-level numeric agreement on the assay; anything less means the
## published object itself has drifted from what hallberg2025.meth.data ships.
test_that("leaf layer: methylation_se matches the committed baseline (verify_against_baseline)", {
  skip_unless_store_available()
  regen <- read_cached_target("methylation_se_regenerated")
  baseline_path <- read_cached_target("baseline_file")
  skip_if_missing(list(regen, baseline_path),
                  c("methylation_se_regenerated", "baseline_file"))

  res <- verify_against_baseline(regen, baseline_path)

  expect_true(res$same_dim)
  expect_true(res$same_rownames_set)
  expect_true(res$same_colnames_set)
  expect_true(res$identical || isTRUE(res$assay_all_equal),
              info = "methylation_se_regenerated no longer matches the committed baseline, even under all.equal()")
})
