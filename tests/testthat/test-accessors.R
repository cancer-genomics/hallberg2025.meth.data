## Card MET-02a ships inst/extdata/methylation_se.rds -- an exact copy of the
## object published from ovarian.subtypes::methylation_se (digest-identical;
## verified separately as part of MET-02a's Done-when). These tests check
## the real data-integrity contract now that a file exists to check,
## mirroring OsSeqExpData/tests/testthat/test-data-integrity.R.

## PHI columns that must never appear in colData -- see repo CLAUDE.md
## "Columns stripped from published data objects".
PHI_COLUMNS <- c("basename", "pgdx_id", "bamfile", "bam_local", "genotype_id",
                  "facet_id", "size")

test_that("methylation_se is exported and is a zero-argument function", {
  expect_true(exists("methylation_se", where = asNamespace("hallberg2025.meth.data"),
                      inherits = FALSE))
  expect_true(is.function(methylation_se))
  expect_length(formals(methylation_se), 0)
})

test_that("methylation_se() loads and has expected structure", {
  se <- methylation_se()
  expect_s4_class(se, "SummarizedExperiment")
  expect_true("beta" %in% SummarizedExperiment::assayNames(se))
  expect_true(nrow(se) > 0)
  expect_true(ncol(se) > 0)
})

test_that("methylation_se() colData carries no PHI columns", {
  se <- methylation_se()
  cols <- colnames(SummarizedExperiment::colData(se))
  present <- base::intersect(PHI_COLUMNS, cols)
  expect_length(present, 0)
})

## ── Object-hash baseline ──────────────────────────────────────────────────────
##
## Hash is a SHA256 digest of the in-memory R object (not file bytes).
## digest::digest() is session-stable for this object type.
##
## IF this hash changes:
##   1. Confirm the change is intentional. MET-02a explicitly forbids
##      regenerating methylation_se from the data-raw/ pipeline to populate
##      this file -- it must stay an exact copy of the committed
##      ovarian.subtypes::methylation_se object until BASE-03 removes that
##      copy. An unexpected change here means the two copies have drifted.
##   2. If intentional (e.g. a later card re-derives the object on purpose),
##      update the baseline:
##        Rscript -e "
##          se <- hallberg2025.meth.data::methylation_se()
##          saveRDS(digest::digest(se, algo = 'sha256'),
##                  'tests/testthat/methylation_se_hash_baseline.rds')
##        "

test_that("methylation_se() object hash matches baseline", {
  baseline_file <- testthat::test_path("methylation_se_hash_baseline.rds")
  skip_if_not(file.exists(baseline_file), "baseline file not found")

  baseline <- readRDS(baseline_file)
  current <- digest::digest(methylation_se(), algo = "sha256")

  expect_identical(
    current, baseline,
    info = paste0(
      "methylation_se() hash changed.\n",
      "  baseline : ", baseline, "\n",
      "  current  : ", current, "\n",
      "See comment above this test for what to do."
    )
  )
})
