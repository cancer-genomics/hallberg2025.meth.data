## DIAGNOSTIC ONLY -- reporting helper for change002_delta_targets.R.
## See code/methylation_provenance_verification_plan.md, "Verification
## safeguards" / Step 0, bullet 3 (two-way delta test).

report_change002_delta <- function(se_lab_tcga_corrected, se_lab_tcga_buggy,
                                    meth_se_corrected, meth_se_buggy) {
  cd_c <- as.data.frame(SummarizedExperiment::colData(se_lab_tcga_corrected))
  cd_b <- as.data.frame(SummarizedExperiment::colData(se_lab_tcga_buggy))
  cn_c <- colnames(se_lab_tcga_corrected)
  cn_b <- colnames(se_lab_tcga_buggy)

  ## Self-consistency: in the corrected object, colData()$lab_id must equal
  ## the physical column name at every position. In the buggy object, this
  ## is exactly where Change 002's bug shows up.
  self_consistent_corrected <- identical(cd_c$lab_id, cn_c)
  buggy_mismatch_ids <- sort(cn_b[cd_b$lab_id != cn_b])

  ## Final methylation_se delta, using explicit shared row+col names (not
  ## positional indexing) so a partial row/col reorder can't produce a
  ## spurious "difference".
  shared_cols <- intersect(colnames(meth_se_corrected), colnames(meth_se_buggy))
  shared_rows <- intersect(rownames(meth_se_corrected), rownames(meth_se_buggy))
  ac <- SummarizedExperiment::assay(meth_se_corrected)[shared_rows, shared_cols, drop = FALSE]
  ab <- SummarizedExperiment::assay(meth_se_buggy)[shared_rows, shared_cols, drop = FALSE]
  differs <- vapply(shared_cols, function(cn) {
    !isTRUE(all.equal(ac[, cn], ab[, cn]))
  }, logical(1))

  changes_md_7 <- sort(c("CGOV166N", "CGOV167N", "CGOV169N",
                          "CGOV477N", "CGOV480N", "CGOV485N", "CGOV486N"))
  final_delta_ids <- sort(shared_cols[differs])

  list(
    n_selab_tcga_cols          = length(cn_c),
    self_consistent_corrected  = self_consistent_corrected,
    n_selab_tcga_mismatch      = length(buggy_mismatch_ids),
    selab_tcga_mismatch_ids    = buggy_mismatch_ids,
    n_final_shared             = length(shared_cols),
    n_final_delta               = sum(differs),
    final_delta_ids             = final_delta_ids,
    changes_md_documented_7    = changes_md_7,
    matches_changes_md_7       = identical(final_delta_ids, changes_md_7)
  )
}
