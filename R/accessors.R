.load <- function(filename) {
  path <- system.file("extdata", filename, package = "hallberg2025.meth.data",
                      mustWork = TRUE)
  readRDS(path)
}

#' Methylation SummarizedExperiment
#'
#' The methylation `SummarizedExperiment` consumed by the manuscript pipeline
#' in `cancer-genomics/hallberg2025`.  All sample identifiers are CG lab IDs.
#'
#' Card `MET-02a` ships the underlying data file (`inst/extdata/methylation_se.rds`)
#' as an exact copy of the object published from `hallberg2025.base::methylation_se`
#' — same object, digest-identical. No consumer has been repointed to this copy
#' yet (`MAN-03`), and the `hallberg2025.base` copy is not removed yet (`BASE-03`);
#' both copies are expected to coexist until that sequence completes.
#'
#' @return A `SummarizedExperiment` object.
#' @export
methylation_se <- function() .load("methylation_se.rds")
