% Hand-written (roxygen2 is not available in this project's renv library --
% see cards/MET-01-make-package-installable.md Notes/log).
% Please edit documentation in R/accessors.R and keep this file in sync.
\name{methylation_se}
\alias{methylation_se}
\title{Methylation SummarizedExperiment}
\usage{
methylation_se()
}
\value{
A \code{SummarizedExperiment} object.
}
\description{
The methylation \code{SummarizedExperiment} consumed by the manuscript
pipeline in \code{cancer-genomics/hallberg2025}.  All sample identifiers are
CG lab IDs.
}
\details{
Card \code{MET-02a} ships the underlying data file
(\code{inst/extdata/methylation_se.rds}) as an exact copy of the object
published from \code{ovarian.subtypes::methylation_se} -- same object,
digest-identical. No consumer has been repointed to this copy yet
(\code{MAN-03}), and the \code{ovarian.subtypes} copy is not removed yet
(\code{BASE-03}); both copies are expected to coexist until that sequence
completes.
}
