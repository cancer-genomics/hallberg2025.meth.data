## ENV-03: environment capture and non-fatal drift reporting for the
## methylation pipeline (hallberg2025.meth.data/data-raw).
##
## This arm has no install_deps.R / no GitHub pins of its own -- everything it
## needs (minfi, the EPIC annotation and manifest packages, SummarizedExperiment,
## the whole CRAN/Bioconductor surface) comes from whatever conda_R/4.5.x
## resolves to on the cluster on a given day. There is no hard gate to mirror
## hallberg2025.seq.data's assert_pins() here, because there is nothing pinned to gate
## against -- capture_environment()/check_environment_drift() are the whole
## story for this arm, and both are diagnostics: they report, they do not
## stop(). This mirrors the *_check / verification targets already in
## _targets.R (see CLAUDE.md's tolerance-regime table) -- a drifted package is
## a reported value, not a build error.
##
## Usage (see _targets.R's `env_check` target):
##
##   source("R/env_check.R")
##   capture_environment()              # list: sessionInfo, Bioconductor
##                                       #   release, explicit package versions
##   check_environment_drift()          # compares the live session against
##                                       #   tests/snapshots/env_hallberg2025_meth_data.rds
##                                       #   and reports (never stop()s)

## The packages ENV-03's Evidence names for this arm, plus BiocManager itself.
methylation_tracked_pkgs <- function() {
    c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
      "IlluminaHumanMethylationEPICmanifest", "SummarizedExperiment",
      "S4Vectors", "BiocManager")
}

pkg_version <- function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

## The full environment snapshot: R version, Bioconductor release, and
## explicit versions for every tracked package. This is what ENV-03 writes to
## tests/snapshots/env_hallberg2025_meth_data.rds.
capture_environment <- function() {
    pkgs <- methylation_tracked_pkgs()
    list(
        pipeline             = "hallberg2025.meth.data",
        captured_at          = Sys.time(),
        r_version            = as.character(getRversion()),
        BiocManager_version  = pkg_version("BiocManager"),
        Bioconductor_release = tryCatch(as.character(BiocManager::version()),
                                         error = function(e) NA_character_),
        sessionInfo          = utils::sessionInfo(),
        package_versions     = setNames(vapply(pkgs, pkg_version, character(1)), pkgs)
    )
}

## Compare a live session (default: the current one) against a saved capture.
## Returns a data.frame of the per-package comparison (invisibly) and, like
## every other *_check target in this DAG, reports rather than aborts.
check_environment_drift <- function(capture_path = file.path("..", "..", "tests",
                                                               "snapshots", "env_hallberg2025_meth_data.rds"),
                                     live = capture_environment()) {
    if (!file.exists(capture_path)) {
        message("No environment capture at ", capture_path, "; nothing to compare against.")
        return(invisible(NULL))
    }
    baseline <- readRDS(capture_path)

    if (!identical(baseline$r_version, live$r_version))
        message("R version drift: ", baseline$r_version, " -> ", live$r_version)
    if (!identical(baseline$Bioconductor_release, live$Bioconductor_release))
        message("Bioconductor release drift: ", baseline$Bioconductor_release,
                " -> ", live$Bioconductor_release)

    pkgs <- union(names(baseline$package_versions), names(live$package_versions))
    b <- baseline$package_versions[pkgs]
    l <- live$package_versions[pkgs]
    same <- mapply(function(x, y) isTRUE(x == y) || (is.na(x) && is.na(y)), b, l)
    out <- data.frame(pkg = pkgs, baseline = unname(b), live = unname(l),
                       drift = !same, stringsAsFactors = FALSE)

    if (any(out$drift)) {
        message("Package version drift vs captured environment:\n",
                paste(sprintf("  %-46s %s -> %s", out$pkg[out$drift],
                               ifelse(is.na(out$baseline[out$drift]), "absent", out$baseline[out$drift]),
                               ifelse(is.na(out$live[out$drift]), "absent", out$live[out$drift])),
                      collapse = "\n"))
    } else {
        message("No package drift vs captured environment (", capture_path, ").")
    }
    invisible(out)
}
