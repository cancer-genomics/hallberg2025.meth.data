## Two-way delta test for CHANGES.md Change 002 (Step 0, bullet 3 of
## code/methylation_provenance_verification_plan.md).
##
## A dedicated targets pipeline, NOT part of the production _targets.R --
## isolates the Change 002 code fix (the two reverted lines in
## build_se_lab_tcga_buggy(), R/historical_buggy_variant.R) from everything
## else, using the SAME frozen frontier and the SAME downstream Step-1
## chain (regenerate_methylation_se_from_obj()) for both variants. Runs in
## a fixed, reproducible package-attach order (tar_option_set below), which
## a prior ad hoc REPL reconstruction of this same comparison did not have --
## that inconsistency in SummarizedExperiment::colData<- behavior across
## script runs is why this exists as a targets pipeline rather than a
## one-off script.
##
## Run from OsMethExpData/data-raw/:
##   Rscript -e "targets::tar_make(script = 'change002_delta_targets.R',
##                                  store = '_targets_change002_delta',
##                                  callr_function = NULL)"
##   Rscript -e "targets::tar_read_raw('delta_report',
##                                     store = '_targets_change002_delta')"
##
## Uses a separate store (_targets_change002_delta) so it never collides
## with the production Step 2 pipeline's _targets/ cache.

library(targets)
library(here)

pkgload::load_all(here("..", "..", "ovarian.subtypes"), quiet = TRUE)

tar_option_set(packages = c("SummarizedExperiment", "S4Vectors",
                             "tibble", "dplyr", "readr", "magrittr"))

tar_source("R/generate.R")
tar_source("R/historical_buggy_variant.R")
tar_source("R/change002_delta_report.R")

list(

  tar_target(bValsselect_file,
             here("..", "..", "extdata", "bValsselect.rds"),
             format = "file"),

  tar_target(combmetadata_file,
             here("..", "..", "extdata", "combmetadata.rds"),
             format = "file"),

  tar_target(se_jhu_file,
             here("..", "..", "output", "methylation.Rmd", "se.rds"),
             format = "file"),

  tar_target(match_table_file,
             here("..", "..", "ovarian.subtypes", "inst", "extdata", "match_table.csv"),
             format = "file"),

  tar_target(signet_file,
             here("..", "..", "ovarian.subtypes", "inst", "extdata", "stomach_muc_signet.csv"),
             format = "file"),

  tar_target(raw_manifest_file,
             here("..", "..", "ovarian.subtypes", "inst", "extdata", "manifest.rds"),
             format = "file"),

  tar_target(raw_manifest, readRDS(raw_manifest_file)),

  tar_target(pkg_manifest,   { data(manifest,   package = "ovarian.subtypes"); manifest   }),
  tar_target(pkg_discordant, { data(discordant, package = "ovarian.subtypes"); discordant }),

  ## -- The two variants under comparison; identical inputs, differ only in
  ##    which build_se_lab_tcga is used. --------------------------------------

  tar_target(se_lab_tcga_corrected,
             ovarian.subtypes:::build_se_lab_tcga(bValsselect_file, combmetadata_file,
                                                   se_jhu_file, raw_manifest)),

  tar_target(se_lab_tcga_buggy,
             build_se_lab_tcga_buggy(bValsselect_file, combmetadata_file,
                                      se_jhu_file, raw_manifest)),

  tar_target(meth_se_corrected,
             regenerate_methylation_se_from_obj(se_lab_tcga_corrected, combmetadata_file,
                                                pkg_manifest, pkg_discordant,
                                                match_table_file, signet_file)),

  tar_target(meth_se_buggy,
             regenerate_methylation_se_from_obj(se_lab_tcga_buggy, combmetadata_file,
                                                pkg_manifest, pkg_discordant,
                                                match_table_file, signet_file)),

  tar_target(delta_report,
             report_change002_delta(se_lab_tcga_corrected, se_lab_tcga_buggy,
                                     meth_se_corrected, meth_se_buggy))

)
