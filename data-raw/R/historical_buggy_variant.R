## DIAGNOSTIC ONLY -- not part of the production pipeline.
##
## Verbatim reconstruction of the PRE-Change-002 build_se_lab_tcga() logic
## (the archived methylation_summarized_experiment.Rmd's positional colData
## assignment), used only by change002_delta_targets.R to verify the Change
## 002 fix against an independent, non-circular comparison. See
## code/methylation_provenance_verification_plan.md, "Verification
## safeguards" / Step 0.
##
## Differs from hallberg2025.base:::build_se_lab_tcga() (hallberg2025.base/R/
## methylation.R) in exactly two respects, both in the "additional JHU"
## block:
##   - no dplyr::arrange(match(lab_id, colnames(se.jhu3))) reordering
##   - no explicit row.names= on the colData<- DataFrame() call
## Everything else is copied verbatim so the only variable between the two
## functions is the fix itself.
##
## NEVER call this outside change002_delta_targets.R.

build_se_lab_tcga_buggy <- function(bValsselect_file, combmetadata_file, se_jhu_file,
                                     manifest) {
  dlevels <- c("Uterine endometrial", "Ovarian endometrioid",
               "Ovarian mucinous", "Colorectal mucinous",
               "Pancreas mucinous", "Stomach mucinous")
  manifest2 <- manifest %>%
    hallberg2025.base:::cancer_names() %>%
    dplyr::select(-tumor_type) %>%
    dplyr::rename(tumor_type = tumor, tumor = tumor.normal) %>%
    dplyr::mutate(study = "JHU", diagnosis = factor(tumor_type, dlevels),
                  tumor = Hmisc::capitalize(tumor))

  lab.and.tcga <- readRDS(bValsselect_file) %>% t()
  metadata <- readRDS(combmetadata_file) %>% tibble::as_tibble()

  df <- tibble::tibble(
    lab_id    = metadata$Sample_Name,
    diagnosis = metadata$Diagnosis,
    tissue    = metadata$Tissue,
    type      = metadata$T.N,
    study     = metadata$batch
  ) %>%
    dplyr::mutate(
      study = ifelse(study %in% 1:2, "JHU", "TCGA"),
      tumor = ifelse(type == "T", "Tumor", "Normal"),
      diagnosis = factor(diagnosis, dlevels)
    ) %>%
    dplyr::select(lab_id, diagnosis, tumor, study)

  se.lab.tcga <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(beta = lab.and.tcga),
    colData = df
  )
  colnames(se.lab.tcga) <- df$lab_id

  se.jhu <- readRDS(se_jhu_file)
  rowindex <- rownames(se.jhu) %in% rownames(se.lab.tcga)
  colindex <- which(!colnames(se.jhu) %in% colnames(se.lab.tcga))
  se.jhu2 <- se.jhu[rowindex, colindex]
  SummarizedExperiment::assays(se.jhu2) <- S4Vectors::SimpleList(
    beta = SummarizedExperiment::assay(se.jhu2) / 1000
  )
  se.jhu3 <- se.jhu2[rownames(se.lab.tcga), ]
  ## BUGGY (pre-fix): no arrange-by-physical-position, no explicit row.names.
  df.addl.jhu <- dplyr::filter(manifest2, lab_id %in% colnames(se.jhu3)) %>%
    dplyr::select(lab_id, diagnosis, tumor, study)
  colnames_before <- colnames(se.jhu3)
  SummarizedExperiment::colData(se.jhu3) <- S4Vectors::DataFrame(df.addl.jhu)
  ## SummarizedExperiment::colData<- on a rowname-less DataFrame has been
  ## observed (2026-07-04) to inconsistently clear colnames(se.jhu3)
  ## depending on package-attach state -- restore them explicitly so this
  ## function's behavior does not depend on that ambiguity. This mirrors
  ## what colnames(se_lab_tcga) <- se_lab_tcga$lab_id (in generate.R's
  ## regenerate_methylation_se_from_obj()) does immediately downstream in
  ## the historical/production chain regardless, so restoring here changes
  ## nothing about the bug itself -- it only makes this function's own
  ## return value deterministic to inspect in isolation.
  colnames(se.jhu3) <- colnames_before

  cbind(se.lab.tcga, se.jhu3)
}
