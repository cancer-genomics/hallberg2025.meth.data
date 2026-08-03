#!/usr/bin/env Rscript
## Step-5 diagnostic: scatterplot of REGENERATED (current, minfi 1.56.0) beta
## values vs the ORIGINAL PUBLISHED beta values, for a random sample of arrays.
##
## Runs fully locally: the regenerated matrices come from the synced
## _targets store; the published matrices are the frozen extdata/bVals*.rds
## (pulled from the cluster -- CG-lab-ID colnames only, no PHI). Both carry
## cleaned CG-ID colnames after the build_batch1_matrices fix, so they align.
##
## Batch 1 is the batch whose funnorm output drifts across minfi versions
## (see the Step-5 diagnosis); set BATCH <- 2 to view the reproducing batch.
##
##   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript \
##     hallberg2025.meth.data/data-raw/plot_beta_concordance.R

suppressMessages({
  library(targets)
  library(ggplot2)
})

## ---- parameters ------------------------------------------------------------
BATCH        <- 1        # 1 (funnorm version-drift) or 2 (reproduces)
N_SAMPLES    <- 5        # number of arrays to draw
SEED         <- 1        # reproducible random draw
N_PROBES_PLOT <- 30000L  # probes subsampled per panel for the scatter (r/maxdiff use ALL probes)

## ---- locate inputs (works under Rscript, source(), or line-by-line) --------
## Walk up from the current dir to the project root: the dir that holds both
## extdata/ and hallberg2025.meth.data/. Robust regardless of the session's getwd().
find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(d, "extdata")) &&
        dir.exists(file.path(d, "hallberg2025.meth.data"))) return(d)
    parent <- dirname(d)
    if (parent == d)
      stop("Could not locate project root (a dir containing extdata/ and hallberg2025.meth.data/) ",
           "above ", start)
    d <- parent
  }
}
PROJ    <- find_root()
data_raw <- file.path(PROJ, "hallberg2025.meth.data", "data-raw")
store    <- file.path(data_raw, "_targets")
extdata  <- file.path(PROJ, "extdata")

regen_target <- sprintf("batch%d_matrices", BATCH)
frozen_file  <- file.path(extdata, if (BATCH == 1) "bVals.rds" else "bVals081820.rds")

## ---- load ------------------------------------------------------------------
message("Reading regenerated ", regen_target, " from ", store)
regen  <- tar_read_raw(regen_target, store = store)$bVals
message("Reading published    ", basename(frozen_file))
pub    <- readRDS(frozen_file)

## align on shared probes + samples
cr <- intersect(rownames(regen), rownames(pub))
cc <- intersect(colnames(regen), colnames(pub))
stopifnot(length(cr) > 0, length(cc) >= N_SAMPLES)
regen <- regen[cr, cc]; pub <- pub[cr, cc]

## ---- random sample draw ----------------------------------------------------
set.seed(SEED)
samp <- sort(sample(cc, N_SAMPLES))
message("Random ", N_SAMPLES, " samples (seed ", SEED, "): ", paste(samp, collapse = ", "))

## ---- per-sample stats on ALL probes ---------------------------------------

samp <- colnames(pub)
stats <- do.call(rbind, lapply(samp, function(s) {
  x <- pub[, s]; y <- regen[, s]
  ok <- is.finite(x) & is.finite(y)
  data.frame(sample = s,
             r        = cor(x[ok], y[ok]),
             max_diff = max(abs(y[ok] - x[ok])),
             frac_gt_0.1 = mean(abs(y[ok] - x[ok]) > 0.1))
}))
print(stats, row.names = FALSE, digits = 4)

## ---- long df for plotting (subsample probes for a legible scatter) ---------
set.seed(SEED + 1)
idx <- if (length(cr) > N_PROBES_PLOT) sort(sample(length(cr), N_PROBES_PLOT)) else seq_along(cr)
dat <- do.call(rbind, lapply(samp, function(s) {
  data.frame(sample = s, published = pub[idx, s], current = regen[idx, s])
}))
dat <- dat[is.finite(dat$published) & is.finite(dat$current), ]

## panel labels with r / max|diff|
labs <- setNames(sprintf("%s\n(r=%.3f, max|Δ|=%.2f)",
                         stats$sample, stats$r, stats$max_diff), stats$sample)

## ---- plot ------------------------------------------------------------------
p <- ggplot(dat, aes(published, current)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  geom_point(alpha = 0.12, size = 0.35, colour = "#225588") +
  facet_wrap(~ sample, nrow = 1, labeller = labeller(sample = labs)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = sprintf("Batch %d beta concordance: current (minfi 1.56.0) vs published", BATCH),
    subtitle = sprintf("%s random arrays (seed %d); %s probes shown / %s total; red = y=x",
                       N_SAMPLES, SEED, format(length(idx), big.mark = ","),
                       format(length(cr), big.mark = ",")),
    x = "Published beta (original)", y = "Current beta (regenerated)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

out <- file.path(data_raw, sprintf("beta_concordance_batch%d.png", BATCH))
ggsave(out, p, width = 3.1 * N_SAMPLES, height = 3.6, dpi = 150)
message("Wrote ", normalizePath(out))
