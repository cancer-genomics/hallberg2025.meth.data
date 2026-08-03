#!/usr/bin/env bash
#SBATCH --job-name=hallberg2025.meth.data
#SBATCH --partition=cancergen,shared
#SBATCH --time=18:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rscharpf@jhu.edu

## Unified methylation provenance-verification pipeline (see _targets.R).
## Runs the whole DAG from the JHU IDATs down to methylation_se in one store.
##
## BEFORE SUBMITTING: set the two IDAT-dir paths at the top of _targets.R.
##
## Invoke from this directory so the relative --output path above and SCRIPT_DIR
## below resolve correctly, e.g.:
##
##   cd /path/to/2025.ovarian.subtypes/hallberg2025.meth.data/data-raw
##   sbatch run_targets.sh
##
## 128G/18h: read.metharray.exp + preprocessFunnorm on the full EPIC IDAT sets are
## the memory peak, and build_orse() holds both full beta+M matrices plus the
## 340 MB orse. Adjust if it OOMs or finishes fast. A warm re-run (only downstream
## targets invalidated) is far cheaper -- targets reuses the cached matrices/orse.
##
## Module init: conda_R/4.5.x is in /jhpce/shared/community/modulefiles; its
## conda/3-24.3.0 dependency is in /jhpce/shared/jhpce/modulefiles. Both paths
## must be added via `module use` before `module load`.

set -euo pipefail

## Under a real `sbatch` submission on this cluster, SLURM stages/executes the batch
## script from a spool directory (e.g. /var/spool/slurm/d/job<id>), so
## ${BASH_SOURCE[0]} resolves there, NOT to this repo -- verified empirically (MET-05).
## $SLURM_SUBMIT_DIR is the one SLURM-provided value that correctly points back to the
## directory sbatch was invoked from. Fall back to BASH_SOURCE for a direct/interactive
## invocation (e.g. `bash run_targets.sh`) where SLURM_SUBMIT_DIR is unset.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

echo "Starting hallberg2025.meth.data unified pipeline"
echo "  Node:    $(hostname)"
echo "  Time:    $(date)"
echo "  Job ID:  $SLURM_JOB_ID"
echo "  Working: $(pwd)"
echo ""

source /usr/share/lmod/lmod/init/bash
module use /jhpce/shared/community/modulefiles
module use /jhpce/shared/jhpce/modulefiles
set +u
module load conda_R/4.5.x
set -u

echo "Rscript: $(which Rscript)"

## Record the pinned minfi / annotation versions alongside this run (minfi is not
## bit-reproducible across versions -- the matrix diffs are only interpretable
## against a known version). tryCatch so a missing informational package can never
## halt the job under `set -e` (cf. job 34217383).
Rscript -e 'for (p in c("minfi","IlluminaHumanMethylationEPICanno.ilm10b2.hg19")) cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED"), "\n")'

## Safeguard A/E: refuse to run against a drifted frozen frontier. Run from the
## project root -- verify_methylation_registry.R uses here(), and data-raw/ has its
## own .here anchor, so here() must resolve at the root, not here.
( cd "$SCRIPT_DIR/../.." && Rscript tests/verify_methylation_registry.R )

Rscript -e "targets::tar_make(callr_function = NULL, reporter = 'verbose')"

echo ""
echo "== Matrix checks (IDATs -> matrices; minfi -> expect within_tol, watch set overlap) =="
for t in check_bVals_b1 check_mVals_b1 check_bVals_b2 check_mVals_b2; do
  echo "--- $t ---"
  Rscript -e "print(targets::tar_read_raw('$t'))"
done

echo ""
echo "== beta correlation vs frozen (assurance; expect all_above = TRUE, min_r > 0.99) =="
for t in cor_bVals_b1 cor_bVals_b2; do
  echo "--- $t ---"
  Rscript -e "x <- targets::tar_read_raw('$t'); cat('n_samples', x\$n_samples, '| min_r', signif(x\$min_r,5), '| median_r', signif(x\$median_r,5), '| all_above', x\$all_above, '\n'); cat('lowest 5:', paste(names(head(x\$r,5)), signif(head(x\$r,5),5), sep='=', collapse=' '), '\n')"
done

echo ""
echo "== check_orse (orse vs frozen se.rds; integer-encoded, tol = 2) =="
Rscript -e "print(targets::tar_read_raw('check_orse'))"

echo ""
echo "== check_se_lab_tcga (JHU beta tol 2e-3, TCGA exact, labels exact) =="
Rscript -e "print(targets::tar_read_raw('check_se_lab_tcga'))"

echo ""
echo "== verification (methylation_se vs baseline) =="
Rscript -e "print(targets::tar_read_raw('verification'))"

echo ""
echo "Unified pipeline complete: $(date)"
echo "Then, on the laptop, confirm the TERMINAL gate: Rscript tests/verify_snapshot.R"
echo "Sync results to laptop with: unison ovarian2025"
