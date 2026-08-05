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
## BEFORE SUBMITTING: set the two IDAT-dir paths at the top of _targets.R (or the
## idat_links/*.path pointer files), and make sure hallberg2025-meth-data-raw.sif
## (built from docker/Dockerfile.meth-data-raw, ENV-17; converted+transferred here,
## ENV-18) sits alongside this script.
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
## ENV-18: the R/package layer now runs inside the pinned hallberg2025-meth-data-raw
## Apptainer image instead of `module load conda_R/4.5.x` -- SLURM submission itself
## is unchanged, per ENV-14's explicit "do not containerize SLURM submission" scope
## note. JHPCE ships the module as `singularity`, not `apptainer` (confirmed
## empirically, ENV-16) -- same docker://-derived image/CLI surface, used here.

set -euo pipefail

## Under a real `sbatch` submission on this cluster, SLURM stages/executes the batch
## script from a spool directory (e.g. /var/spool/slurm/d/job<id>), so
## ${BASH_SOURCE[0]} resolves there, NOT to this repo -- verified empirically (MET-05).
## $SLURM_SUBMIT_DIR is the one SLURM-provided value that correctly points back to the
## directory sbatch was invoked from. Fall back to BASH_SOURCE for a direct/interactive
## invocation (e.g. `bash run_targets.sh`) where SLURM_SUBMIT_DIR is unset.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "Starting hallberg2025.meth.data unified pipeline"
echo "  Node:    $(hostname)"
echo "  Time:    $(date)"
echo "  Job ID:  $SLURM_JOB_ID"
echo "  Working: $(pwd)"
echo ""

## HOSTNAME is a bash builtin, not exported to child processes by default in a
## non-login batch-script shell -- and the singularity/3.11.4 modulefile itself
## calls os.getenv("HOSTNAME") to gate "only loadable on a compute node", which
## throws a lua error ("bad argument #1 to 'match' (string expected, got nil)")
## if it's unset. Confirmed empirically (ENV-18) from a bare `srun ... bash -c`
## shell -- export it before `module load` so this always resolves regardless of
## the invoking shell.
export HOSTNAME="$(hostname)"

source /usr/share/lmod/lmod/init/bash
module use /jhpce/shared/community/modulefiles
module use /jhpce/shared/jhpce/modulefiles
set +u
module load singularity/3.11.4
set -u

IMAGE="$SCRIPT_DIR/hallberg2025-meth-data-raw.sif"
## The two IDAT dirs (idat_links/*.path, MET-05/MET-06's env-var indirection) share
## this parent -- bind it once. singularity does NOT auto-bind cluster storage
## mounts like /dcs07 or /dcs11 (confirmed empirically, ENV-18 -- a container without
## explicit --bind can't even see them), so both this and PROJECT_ROOT (covering the
## repo-relative `../..` reads: hallberg2025.base, tests/) must be passed explicitly.
IDAT_ROOT=/dcs07/scharpf/data/data-warehouse/methyl-epic
## The batch2 frozen-frontier symlinks in extdata/ (bVals081820.rds, meth_081720,
## etc. -- see this pipeline's own docs, "Batch2 history") resolve to absolute
## /dcl01/scharpf1/... targets. Same missing-auto-bind issue as /dcs07/dcs11 above,
## confirmed empirically (ENV-18): without this, targets sees them as missing files
## even though they resolve fine on the bare host.
DCL01_ROOT=/dcl01/scharpf1

## data-raw/.Rprofile re-sources the project root's .Rprofile (see its own header),
## which activates uvr's host library (.uvr/library) -- exactly the "global dotfile
## clobbers .libPaths()" landmine already known elsewhere in this project, but hit
## here for a new reason: R sources .Rprofile from cwd regardless of container
## boundaries, since the bind-mounted host filesystem makes it visible inside the
## container too. Left unset, the container's own baked-in packages disappear behind
## uvr's (non-functional-from-inside-the-container) library and every Rscript call
## below fails with "object 'minfi' not found" (confirmed empirically, ENV-18).
## SINGULARITYENV_ is how singularity passes a var into the container's environment.
export SINGULARITYENV_R_PROFILE_USER=/dev/null

sing() {
  singularity exec --bind "$PROJECT_ROOT" --bind "$IDAT_ROOT" --bind "$DCL01_ROOT" "$IMAGE" "$@"
}

echo "Image:   $IMAGE"
echo "Rscript: $(sing which Rscript)"

## Record the pinned minfi / annotation versions alongside this run (minfi is not
## bit-reproducible across versions -- the matrix diffs are only interpretable
## against a known version). tryCatch so a missing informational package can never
## halt the job under `set -e` (cf. job 34217383).
sing Rscript -e 'for (p in c("minfi","IlluminaHumanMethylationEPICanno.ilm10b2.hg19")) cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED"), "\n")'

## Safeguard A/E: refuse to run against a drifted frozen frontier. Run from the
## project root -- verify_methylation_registry.R uses here(), and data-raw/ has its
## own .here anchor, so here() must resolve at the root, not here.
( cd "$PROJECT_ROOT" && sing Rscript tests/verify_methylation_registry.R )

sing Rscript -e "targets::tar_make(callr_function = NULL, reporter = 'verbose')"

echo ""
echo "== Matrix checks (IDATs -> matrices; minfi -> expect within_tol, watch set overlap) =="
for t in check_bVals_b1 check_mVals_b1 check_bVals_b2 check_mVals_b2; do
  echo "--- $t ---"
  sing Rscript -e "print(targets::tar_read_raw('$t'))"
done

echo ""
echo "== beta correlation vs frozen (assurance; expect all_above = TRUE, min_r > 0.99) =="
for t in cor_bVals_b1 cor_bVals_b2; do
  echo "--- $t ---"
  sing Rscript -e "x <- targets::tar_read_raw('$t'); cat('n_samples', x\$n_samples, '| min_r', signif(x\$min_r,5), '| median_r', signif(x\$median_r,5), '| all_above', x\$all_above, '\n'); cat('lowest 5:', paste(names(head(x\$r,5)), signif(head(x\$r,5),5), sep='=', collapse=' '), '\n')"
done

echo ""
echo "== check_orse (orse vs frozen se.rds; integer-encoded, tol = 2) =="
sing Rscript -e "print(targets::tar_read_raw('check_orse'))"

echo ""
echo "== check_se_lab_tcga (JHU beta tol 2e-3, TCGA exact, labels exact) =="
sing Rscript -e "print(targets::tar_read_raw('check_se_lab_tcga'))"

echo ""
echo "== verification (methylation_se vs baseline) =="
sing Rscript -e "print(targets::tar_read_raw('verification'))"

echo ""
echo "Unified pipeline complete: $(date)"
echo "Then, on the laptop, confirm the TERMINAL gate: Rscript tests/verify_snapshot.R"
echo "Sync results to laptop with: unison ovarian2025-login"
