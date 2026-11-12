#!/bin/bash -l
# One Slurm task: preparation, one chain, or final postprocessing.
set -euo pipefail
: "${SLURM_JOB_ID:?This script runs inside an allocation}"
: "${CONUS_REPO:?}" "${CONUS_RUN_DIR:?}" "${CONUS_STAGE:?}"
cd "$CONUS_REPO"
module purge
module load cpu
module load "$CONUS_GCC_MODULE"
module load "$CONUS_MATLAB_MODULE"
module list
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CONUS_CHAIN_ID=${SLURM_ARRAY_TASK_ID:-0}
export CONUS_RUNTIME_SCRATCH="/scratch/$USER/job_${SLURM_JOB_ID}/conus_${CONUS_STAGE}_${CONUS_CHAIN_ID}"
mkdir -p "$CONUS_RUNTIME_SCRATCH" "$CONUS_RUN_DIR/logs"
export TMPDIR="$CONUS_RUNTIME_SCRATCH"
# Atomic directory lock protects a chain's checkpoint against duplicate jobs.
lock_name=$CONUS_STAGE
if [[ $CONUS_STAGE == benchmark || $CONUS_STAGE == chains ]]; then lock_name="chain_${CONUS_CHAIN_ID}"; fi
lock="$CONUS_RUN_DIR/.lock_${lock_name}"
mkdir "$lock" || { echo "Run/chain already locked: $lock. Check its job.txt and squeue before removing a stale lock." >&2; exit 73; }
printf '%s %s %s\n' "$SLURM_JOB_ID" "$(hostname)" "$(date -Iseconds)" > "$lock/job.txt"
trap 'rm -f "$lock/job.txt"; rmdir "$lock" 2>/dev/null || true' EXIT
printf 'Stage=%s chain=%s job=%s host=%s start=%s\nRepo=%s\nResults=%s\nScratch=%s\n' \
    "$CONUS_STAGE" "$CONUS_CHAIN_ID" "$SLURM_JOB_ID" "$(hostname)" "$(date -Iseconds)" "$CONUS_REPO" "$CONUS_RUN_DIR" "$CONUS_RUNTIME_SCRATCH"
# Build Linux HVf only when absent. A stale existing binary must be rebuilt
# deliberately with all users of this checkout stopped; never replace it here.
if [[ $CONUS_STAGE == prepare && ! -f $CONUS_REPO/native/HVf_VTI/bin/HVf_VTI ]]; then
    build_lock="$CONUS_REPO/native/HVf_VTI/.build_lock"
    mkdir "$build_lock" || { echo "Build already locked: $build_lock" >&2; exit 73; }
    if sh native/HVf_VTI/build.sh; then rmdir "$build_lock"; else rmdir "$build_lock"; exit 1; fi
fi
for binary in surf96 hrftn96 tprep96 tdisp96 tregn96 tlegn96 tdpegn96 trftn96; do
    [[ -x $HBI_CPS_BIN/$binary ]] || { echo "Missing executable: $HBI_CPS_BIN/$binary" >&2; exit 69; }
done
[[ -x $HBI_HVF_BIN/HVf ]] || { echo "Missing isotropic validation executable: $HBI_HVF_BIN/HVf" >&2; exit 69; }
# A nonzero MATLAB exit propagates to Slurm; no background jobs hide failures.
srun --unbuffered --ntasks=1 matlab -singleCompThread -batch "run(getenv('CONUS_DRIVER'));"
printf 'Stage finished: %s\n' "$(date -Iseconds)"
