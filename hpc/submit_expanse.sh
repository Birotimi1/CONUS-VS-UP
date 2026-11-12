#!/usr/bin/env bash
# Submit one explicit stage. This file does not SSH or submit automatically.
set -euo pipefail
stage=${1:?Usage: submit_expanse.sh prepare|benchmark|chains|postprocess RUN_ID [--dry-run]}
run_id=${2:?Supply a unique run ID, reused for benchmark/chains/postprocess}
[[ $run_id =~ ^[A-Za-z0-9_-]+$ ]] || { echo 'Use letters, digits, underscore or hyphen in RUN_ID.' >&2; exit 64; }
case "$stage" in prepare|benchmark|chains|postprocess) ;; *) echo 'Unknown stage' >&2; exit 64;; esac
repo=$(cd "$(dirname "$0")/.." && pwd)
export CONUS_REPO="$repo" CONUS_STAGE="$stage"
export CONUS_RUN_DIR="${CONUS_RESULTS_ROOT:-$repo/DATA/transdimensional}/$run_id"
export HBI_CPS_BIN="${HBI_CPS_BIN:-/expanse/lustre/projects/syu127/brotimi/Inversion/bin}"
export HBI_HVF_BIN="${HBI_HVF_BIN:-$HBI_CPS_BIN}"
export CONUS_MATLAB_MODULE="${CONUS_MATLAB_MODULE:-matlab/2023b}"
export CONUS_GCC_MODULE="${CONUS_GCC_MODULE:-gcc}"
export CONUS_CHAIN_SECONDS="${CONUS_CHAIN_SECONDS:-158400}" # 44h; checkpoint between iterations
export CONUS_DRIVER="$repo/DATA/run_expanse_250km.m"
args=(--parsable --job-name="conus_${stage}" --account="${CONUS_ACCOUNT:-syu127}"
      --partition=shared --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G
      --time=48:00:00 --constraint=lustre --export=ALL
      --output="$CONUS_RUN_DIR/logs/${stage}_%A_%a.out"
      --error="$CONUS_RUN_DIR/logs/${stage}_%A_%a.err")
if [[ $stage == benchmark || $stage == chains ]]; then args+=(--array=1-2%2); fi
if [[ -n ${CONUS_AFTER_JOB:-} ]]; then
    [[ $CONUS_AFTER_JOB =~ ^[0-9]+$ ]] || { echo 'CONUS_AFTER_JOB must be a job ID.' >&2; exit 64; }
    args+=(--dependency="afterok:$CONUS_AFTER_JOB")
fi
if [[ ${3:-} == --dry-run ]]; then
    printf 'Run directory: %s\n' "$CONUS_RUN_DIR"
    printf 'sbatch'; printf ' %q' "${args[@]}" "$repo/hpc/expanse_job.sh"; printf '\n'
    exit 0
fi
[[ $# -eq 2 ]] || { echo 'Unexpected argument' >&2; exit 64; }
command -v sbatch >/dev/null || { echo 'Run this script on Expanse after authenticating.' >&2; exit 69; }
mkdir -p "$CONUS_RUN_DIR/logs"
job=$(sbatch "${args[@]}" "$repo/hpc/expanse_job.sh")
printf '%s stage=%s job=%s run=%s\n' "$(date -Iseconds)" "$stage" "$job" "$CONUS_RUN_DIR" | tee -a "$CONUS_RUN_DIR/logs/submissions.log"
