# Expanse: 250 km transdimensional synthetic run

Prepared locally for Birotimi's `syu127` allocation. **No Expanse job has been
submitted or remote installation verified.** The user requested local scripts
only after the SSH connection required renewed authentication.

The defaults follow the old HBI installation:

- Project/account: `syu127`; user: `brotimi`.
- CPS and original isotropic HVf: `/expanse/lustre/projects/syu127/brotimi/Inversion/bin`.
- MATLAB module: `matlab/2023b`; compiler: GNU `gcc` module / `gfortran`.
- Suggested new checkout: `/expanse/lustre/projects/syu127/brotimi/CONUS-VS-UP`.
- Actual repository path is derived from the submission script, not hard-coded.

Verify these on Expanse (`module avail matlab`, `module avail gcc`, allocation
availability). The previous HBI batch generator is `Codes/Expanse_codes/30k_tests`
in `HVSR-DC-JI-Inv`: it allocated a full 128-core node and backgrounded 20 MATLAB
processes. Here we use two independent one-core array tasks on `shared`, with
8 GB each and `-singleCompThread`. Extra cores do not parallelize the current
serial Fortran solver. No GPU is requested: these solvers are CPU code.

SDSC currently documents a 48-hour shared-partition limit, node-local scratch,
and `--constraint=lustre` for jobs accessing the project filesystem:
[SDSC user guide](https://www.sdsc.edu/systems/expanse/user_guide.html).
Slurm/module behavior still needs verification on the actual account.

## Install source, retaining the old HBI project

On Expanse, clone this branch into a separate new checkout if it does not already
exist. Do not copy macOS binaries or point MATLAB at the old HBI engine.

```bash
cd /expanse/lustre/projects/syu127/brotimi
git clone --branch codex/thin-five-observable-synthetic \
  https://github.com/Birotimi1/CONUS-VS-UP.git CONUS-VS-UP
cd CONUS-VS-UP
export HBI_CPS_BIN=/expanse/lustre/projects/syu127/brotimi/Inversion/bin
export HBI_HVF_BIN="$HBI_CPS_BIN"
export CONUS_ACCOUNT=syu127
export CONUS_MATLAB_MODULE=matlab/2023b
export CONUS_GCC_MODULE=gcc
```

The original isotropic HVf remains necessary for independent isotropic-limit
validation. The VTI workflow additionally needs Linux `tprep96`, `tdisp96`,
`tregn96`, `tlegn96`, `tdpegn96`, `trftn96`, plus `surf96` and `hrftn96`.
A minimal old `surf96` installation may not contain the TI/RF programs. Missing
executables stop the job; there is no isotropic fallback for an anisotropic Earth.

The preparation job compiles `native/HVf_VTI` with gfortran only if no executable
exists. Existing source/executable hash mismatches stop the run. Rebuild a stale
binary manually with `sh native/HVf_VTI/build.sh` only when no jobs use that
checkout. Build manifests support either Linux `sha256sum` or macOS `shasum`.

## Four explicit stages

Choose a unique run ID. The shell script supports `--dry-run` to show the exact
submission without calling Slurm or making result directories.

```bash
bash hpc/submit_expanse.sh prepare ECM1_250KM_VTI_001 --dry-run
bash hpc/submit_expanse.sh prepare ECM1_250KM_VTI_001
```

1. **prepare** verifies the algorithm, validates the Linux native build, generates
   observations from the ECM1/AK135 truth, checks sensitivity/RF controls, and
   writes an atomic `prepared.mat`. It has no production MCMC. The prescribed
   settings are 2 chains, 15,000 steps, 8,000 burn-in, linked crustal VTI,
   monotonic sediments and allowed crust/mantle LVZs. HVSR is 0.1–25 Hz with
   `apsv=ash=1e-3`; the other period bands extend through 100 seconds.
   Completed gate/data groups are reused only under the saved exact source,
   parameter and binary identity. Restarting preparation never silently adopts
   a changed model or executable. Individual partially completed gate groups
   may still be recomputed.
2. **benchmark** runs the two independent dispersed initial Earths, one per array
   task, and saves iteration-zero checkpoints. This is the relevant cost check
   before committing to 30,000 iterations; the true Earth can be much cheaper
   than a thick-sediment prior draw. Existing checkpoints are reused.
3. **chains** continues each independent checkpoint. Checkpoints are written
   after every completed iteration. After 44 hours, a chain stops at an iteration
   boundary. Re-submit the same stage/run ID until both reach 15,000. No Slurm
   automatic requeue or extra compute submission happens behind the user's back.
4. **postprocess** requires both completed `chain.mat` files. It assembles retained
   models, figures and diagnostics, and recomputes the final native likelihood.
   Partial chains are not reported as finished results.

Submit the next stage after inspecting the previous one, or explicitly set a
Slurm dependency for the preparation-to-benchmark transition:

```bash
# Replace 12345 with the successful/submitted preparation job ID.
CONUS_AFTER_JOB=12345 bash hpc/submit_expanse.sh benchmark ECM1_250KM_VTI_001
# After reviewing initialization cost:
bash hpc/submit_expanse.sh chains ECM1_250KM_VTI_001
# Repeat chains with the same RUN_ID when its wall-time chunk ends.
# Only after both chains have all iterations:
bash hpc/submit_expanse.sh postprocess ECM1_250KM_VTI_001
```

`afterok` means a Slurm stage exited successfully; a time-budgeted chain can exit
successfully while still partial. Check iteration counts before postprocessing.
The 44-hour budget is checked between iterations, not inside one. If Slurm kills
an unfinished iteration, the last completed checkpoint remains valid; resumption
replays from its saved RNG state. A hard-killed job can leave a lock directory:
inspect its `job.txt` and `squeue` before manually removing a stale lock.

A full-band true-Earth forward previously took about 17 minutes locally, and some
prior draws are more expensive. Running two chains concurrently does not reduce
individual forward cost. We have not established that 15,000 iterations will
finish within a practical allocation. The benchmark stage deliberately exposes
that uncertainty without reducing mode completeness, damping accuracy or periods.

## Paths, progress and restart

All run outputs stay under `<repo>/DATA/transdimensional/<RUN_ID>` by default.
Override `CONUS_RESULTS_ROOT` with another absolute persistent project path if
needed. `logs/` contains Slurm stdout/stderr, MATLAB diaries and submission records.
Each `chain_01`/`chain_02` contains `initial_model.mat`, `progress.txt`,
`checkpoint.mat`, and eventually `chain.mat`.

Native scratch is `/scratch/$USER/job_$SLURM_JOB_ID/...`. It changes between jobs
without changing the stored scientific configuration; the saved parameter fallback
is `/tmp`. Fatal native evidence is copied into the persistent run folder before
scratch cleanup. Slurm scratch itself is temporary, not an archive.

The monitor prints every 60 seconds during a native call. The configured default
native-call timeout is 7,200 seconds. A timeout **aborts** with
`thbi:native_timeout`; it is not a statistical rejection. Inspect the preserved
model/logs. Altering model/solver settings requires a fresh run ID; do not mix
new source or binaries into an existing checkpoint. `prepared.mat` also detects
changes before starting a new array task. Only runtime scratch and the 44-hour
job budget are excluded from checkpoint identity.

Useful commands on Expanse:

```bash
squeue -u "$USER"
tail -f DATA/transdimensional/ECM1_250KM_VTI_001/logs/chains_*.out
cat DATA/transdimensional/ECM1_250KM_VTI_001/chain_01/progress.txt
sacct -j JOB_ID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

The local stopped run is preserved; these source/runtime changes require a new
prepared run, not continuation of old pre-change checkpoints. The independent
VTI body-wave benchmark and deferred forward-consistency question remain open.
