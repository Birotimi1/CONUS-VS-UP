# Damping, progress, and Expanse preparation

The transdimensional `reference100_vti` profile already uses `apsv=ash=1e-3`,
matching the older production HBI code. `THBI/b3_FORWARD_MODEL.m` passes these
values to HVf_VTI for HVSR. They are numerical complex-frequency stabilization
parameters of the P-SV/SH body-wave integral, applied through all layers,
including sediments and mantle. They are not a sediment-only material Q model.
Ellipticity uses only the fundamental Rayleigh surface wave (`nks=0`) and retains
zero body-wave damping, as in the old workflow. We did not add duplicate damping
or change the physical/sampled Earth.

The native parser maps `-apsv` and `-ash` to `PSVDAMP` and `SHDAMP`;
`native/HVf_VTI/BW_INTEGRALS_DAMPED.f90` forms the two complex frequencies, then
propagates each through every layer. The configured 1e-3 values override the
upstream 1e-5 defaults. Diagnostic startup messages now expose these settings.

## Why the old run appeared silent

The user stopped the run. Its `chain_01` folder in
`DATA/transdimensional/ECM1_19178_AK135_VTI_HK_CCP_20261111_120701` contains no
checkpoint. The previous live inspection showed native operators advancing;
absence of printed iterations did not establish a process deadlock. Its saved
parameters used quiet forwards and 100-step reporting/checkpoints. The updated
local launcher already had per-step reports, but those edits did not affect an
existing MATLAB run.

Damping was not absent. Full modal sums, strict search tolerance, numerical
sublayers and broad prior-drawn models can be expensive. No scientific tolerance,
mode count, sediment constraint or prior has been relaxed to make logging work.

## Changes

- `run_one_chain` announces initialization and forward stage, writes a persistent
  `progress.txt`, preserves the proposed initial Earth, and checkpoints iteration
  zero after its first successful likelihood. `stop_iter=0` supports an honest
  initialization benchmark. A runtime wall budget stops only between iterations,
  retaining exact RNG/state continuation; operational scratch can vary per job.
- `run_native_command` and `scripts/native_watchdog.sh` echo a heartbeat while
  the external program runs. An optional timeout/interruption aborts the job,
  never rejects a model on the basis of computational cost. The numerical
  command, its stdout data files and native stderr logs are preserved.
- Native preflight and RF utility calls also use the monitor. MATLAB's
  [system echo option](https://www.mathworks.com/help/matlab/ref/system.html)
  exposes these messages while the command is running.
- Added the [Expanse stage launcher and submission instructions](../hpc/README.md).
  Preparation runs once; separate array tasks benchmark/start each chain and
  resume completed iterations. Full scientific code blocks remain visible in
  `DATA/run_expanse_250km.m`.
- GNU/Linux and macOS SHA-256 tools are both supported in the native build guard.

Expanse was not accessed beyond one failed noninteractive authentication attempt.
The user then requested scripts locally only; no remote transfer, build or
submission was performed. Site defaults were checked against the
[SDSC user guide](https://www.sdsc.edu/systems/expanse/user_guide.html), but the
actual account, modules, CPU time and complete Linux CPS installation still
require on-cluster verification.

## Executed checks and limits (reviewed 2026-11-12T02:53:35-05:00)

The following MATLAB checks ran on November 11, before the usage interruption;
their existing outputs were reviewed and preserved on November 12. They were
not repeated during the checkpoint task.

| Check | Observed result |
|---|---|
| MATLAB syntax audit | 136 stage/helper/launcher files checked; no syntax errors. The subsequently added native-runtime test also executed successfully. |
| Analytic algorithm suite | RJ round trips, 1,080 two-stage paths, exact 36-state transition matrix, reference/VTI/RF likelihood checks and retained-model processing passed. Maximum detailed-balance error was 1.7347e-18. |
| Exact chain continuation | Continuous versus iteration-zero checkpoint, wall-budget pause, iteration-91 checkpoint and resume agreed; all 11 reference moves were exercised. |
| Native process handling | Live heartbeat, unchanged nonzero exit status, fatal timeout and descendant cleanup passed. Injected timeout left the iteration-zero checkpoint intact and did not count as a rejected model. |
| Native damping fixture | Original HVf and HVf_VTI in the isotropic limit agreed at both damping values (0.001 and 0.01), maximum relative difference 1.9915e-6. Changing damping changed H/V, maximum absolute difference 0.070591. This is a two-layer, three-frequency fixture, not an independent anisotropic body-wave benchmark. |
| Public forward wrapper | Four surface observables, native Ps synthetic waveform and Sp CCP returned finite responses from the 250-km Earth through the monitor. The fixture used three surface periods (5, 20, 80 s) and reduced mode/quadrature settings only within that test; production settings were unchanged. Ps H–kappa was covered by the separate analytic/reference tests. |
| Shell preparation | Shell syntax checks and a mocked `sbatch` check passed for all four submission stages, array/dependency arguments, paths containing spaces and invalid run IDs. No scheduler job was submitted. |
| Local native rebuild | The portable hash-tool build completed; source/executable validation passed. No scientific Fortran routine was changed in this update. |

Preserved evidence: `DATA/transdimensional/RUNTIME_EXPANSE_AUDIT_20261112_025129/`. It includes the original MATLAB logs,
audit scripts, damping report and small native outputs. These generated files
are ignored by Git and remain on the shared drive; source, tests and this
report are versioned.

The preceding November 10 full numerical gates remain documented in
[the native audit](THBI_NATIVE_AUDIT_2026-11-10.md); they were not rerun for the
process-monitor changes. Passing these checks does not prove production
convergence or practical 15,000-step runtime. Remote Linux build/module/license
checks, initial-model timings on Expanse and the full recovery experiment
remain outstanding. Independent anisotropic body-wave validation and the
user-deferred **? OPEN** forward-consistency question also remain outstanding.

**November 12 follow-up:** the independent body-wave implementation and additional
production-grid body quadrature checks now pass on the specified fixtures; see
[the body-wave audit](THBI_BODY_WAVE_VALIDATION_2026-11-12.md). The earlier
results and limitations above retain their original dates. Expanse work is
now left to the user.
