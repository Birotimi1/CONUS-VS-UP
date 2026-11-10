# Native workflow audit — ongoing, not production ready

This records the recovered November 9 audit and the November 10 continuation.
MATLAB R2022a Update4, MACI64. The user authorized native validation and small
execution checks; the production inversion remains user-run.

## Completed evidence

- MATLAB prior/material/RF helpers, 121 VTI move pairs, 1,080 layered DR paths,
  exact 36-state stationarity (6.94e-18), balance (1.73e-18), 6,000 independent
  continuous-prior transitions (largest paired drift 2.033 SE) passed.
- The 240-step analytic-forward fixture passed all eleven move families,
  exact continuation, independent starts, response caching and reconstruction.
  The November 10 rerun additionally verifies a scalar checkpoint spec with
  the entire source cell and a subsequently attached native signature.
- Native isotropic-limit and nonzero-anisotropy comparisons passed. Independent
  CPS/HVf fundamental phase/group/ellipticity relative differences were
  7.57e-6, 8.20e-6 and 5.54e-6. Xi=.8 and1.2 plus a mantle LVZ also passed;
  largest modal discrepancy 9.18e-6. This covers fundamental modes only.
- Full configured synthetic curves/RF products were generated in52.6s under
  the earlier unconverged 10-mode configuration. These are diagnostic outputs,
  not validated production data. Full-grid split/merge/deep-Vs/kappa/xi
  sensitivity checks completed; round-trip differences were zero. RF controls
  passed homogeneous, bilayer, deep conversion, duration and water-level tests.
- Corrected CPS `tdpegn96` calls to include required `-C`/`-U` selectors;
  corrected the native parity test's merged precision-config lookup; added
  pointwise isotropic-limit checks so a large resonance cannot mask errors.
- Fixed `run_one_chain`'s cell-valued struct constructor: unwrapped
  `source_snapshot()` created an Nx2 struct array. Native signatures now attach
  to a scalar spec. The native chain1 completed12steps, including stop/resume
  at5, after this fix. Its second chain stopped at a native Love-search error;
  final c1–c5/plot checks have not completed in that execution audit.
- Fixed the Fortran cancellation-mask offset, made failed dispersion return
  nonzero before H/V output, and reject legacy exit-zero error reports in MATLAB.
  The observed Love failure now propagates and preserves the failed model/files.
- Added mandatory HVSR convergence gates to both active synthetic launchers.
  Do not remove/loosen them to force a production run.

## Unresolved numerical issue

At20Hz, the same true Earth gives HVSR28.0807 (10modes/1000nodes),26.1876
(40/4000),15.2288 (80/4000),14.3433 (160/4000). The320-mode result equals160
because the old search reports only134modes; CPS independently finds many
more. Agreement between160/320 therefore does not establish convergence.
Tighter root precision can also cancel a search while the old executable
returns a finite partial curve and exit0. This is now an explicit failure.

An isolated double-precision build alone still misses roots. An isolated
finer-search build matches CPS's first320 Rayleigh modes at19/20/25Hz to
about1.5e-5 relative velocity, but terminates at higher modes; Love root
comparisons remain under investigation. These diagnostic builds live in
ignored experiment directories and are not the production executable.

The two12-step smoke chains use representative curve coordinates and full
RF products. They test execution only and cannot demonstrate recovery,
mixing, uncertainty calibration or resolution through250km. The current
native failure must be resolved before those checks can complete.

**? OPEN:** dispersion-only extension/shared-Earth consistency comparison
remains deferred by the user.

## Preserved output paths

All following directories are below `DATA/transdimensional/` (ignored by Git):

- `NATIVE_AUDIT_20261109_182057` — original missing CPS-selector failure.
- `NATIVE_AUDIT_20261109_182215` — corrected native fundamental comparison.
- `VTI_ENVELOPE_20261109_182626` — xi bounds and mantle LVZ.
- `WORKFLOW_VALIDATION_20261109_182418` — full grids, sensitivity and RF
  controls; original scalar-checkpoint failure.
- `HVSR_CONVERGENCE_20261109_182743` —10/20/40/80/160/320-mode diagnostics.
- `MODAL_AUDIT_20261110_043915` — independent CPS modes, double precision,
  fine-search diagnostics; inspect each report before using any finite output.
- `NATIVE_SMOKE_20261110_044432` — corrected checkpoint, native chain1 and
  preserved chain2 Love-search failures. Incomplete execution audit.

Local logs include `/private/tmp/conus_full_workflow_audit.log`,
`conus_checkpoint_fix_test.log`, `conus_native_smoke_resume.log` and the
November9 core/native/envelope/convergence logs. Shared-drive outputs are
not part of the Git backup. No production inversion has been launched.

## Continuation status — 2026-11-10T10:26:36-05:00

The earlier missing-mode issue is resolved on the truth and preserved Love-cutoff/LVZ fixtures: complete modal counts agree with independent CPS, and search/root/quadrature refinement passed. See the daily entry above for measured errors and output paths. CPS branch artefacts and period-grid dependence required additional wrapper checks; isolated-period retry validation and the final native two-chain audit are still in progress. All three within-zone monotonicity flags are now false at the user's request. Earlier results and failures above remain historical evidence.
