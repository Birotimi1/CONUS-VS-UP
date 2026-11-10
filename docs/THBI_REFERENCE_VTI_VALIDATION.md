# Reference/VTI implementation validation

This is the historical implementation-only record. Native execution began
later on November9; current results and unresolved HVSR failures are in
[the November10 audit](THBI_NATIVE_AUDIT_2026-11-10.md).

Recorded 2026-11-09T13:32:49-05:00. Branch `codex/thin-five-observable-synthetic`.
These are implementation checks, not native physics or recovery results.

## Executed

- Octave `test_thbi_reference100()` after HK interpolation/Sp DOF changes:
  passed 100 full delayed-rejection paths, true-model/material checks, LVZ
  support, HK analytic phases, Sp sign/migration, cross-convolution, RF
  likelihood and rank/ESS controls.
- `test_thbi_vti()`: passed THBI-average tensor identities, independent
  Christoffel eigenvalue residuals, isotropic root limit, VTI HK arrivals,
  exact layer/refinement migration, CPS model96 read-back and scale-invariant
  observed RF DOF. Final rerun also checked legacy isotropic entry points
  reject VTI before native calls or file writes.
- `test_thbi_reference100('reference100_vti')`: passed all **121** move-family
  pairs, physical inverse paths, proposal densities and full DR path flux.
- `test_thbi_reference_chain('reference100_vti')`: passed the 240-step injected
  analytic-forward fixture, all eleven families attempted, stage-two accepts,
  per-state fresh likelihood/cache checks, xi/kappa storage and reconstruction,
  independent chains and exact continuation after iteration91. Chain1 ended
  K[1,2,5], logL−58.87, acceptance0.66, 192 forward-fixture calls; chain2
  K[3,1,3], logL−57.68, acceptance0.70, 197 calls. Noise-only fixture used one
  forward call. These numbers describe a **toy fixture**, not Earth recovery.
- Regressions: `test_thbi_rj`, `test_thbi_extended100`, `test_thbi_layered_dr`
  (1,080 paths), `test_thbi_dr_stationarity('layered_path')`, and
  `test_thbi_chain('extended100','layered_path',false)` passed.
- Exact discrete random-DR fixture: 36 states, 648 second-stage paths;
  stationarity error1.3878e−17 and detailed-balance error1.7347e−18. The naive
  retry negative control remains nonstationary (8.2046e−4).
- GNU Fortran16.1.0 compiled `native/HVf_VTI/HV.f90` successfully. The compiled
  source/executable manifests and hashes of seven installed CPS programs
  were verified with `check_native_vti_build`. No scientific native executable
  was run by these build/hash checks.
- MATLABR2022a standalone `mlint`: 112 MATLAB files, exit0, no syntax/unreadable
  file errors. A subsequent incremental pass checked the last six edits,
  exit0. Remaining messages concern unused suppressions, dlmread deprecation,
  a legacy unreachable debug-plot block and minor style. Logs:
  `/private/tmp/conus_vti_mlint_final.log`,
  `/private/tmp/conus_vti_mlint_incremental.log`.

The first static invocation with relative paths could not read its inputs;
that invocation was not counted as a syntax check. Absolute-path passes above
succeeded. Octave skipped MATLAB-only histcounts-dependent summary/graphics
branches; those remain for the user's MATLAB run.

## Not executed

- Native HVf/CPS isotropic-limit and nonzero-anisotropy parity.
- Native true-model curves, waveform/RF-product checks or deep sensitivity.
- Full MATLAB figures or a physical two-chain inversion.

User-run `test_thbi_vti_native` compares all six response blocks in the
isotropic limit and requires anisotropy response. It then independently
compares anisotropic HVf fundamental phase/group/ellipticity with CPS.
`test_thbi_reference_native` checks native split/merge reversibility and
responses to deep Vs, kappa and xi perturbations. The DATA launchers run
these gates before inversion and preserve logs/outputs on failure.

No new scientific inversion-output directory was generated in this session.
Temporary toy fixtures cleaned themselves up. Native build products remain
ignored in `native/HVf_VTI/bin/` and `build/`. The official CPS source archive
used for schema inspection is `/private/tmp/conus_CPS330_source.tgz` (not Git).

**? OPEN:** dispersion-only extension/shared-Earth comparison is still
explicitly deferred. The new physics implementation needs the user's native
validation; no claim of shallow recovery or resolution through250km is made.
