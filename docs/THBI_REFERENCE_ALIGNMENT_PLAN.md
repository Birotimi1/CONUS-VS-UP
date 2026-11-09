# THBI reference alignment — implemented; native validation pending

User requested on 2026-11-09: resolve the comparative audit differences while
retaining layered Earth, CPS/HVf and corrected random transdimensional DR.
Explicit RF choice: Ps H–κ and Sp CCP. User deferred the dispersion-only
extension/shared-Earth question; mark it **? OPEN**, do not change it.

- [x] Verify reference permits crustal/mantle LVZs.
- [x] Locate old synthetic provenance: manually assigned 70-km Earth, not an ECM1 extraction.
- [x] Download official ECM1 crust and sediment archives; shared ECM1 directory is empty.
- [x] Extract a genuine thin-sediment continental column with source rows/hashes/units.
- [x] User selected AK135 reference mantle; explicitly labelled and extracted.
- [x] Adapt reference Ps H–κ stack processing/energy residual and Sp P–SV, migration, parent-pulse, taper, depth weighting and cross-convolution residual to CPS.
- [x] Add sampled crustal Vp/Vs, reference material relations and layered near-Moho condition; keep LVZ support.
- [x] Broaden sediment/Moho bounds appropriately and use dispersed starts.
- [x] Add reversible coupled discontinuity proposals with audited densities.
- [x] Add likelihood selection for reference RF residuals/absolute sigma and degrees of freedom without changing the tested RJ/DR correction.
- [x] Update state storage/reconstruction, restart schema, summaries, RF-specific plots, and launchers.
- [x] Add modern convergence diagnostics and native preflight/sensitivity scripts for user-run MATLAB.
- [x] Run solver-free tests and MATLAB static analysis; do not launch production/native inversions.
- [x] Update documentation/provenance/daily log.

Checkpoint/push outcomes are recorded in Git history and the daily log;
publication is verified after each local commit, not asserted in advance.

Do not claim exact anisotropic THBI equivalence with isotropic CPS/HVf.
Do not restore legacy uncorrected retries, kernel-reset machinery or posterior
culling as if they were requirements of the layered target. Keep old profiles
readable/reproducible where practical; the updated launcher selects the new profile.

## Scope update received during implementation

User explicitly selected **extend anisotropy before testing** on2026-11-09.
The new isotropic reference100 implementation is an intermediate checkpoint,
NOT the final requested test configuration. Do not launch its native inversion.
The current native launchers select reference100_vti and require a native
comparison gate before inversion. The extension is implemented and compiled;
scientific native validation is still user-run. See THBI_WORKFLOW.md.

- [x] Package/extend the available HVf Fortran source for VTI (full PSV and SH).
- [x] Use CPS tprep96/tdisp96/tregn96 and trftn96 with the same elastic tensor.
- [x] Adapt HK kinematics and Sp migration to VTI slownesses; preserve phase signs.
- [x] Sample/reference anisotropy consistently with THBI, including linked P
      anisotropy and stability tests; update RJ storage, priors and proposals.
- [x] Supply user-run native isotropic-limit and anisotropic parity tests.

CPS executables exist in the shared AGU_models/bin_v3.30 directory. HVf source
exists under AGU_models/SURF96_Rayl_1Dinv/src/HVf/HV-DFA_source_FORTRAN_V15.
THBI example uses xi in[0.8,1.2] in crust, mantle xi fixed1, phi=1-(xi-1), eta=1.
Its optional Sp DOF43 override is specific to Hopper CCP data; the synthetic
loader now uses the underlying Silver-Chan spectral estimate instead.
