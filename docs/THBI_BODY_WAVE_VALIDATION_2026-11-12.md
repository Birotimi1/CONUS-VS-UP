# Independent anisotropic body-wave validation — November 12, 2026

The independent implementation comparison passes for eleven Earth models.
It checks each native **horizontal P–SV, horizontal SH and vertical P–SV**
component separately, so agreement cannot result from cancellation in H/V.
The maximum relative native/reference difference is **1.3844e-8**, and the
250-km ECM1/AK135 truth differs by at most **8.6557e-10**.

Run the local MATLAB script `DATA/validate_thbi_body_waves.m`. No inversion or
Expanse operation is involved. The active local forward-only, inversion and
whole-workflow validation launchers now require the independent body gate.
The existing Expanse scripts were left untouched at the user's request.

## Independent construction and its own checks

`functions/reference_vti_compliance.m` derives the forced response from
physical displacements and tractions, rather than importing HVf's analytic
polarization vectors, transformed vertical coordinate or Wang normalization.
For `exp(i*k*x - i*w*t)` with z positive downward, its state is
`y = [ux, uz, txz, tzz]`. Elastic constitutive relations and momentum balance give:

```text
txz = L (ux' + i k uz)       tzz = F i k ux + C uz'
y' = M y

M = [ 0                       -i k         1/L       0
     -i k F/C                  0           0         1/C
     (A-F²/C)k²-rho w²         0           0        -i k F/C
      0                       -rho w²     -i k       0 ]

SH: [uy, tyz]' = [ 0                  1/L ] [uy, tyz]
                 [ N k²-rho w²        0  ]
```

The equations above are derived directly for this benchmark. A generic MATLAB
`eig` calculation selects the decaying half-space subspace at positive
imaginary frequency. Impedance is traction/displacement: `Z = Td / Ud`.
For a finite layer above bottom impedance Zb:

```text
Rbottom = (Tu - Zb Uu)^(-1) (Zb Ud - Td)
Rtop    = exp(-lambda_up h) Rbottom exp(lambda_down h)
Ztop    = (Td + Tu Rtop) / (Ud + Uu Rtop)
G       = Ztop^(-1)
```

Both exponentials decay; large growing transfer matrices are avoided. The
implementation scales the ODE before diagonalization and restores SI
compliance units afterward. Its traction normal points downward, so its
spectral diagonal has the same negative-imaginary convention as `GBW.dat`.
Azimuthal integration gives factors `1/(4*pi)` for each horizontal component
and `1/(2*pi)` for the vertical component in the integral of `k*Im(G)`.

The reference is checked independently before comparison with native output:

- Isotropic P/SV potential solution: relative norm error 1.61e-15.
- Isotropic and VTI SH half-space closed forms: errors below 5e-16.
- VTI normal-incidence P/SV impedances: error 3.30e-16.
- A separate direct `expm` layer recurrence on a thin anisotropic stack:
  errors below 1.6e-15. The dissipative compliance has the expected passive sign.
- Closed-form integrated VTI SH half-space: the documented rectangle mesh
  differs by 6.7543e-4 at 12,000 intervals. This tests normalization and endpoint
  quadrature independently, rather than comparing two copies of one sum.

The separation of body/surface Green contributions follows the original
[HV-DFA formulation](https://arxiv.org/abs/1604.06406). The independent physical
ODE and bounded impedance calculation above were implemented for this audit;
they are not an external package being presented as a third-party certification.
The existing independent surface/RF comparisons use the documented
[CPS TI programs](https://www.eas.slu.edu/eqc/eqc_cps/CPS/CPS330/cps330o.pdf).

## Native comparison

| Earth | Maximum component relative difference |
|---|---:|
| isotropic_halfspace | 1.868e-15 |
| vti_halfspace | 3.797e-15 |
| general_vti_layers | 1.384e-08 |
| ecm1_ak135_truth | 8.656e-10 |
| crust_xi_min | 8.887e-10 |
| crust_xi_max | 8.434e-10 |
| mantle_lvz | 8.689e-10 |
| prior_draw_01 | 3.357e-10 |
| prior_draw_02 | 8.557e-10 |
| prior_draw_03 | 1.068e-08 |
| prior_draw_04 | 2.805e-09 |

The truth uses the production 4,000-interval body mesh at 0.1, 0.2, 0.5, 1, 5,
18 and 25 Hz. Other comparison fixtures use 512 intervals at 0.1, 1 and 25 Hz
only to isolate propagation/normalization agreement on identical meshes; this
does not change production settings or establish their quadrature convergence.
The general VTI cases include anisotropic half-spaces and independent phi/eta.
The actual inversion continues to sample linked **crustal** anisotropy only;
sediment and mantle remain isotropic, and sediment Vs remains nondecreasing.

The native inputs use the same maximum-2-km numerical subdivision as b3;
the independent reference retains the unsplit physical layers. Their agreement
therefore also checks that numerical subdivision preserves the physical response.
An initial direct native call with thick unsplit 250-km layers produced nonfinite
PSV terms at 18/25 Hz and failed loudly. That invocation bypassed the production
wrapper. The corrected benchmark uses the production input contract; this
limitation remains documented and the numerical subdivision must be retained.
No scientific Fortran, damping, prior or acceptance equation was altered.
A source-recording filename typo in the new launcher test was also corrected;
that failed invocation stopped before numerical testing. Both failure records
remain preserved, rather than being reported as passes.

## Additional production-grid accuracy checks

Fresh body-only native calculations at 4,000 and 8,000 integration intervals
passed for the truth, both crustal xi bounds, the mantle LVZ and four valid
prior draws. The maximum component change was **0.00101036 (0.101%)**, below the
2% check threshold. These checks isolate quadrature from the earlier
same-grid physics comparison; neither mode truncation nor H/V cancellation
can make this test pass. They are implemented in
`THBI/tests/test_thbi_body_quadrature.m` and required by the local launchers.

## Existing workflow evidence reviewed

The saved November 10 full-workflow report records `passed=true` and
`execution_passed=true`. Its surface-mode comparisons, Ps/Sp controls, two
six-step native chains, exact continuation, c1–c5 and fresh final likelihood
checks remain the evidence for those unchanged scientific stages. The November
11 analytic/runtime audit additionally passed the new iteration-zero checkpoint,
wall-budget pause, native timeout behavior and live solver wrapper checks.
Those successful scientific tests were not rerun just to produce a new date.

All Fortran scientific-source hashes and all eight CPS executable hashes match
the preserved modal comparison evidence. Reading the saved 4,000/8,000-point
truth body contributions separately gives a maximum component change of
2.8523e-4 (0.02853%) over the active band. The earlier complete-H/V refinement
check was therefore not concealing a large body-component cancellation.

Evidence is under `DATA/transdimensional/BODY_REFERENCE_20261112_031310/`:
`body_reference_report.mat`, `validation.log`, `validation_inventory.json`,
`static_checks_final.log` (142 MATLAB files), `quadrature/body_quadrature_report.mat`, and the per-model input/native-output directories. The MAT
report includes the exact reference/test source text. Earlier failed and
seven-case intermediate tests remain in the neighboring timestamped folders.
Generated evidence is preserved on the shared drive and excluded from Git.

The completed implementation/forward checks do **not** establish posterior
mixing, recovery of the shallow interfaces, calibrated uncertainty, resolution
at 250 km, or practical full-production runtime. Those require the user-run
recovery experiment and its diagnostics. The explicitly user-deferred
**? OPEN — dispersion-only extension/shared-Earth reference comparison** remains
open. Expanse installation/submission is outside this validation task.
