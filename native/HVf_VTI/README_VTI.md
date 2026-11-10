# HVf VTI source extension

This folder contains a modified copy of the project's damped HV-DFA V15
Fortran source. Original authors: Antonio García-Jerez, José Piña-Flores and
the HV-INV team; the modal search credits Marc Wathelet/Geopsy. The copied
code is GPL-3.0; retain `LICENSE.txt`, author headers and `UPSTREAM.json`.
Original source hashes identify the exact local source used. Upstream:
<https://github.com/agarcia-jerez/HV-DFA>.

**Status (November10):** complete modal sums now pass search/root/quadrature
refinement and independent CPS surface-wave root/Green-amplitude checks on the
truth and preserved LVZ/Love-cutoff fixture. The final native workflow audit
is still in progress. See `docs/THBI_NATIVE_AUDIT_2026-11-10.md`. Native
failures return nonzero before H/V output; a finite partial sum is not valid.

## Build and model convention

From the repository root:

```sh
sh native/HVf_VTI/build.sh
```

Requires gfortran and shasum. The build uses double-precision default REAL. Builds `bin/HVf_VTI` and records source/binary
SHA256 manifests. Generated binaries/build files are ignored by Git.
`check_native_vti_build` detects stale sources or replaced executables and
records actual CPS executable hashes. Chain restarts also compare these.

Use `-vti` **before** `-f model.txt`. The model starts with its row count,
then rows `h Vp_average Vs_average rho xi phi eta`. The first four columns
are SI units (m, m/s, m/s, kg/m³); the last three are dimensionless. Final
thickness is zero. Without `-vti`, four-column isotropic input remains valid.
The frequency/mode/damping flags retain the upstream interface. The added
`-dphase` bounds accumulated vertical phase change during modal search
(default pi/32); halving it provides a search-resolution check. `-complete`
requires the search to terminate below both mode capacities or exits with an
error. `modal_counts.dat` records actual Rayleigh/Love counts by frequency.
The active profile uses4096as capacity,4000body nodes,1e-8percent root
precision and a pi/512phase-step bound. Identical numerical sublayers share
one phase-bound evaluation; this preserves the same physical integral. These are not4096physical layers or an assertion that4096modes
exist at every frequency. Physical absence of a Love fundamental is allowed.
Previous10-mode synthetics were not converged and must be regenerated.

The active inversion mirrors THBI's example: isotropic sediment and mantle,
one crustal xi in [0.8,1.2], phi=2−xi, eta=1. The synthetic crust uses an
**assigned** xi=1.06. This is not anisotropy measured by ECM1 or AK135.
The tensor routines accept more general stable VTI tensors, but that broader
parameter space is not established by this initial test configuration.

## Operators and derivation

The THBI parameter averages satisfy
`Vs²=(2 Vsv²+Vsh²)/3`, `Vp²=(Vpv²+4 Vph²)/5`.
These are the reference's velocity convention, not a general TI Voigt
isotropic projection. With `xi=(Vsh/Vsv)²`, `phi=(Vpv/Vph)²`:

```
Vsv² = 3 Vs²/(xi+2)      Vsh² = xi Vsv²
Vph² = 5 Vp²/(phi+4)     Vpv² = phi Vph²
A = rho Vph²             C = rho Vpv²
L = rho Vsv²             N = rho Vsh²
F = eta (A−2L)
```

Require L>0, N>0, C>0, A>N and (A−N)C>F². MATLAB and Fortran both check
stability. The matrix's Vs/Vp columns remain the THBI parameter averages;
they are never substituted for directional velocities in a VTI propagator.

For horizontal wavenumber k and depth factor exp(−gamma*z), the P/SV roots
solve

```
(L gamma² − A k² + rho omega²)
  (C gamma² − L k² + rho omega²) + (F+L)² k² gamma² = 0.
```

`vti_elastic.f90` constructs both P/SV eigenvectors and their inverse in
the upstream transformed state `[u,v,L(u'−k v),F k u+C v']`. Outgoing roots
retain decay into the half-space. SH uses
`gamma_SH²=(N k²−rho omega²)/L`. Upstream Wang normalization and analytic
depth integrations remain in place.

Rayleigh group velocity uses the same stationary-energy identity as HV-DFA.
Let I0 be the density-weighted displacement norm,
I1=integral(A|u|²+L|v|²), and ID2=integral(L|u'|²+C|v'|²).
Write the elastic energy as k²I1+k*Icross+ID2=omega²I0. Differentiating
the stationary eigenvalue and eliminating Icross gives

`vg = (k² I1 + omega² I0 − ID2)/(2 omega k I0)`.

Thus F enters the eigenfunctions/roots even though it disappears from this
particular expression for group velocity. Love's horizontal energy integral
uses N; its vertical traction uses L. This distinction is essential.

Body-wave integrations use separate PSV and SH meshes ending at the
half-space shear branch points omega/Vsv and omega/Vsh. The original damping
and rectangle endpoint convention (including k=0.1*dk at the near-zero
endpoint) are retained. No integration beyond the radiation branch is added.
The active profile has an isotropic mantle/half-space. Numerical refinement passed on the documented fixtures; this is not an
independent validation of the anisotropic body-wave physics.

## Modified files

| File | Change |
|---|---|
| `vti_elastic.f90` | New tensor conversion, Christoffel roots, eigenvectors, inverse, homogeneous Rayleigh bound |
| `HV.f90` | Initialize tensor; branch-specific modal bounds |
| `read_command_line.f03` | Seven-column VTI input; phase-step and complete-sum flags |
| `Dispersion.f90` / `RootSolverTemplate.f90` | Adaptive phase-bounded search, refined root tolerance, cancellation-mask repair and physical Love cutoffs |
| `aux_procedures.f90` | Homogeneous bound uses VTI secular equation; original isotropic matrix routines retained for comparison |
| `y.f90` | VTI Rayleigh/SH propagation; normalized determinant phase prevents rank loss in deep evanescent layers |
| `WangMethod_Rayleigh_new.f90` | VTI eigenfunctions; A/C/L energy integrals |
| `WangMethod_Love.f90` | SH root and N-weighted horizontal energy |
| `BW_INTEGRALS_DAMPED.f90` | VTI P/SV and SH body-wave Green functions |

The other copied source files retain their algorithms apart from encoding
and line-ending normalization; `build.sh` additionally enforces double precision. No original shared-drive source was edited.

CPS uses the same constants via model96 directional velocities and
`VPF=sqrt(F/rho)`. Primary references: [CPS VTI tests](https://www.eas.slu.edu/eqc/eqc_cps/VTI/index.html)
and [CPS surface-wave program sequence](https://rbherrmann.github.io/ComputerProgramsSeismology/TUTORIAL/RDOC/index.html).
