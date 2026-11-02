# Thin-model forward validation

Executed on 2026-11-02 using GNU Octave 11.3.0 and the existing macOS ARM64
HVf/CPS binaries on the shared filer. Numerical execution completed for all
five observables; the Python figure was rendered from those saved outputs.
The MATLAB execution path has not been separately run in this environment.

Result directory (relative to the repository):
`DATA/synthetics/ECM1_thin_20261102_085356/`.
Generated results are retained locally and ignored by Git. Reproduce them
with the commands in the main README.

| Check | Result |
|---|---|
| Original thin h and Vs values | Preserved exactly |
| Physical geometry | 9 finite layers to 70 km, then half-space; Moho at 40 km |
| Numerical subdivision | 40 rows; material properties preserved on both sides of every physical interface |
| Actual native solver input files | Identical refined Earth within text serialization precision for all five solvers |
| HVSR | 80/80 finite samples |
| Ellipticity | 50/50 finite samples |
| Phase velocity | 35/35 finite samples |
| Group velocity | 40/40 finite samples |
| Receiver function | 4096/4096 finite samples |
| CSV versus MAT exports | Agree within text serialization precision |
| Original thin-run comparison | All five clean arrays match exactly (maximum absolute difference zero) |
| Noise | None; obs=clean and sig=0 |

The reference comparison used the existing HBI result at
`../joint_forward_test/ECM1_thin/ECM1_20261018_230849/ecm1_thin/`.
The new run used the imported forward calculators, not the saved reference
arrays. The original run was read only after the new calculations completed.

The independent receiver-function checks passed:

- A homogeneous half-space has negligible late conversion energy.
- A 35 km uniform crust produces a positive Moho Ps peak at **4.350 s**, versus
  **4.349 s** from the independently calculated ray-theory travel-time lag.
- Splitting the same homogeneous layer preserves its receiver function.
- Numerical layer refinement preserves material properties and the waveform.
- An invalid ray parameter of 6 s/km is rejected.

For the thin multilayer model, the primary Moho Ps reference is **4.75746 s**.
That reference is not an interpretation of every waveform peak; the full
receiver function includes other conversions and reverberations.

`solver_provenance.json` in the run directory records the resolved binary
paths and SHA256 hashes. Native input files and logs are retained under the
five observable subdirectories. This is a forward-consistency test, not a
demonstration of uniqueness, data fitting, or inversion recovery.
