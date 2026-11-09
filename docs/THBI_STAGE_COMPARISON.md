# Main-body stage comparison with the local THBI checkout

Reference read: `/Users/birotimi/Downloads/THBI_MCMC`.
Working code: `CONUS-VS-UP/HBI` on the shared filer.

The important distinction is each file's role and the work visible in its
main body, not its line count. The earlier version put much of b3 below a
generic forwarding call and assigned different roles to a1/a2. These stages
now expose the equivalent fixed-layer work directly.

| THBI stage | What is in its main body | Our corresponding main-body implementation |
|---|---|---|
| `a0_STARTUP_BAYES` | Project, MATLAB dependencies, native solver paths, execution/data directories | `a0_STARTUP_BAYES`: project and MATLAB paths, CPS/HVf discovery, data/results/scratch directories; optional `paths` output |
| `a1_TEST_CONDITIONS` | Explicit checks on allowed velocities, interfaces, gradients and spline parameters | `a1_TEST_CONDITIONS`: existing fixed-layer parameter counts, finite values, thickness/Vs bounds, monotone Vs and log-noise support; no unrelated spline constraints added |
| `a2_LOAD_DATA` | Synthetic/observed-data branches, generation/loading, processing and data saving | `a2_LOAD_DATA`: prescribed thin h/Vs and material equations, five forward datasets, RF time selection, errors/covariance, noise realization and saving; this initial workflow remains synthetic-only |
| `b1_INITIATE_MODEL` | Initial sediment/crust/mantle parameters, spline construction, noise parameters, assembled model | `b1_INITIATE_MODEL`: reference theta, independently drawn Gaussian perturbations, conditional reflection, Vp/density equations, physical matrix and five noise multipliers; no b2 call is used to hide initialization |
| `b3_FORWARD_MODEL_BW` / `b3_FORWARD_MODEL_SW_precise` | Data-type branches, wave/model preparation, native solver calls, waveform processing and prediction assignment | `b3_FORWARD_MODEL`: shared model preparation, HVSR/ellipticity mode settings and HVf command/output code, CPS phase/group calculations, native RF command/header checks and arrival references, then fixed-grid prediction assignment |

The main b3 body now contains the calculations previously inside
`forward_five_observables`, `hv_prediction`, `dispersion_prediction`, and
`forward_receiver_function_cps`. Those forwarding routines are no longer
part of the HBI b3 file. Small path, layer-validation/refinement, model-file
and SAC utilities are in the repository-level `functions/` folder as
inspectable implementations, shared with the zoned THBI stages.

THBI's MINEOS, PropMat/telewavesim, spline and transdimensional physics are not
substituted for our validated fixed-layer HVf/CPS configuration. The current
run continues to use the same h/Vs values, priors, observation grids, solver
flags, noise model and hierarchical likelihood. Actual wave-propagation
physics remains inside the native executable programs.

`DATA/run_ecm1_thin_hbi.m` visibly executes the sequence: paths, parameters,
run directory, a2 data preparation, independent chains, c1 summaries, c2
plots, and final saves. To edit the inversion's thin synthetic truth, edit
its h/Vs block in `a2_LOAD_DATA.m`. Existing standalone clean-forward scripts
use `DATA/generate_ecm1_synthetics_five_observables.m` and helpers in `functions/`.

Verification for this reorganization was static: the moved HVf, CPS
dispersion and RF blocks were compared against the pre-edit source, allowing
only the needed variable/scope changes. MATLAB, Octave, native solvers and
inversion were not run; the user will launch the MATLAB script.
