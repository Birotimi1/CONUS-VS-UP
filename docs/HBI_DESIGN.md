# Fixed-layer five-observable HBI

This is the first hierarchical inversion stage of CONUS-VS-UP. It fits
diffuse-field HVSR, fundamental Rayleigh ellipticity, Rayleigh phase velocity,
Rayleigh group velocity, and an incident-P receiver function with one flat,
isotropic layered Earth. It does not change the number of physical layers.

## Sources and coding conventions

The implementation uses numbered MATLAB functions and `par`, `model`,
`data`, `predata`, `misfit`, and `datahparm` structures, following Zach Eilon
and Brennan Brunsvik's THBI_MCMC organization and the existing HBI workflow.
It has no dependency on the original MATLAB source folders at runtime.
Native HVf and CPS programs remain external, with configurable paths.

Reference THBI checkout: `/Users/birotimi/Downloads/THBI_MCMC`, origin
`https://github.com/brennanbrunsvik/THBI_MCMC.git`, HEAD
`585ee81` (the RF P-incidence-angle correction). The checkout has pre-existing
local changes; it was read without modification. THBI is GPLv3; this
repository retains its GPLv3 license. The b1/b2/b4/b5/b6/b9 separation and
data-type-specific hierarchical formulation are adapted here; the spline,
anisotropy, transdimensional proposals, MINEOS kernels, and PropMat machinery
are not imported into a fixed-layer CPS calculation.

HBI reference: `/Users/birotimi/Downloads/hbi_all_codes/Codes/Mod5_HBI/functions_hbi_new`.
The physical bounds and proposal steps come from `c1_define_hbi_prior.m`;
the reflected coordinate proposal from `c5_perturb_theta.m`; Gaussian error
scaling from `c3_compute_likelihood_damp.m`; and the lognormal noise prior
from `Functions/compute_log_prior_sigma.m`. The reflection helper is copied
and simplified. The original Vp/density conversion and order of applying
floors are preserved. See [forward provenance](PROVENANCE.md) for the copied
forward helpers and their original hashes.

## Function map

| File in HBI | Purpose |
|---|---|
| `a0_STARTUP_BAYES` | Explicit project, MATLAB, solver, input/output, and scratch paths; returns `paths` |
| `a1_DEFINE_PARAMETERS` | Full priors, solver settings, grids, noise, run lengths |
| `a1_TEST_CONDITIONS` | Parameter counts, finite values, thickness/Vs bounds, monotonicity, noise bounds |
| `a2_LOAD_DATA` | Thin true model, five synthetic datasets, errors/covariance, seeded noise, saved solver inputs |
| `a2_PREPARE_DATA` | Standalone covariance preparation for independently supplied data and tests |
| `a3_CALC_LOG_PRIOR` | Apply a1's support constraints and evaluate the Gaussian log-noise prior |
| `b1_INITIATE_MODEL` | Full Vp/density equations, model assembly, and independent starting-model jitter |
| `b2_PERTURB_MODEL` | Propose one thickness, Vs, or log noise multiplier |
| `b3_FORWARD_MODEL` | Full MATLAB solver interfaces, layer refinement, model writers, native commands, SAC/dispersion readers, and grid checks, all in this file |
| `b4_CALC_MISFIT` | Whiten residuals using fixed base covariance |
| `b5_CALC_LIKELIHOOD` | Normalized Gaussian likelihood for all five datasets |
| `b6_IFACCEPT` | Log-space Metropolis-Hastings decision |
| `b9_SAVE_RESULT` | Record each state, including rejected/repeated states |
| `run_one_chain` | Sampling, burn-in tuning, caching, checkpoint/resume |
| `run_all_chains` | Serial or MATLAB process-worker chains |
| `run_synthetic_hbi` / `resume_hbi` | New run / exact continuation |
| `c1_PROCESS_RESULTS` | Retained-sample summaries and classical split R-hat |
| `c2_PLOT_RESULTS` | MATLAB fit, shallow-profile, and chain-diagnostic figures |

### Editing and debugging the numbered stages

The numbered files contain their actual implementations, with `%%` sections
for navigation in the MATLAB Editor. They remain functions with explicit
inputs/outputs so each chain has its own workspace. In a0, a1, a2, b1 and b3,
the main bodies expose the stage's work directly. b1 contains the actual
initialization/reflection loop. b3 contains the HVf/CPS command construction,
execution, output checks and RF processing directly in its main body; it no
longer calls a generic `forward_five_observables` function. Low-level file,
path and layer utilities, quantiles/R-hat, and profile/report writers are in
the repository-level `functions/` folder, shared with THBI. c1 and c2 retain
their main summary and plotting blocks. Set Editor breakpoints in the stages
or helpers, or use `dbstop if error` before running the separate launcher.

Edit settings in `a1_DEFINE_PARAMETERS.m` and material-property equations in
`b1_INITIATE_MODEL.m`. Calling b1 with only `par` creates a jittered chain
start; supplying `theta, eta` builds that exact model without initialization
jitter. b2 uses this second form, keeping the material equations in one
place. `HBI/parms/bayes_inv_parms` and the repository helpers
`functions/make_layer_model` and `functions/calc_log_prior` remain compatibility
entry points; the HBI workflow directly uses the numbered implementations.
The latter two check that the HBI stages are active and never alter the path.

The true h/Vs arrays for the HBI synthetic are directly in `a2_LOAD_DATA`.
The earlier `z1_SYNTH_DATA` name is retained only for existing callers/tests.
Both `a2_LOAD_DATA` and the chain's structural proposals use the complete
forward implementation in `b3_FORWARD_MODEL.m`. Editing its HVSR, ellipticity,
dispersion, or RF sections therefore changes both generation and inversion.
HVf, surf96, and hrftn96 are still external compiled physics solvers; their
MATLAB model writing, command construction, parsing, and validation code is
now visible within b3. The earlier standalone clean-forward demonstration
is retained in `functions/forward_five_observables.m` and is not called by HBI's b3.
The separate DATA run script now shows output setup, a2, the chains, c1, c2,
and saves as sequential blocks rather than delegating the entire workflow to
`run_synthetic_hbi`. The latter remains available for programmatic callers.
See [direct comparison with THBI](THBI_STAGE_COMPARISON.md).

## Model and priors

There are nine finite-layer thicknesses and ten Vs values, including a
half-space. Numerical subdivision to at most 2 km introduces no additional
parameters. Vp and density are deterministic empirical functions of Vs;
they are not independently inferred in this test. Monotonic Vs and uniform
physical bounds match the existing HBI assumptions. Interface 7 is labelled
the Moho for this synthetic configuration; its depth is the sum of the first
seven thicknesses, not an additional parameter. The lower boundary depth is
also allowed to vary.

The specified thin truth has 10 m and 40 m soft layers, deeper shallow
interfaces at 130 m and 1 km, a 40 km Moho, and a half-space below 70 km.
The third-layer truth thickness (0.08 km) lies on the inherited prior's lower
bound. Posterior asymmetry there is expected; the model does not have a
symmetrical recovery interval around that true value.

Five independent hierarchical multipliers gamma_d scale the stated base
errors. We sample eta_d = ln(gamma_d), with a truncated normal prior of mean
zero and standard deviation 2 in **natural-log units**, bounded by
ln(0.1) and ln(50). Equivalently, gamma has the original HBI bounded
lognormal prior. Its median is 1 before truncation; 1 is not the mode of a
lognormal density with width 2. The Jacobian cancels the `-ln(gamma)` term
when writing the target in eta coordinates. This is tested explicitly.

## Likelihood

For each data type d, r_d = observed_d - predicted_d and
C_d = exp(2 eta_d) C0_d. With C0_d = L0_d L0_d', the cached energy is
E2_d = ||L0_d \ r_d||^2. We use

```
log L_d = -0.5 [N_d log(2 pi) + log det(C0_d)
                    + 2 N_d eta_d + exp(-2 eta_d) E2_d]
log posterior(theta, eta) = sum_d log L_d + log prior(theta, eta).
```

There is no ad hoc division by the number of samples and no silently
discarded part of a failed prediction. Covariance determinants and noise
normalization remain in the target. No proposed model is allowed to change
the set of observations used in its likelihood. A noise-only proposal
reuses predictions and whitened energies; it does not call any solver.

The RF is the native signed CPS R/Z transfer response, with p = 0.06 s/km,
Gaussian alpha = 2.5, dt = 0.025 s, and a 4096-point calculation. We select
every eighth native sample in -2 to 25 s (136 points at 0.2 s spacing).
The **added synthetic RF noise**, not the Earth response, has prescribed
exponential correlation exp(-|t_i-t_j|/0.4 s). The same covariance generates
the noise and whitens the residuals. The other four datasets have independent
Gaussian errors in this controlled test. Cross-data-type covariance is zero
by construction. This is not a model of the dependence of real HVSR,
ellipticity, or deconvolved RF observations; that must be assessed with real
data before applying this likelihood to a station.

Base standard deviations are max(5% of HVSR, 0.05), max(3% of ellipticity,
0.015), 0.03 km/s for phase, 0.04 km/s for group, and 0.015 RF amplitude.
The true multipliers are all 1. Synthetic errors are added without clipping.
The forward-only example's zero-sigma arrays are deliberately not reused.

THBI's body-wave cross-convolution misfit is not the same observable as an
already formed R/Z receiver function. Its RF processing/misfit is therefore
not transplanted into this direct synthetic waveform likelihood.

## Sampling and reproducibility

Random-scan reflected Gaussian proposals update one coordinate at a time.
The mixture probabilities are state independent. Conditional Vs bounds use
the unchanged neighboring velocities, so forward/reverse proposals have the
same reflection interval. The MH ratio is the posterior ratio at T = 1.
Proposal steps adapt only during burn-in and are fixed for retained states.
This version uses single-stage proposals; legacy cyclic reheating and
delayed-rejection heuristics are not part of the target sampler.

All post-burn-in states, including repeats, contribute to the summaries.
The plotted "best visited fit" minimizes the sum of **base-covariance**
energies among visited states including burn-in. It is not a posterior mean,
MAP estimate, or evidence of uncertainty coverage. Its score cannot be
improved merely by inflating gamma. Retained intervals and classical split
R-hat are shown for inspection; a short pilot is not a convergence study.
Production interpretation also requires effective sample size, more dispersed
starts, longer chains, sensitivity to priors, and repeated noise realizations.

Each chain gets a separate local temporary directory and a separate output
folder on the filer. CPS `surf96` caches its previous model in `tmpsrfi.*`;
these files are cleared before every new phase/group calculation. A native
true -> perturbed -> true regression test guards against stale predictions.
No forward binaries, original source directories, or user inputs are changed.

Checkpoints save the current state, cached predictions, RNG state, proposal
tuning, chain histories, configuration, and fixed data. Writes use a temporary
file followed by rename. Resume rejects changed configuration/data and is
tested against an uninterrupted run. Use the same code revision and binary
versions when resuming. Known forward failures reject the proposal and are
counted; ten consecutive failures stop the run and preserve solver files.
Configuration/programming errors stop immediately. Process workers are
required for parallel execution because the external solvers change cwd.

## Running in MATLAB

From the repository root:

```matlab
run('DATA/run_ecm1_thin_hbi.m')  % 2 x 240 steps, first 120 burn-in
```

The default is serial and needs no parallel toolbox. To use two process
workers or launch a longer run:

```matlab
addpath('HBI'); a0_STARTUP_BAYES;
run_synthetic_hbi('pilot', struct('parallel',true));
% After inspecting the pilot, choose production settings explicitly:
run_synthetic_hbi('production', struct('parallel',true)); % 4 x 30,000
% Resume an interrupted run without changing settings:
resume_hbi('/absolute/path/to/DATA/inversions/ECM1_thin_HBI_...');
```

The separate launcher exposes `par` for editing before any computation. It
opens the three figures in MATLAB as well as saving them. Set
`par.res.show_figures = false` for a noninteractive run. After completing an
inversion, `run('DATA/plot_ecm1_thin_hbi.m')` reopens the saved plots without
recalculating the forward models or starting MCMC.

Outputs under `DATA/inversions/ECM1_thin_HBI_<profile>_<timestamp>/` include
`inputs.mat`, per-chain `checkpoint.mat` and `chain.mat`, `results.mat`,
`truth_forward/`, `joint_fit.png/pdf`, `diagnostics.png/pdf`,
`shallow_model.png`, and `inversion_report.txt`. Results are ignored by Git.
The launcher and parameter function are tracked. Originals in Downloads
remain untouched. `HBI_CPS_BIN`, `HBI_HVF_BIN`, and `HBI_SHARED_PROJECT`
configure external solvers as in the forward-only example.

Run checks with `addpath('HBI','tests'); a0_STARTUP_BAYES; test_hbi_core;
test_hbi_forward_reuse;`. The core checks use an analytic cheap forward map;
the latter executes all five native forward solvers in MATLAB.

This stage is CPU based. It establishes the layered model, likelihood, and
data interfaces needed before transdimensional moves or GPU solvers can be
introduced and validated.
