# Zoned transdimensional joint inversion

The current second-stage split/merge algorithm and MATLAB check script are
documented in [THBI_RANDOM_SECOND_MOVE.md](../docs/THBI_RANDOM_SECOND_MOVE.md).


[Complete workflow and script-by-script debugging map](../docs/THBI_WORKFLOW.md) — updated 2026-11-08.

First layered-Earth RJMCMC implementation for HVSR, ellipticity, Rayleigh
phase and group velocity, and P receiver functions. The numbered stages
follow the organization of Zach Eilon / Brennan Brunsvik's local
`THBI_MCMC` checkout. Their main bodies contain the actual stage calculations;
low-level utilities are in the repository-level [`functions/`](../functions/README.md) folder,
shared with HBI and the standalone thin forward example.

**Status: the user-run two-chain 10k synthetic has saved completed results
and successful saved-state/final-forward checks. Its reports show native
proposal failures and substantial between-chain differences; convergence
and recovery are not established.** The new chain-colored plotting script
has been statically checked but not executed by the agent. The existing
fixed-layer `HBI/` workflow is retained separately.
See the [comparison review and remaining validation work](../docs/THBI_REVIEW.md).
The [dependency audit](DEPENDENCIES.md) lists the packaged helpers and
external MATLAB/native requirements.
The [complete function-suite applicability inventory](../docs/THBI_FUNCTION_SUITE.md)
covers all 127 reference helper files. The earlier dependency audit alone
did not establish full feature parity. The project now contains 48 helpers,
including a dedicated `functions/plotting/` folder.

## Current 250-km / 100-s Ps–Sp test (2026-11-08)

The `extended100` profile preserves the thin shallow truth and extends the
variable mantle domain to 250 km. HVSR reaches 0.01 Hz; ellipticity, phase
and group velocity reach 100 s. Ps and Sp passbands reach 100 s, with
819.2-s native records and a separate 60–100-s diagnostic band. Six
likelihood blocks represent five observable families; filtered RF views
are never counted as independent observations.

```matlab
run('DATA/check_thbi_layered_dr.m')        % kernel tests; no native solvers
run('DATA/check_thbi_extended100.m')        % native forward preflight and figures
run('DATA/run_ecm1_thin_thbi_100s_ps_sp.m')  % preflight, then 2 x 15k / burn-in 8k
```

The second-stage kernel permits physical-layer births and deaths from the
rejected Earth, with the full delayed-rejection correction. Following the
user's clarification, the second move is random under the ordinary move/zone
mixture (`layered_path`), without a same-zone preference. Start a fresh run
(checkpoint schema v6); old checkpoints contain a different configuration.
The 60-s profile remains available for comparison.
The current 100-s launcher overrides the profile defaults to 15,000 iterations
and 8,000 burn-in per chain (7,000 retained states each).

**Validation:** numerical kernel, exact balance, continuous-prior ensemble,
RF filtering and toy-chain/restart tests have run in Octave. MATLAB could
not obtain its license (-96). Native CPS/HVf, MATLAB-specific summaries and
the new full inversion have not been executed by the agent. Short prior-only
count diagnostics exposed slow/uneven sampling and must not be interpreted
as proof of convergence. See the current validation record in
[the earlier DR derivation and tests](../docs/THBI_TARGETED_DR.md) and
[the RF/data comparison and 250-km settings](../docs/THBI_EXTENDED100_REVIEW.md).
The constant deep true mantle contains no new discontinuity; 100-s coverage
alone does not establish resolution at 250 km.


## MATLAB entry points

From the `CONUS-VS-UP` repository root, run these yourself when ready:

```matlab
% Algebraic RJ, likelihood, caching and checkpoint checks; no native solvers.
run('DATA/check_thbi_implementation.m')

% Known thin Earth -> five noisy synthetics -> two RJ chains -> saved plots.
run('DATA/run_ecm1_thin_thbi.m')

% Extended test: two 10,000-step chains, burn-in 4,000, checks and plots.
run('DATA/run_ecm1_thin_thbi_2chains_10k.m')

% Replot a completed run without sampling or native calculations.
run('DATA/plot_ecm1_thin_thbi.m')

% Accepted Vs models by chain, exact physical steps, depth zooms and final models.
run('DATA/plot_ecm1_thin_thbi_colored_chains.m')

% Rebuild posterior/model-suite/diagnostics from a completed run, without MCMC.
run('DATA/postprocess_ecm1_thin_thbi.m')
```

In `DATA/check_thbi_implementation.m`, the optional prior-count experiment
and CPS/HVf check are disabled initially. Set `run_prior_check=true` to
sample the known prior without likelihoods; set `run_native_check=true` to
check all five responses before a layer split, after the split, and after
the inverse merge. The default cheap chain check is an implementation test,
not a geophysical inversion.

Edit `parms/bayes_inv_parms.m` for model bounds, move probabilities and
steps, data errors, solver settings, and chain lengths. The pilot is two
serial chains of 1,000 iterations with 500 discarded as burn-in. These
lengths are debugging defaults, not a convergence recommendation. The
`production` profile provides four 30,000-step chains, burn-in 10,000.
Select the profile in the DATA launcher. `par.inv.parallel=true` requires
MATLAB process workers and Parallel Computing Toolbox; thread workers are
rejected because native calculations change working directories.

The dedicated `DATA/run_ecm1_thin_thbi_2chains_10k.m` launcher sets the requested
2 x 10,000 iterations and 4,000 burn-in explicitly in its main body. It retains
6,000 states per chain (12,000 total, including rejected repeats), saves 600
cached response draws per chain, and checkpoints/prints progress every 100
iterations. It runs serially by default with independent chain seeds. Its
settings do not change the shorter pilot or the production profile.

Before native sampling, this launcher runs the four short implementation
checks and a CPS/HVf split/merge response check using the actual synthetic
data grids. After sampling, it checks all stored physical states, accepted
count transitions, complete response histories and retained-state totals;
builds the posterior/model suite; recomputes the final retained Earth's five
responses and compares its residual energies/likelihood with the saved chain;
and creates the diagnostic figures. Native preflight and final-forward checks
are enabled for this launcher. Set `run_implementation_checks=false` only to
skip the four existing algebraic/toy checks; it does not skip native checks.

Results go to `DATA/transdimensional/ECM1_thin_THBI_2chains_10k_burn4k_<timestamp>/`.
Inspect `test_report.txt`, `inversion_report.txt`, `joint_fit.png`,
`layer_counts.png` and `diagnostics.png` first. The launcher reports solver
failures and lack of retained birth/death acceptance as diagnostic issues;
execution completion does not establish convergence or successful recovery.
Each run keeps its own launcher copy, inputs, native preflight records, chain
checkpoints, posterior/model suite and final forward results. The existing
`resume_thbi(result_dir)` can continue interrupted chains with unchanged code
and settings; it performs ordinary postprocessing but does not execute this
launcher's additional validation/report blocks or fresh c5 check.

For the new plots adapted from our shared `Figure_codes` scripts, use
`DATA/plot_ecm1_thin_thbi_colored_chains.m`. It shows every accepted
post-burn-in model by default, uses all retained states for uncertainty
summaries, and saves PNG/PDF figures under the run's `Figures_THBI_Chains/`
folder. It also distinguishes the best retained model from the last chain
state. See [plot selection and output guide](../docs/THBI_CHAIN_PLOTS.md).

For breakpoints, useful variables are `model`, `ptb`, `trial`, `state`,
`misfit`, and `par`. In MATLAB, `dbstop if error` stops on unhandled errors.
The optional `request` argument of b2 forces a move for debugging; normal
sampling does not pass it.

Do not use `addpath(genpath(repo))` to mix this folder, fixed-layer HBI, and
the original THBI checkout. Run the DATA launcher from the repository root;
a0 adds this implementation's directories and checks stage resolution.

## Main stages and separate functions

| File | Calculations visible in the main body |
|---|---|
| `a0_STARTUP_BAYES.m` | Project, parameter, function, solver and output paths; path-collision checks |
| `parms/bayes_inv_parms.m` | All run, model, proposal, noise, synthetic and output settings |
| `a1_TEST_CONDITIONS.m` | Zone boundaries, layer counts, thickness fractions, velocities, physical jumps and noise support |
| `a2_LOAD_DATA.m` | Thin truth, Vp/density, configured forward predictions, observation selection, covariance, seeded noise and native-file saving |
| `a2_PREPARE_DATA.m` | Fixed covariance validation and Cholesky whitening for supplied data structures |
| `a3_CALC_LOG_PRIOR.m` | Count mass, dimension-dependent interface/velocity normalization, boundary and log-noise densities |
| `a3_BUILD_EMPIRICAL_PRIOR.m` | Independent conditioned-prior models, depth profiles and count distributions |
| `b1_INITIATE_MODEL.m` | Reference/prior initialization, zone parameters, independent jitter and data errors |
| `b2_PERTURB_MODEL.m` | Velocity, interface, boundary, birth, death and noise proposals; split/merge equations and forward/reverse densities |
| `b3_FORWARD_MODEL.m` | Shared Earth, HVf input/commands/parsing, CPS phase/group blocks, CPS RF calculation and fixed-grid predictions |
| `b4_CALC_MISFIT.m` | Whitened residual energies and amplitude RMS |
| `b5_CALC_LIKELIHOOD.m` | Normalized Gaussian likelihoods with sampled error multipliers |
| `b6_IFACCEPT.m` | Log likelihood, prior, proposal and Jacobian acceptance terms |
| `b9_SAVE_RESULT.m` | Variable-size states, rejected repeats, proposal diagnostics and a fixed subsample of cached five-observable responses |
| `run_one_chain.m` | RNG/checkpoint handling and full propose/forward/accept/store loop |
| `run_all_chains.m` | Serial or independent process execution |
| `c1_PROCESS_RESULTS.m` | Common-depth velocities, count and response quantiles, birth/death acceptance by zone and mixing by depth |
| `c2_PLOT_RESULTS.m` | Five-data fit bands, prior/posterior counts, chain histories, shallow profiles and velocity mixing |
| `c2_BUILD_POSTERIOR.m` | Every retained scalar/zone/count/noise sample, preserving repeated states and chain/iteration identity |
| `c3_BUILD_MODEL_SUITE.m` | Piecewise-constant Vs/Vp/density from every retained Earth on the common depth grid |
| `c4_FINAL_MODEL.m` | Pointwise material quantiles and a separately selected actual retained Earth |
| `c5_FINAL_FORWARD_MODEL.m` | Optional fresh all-response calculation and likelihood of that retained Earth |

`../functions/make_mod_from_parms.m` assembles physical layers and derives Vp
and density. `../functions/model_from_chain.m` reconstructs a saved state.
Other utilities handle reflection, quantiles, classical split R-hat,
checkpoint writes, source snapshots, validation and reports.
The same repository `functions/` folder contains shell quoting, path resolution, model
validation/subdivision, CPS model writing and SAC reading. The actual
native-solver commands remain in b3.

Both a0 startup stages select the common helper folder and check each helper's
resolved path. The old project-local helper directories have been removed;
neither Downloads checkout is needed for MATLAB function resolution. See the
[file-by-file guide and migration record](../functions/README.md).

## Earth and inferred quantities

The sediment, crust and mantle zones have independent variable numbers of
piecewise-constant layers. For a zone with K finite layers, infer K Vs
values and K-1 ordered internal interface fractions. Infer the sediment
base, Moho, half-space Vs, and five log-noise multipliers as well. The
half-space starts at the fixed 70 km base in this first thin test.

| Setting | Sediment | Crust | Mantle |
|---|---:|---:|---:|
| Allowed finite layer counts | 1-8 | 1-10 | 1-8 |
| Vs bounds (km/s) | 0.12-3.20 | 2.40-4.30 | 3.60-5.00 |
| Minimum fraction of zone thickness per layer | 0.002 | 0.001 | 0.001 |
| Monotonic Vs within zone | Yes | No | No |

The sediment base ranges from 0.05 to 2 km and Moho from 20 to 60 km.
Half-space Vs ranges from 4.1 to 5.1 km/s. Positive basement/Moho jumps and
a maximum 25% Moho Vs increase relative to the mean velocity across the
interface, plus a 0.8 km/s absolute jump limit, condition the prior. These
match the active reference THBI checks (its parameter-file comment says
30%, but its actual a1 threshold is 25%). The base count mass
is proportional to 1/K; after cross-zone conditions, its marginal need
not still be exactly 1/K. The empirical prior stage estimates that
conditioned distribution for comparison with posterior counts.

The truth is the existing thin ECM1 model: 9 finite layers to 70 km and
one half-space. Its zone counts are [3, 4, 2], sediment base 0.13 km, and
Moho 40 km. The 0.05 km soft cap refers to the first two layers. The
reference starting counts are [3, 3, 2], with a 0.25 km sediment base and
36 km Moho; the chain never receives the truth object as an argument.
Before calling a solver, a2 verifies the truth's zoned representation and
structural prior support. A changed bound that excludes the truth stops the
synthetic experiment with an explanatory error.

Vp and density follow the existing HBI empirical equations and floors.
One identical flat, isotropic, elastic matrix is supplied to every solver.
HVf computes diffuse-field HVSR and fundamental Rayleigh ellipticity;
CPS `surf96` computes dispersion and `hrftn96` computes the P receiver
function. Numerical subdivision to at most 2 km does not add inferred
interfaces or change K. The current code remains CPU based.

The five hierarchical variables scale standard deviations:
`C_d = exp(2*eta_d)*C0_d`. RF errors use fixed exponential time correlation
(0.4 s); other synthetic errors are independent. Data types are independent
in this controlled test. No model-dependent masks or dropped samples are
used. Noise-only moves reuse structural predictions. Full equations and
split/merge derivation are in
[the transdimensional design note](../docs/TRANSDIMENSIONAL_DESIGN.md).

## Results and continuation

The launcher creates `DATA/transdimensional/ECM1_thin_THBI_<profile>_<time>/`:

- `inputs.mat`, `truth_forward/`: truth, data, covariance, solver settings
  and native truth inputs/outputs.
- `prior.mat`: independent conditioned-prior samples summarized on the
  depth/count grids (set `par.res.nprior=0` to skip).
- `chain_XX/checkpoint.mat`, `chain_XX/chain.mat`: full variable-size
  histories, current predictions, sampled response arrays, RNG and
  source/configuration snapshots. `par.res.prediction_stride=10` stores
  cached responses every tenth retained state, starting at burn-in + 1;
  zero disables that storage. Rejected repeats are included. No extra
  solver calls are needed, and all structural states remain saved.
- `results.mat`: chains and common-depth/count summaries.
- `postprocessing.mat`: retained posterior quantities, the full Vs/Vp/density
  profile suite and final summary, saved with MATLAB v7.3. Storage grows with
  retained states times depth points.
- `joint_fit.png/.pdf`, `layer_counts.png/.pdf`, `diagnostics.png`,
  `shallow_model.png`, `velocity_mixing.png`, `inversion_report.txt`.
- Extended figures: full/shallow posterior heatmaps, model suite, scalar
  tradeoffs, physical parameterization, proposal acceptance, interface trends,
  prior/posterior comparison, per-chain progress and signed RF wiggles.

The final summary contains both pointwise medians and a representative
physical model. The latter is the retained state with minimum base-covariance
residual energy, excluding burn-in; ties use the first occurrence. It is not
a cross-dimension MAP. The optional c5 stage forwards this actual Earth, not
the pointwise median profile. Set `run_final_forward=true` in the run or
postprocessing script to enable that additional native calculation; that
legacy pilot switch defaults to false. Both dedicated 10k launchers always
perform their final forward check. Saved-response intervals still represent the sampled distribution.

New runs use the retained model with minimum base-covariance residual
energy for both the fit figure and final forward calculation. Its exact cached
response is stored even off the response-saving stride. Legacy completed
runs without that cache retain explicitly labeled visited-fit curves, possibly
including burn-in. Neither selection is a cross-dimension MAP or the response
of a posterior mean. Quantiles pool retained states; interpret
them only after checking mixing, birth/death acceptance and independent
chains. Layer counts exclude the half-space and numerical subdivisions.
Response bands summarize the retained model responses without adding
observational noise; they are not predictive intervals for new observations.
Classical split R-hat is supplied for Vs at every output depth, counts,
boundaries, half-space Vs and noise. Infinite values mean distinct constant split-chain segments with zero within
variance; identical constants have undefined R-hat (NaN). Rank-normalized R-hat and ESS are not
implemented in this first version.

To resume an interrupted run with unchanged source, settings and data:

```matlab
addpath('THBI'); a0_STARTUP_BAYES();
resume_thbi('/absolute/path/to/the/original/run/folder')
```

Resume completes the originally requested number of iterations. It does
not extend a finished chain or accept edited settings. Change the code or
prior in a new output directory. Per-chain native scratch uses local
temporary storage. Completed truth files are retained. Failed truth,
initial-model, first native-proposal and fatal-proposal cases preserve their
Earth, data, parameters, exception report and native files before cleanup.
Unexpected MATLAB errors propagate immediately. Known native output
failures are recorded as rejected states, with abort after the configured
failure limit; inspect them before interpreting a posterior.

## Source relationship

Reference: `https://github.com/brennanbrunsvik/THBI_MCMC`, local commit
`585ee81f55961f31cf018163eabdbf972f5d0652` in
`/Users/birotimi/Downloads/THBI_MCMC` (GPL-3.0 checkout).
Its numbered organization, zoned parameter structures, physical-condition
stage, proposal categories and empirical-prior/posterior workflow guide
this implementation. This is not a complete port of its spline algorithm.

The working five-observable forward and likelihood blocks were copied
from this repository's `HBI/`, with small forward utilities moved into
their own folder. The split/merge kernel, normalized layered prior,
variable-size storage and depth-based processing were written for this
parameterization. MINEOS, spline libraries and the original RF engines
are not dependencies. Original source hashes and reference-file hashes
are recorded in `IMPORTS.json`; neither Downloads checkout was edited.
