# CONUS-VS-UP

> **Updated 2026-11-10:** the active `reference100_vti` test uses a published
> thin ECM1 column, labelled AK135 mantle to 250 km, linked crustal anisotropy,
> Ps H–κ and Sp CCP. Native launchers require the new CPS/HVf comparison gate.
> Algorithm checks, native numerical gates, two short native chains,
> checkpoint resume and postprocessing passed. The full recovery experiment
> remains user-run. Sediment Vs is monotone; crust/mantle LVZs are allowed.
> See the [current workflow](docs/THBI_WORKFLOW.md) and
> [completed execution audit](docs/THBI_NATIVE_AUDIT_2026-11-10.md).
> Older extended100/constant-mantle results below are historical.


The current second-stage split/merge algorithm and MATLAB check script are
documented in [THBI_RANDOM_SECOND_MOVE.md](docs/THBI_RANDOM_SECOND_MOVE.md).


[Complete workflow and script-by-script debugging map](docs/THBI_WORKFLOW.md) — updated 2026-11-12.

[Expanse preparation and submission scripts](hpc/README.md) — prepared locally;
remote installation and production performance still need verification.
[Damping and live-progress validation](docs/THBI_DAMPING_AND_EXPANSE_2026-11-11.md)
records the existing HBI damping, initial checkpoints and fatal timeout handling.

[Daily work logs](docs/daily/) and [local commit / GitHub preservation procedure](docs/DAILY_WORKFLOW.md).

Shear-wave velocity modelling of the continental United States using joint
HVSR, ellipticity, surface-wave dispersion, and receiver functions.

The **known thin layered Earth** is used for a clean forward calculation,
a fixed-layer hierarchical Bayesian inversion, and a new zoned
transdimensional implementation for all five observables. The new `THBI/`
code has completed a user-run two-chain thin synthetic test; the saved
diagnostics show that convergence and recovery need further work.
GPU acceleration remains future work.

## Current 250-km / 100-s Ps–Sp test (2026-11-10)

The `reference100_vti` profile uses ECM1 cell 19178 (10-m sediment),
labelled AK135 mantle to 250 km and assigned crustal xi=1.06. HVSR spans
0.1–25 Hz; ellipticity, phase/group velocity and Ps/Sp passbands reach
100 s. Ps uses H–kappa and Sp uses CCP products. Six likelihood blocks
represent five observable families; filtered RF views are never counted
as independent observations.

```matlab
run('DATA/generate_ecm1_ak135_five_observables.m') % known Earth, native gates and curves only
run('DATA/run_ecm1_thin_thbi_100s_ps_sp.m')  % preflight, then 2 x 15k / burn-in 8k
```

The second-stage kernel permits physical-layer births and deaths from the
rejected Earth, with the full delayed-rejection correction. Following the
user's clarification, the second move is random under the ordinary move/zone
mixture (`layered_path`), without a same-zone preference. Start a fresh run
(checkpoint schema v8); old waveform/isotropic checkpoints are incompatible.
The older profiles remain available for historical comparison.
The current 100-s launcher overrides the profile defaults to 15,000 iterations
and 8,000 burn-in per chain (7,000 retained states each).

**Validation:** MATLAB R2022a completed the algorithm suite, independent
CPS/HVf surface-wave comparisons, complete HVSR sums/refinement, RF controls,
full synthetic generation, two six-step native chains, checkpoint resume,
c1–c5 and plots. The short chains test execution only; their fits do not
demonstrate recovery. The audit's reduced execution grid does not replace
its full-band numerical gates. See the
[measured results and remaining limits](docs/THBI_NATIVE_AUDIT_2026-11-10.md).
Neither 100-s coverage nor passing execution proves 250-km resolution.
The user-deferred shared-Earth comparison remains **? OPEN**.


## New transdimensional workflow

The new `THBI/` folder follows Zach/Brennan's numbered stages with explicit
calculation blocks, a separate `parms/` directory and shared utilities in the
repository-level [`functions/`](functions/README.md) folder. Sediment, crust and mantle layer counts change through reversible
split/merge moves. CPS and HVf calculate all five responses from one Earth.

Run manually in MATLAB from the repository root:

```matlab
run('DATA/check_thbi_implementation.m')  % cheap checks; no native solvers by default
run('DATA/run_ecm1_thin_thbi.m')         % synthetic RJ pilot and plots
```

For the two-chain, 10,000-iteration test with 4,000 burn-in, run:

```matlab
run('DATA/run_ecm1_thin_thbi_2chains_10k.m')
```

This separate launcher includes implementation and native split/merge checks,
12,000 retained states across both chains, saved-state checks, a fresh final
forward calculation, and plots. It checkpoints every 100 iterations and runs
the two chains serially by default. Read its `test_report.txt` and the mixing
diagnostics before interpreting recovery. It has not been executed by the agent.

Edit `THBI/parms/bayes_inv_parms.m` for priors and run settings. Results go
under `DATA/transdimensional/`. Use `DATA/plot_ecm1_thin_thbi.m` to replot a
completed run. See the [stage-by-stage guide](THBI/README.md) and
[target/proposal equations](docs/TRANSDIMENSIONAL_DESIGN.md). The fixed-layer
HBI workflow below remains available for comparison.

For accepted Vs models colored by chain, exact variable-layer steps, shallow
zooms and per-chain final models, run
`DATA/plot_ecm1_thin_thbi_colored_chains.m`. It reads the completed results and
saves figures in that run's `Figures_THBI_Chains/` folder without sampling or
forward calculations. See the [plotting guide](docs/THBI_CHAIN_PLOTS.md).
The [comparison review](docs/THBI_REVIEW.md) records corrected gaps,
intentional differences from both source workflows and pending MATLAB checks.
The [full function-suite comparison](docs/THBI_FUNCTION_SUITE.md) accounts for
all 127 reference helpers and identifies the applicable adaptations and
unimplemented extensions. The synthetic run now saves a separate posterior,
Vs/Vp/density model suite and physical final-model selection, with extended
diagnostics. Rebuild these for an existing completed run with
`DATA/postprocess_ecm1_thin_thbi.m`.

## Run the hierarchical synthetic inversion

In MATLAB, from the repository root:

```matlab
run('DATA/run_ecm1_thin_hbi.m')
```

This runs a short two-chain pilot and writes fit, velocity, and chain-diagnostic
figures under `DATA/inversions/`. A pilot tests the implementation and fitting
behavior; it does not establish posterior convergence or unique recovery.
The new `HBI/` folder follows THBI's numbered MATLAB-function organization
and retains our fixed-layer priors and forward calculations.
See [function map, equations, source attribution, and run settings](docs/HBI_DESIGN.md).

The launcher opens figures in MATLAB. To replot a completed run without
sampling again, run `DATA/plot_ecm1_thin_hbi.m`.

The full stage implementations are in `HBI/a0_...` through `HBI/c2_...`,
with labeled sections for editing and breakpoints. Edit run/model/noise
settings in `HBI/a1_DEFINE_PARAMETERS.m`; Vp/density equations in
`HBI/b1_INITIATE_MODEL.m`; and the complete MATLAB interfaces to all five
forward solvers directly in the main body of `HBI/b3_FORWARD_MODEL.m`.
`HBI/a1_TEST_CONDITIONS.m` contains the constraint blocks, and
`HBI/a2_LOAD_DATA.m` contains the HBI thin truth, synthetic data and noise
preparation. The separate DATA launcher shows the full sequence through
plotting. See the [comparison with THBI's stage roles](docs/THBI_STAGE_COMPARISON.md).

To use MATLAB process workers: `addpath('HBI'); a0_STARTUP_BAYES;
run_synthetic_hbi('pilot',struct('parallel',true));`. Change the profile to
`'production'` for four 30,000-step chains after reviewing the pilot.

## Run the thin synthetic test

In MATLAB, open the repository and run:

```matlab
run('DATA/generate_ecm1_synthetics_five_observables.m')
```

Or from the repository root with GNU Octave:

```sh
octave --no-gui --quiet --eval "run('DATA/generate_ecm1_synthetics_five_observables.m')"
python3 scripts/plot_synthetic_thin.py
```

The plotter reads the latest completed run. For an earlier run, pass its path:

```sh
python3 scripts/plot_synthetic_thin.py DATA/synthetics/ECM1_thin_YYYYmmdd_HHMMSS
```

The six-panel PNG/PDF shows the prescribed Vs profile (including a 200 m
shallow inset), HVSR, ellipticity, phase velocity, group velocity, and P
receiver function. Plotting needs the packages in `requirements-plot.txt`;
it does not rerun or replace any forward calculation.

## Prescribed truth

This is the same thin ECM1 synthetic model used in the existing HBI test.
It is not a measured or recovered real-site structure. Thicknesses and Vs are
assigned; Vp and density use the original HBI Brocher polynomials and floors.
All values, including derived Vp/density, are exported in `true_model.csv`.

| Depth interval (km) | Thickness (km) | Vs (km/s) |
|---|---:|---:|
| 0-0.010 | 0.010 | 0.20 |
| 0.010-0.050 | 0.040 | 0.55 |
| 0.050-0.130 | 0.080 | 1.50 |
| 0.130-1.000 | 0.870 | 2.80 |
| 1-9 | 8 | 3.55 |
| 9-23 | 14 | 3.75 |
| 23-40 | 17 | 3.98 |
| 40-48 | 8 | 4.30 |
| 48-70 | 22 | 4.45 |
| Below 70 (half-space) | infinite | 4.54 |

The 50 m soft-cap label refers to the first two layers; other shallow velocity
transitions occur at 130 m and 1 km. The Moho is at 40 km.

## Five observables and shared physics

| Observable | Sampling | Forward calculator |
|---|---|---|
| Diffuse-field HVSR | 80 frequencies, 0.2-25 Hz | HVf: 10 Rayleigh + 10 Love modes, 1000 body-wave integration points |
| Fundamental Rayleigh ellipticity | 50 periods, 3-25 s | HVf: one Rayleigh mode, no Love or body-wave contribution |
| Fundamental Rayleigh phase velocity | 35 periods, 9-35 s | CPS `surf96` |
| Fundamental Rayleigh group velocity | 40 periods, 3-25 s | CPS `surf96` |
| Incident-P receiver function | 4096 samples, dt=0.025 s | CPS `hrftn96` |

Every calculator receives the **identical flat, isotropic, elastic layered
Earth**, including the half-space. Thick homogeneous layers are subdivided to
at most 2 km for numerical stability in HVf. This makes 40 numerical rows from
the 10 physical rows and adds no physical interfaces or parameters. Both
matrices are saved. CPS Qp=Qs=0 denotes infinite material Q; HVf body-wave
integration damping is `apsv=ash=0.001`, not inferred physical attenuation.

The receiver function uses p=0.06 s/km, Gaussian alpha=2.5, 10 s pre-arrival,
and doubled internal FFT length (`-2`). Its native signed R/Z transfer
amplitude is retained, including conversions and reverberations. The Moho Ps
marker is an independently calculated primary ray-theory lag, not a picked
waveform peak. The displayed window is -2 to 25 s; the full trace is saved.
No observed data are used. The forward-only example adds no noise;
the HBI and THBI inversion launchers add seeded synthetic noise.

## Executables and paths

Use a working CPS directory containing `surf96`, its companion programs, and
`hrftn96`, plus the HV-DFA `HVf` executable. On the shared filer, the launcher
finds the existing sibling `SharedData/HVSR-Joint-Inversion/AGU_models` binary
directories. On another machine set:

```sh
export HBI_CPS_BIN=/path/to/cps/bin
export HBI_HVF_BIN=/path/to/hvf/bin
```

`HBI_SHARED_PROJECT` overrides the sibling-project location.
`HBI_SYNTH_FIVE_OUTPUT` overrides the result parent directory. No executable
is copied into Git; each run records the resolved solver paths. The helpers
are MATLAB/Octave compatible; local validation uses Octave and Python.

## Outputs and checks

Each run creates `DATA/synthetics/ECM1_thin_<timestamp>/` containing:

- `synthetic_thin.mat`: truth, five clean observable arrays, configuration,
  shared input matrices, and forward status.
- `hvsr.csv`, `ellip.csv`, `cph.csv`, `ugr.csv`, `rf.csv`: calculated curves.
- `true_model.csv`, `solver_model.csv`: physical and numerical layers.
- `five_observables.png` / `.pdf`: generated by the Python plotting command.
- `manifest.json`, `forward_report.txt`, `validation.json`: settings and checks.
- Five solver subdirectories: exact inputs, native outputs, and logs.

`DATA/synthetics/LATEST_RUN.txt` points to the latest successful calculation.
Generated results are ignored by Git; code and documentation are versioned.
A failed or incomplete observable aborts generation and preserves diagnostics
rather than silently dropping samples. All five stored input matrices are
checked for equality. Clean `.sig` arrays are zero and must not be used as
likelihood standard deviations.

`tests/test_receiver_function_cps.m` checks a homogeneous medium's late
energy, two-layer Moho Ps lag, invariance to layer subdivision, and rejection
of invalid ray-parameter units. Add `functions` and `tests` to MATLAB/Octave's
path, then call `test_receiver_function_cps(cps_bin, new_output_directory)`.

See [import provenance](docs/PROVENANCE.md). Successful forward generation
demonstrates consistency of this configuration, not uniqueness or inversion
recovery.
