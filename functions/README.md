# Project MATLAB functions

[Complete workflow and script-by-script debugging map](../docs/THBI_WORKFLOW.md) — updated 2026-11-12.

This is the single helper folder for the fixed-layer HBI, zoned layered
THBI, and standalone thin synthetic workflows. The initial 30 files were
inspected individually during consolidation on 2026-11-04. The subsequent
full-suite comparison adds 12 helpers, including nine in `plotting/`, for
**42 MATLAB helper files** in that historical audit. The 2026-11-07 update adds
five helpers for **47 total**. The 2026-11-08 layered DR update adds one for
**48 helpers at that checkpoint**; see [the layered DR derivation](../docs/THBI_LAYERED_DR.md). Later reference/VTI helpers are also in this folder. There is no MATLAB helper dependency on either
Downloads source checkout. See the [127-file reference comparison](../docs/THBI_FUNCTION_SUITE.md)
for adaptations, covered roles and unimplemented extensions.

The numbered stages remain in `HBI/` and `THBI/`, with their main calculation
blocks visible for debugging. Those folders contain different algorithms
with some identical stage names; use the relevant DATA launcher from the
repository root. Each a0 adds this helper folder and checks MATLAB's stage
and helper resolution. Do not add every repository/source directory with
`genpath`.

## File-by-file guide

| Function | What it does | Used by |
|---|---|---|
| `birth_contrast_logpdf` | Exact narrow/broad Gaussian mixture density shared by birth and reverse death | THBI RJ |
| `delayed_rejection_layered` | Full two-leg layered RJ path acceptance, Q2 targeting and both rejection factors | THBI extended60/extended100 sampler |
| `delayed_rejection_damp` | Fixed-coordinate second-stage acceptance for the optional comparison mode | THBI sampler |
| `receiver_function_bands` | Explicit RF period passband and diagnostic views, without independent-data duplication | THBI RF |
| `waterlevel_receiver_function` | Signed regularized daughter/parent deconvolution; Sp uses Z/R | THBI Sp |
| `plotting/plot_RF_BANDS` | Full-band and short/intermediate/30-60-s RF views | Extended60 launchers |
| `make_mod_from_parms` | Converts zoned fractions/Vs to physical layers, Vp/density, and model metadata | THBI |
| `model_from_chain` | Reconstructs a variable-layer model from one saved chain state | THBI results/tests |
| `validate_configuration` | Checks zoned priors, moves, noise, run lengths and output grids | THBI |
| `reference_vti_compliance` | Independent displacement/traction equations and bounded impedance; analytic/expm controls; not used as the inversion forward solver | Native body-wave validation |
| `reference_vti_body_integrals` | Integrates independent compliance into separate PSV-horizontal, SH-horizontal and PSV-vertical body contributions | Native body-wave validation |
| `reflect_to_interval` | Reflects a Gaussian step into fixed conditional bounds | HBI and THBI proposals |
| `sample_quantiles` | Interpolates empirical quantiles along sample rows | HBI and THBI summaries |
| `split_rhat` | Computes classical split R-hat on fixed-size metrics | HBI and THBI summaries |
| `save_checkpoint` | Writes a temporary MAT checkpoint, then renames it | Both samplers |
| `source_snapshot` | Saves exact THBI and shared-helper source text for resume checks | THBI sampler |
| `preserve_forward_failure` | Saves failed model, data, settings, RNG, error and native files | THBI |
| `remove_scratch` | Removes the completed chain's temporary solver directory | HBI cleanup callback |
| `plot_layer_model` | Draws a stepwise Vs profile including the half-space | HBI and THBI figures |
| `write_hbi_report` | Writes fixed-layer fit, noise, parameter and mixing summaries | HBI figures |
| `write_thbi_report` | Writes variable-layer counts, fit, failure, move and mixing summaries | THBI figures |
| `find_bin` | Resolves configured, environment, PATH or shared native executable directories | All forward workflows |
| `merge_defaults` | Fills missing configuration fields without overwriting supplied values | All forward workflows |
| `quote_shell_path` | Quotes native paths for the POSIX shell | All forward workflows |
| `validate_layer_model` | Checks four-column isotropic / seven-column VTI matrices, tensor stability and CPS layer limit | All forward workflows |
| `refine_layer_model` | Merges identical material and subdivides thick layers for numerical evaluation | All forward workflows |
| `write_elastic_cps_model` | Writes flat elastic CPS model96 in isotropic or VTI form | All forward workflows |
| `read_cps_sac` | Reads little/big-endian SAC traces and RF header settings | All RF workflows |
| `forward_five_observables` | Coordinates the standalone clean five-response demonstration | Standalone synthetic launcher |
| `hv_prediction` | Writes HVf input, computes HVSR/ellipticity and maps its frequency output | Standalone forward demonstration |
| `dispersion_prediction` | Clears CPS cache, writes inputs and maps phase/group output | Standalone forward demonstration |
| `forward_receiver_function_cps` | Calls CPS hrftn96, checks SAC settings and computes reference Ps delays | Standalone demonstration/RF checks |
| `ecm1_thin_true_model` | Defines the prescribed thin ECM1 truth and empirical material properties | Standalone demonstration/HBI checks |
| `run_synthetic_thin` | Generates and saves five clean observables and a manifest | Standalone synthetic launcher |
| `write_synthetic_csv` | Writes a header and numeric forward-example CSV data | Standalone synthetic output |
| `write_synthetic_text` | Writes forward-example JSON and latest-run pointer text | Standalone synthetic output |
| `make_layer_model` | Compatibility entry point to the full fixed-layer HBI b1 implementation | Existing HBI callers/tests |
| `calc_log_prior` | Compatibility entry point to the full fixed-layer HBI a3 implementation | Existing HBI callers/tests |
| `count_N_layers_in_mod` | Counts finite physical layers by zone and distinguishes the half-space | THBI posterior and parameterization |
| `collate_allmodels_perchain` | Collects every post-burn-in state with chain/iteration identity, including repeats | THBI posterior |
| `sample_layer_model` | Evaluates Vp/Vs/density at specified depths with exact step interfaces | THBI model suite |
| `plotting/plot_HEATMAP_ALLMODELS` | Plots all-state Vs/Vp distributions at full or shallow depth | Extended THBI figures |
| `plotting/plot_SUITE_of_MODELS` | Draws a deterministic subset of retained profiles by chain | Extended THBI figures |
| `plotting/plot_TRADEOFF` | Plots scalar structural/noise tradeoffs across variable layer counts | Extended THBI figures |
| `plotting/plot_PARAMETERISATION` | Shows the physical zoned layers and half-space explicitly | Extended THBI figures |
| `plotting/plot_PERTURBATIONS` | Separates first-stage and second-stage attempted/accepted proposals | Extended THBI figures |
| `plotting/plot_PRIORvsPOSTERIOR` | Compares conditioned count masses and velocity envelopes | Extended THBI figures |
| `plotting/plot_LAYER_TRENDS` | Shows physical interface and count histories | Extended THBI figures |
| `plotting/plot_progress_chain` | Inspects complete/partial saved-chain progress without sampling | Extended THBI figures/manual inspection |
| `plotting/plot_RF` | Draws signed RF wiggles on explicit axes | Extended THBI RF figures |

The two compatibility entry points check that HBI is active. They never
change MATLAB's search path or select a different inversion algorithm.
The inversion b3 stages retain the HVf, dispersion and RF calculation blocks
in their main bodies; they use the common layer/IO helpers, not the standalone
demonstration coordinator. The synthetic truth used by each inversion stays
visible in that workflow's a2 stage.

## Migration and checks

`MIGRATION.json` is the historical record of the first 30-function relocation:
previous file/local-function locations and hashes at that time. The next full-suite snapshot is `THBI_SUITE_AUDIT.json`; the 60-s/DR update is
recorded in `../docs/THBI_NEXT_RUN_REVIEW.md` and its static-check manifest. Exact duplicates were
removed from the old project `HBI/functions`, `THBI/functions`, and `src`
directories. The source checkouts in Downloads were not changed. The earlier
import record in `THBI/IMPORTS.json` remains historical provenance.

The source guard now snapshots `THBI/**/*.m` and `functions/**/*.m`, with paths
relative to the repository. It does not scan generated DATA output. Old THBI
checkpoints have a different source snapshot and will be rejected for resume;
start a new run with this source layout.

Verification is static: source-body comparison, standalone MATLAB Code
Analyzer and call inventory. No MATLAB test, native solver or inversion was
executed during this consolidation. See the [dependency audit](../THBI/DEPENDENCIES.md).

CPS and HVf remain external executable installations, with their compiled
runtime libraries. They are not MATLAB functions and have not been copied
into this folder. Their paths are set in a0 / `par.forc` or through the
existing environment overrides.


## November 9 reference/VTI additions

The active stage map is [THBI_WORKFLOW.md](../docs/THBI_WORKFLOW.md).
`vti_elastic` and `vti_vertical_slowness` define the shared tensor/kinematics;
`check_native_vti_build` checks source/executable provenance. `layered_conversion_time`
uses those same VTI roots. `ccp_depth_window`, `ccp_depth_weights`,
`ccp_cross_convolution_misfit` and `rf_effective_dof` support the Sp product;
`Rotate_XZ_to_PSV` retains explicit Z-down conventions. Reference material
relations, `mcmc_rank_diagnostics`, `model_fit_score` and RF/anisotropy plots
are packaged here. See [import provenance](../docs/THBI_REFERENCE_IMPORTS.json).

The legacy `forward_five_observables`, `hv_prediction` and
`dispersion_prediction` entry points explicitly reject VTI input. The active
VTI native calculation blocks are in `THBI/b3_FORWARD_MODEL.m`.
