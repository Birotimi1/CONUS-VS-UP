# Imported forward-solver provenance

These helpers were copied from the existing, previously exercised local HBI forward test.
External HVf/CPS executables are used in place and are not redistributed.

Source directory: `/Users/birotimi/Downloads/hbi_all_codes/AGU_models/five_observable_test`

| File | Original SHA256 |
|---|---|
| forward_five_observables.m | `bab0ac3e0f80746a19bded30e01055cabc8653ccb408a003a7f492f70b043a44` |
| forward_receiver_function_cps.m | `4e53f30c6bf80ddcba5544b55b9f0ede65f6b1cb142e7d24a94c0d733a76d62c` |
| quote_shell_path.m | `2076888706200497168eadb5ea542f4d88ecb315e195d36fb27ad11f17967e89` |
| read_cps_sac.m | `231cf17a378c832193120e1e024425e1b1cd384154bd95b5681244289bec8c82` |
| refine_layer_model.m | `8fae0d699e16d8b9d6627633aba3d5521fbc2103fa9fd5a12b19e74cfb448a18` |
| validate_layer_model.m | `c0e4aa47606a06e0ece6a56e7bc614620b3651c5a8b0569ccc576d3753db689f` |
| write_elastic_cps_model.m | `364c3e965ca28304acfd660fa885875ef07adc5fb8ab91a9ff1a6e7daa085be2` |

The prescribed ECM1 thin layer thicknesses, Vs values, and Brocher conversion equations are retained from `generate_ecm1_five_observables.m` in that source directory.
Only the thin case is included here. No observed dataset, prior/noise estimation, or MCMC is required by this forward-only example.
The source/refined layer matrices and exact solver inputs/logs are written for each new run.

For the hierarchical workflow, the MATLAB forward calculation blocks are
included directly in `HBI/b3_FORWARD_MODEL.m`, with section headings for
inspection and debugging. Low-level helpers are in repository `functions/`.
Synthetic generation and HBI
inversion both call that numbered stage. This relocation preserves the
solver commands, equations, units, and parsing; its default repository-path
calculation accounts for the different file location. The prior, parameter,
and model-construction implementations now reside in HBI's numbered a1, a3,
and b1 stages; older helper names remain compatibility entry points.

## Zoned layered THBI implementation (2026-11-03)

`THBI/` copies the existing HBI a0/a2/b3/b4/b5 calculation blocks and adapts
the run infrastructure. Forward utilities extracted from HBI b3 were initially
in `THBI/functions/forward/` and now reside in repository `functions/`;
numerical commands, model serialization, RF
settings and grid mapping are retained. Checkpoint and plotting utilities
also originate in the existing HBI implementation. Original file hashes
and local reference-file hashes are recorded in `THBI/IMPORTS.json`.

The architectural reference is Zach Eilon / Brennan Brunsvik's
[`THBI_MCMC`](https://github.com/brennanbrunsvik/THBI_MCMC), read locally at
`/Users/birotimi/Downloads/THBI_MCMC`, commit
`585ee81f55961f31cf018163eabdbf972f5d0652`. Its GPL-3.0 license is consistent
with the existing repository license. Its numbered stage organization,
three-zone parameter structures and proposal categories guide the setup.
No original MINEOS/spline/RF-engine dependency was imported. This is a new
layered split/merge formulation, not a verbatim port of the reference's
spline-knot sampler. The layer prior, RJ proposal densities, variable-size
history and common-depth result processing are implemented explicitly.

Neither Downloads source checkout was modified. Native executables remain
external and are used in place. New THBI runtime validation is pending;
see `docs/VALIDATION.md`.

The subsequent comparison review adds hashes of the original HBI production
files read, under `reviewed_original_hbi` in the import manifest. It changes
THBI error propagation and diagnostics while retaining native numerical
commands/parsing; it also aligns the active Moho tests with the reference
code and adds sampled-response/depth-mixing summaries. See `THBI_REVIEW.md`.

## Common project functions (2026-11-04)

The initial 30 shared MATLAB helpers and standalone synthetic functions were
consolidated in `functions/`. Their original project paths, local-function names, pre-move
hashes and final hashes are recorded in `functions/MIGRATION.json`.
Duplicate implementations were consolidated. The HBI/THBI numbered stages
retain their main calculation blocks. The original import hashes above and
in `THBI/IMPORTS.json` are historical and have not been rewritten.

The only non-relocation edits in the helpers are the default repository-path
calculation, explicit HBI compatibility-path guards, and source-snapshot
scope. Startup/worker paths and helper call names were updated accordingly.
Neither Downloads source checkout nor any native executable was modified.
No numerical execution was performed for this consolidation.

## Full function-suite comparison and adaptations (2026-11-04)

The user's full reference listing exposed missing retained-model processing
and plotting support beyond the pilot's dependency closure. All 127 original
helper files are now accounted for in `docs/THBI_FUNCTION_SUITE.md`. Original
implementations for the adapted helpers and c2-c5 stages were read in full;
the broader inventory identifies scope differences rather than claiming all
specialty code has been validated line by line.

Twelve reference helper roles were adapted into the layered schema (the two
heatmaps share one implementation), with one new step-evaluation helper.
The project has 42 helpers, nine in `functions/plotting/`, and four new
main-body c2-c5 stages. These are rewritten counterparts, not verbatim/API-
compatible copies of spline/anisotropic structures. Sources, disposition,
hashes and validation metadata are in `functions/THBI_SUITE_AUDIT.json`.
`functions/MIGRATION.json` remains the historical pre-extension relocation
record. The original GPL-3.0 architectural attribution and license remain.
No Downloads source, sampler proposal/acceptance equation or native executable
was changed for these additions. No numerical test or plotting script was run.

## Accepted-chain synthetic plotting adaptation (2026-11-06)

`DATA/plot_ecm1_thin_thbi_colored_chains.m` adapts the chain palette, Vs ensemble
views, depth zooms, interface side panels and start/final comparison roles from
our shared project `Figure_codes/plot_model_unified_colored_chains.m` and
`Figure_codes/plot_model_unified_all_2.m`. Their relevant model-selection,
summary and plotting blocks were inspected; the original files are unchanged.

The fixed-layer `theta` representation is replaced by per-state zoned model
reconstruction. Physical step curves retain every interface. Posterior depth
summaries preserve all retained repeats with equal weights, without the
reference plotting scripts' extra likelihood weights/percentile filters.
Best retained models minimize base-covariance residual energy and are not
called MAP. Five-observable overlays use the actual saved response iteration
IDs; no new predictions are generated. The new script and documentation do
not change the sampler, solver or checkpoint source files.

## 2026-11-07: layered RJ/DR and 60-s Ps/Sp extension

The original Downloads HBI delayed-rejection file and THBI spline-birth,
spline-death and delay_reject files were inspected again. They remain
unchanged. The new `functions/delayed_rejection_damp.m` is a reimplementation
of the second-try role using the full Green-Mira acceptance, not a verbatim
copy of the original likelihood-only retry. The layer split/merge remains
an independently derived reversible transformation with a unit Jacobian.
Existing THBI attribution/license terms remain in force.

CPS `hrftn96.f` from the official `NP330.Oct-12-2026.tgz` distribution was read
for component flags, sampling, filtering and signs. No new CPS source or
binary was installed into this repository. Sp uses the existing CPS SV
components followed by project-owned water-level deconvolution. See
[the processing details and primary-source links](THBI_NEXT_RUN_REVIEW.md).


## Targeted layered DR and 100-s extension (2026-11-08)

The new delayed-rejection equation is derived for our physical split/merge
maps. It is not a copied spline-knot factor or a claimed equation from the
Brennan supporting information. Q2's same-zone preference is an explicit
local choice (weight 0 recovers random second-zone selection). Sources,
reverse-path density/Jacobian terms and numerical tests are recorded in
[THBI_TARGETED_DR.md](THBI_TARGETED_DR.md). The Ps H–kappa–xi and Sp CCP
processing in the paper are distinguished from our CPS waveform synthetics
in [THBI_EXTENDED100_REVIEW.md](THBI_EXTENDED100_REVIEW.md). Original source
checkouts in Downloads were read, not edited.


User clarification on 2026-11-08 selected the paper's random second move.
`extended60` and `extended100` now use `layered_path`, no same-zone preference.
The full reverse-path acceptance correction is retained. The optional targeted
mode and its historical tests remain available for explicit comparisons.
See [THBI_RANDOM_SECOND_MOVE.md](THBI_RANDOM_SECOND_MOVE.md).
