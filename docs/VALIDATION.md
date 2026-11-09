# Thin-model forward validation

Publication review (2026-11-08): the latest local 100-s launcher settings are
**two 15,000-iteration chains with 8,000 burn-in each**. Its run label and
console/saved-report counts now follow those settings. The actual settings
block passed configuration evaluation in Octave, and Code Analyzer reported
no launcher diagnostics. Earlier numerical manifests precede this launcher-only
edit; their historical hashes are preserved. See [launcher validation](THBI_LAUNCHER_VALIDATION.json).


Current configuration after user clarification: **random second move from the rejected Earth**.
See [THBI_RANDOM_SECOND_MOVE.md](THBI_RANDOM_SECOND_MOVE.md). The earlier targeted
validation below remains a historical record, not a claim that those runs used
the new default.


## Current targeted DR and extended100 validation — 2026-11-08

The current numerical/static record is [THBI_TARGETED_VALIDATION.json](THBI_TARGETED_VALIDATION.json),
with derivation and limitations in [THBI_TARGETED_DR.md](THBI_TARGETED_DR.md).
The solver-free suite passed in Octave, including the new 250-km/100-s
profile. A 150k controlled prior diagnostic passed after two failed 30k
checks; these do not establish production mixing. MATLAB license checkout
failed (-96); no native CPS/HVf or full inversion was executed by the agent.
The 122-file static scan has no syntax errors. See
[the new native preflight and data comparison](THBI_EXTENDED100_REVIEW.md).
The records below describe earlier stages and remain historical evidence.


## Layered second-stage RJ update (2026-11-08)

See [THBI_LAYERED_DR.md](THBI_LAYERED_DR.md) for the path derivation, affected
files and user-run checks. No numerical checks have been run by the agent.
The current static record is `THBI_LAYERED_DR_STATIC.json`; previous
`THBI_WORKFLOW_*` and `THBI_EXTENDED60_*` records are historical.


## Whole-workflow review (2026-11-08)

See [the workflow map](THBI_WORKFLOW.md) and [new audit findings](THBI_WORKFLOW_AUDIT.md).
Corrections align retained-fit plots with final-model selection, honour optional
output settings, distinguish undefined constant-chain R-hat, and validate
numerical/model/header inputs. Regression scripts were extended but not run.
Static verification is recorded in `THBI_WORKFLOW_STATIC.json` and
`THBI_WORKFLOW_MLINT.txt`. The earlier static manifests below remain historical.


## 60-s Ps/Sp and delayed-rejection update (2026-11-07)

See [the new reference comparison and test gates](THBI_NEXT_RUN_REVIEW.md).
Verification for this update is source review and standalone MATLAB Code
Analyzer only. New numerical/native tests are written but **not executed**:
`test_thbi_extended60`, `test_thbi_receiver_functions`, and extended60 modes
of the existing chain and failure tests. The DATA check script plots the
extended true observables and bands without production MCMC; the separate
run script performs the requested 2 x 10k / 4k burn-in test. Static results
and source hashes are recorded in `THBI_EXTENDED60_STATIC.json`.


## Saved user-run results and new chain plots (2026-11-06)

Read-only inspection found the completed user-run experiment
`DATA/transdimensional/ECM1_thin_THBI_2chains_10k_burn4k_20261105_003201/`.
Its saved reports record two complete 10,000-step chains, 4,000 burn-in,
successful state/count and native split/merge checks, and exact agreement
between the final recomputed and stored residual energies/likelihood.
The reports also record three native proposal failures in chain 1 (two after
burn-in), none in chain 2, and substantial between-chain discrepancies in the
split R-hat diagnostics. Completion does not establish convergence/recovery.
These are inspected user-run records, not new agent-executed tests.

`DATA/plot_ecm1_thin_thbi_colored_chains.m` now adapts the supplied plotting
scripts for this storage format. Standalone Code Analyzer inspection found no
syntax errors; its two unreachable-code notices correspond to the disabled
hidden-figure and optional `.fig` branches under the default display settings.
Its project calls resolve to startup, model reconstruction/conditions,
step sampling and quantiles. It makes no sampler, native-forward or shell
calls. The saved MAT metadata agrees with the expected chain and prediction
array layout. MATLAB execution and visual verification of the new figures
were left to the user. See [plotting definitions and outputs](THBI_CHAIN_PLOTS.md).

## Extended transdimensional test launcher (2026-11-05)

`DATA/run_ecm1_thin_thbi_2chains_10k.m` is prepared for manual MATLAB execution:
two serial chains of 10,000 iterations, 4,000 burn-in, 6,000 retained states per
chain, and 100-step checkpoints. Its main body exposes the settings, native
split/merge experiment, saved-state/count checks, posterior/model-suite stages,
fresh final-model forward comparison, plotting and diagnostic report.

The launcher invokes the four existing short implementation checks before
generating the five synthetic datasets. Native split/merge comparisons use
the actual run grids and a reused solver directory; the final check compares
recomputed residual energies and log likelihood with the selected retained
chain state (relative tolerance 1e-6, denominator 1+absolute saved value).
Solver-failure totals and accepted post-burn-in births/deaths are reported
separately from execution completion. No R-hat or recovery threshold is
treated as an automatic guarantee of convergence.

This addition changes no sampler, prior, likelihood or solver calculations.
Only standalone static source analysis is performed by the agent; MATLAB,
native checks, inversion, figures and recovery assessment remain for the user.

## Previously executed thin forward example

Executed on 2026-11-02 using GNU Octave 11.3.0 and the existing macOS ARM64
HVf/CPS binaries on the shared filer. Numerical execution completed for all
five observables; the Python figure was rendered from those saved outputs.
That initial forward example was run in Octave; subsequent MATLAB checks
are described below.

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

## HBI implementation checks (MATLAB R2022a)

Before the instruction to leave execution to the user, the following checks
completed in MATLAB 9.12.0.2009381 (R2022a Update 4), MACI64:

- `test_hbi_core`: analytic correlated Gaussian likelihood, noise
  normalization, transformed lognormal prior, invalid data/covariance,
  prediction reuse for noise-only proposals, exact interrupted/resumed
  histories, independent chain seeds/files, and storage of rejected states.
  These use a cheap analytic forward map, not a geophysical inversion.
- `test_hbi_forward_reuse`: true -> perturbed -> true using all five native
  solvers in the same scratch directory. This exposed and corrected CPS
  `surf96` reuse of its `tmpsrfi.*` model cache. After the correction all five
  responses changed under perturbation and returned to their original values.
  Synthetic generation and inversion predictions agreed to 1e-10.

No five-observable MCMC inversion was launched. The final run launcher,
result processing, MATLAB figure display/export, and production recovery
remain for the user to exercise in MATLAB. No inversion results, posterior
convergence, or successful recovery are claimed. The code is ready to launch
with `DATA/run_ecm1_thin_hbi.m`; it does not launch itself on opening a file.

## Numbered-stage refactor

The user subsequently reported a successful run and requested that the
numbered stage files expose the full code for debugging. Settings and prior
equations were moved to a1/a3, model construction to b1, and the full MATLAB
forward interfaces and their helper bodies into b3. Synthetic generation
now calls the same b3 implementation as inversion. Source comparisons
confirmed that the moved forward computations were retained (with the
default repository path adjusted for the new location), and that misfit,
likelihood, saved-state, summary, and plotting calculations were unchanged.
No MATLAB, Octave, native solver, or inversion was launched for this refactor;
runtime verification is left to the user as requested.

After a further direct comparison with THBI's a0/a1/a2/b1/b3, a1 condition
checks and a2 data preparation were given their corresponding stage roles.
The b1 initialization loop and b3 HVf/CPS/RF blocks now execute in their main
function bodies. The DATA launcher exposes the full run sequence. Static
comparison again confirmed preservation of the moved native solver command,
parsing and numerical blocks with scope/variable renaming accounted for.
No MATLAB or numerical solver execution was performed for these edits.

## Transdimensional implementation (2026-11-03)

The new `THBI/` stages, separate functions, tests and DATA launchers were
written without launching MATLAB, Octave, CPS, HVf or an inversion. The
standalone R2022a `mlint` Code Analyzer was used to inspect MATLAB source;
it does not execute the scripts. Source checks compare the copied forward
calculation, covariance and likelihood blocks with the existing HBI files.

Prepared checks, **not yet executed**:

| Check | Intended coverage |
|---|---|
| `test_thbi_rj` | Split/merge inverses, absolute proposal/prior ratios, unequal move probabilities, restricted split support, numerical Jacobian, count bounds, ordered-uniform normalization and Moho conventions |
| `test_thbi_chain` | Correlated Gaussian target, noise/cache, rejected states, response storage and summaries, depth diagnostics, variable-size histories, checkpoint continuation, seeds and invalid input |
| `test_thbi_failures` | Deliberate initial/native/programming failures; immediate exception propagation and preservation of failed models and native scratch |
| `test_thbi_prior_counts` | Optional prior-only count frequencies against analytic masses after disabling cross-zone conditioning |
| `test_thbi_forward` | Optional five-native-solver responses before/after a birth, inverse merge, common input Earth and cache reuse |
| `DATA/run_ecm1_thin_thbi.m` | User-run thin synthetic pilot, result processing and MATLAB plots |

The prior-count comparison is stochastic and needs adequate mixing; it is
not a proof of detailed balance. The default cross-zone-conditioned prior
does not generally have exactly the base 1/K count marginal, which is why
the analytic test disables those conditions and the main workflow estimates
its conditioned prior independently.

No successful RJ runtime, synthetic recovery, posterior convergence,
uncertainty coverage or GPU speedup is claimed. The existing fixed-layer
HBI validation above does not validate the new transdimensional kernel.

The subsequent source comparison is recorded in [THBI_REVIEW.md](THBI_REVIEW.md).
It corrected the Moho convention, exception handling and diagnostic
preservation; added sampled-response bands and depth-based mixing; and
expanded the unexecuted MATLAB checks. The new b3 retains native commands
and numerical parsing, with stricter exception handling. The fixed-layer
HBI and Downloads source files were left unchanged. The standalone Code
Analyzer was rerun after these edits; no MATLAB test or solver was executed.

## Shared helper consolidation (2026-11-04)

Thirty individually reviewed helper files now live in repository `functions/`.
The old project HBI/THBI helper folders and standalone `src` helper copies
were consolidated; the Downloads sources were untouched. Main calculations
remain in the numbered stages. Both startups and process-worker branches
select the new shared folder. Compatibility helpers check the active HBI
stage paths without switching algorithms.

Static verification covered 84 MATLAB files with the standalone R2022a Code
Analyzer and call listing. No syntax errors or unresolved project function
calls were identified. Notices concern style/redundant suppressions and the
disabled optional test branches. Source comparison covered all 30 helpers,
six affected stage bodies and all seven extracted HBI b3 utilities. The
source/migration hashes and per-file roles are recorded in
`functions/MIGRATION.json` and `functions/README.md`.

THBI checkpoints now track shared-helper source as well as the THBI source.
Pre-consolidation checkpoints have a different source snapshot and are
rejected for continuation. No numerical calculation, MATLAB test, native
solver or inversion was launched for these changes. The earlier successful
forward/HBI checks do not validate execution after this consolidation.

## Full reference-suite comparison and result stages (2026-11-04)

The 127-file reference inventory distinguishes dependency closure from feature
parity. Twelve helper roles were adapted, with an additional exact step-layer
evaluation helper; the two original heatmap files share a single new helper.
The project now has 42 helpers and separate c2 posterior, c3 suite, c4 final
summary and c5 final-forward stages. See `THBI_FUNCTION_SUITE.md` for every
reference file's disposition and the actual remaining scope.

`test_thbi_results` is added to the manual check launcher. It checks single-
chain burn-in removal, preservation of rejected-repeat weights, exact layer
counts, values on either side of an interface, half-space evaluation and
selection of an actual retained final model. It has not been run. The
optional c5 calculation is disabled by default in the user launchers.

Standalone Code Analyzer and call inventory cover the expanded source tree;
the current counts/findings are recorded in `functions/THBI_SUITE_AUDIT.json`.
No MATLAB/Octave execution, native calculation, graphical rendering, posterior
recovery test, or GPU work was performed. Static source inspection does not
establish runtime correctness or convergence.
