# Whole-workflow source audit — 2026-11-08

> Historical review. Current settings and executed tests are in [THBI_TARGETED_DR.md](THBI_TARGETED_DR.md) and [THBI_EXTENDED100_REVIEW.md](THBI_EXTENDED100_REVIEW.md).


**Subsequent update (2026-11-08):** extended60 now uses two full layered
proposal stages, including birth/death in stage two. See
[THBI_LAYERED_DR.md](THBI_LAYERED_DR.md). Fixed-scalar-only descriptions
and earlier static manifests below describe the previous implementation.


Scope: the active `DATA/run_ecm1_thin_thbi_60s_ps_sp.m` path, all numbered
THBI stages, the MCMC/DR controller, model/forward/IO/statistical helpers,
checkpointing, retained-state processing, plots, preflight and test scripts.
The separate fixed-layer HBI is included in static analysis because it
shares helpers. No original Downloads checkout was edited.

**Verification is source review and standalone MATLAB Code Analyzer.**
No MATLAB/Octave session, native solver, numerical test, MCMC or plot was
executed by the agent. The regression checks below are code written for
the user to run; they are not reported as passed numerical results.

The [workflow map](THBI_WORKFLOW.md) identifies every active stage, inputs,
outputs, data keys, saved files, control-flow branches and modification points.

## Corrections from this review

| Finding | Why it matters | Correction / regression coverage |
|---|---|---|
| c1/joint-fit plots used the best visited state, including burn-in, while c4/c5 used the best retained state. | A “final fit” figure and the fresh final forward calculation could describe different Earths. | `run_one_chain` now keeps exact `bestfit_retained` predictions and its iteration, independent of response stride. c1 and c4 use the same minimum-base-E2 selection and tie ordering; launchers/plots assert matching chain/iteration. `test_thbi_results` explicitly gives burn-in a better fit than all retained states and checks that it is excluded. |
| The dedicated 10k launchers unconditionally generated the empirical prior. | Setting documented `nprior=0` caused an assertion instead of skipping the comparison stage. | Both launchers now guard the stage; configuration validates a nonnegative integer. The prior definition in a3_CALC_LOG_PRIOR still applies when comparison draws are disabled. |
| The dedicated launchers assumed `prediction_stride>0` in post-run checks. | Disabling cached response histories could break validation after a costly completed run. | Zero-stride runs now expect empty prediction histories and skip nonexistent prediction fields. A toy-chain regression checks empty response summaries while retaining the exact best retained fit. |
| All zero-within-variance R-hat values were marked Inf. | Identical constant chains have 0/0, so Inf was an incorrect diagnostic classification. | Identical constant split chains return NaN; differing constant split chains return Inf. Both cases are tested. This does not claim either case establishes convergence. |
| Layer validation checked Vp squared but not positive Vp. | A negative Vp could pass the bulk-modulus inequality and reach native solvers. | Explicit positive Vp check; regression passes a negative-Vp half-space and expects rejection. |
| Finite prediction samples could still overflow the whitened residual energy; NaN proposal factors could reach acceptance arithmetic. | Numerical problems must not silently be converted into a likelihood rejection or acceptance. | b4/b5 check representable energies/likelihoods, both acceptance functions reject undefined ratios by raising errors, and the DR Gaussian ratio uses a factored difference. Error-path tests cover nonfinite energies, overflow and NaN MH factors. |
| SAC timing/filter/ray-parameter header NaNs were not explicitly screened. | NaN can evade `abs(header-request)>tolerance`, then cause a later grid/indexing failure. | Reader validates finite required headers. A synthetic malformed SAC file regression expects `fiveobs:sac`. |
| Fractional stop iterations/invalid chain IDs and invalid synthetic seeds were not caught at entry. | They could produce incomplete histories, unintended seeds, or fail only after preparation. | Explicit integer/range checks before chain setup/native work; toy tests exercise invalid values. |
| Postprocessing rebuilt c1 in memory without updating the saved summary used by replotting. | Later replotting could reuse stale summary metadata. | The user-run postprocessing script now saves its rebuilt summary back to `results.mat`. |

The model priors, six data blocks, period bands, geological constraints and
requested two-chain/10k/4k settings were not retuned in this review. The new
guards stop invalid numerical states; they are not evidence that the previous
default run encountered those particular failures.

## Selection and compatibility details

`chain.bestfit` remains the historical all-visited diagnostic. New chains also
store `bestfit_iteration`, `bestfit_retained`, and
`bestfit_retained_iteration`. The retained cache holds an actual state plus
its existing predictions/misfit. No additional native call is needed, and
it is available even when `prediction_stride=0`.

c1 chooses the minimum sum of base-covariance residual energies across all
retained states. Ties select the earliest retained iteration in the earliest
chain, matching c4's concatenated-state ordering. The figure's cached fit
and c5's fresh forward calculation therefore refer to the same Earth for
new-format runs. This is not MAP or a response computed from median Vs.

Old completed runs lack the new cache. Replotting remains possible, but their
available visited-fit curve is labeled **legacy**, possibly including burn-in.
The code does not invent missing retained predictions or silently recompute
them. Checkpoint schema is now `zoned_splitmerge_dr_v4`; old source snapshots
cannot be resumed under this edited source. Keep old results and start fresh.

## Algorithm relationships checked in source

- The true structure is declared in a2; its predictions are calculated by the
  same b3 used for trial and final models. True parameters are not supplied
  to the chain's proposal/acceptance functions.
- a2's inline covariance preparation matches a2_PREPARE_DATA. The latter is
  used for supplied data/toy tests, not called by the synthetic generator.
- b1 and b2 assemble physical layers with `make_mod_from_parms`; counts are
  independent of b3's numerical layer subdivision.
- a1 conditions, normalized a3 densities, reversible b2 split/merge proposal
  densities, unit Jacobian and b6 acceptance are expressed in the same
  interface-fraction/Vs/noise coordinates. The mixture density is shared by
  birth and reverse death.
- DR keeps the same scalar coordinate and K, starts its smaller proposal at
  the accepted state, includes both rejection factors, and does not apply
  to births/deaths or a failed first native solve.
- Noise-only moves reuse structural predictions/base E2 while recomputing
  the full normalized likelihood. Gamma therefore cannot be updated by a
  residual-only formula that omits the covariance normalization.
- b9 retains repeated states; c1/c2/c3 use their full iteration weight. Stored
  response samples are selected by iteration, not accepted-only thinning.
- Source/configuration/data plus RNG are guarded on resume. Native binaries
  are external and not hashed inside that checkpoint guard; keep them fixed
  for a run and its continuation.
- Fresh final-model E2/logL are compared with the exact saved selected row.
  Replotting and ordinary postprocessing do not execute native solvers unless
  the explicitly named optional final-forward switch is enabled.

## Numerical checks to run in MATLAB

From the repository root:

```matlab
run('DATA/check_thbi_extended60.m')
```

This already calls the revised `test_thbi_results`, `test_thbi_chain` (legacy
and extended60), `test_thbi_extended60`, RJ and failure tests before native
RF controls and full six-block split/merge checks. It makes truth/RF-band
figures and stops before production inversion.

An optional longer prior-only diagnostic now accepts either profile:

```matlab
addpath('THBI'); a0_STARTUP_BAYES();
test_thbi_prior_counts('extended60')
```

It uses no CPS/HVf, disables cross-zone conditioning to obtain an analytic
1/K count target, and enables the profile's DR/mixture machinery. It changes
count limits and contrast scales for this controlled experiment. It is a
stochastic diagnostic, not a guarantee that production chains mix.

After inspecting the preflight outputs, the requested inversion launcher is:

```matlab
run('DATA/run_ecm1_thin_thbi_60s_ps_sp.m')
```

## Remaining limits, not silently treated as solved

- Native Ps/Sp sign/timing/amplitude and stability tests are written but have
  not been executed by the agent. Neither have the new toy/regression tests.
- The current data stage is synthetic-only. Observed station ingestion,
  event/stack metadata and matching real RF processing are not wired up.
- Vp and density are deterministic functions of Vs; Q and anisotropy are not
  inverted. RF tradeoffs with those quantities are therefore restricted.
- Cross-data covariance is assumed zero; RF covariance shape is fixed and
  only its scale is inferred. Those assumptions need reassessment for real
  data, especially related dispersion/ellipticity measurements and RF bands.
- These are consistent-forward synthetics, not an independent check of the
  entire forward physics. HVSR mode-count/refinement convergence and model
  mismatch sensitivity have not been quantified by this source audit.
- The old two-chain run did not establish mixing/recovery. Fixed burn-in,
  classical split R-hat and a successful script exit are insufficient to
  establish uncertainty calibration. Rank-normalized R-hat/ESS are absent.
- Both long RF bands and a deeper variable model domain can add information;
  they do not guarantee recovery of a particular shallow interface.

The final standalone Code Analyzer scan covered **113 MATLAB files**,
including **47 shared helpers**, with no syntax errors. Remaining messages
are redundant suppression comments, branches controlled by disabled optional
switches, one indentation warning, one find-index performance suggestion and
one unnecessary comma in the legacy HBI tests. All 113 recorded source
hashes match the reviewed files; the workflow documentation links resolve.

Static diagnostics and source hashes are saved in
[THBI_WORKFLOW_STATIC.json](THBI_WORKFLOW_STATIC.json) and
[THBI_WORKFLOW_MLINT.txt](THBI_WORKFLOW_MLINT.txt). The Code Analyzer call
inventory is [THBI_WORKFLOW_CALLS.txt](THBI_WORKFLOW_CALLS.txt); it is a
static inventory, not an executed call trace. The earlier
`THBI_EXTENDED60_*` files remain the historical prior audit.
