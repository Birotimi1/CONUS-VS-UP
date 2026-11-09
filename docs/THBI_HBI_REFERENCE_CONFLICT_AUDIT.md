# Current layered THBI versus old HBI and Zach/Brennan THBI

Source review: 2026-11-08–09. Baseline current commit: `70e80c7`.
This is a comparison of the local source, not a new native inversion or a claim
that any previous chain converged. No inversion code, priors or run settings
were changed during this audit.

The current engine is a distinct layered RJ-HBI algorithm. It preserves the
old HBI's material conversion and hierarchical error multipliers, and adopts
the reference THBI's staged organization and random second proposal from a
rejected model. Its forward operators, structural prior and posterior target
are not identical to either predecessor. Several differences are necessary
corrections; others remain scientific design choices or validation gaps.

## Source scope

- Current: `/Volumes/AS-Filer/EES/jbrussel/SharedData/US_TA/CONUS-VS-UP`.
- Old HBI: `/Users/birotimi/Downloads/hbi_all_codes`, especially
  `Codes/Mod5_HBI/functions_hbi_new` and `Codes/Mod1_forward_models`.
- Reference THBI: `/Users/birotimi/Downloads/THBI_MCMC`, using its main stages,
  `functions/master_subfuncs/delay_reject.m`, and
  `example_github/parms/bayes_inv_parms.m`. Other reference parameter sets can
  select different data and settings; the example is not every published run.

The companion [evidence record](THBI_HBI_REFERENCE_CONFLICT_EVIDENCE.json)
records inspected-source hashes and solver-free numerical checks.

## Workflow and edit locations

| Stage | Current entry / calculation | Comparison and consequence |
|---|---|---|
| Run settings | `DATA/run_ecm1_thin_thbi_100s_ps_sp.m:15` → `THBI/parms/bayes_inv_parms.m` | Launcher actually overrides profile defaults to **2 × 15,000 iterations, burn-in 8,000**. Profile defaults alone say 10,000/4,000. |
| Paths | `THBI/a0_STARTUP_BAYES.m:15` | Packaged stages/helpers take precedence and are checked with `which`; native CPS/HVf executables remain external dependencies. |
| Data | `THBI/a2_LOAD_DATA.m:8` | Builds the known Earth, calls b3, then adds seeded noise. This is synthetic generation, not the old station-datapack loader or THBI's real-data adapter. |
| Data covariance | `THBI/a2_LOAD_DATA.m:100`; `THBI/a2_PREPARE_DATA.m:6` | Coordinates, errors and whitening are fixed before sampling. The synthetic loader contains the block inline; keep it consistent with the supplied-data helper. |
| Prior and conditions | `THBI/a1_TEST_CONDITIONS.m`; `THBI/a3_CALC_LOG_PRIOR.m` | Explicit zoned support plus normalized, dimension-dependent prior density. Feasibility alone is insufficient for RJ acceptance. |
| Prior diagnostic | `THBI/a3_BUILD_EMPIRICAL_PRIOR.m` | Independent draws from the conditioned joint prior. Reference THBI uses prior MCMC; old HBI draws and sorts fixed-layer velocities. None means learning the prior from the synthetic observations. |
| Chain start | `THBI/b1_INITIATE_MODEL.m:5` | Default is a small jitter around one reference Earth and the same counts `[3 3 2]` in both chains. Independent seeds do not make these dispersed prior starts. |
| Proposal | `THBI/b2_PERTURB_MODEL.m` | Velocity, interface, boundary, split, merge, or log-noise. Physical layers replace spline coefficients/knots. |
| Model assembly | `functions/make_mod_from_parms.m` | Zone fractions → physical thicknesses; Vs → Vp/density. Numerical refinement does not change the sampled layer count. |
| Forward | `THBI/b3_FORWARD_MODEL.m` | HVf: HVSR and fundamental Rayleigh ellipticity; CPS surf96: fundamental Rayleigh phase/group; CPS hrftn96 and Sp deconvolution: RFs. Same refined Earth for every operator. |
| Misfit / likelihood | `THBI/b4_CALC_MISFIT.m`; `THBI/b5_CALC_LIKELIHOOD.m` | Fixed-covariance residual energies followed by a normalized Gaussian for each of six blocks. |
| First acceptance | `THBI/b6_IFACCEPT.m:15` | Log likelihood + log prior + reverse/forward proposal + Jacobian, at temperature one. |
| Second stage | `THBI/run_one_chain.m:141`; `functions/delayed_rejection_layered.m` | After a valid first rejection, draw randomly from that rejected Earth and apply the full reverse-path correction. Birth/death are eligible in both stages. |
| Store / restart | `THBI/b9_SAVE_RESULT.m`; `THBI/run_one_chain.m:236` | Store all iterations including repeats. Checkpoints preserve state/RNG and reject changed source/settings/data. |
| Posterior | `THBI/c1_PROCESS_RESULTS.m`; `c2_BUILD_POSTERIOR.m`; `c3_BUILD_MODEL_SUITE.m` | Pool post-burn-in states on common depths, not by layer-slot index. No likelihood-based culling. |
| Final fit / plots | `THBI/c4_FINAL_MODEL.m`; `c5_FINAL_FORWARD_MODEL.m`; `c2_PLOT_RESULTS.m` | Recompute one actual retained Earth; depth medians remain summaries. The representative minimizes base residual energy, not posterior density. |

## 1. RFs are different observables from the reference example

**Scientific mismatch if reference data or results are treated as interchangeable.**

The reference example enables `RF_Sp_ccp` and `HKstack_P`, alongside Rayleigh
and Love phase velocity. `SW_HV` is commented out in that particular example
(`example_github/parms/bayes_inv_parms.m:18`), although its engine supports it.

Its Sp path uses a depth-domain CCP comparison, P–SV processing, a simplified
parent pulse, and depth weighting that suppresses the shallow part relative to
the deeper part. Its b4 uses normalized cross-convolution error. Its Ps path
uses H–κ stack energy error, rather than directly fitting our R/Z trace.
See reference parameters:161–199 and `b4_CALC_MISFIT.m:33,119,191`.

Current Ps is a CPS R/Z transfer response; current Sp is regularized signed
Z/R from CPS SV components, retaining negative conversion lags. Both are
fitted as time-domain traces (`THBI/b3_FORWARD_MODEL.m:179,243`;
`THBI/b4_CALC_MISFIT.m:8`). There is no CCP migration, H–κ likelihood,
simple-parent substitution or reference Sp depth weighting in this workflow.

This is coherent for our matched synthetic experiment. It does **not**
reproduce their RF likelihood. Before real RF inversion, choose the measured
product and implement matching ray parameters, components, signs,
deconvolution, filters, sampling and uncertainty. Do not feed a CCP/H–κ
product into the existing waveform block. We do not need to adopt those
products merely because the reference uses them.

## 2. Vp/Vs and material assumptions follow old HBI, not reference THBI

**Intentional simplification with a real-data limitation.**

Old HBI derives Vp and density from Vs using Brocher polynomials and floors
(`c3_compute_likelihood_damp.m:328`). Current assembly and synthetic truth
retain those equations and their floor order (`functions/make_mod_from_parms.m:15`;
`THBI/a2_LOAD_DATA.m:20`). Vp/Vs is therefore depth-dependent but deterministic;
it is not an independently sampled parameter.

Reference THBI samples crustal Vp/Vs and crustal anisotropy in this example,
and uses separate mantle material relations. See its
`functions/make_mod_from_parms.m:42,58` and example parameters:63–86.
Its proposals can also perturb Moho depth with Vp/Vs (`b2_PERTURB_MODEL.m:424`).

Our synthetic generates and inverts with the same deterministic relations,
so this is not an internal truth/forward mismatch. It leaves the experiment
unable to demonstrate robustness to a different true Vp/Vs, density relation
or anisotropy. Such mismatch could be absorbed by layer depths, Vs or gamma
in real data. Test material-model mismatch before interpreting those results;
add a sampled Vp/Vs parameter only with an explicit prior and corresponding
updates to every forward operator.

## 3. Shallow structure has a substantially different prior and proposal geometry

**Important for the failed shallow recovery; cause not established without chain diagnostics.**

Old HBI has nine indexed thicknesses plus ten indexed velocities, specific
bounds for each, and global Vs monotonicity (`c1_define_hbi_prior.m:8–27`).
Current K varies in all three zones, with shared Vs bounds within each zone
and monotonicity only in sediments (`THBI/parms/bayes_inv_parms.m:45–70`).
Reference THBI uses a sediment gradient and variable spline knots in the
crust/mantle; its knot count is not our physical-layer count.

Consequently, losing the old fixed-layer depth and velocity restrictions
opens many more possible shallow structures. Correct birth/death balance
does not guarantee rapid exploration of that larger posterior. The true
interfaces at **0.01, 0.05, 0.13 and 1 km are inside current prior support**;
they are not categorically prohibited by the layered representation.

At the true zone boundaries (0.13, 40 and 250 km), the minimum layer
thicknesses are 0.26 m, 39.87 m and 210 m. Interface proposal standard
deviations translate to 5.2 m, 159.48 m and 2.1 km. These are consequences of
fractional coordinates, not a universal physical resolution. Moving a zone
boundary also moves internal interfaces in physical depth because their
fractions remain fixed (`THBI/b2_PERTURB_MODEL.m:102–114`).

Old HBI's second attempt used adjacent-layer compensation. Reference THBI
also has discontinuity mean/jump moves. Current scalar moves plus split/merge
do not include all these coupled moves. Per the user's choice, the active
second attempt is random, not the old compensating proposal.

The current split/merge preserves thickness-weighted **Vs**, not travel time
or predicted data. For example, merging the true first two layers changes
their vertical S transit time from 0.12273 s to 0.10417 s (−15.12%). That is
not a reversibility error, but it illustrates why deleting a strong shallow
interface can be hard to accept. Inspect per-zone stage-one/stage-two
acceptance, K transitions, gamma and depth-dependent mixing before selecting
a different reversible proposal or longer run.

Two additional prior differences matter beyond this thin test:

- Current sediment base is restricted to 0.05–2 km, so zero-thickness
  sediments and deeper sediment-zone bases are unavailable. Reference example
  permits 0–5 km. Current settings are not yet a general CONUS prior.
- Reference `a1_TEST_CONDITIONS.m:86` has an active near-Moho gradient
  condition. Current a1 retains positive/limited Moho jumps but has no
  equivalent finite-depth constraint on adjacent layered structure. A spline
  derivative rule cannot simply be copied onto discontinuous layers; define
  a physical-depth condition if this restriction is scientifically wanted.

## 4. The 250-km truth does not test recovery of a deep discontinuity

**Mismatch between test scope and a claim of deep structural recovery.**

Current truth extends the original 4.54 km/s half-space unchanged from
70 to 250 km and uses the same material below 250 km
(`THBI/a2_LOAD_DATA.m:15–18`). There is no true interface contrast at the new
base or anywhere below 70 km. `refine_layer_model` correctly merges the
redundant finite segment into its identical half-space; it does not truncate
a distinct feature. The inversion itself still has variable structure to
250 km and a separately sampled half-space velocity.

The present test exercises the larger search domain and longer-period
operators. Recovering a homogeneous deep mantle would be useful, but would
not establish recovery of a LAB, low-velocity zone or discontinuity near
250 km. Add a separate known deep-contrast truth, retaining the thin shallow
layers, and evaluate depth sensitivity and recovery. A longest period of
100 s is not by itself a depth-resolution test. The current forward model
also fixes the transition to its half-space at 250 km, so structure beneath
that boundary is not freely parameterized.

## 5. Hierarchical noise is consistent with old HBI, but joint weighting differs from THBI

**Same old-HBI gamma prior in new coordinates; different reference likelihood.**

Old HBI samples gamma through a log proposal, evaluates a lognormal density
including `-log(gamma)`, and adds the proposal Jacobian to acceptance.
Current code directly samples eta = ln(gamma) with a truncated normal.
The change-of-coordinate terms cancel correctly; adding another `-eta` or
Jacobian to current acceptance would double count it. Both use log-width 2,
center gamma=1 and bounds 0.1–50. See old `Functions/compute_log_prior_sigma.m:41`
and `c7_run_mcmc_cyclic_damp.m:217`, versus current `a3_CALC_LOG_PRIOR.m:38`.
The numerical comparison verified matching prior differences.

Current b5 uses the complete Gaussian with C = gamma² C0, including its
normalization. The six blocks contain 104 HVSR, 88 ellipticity, 68 phase,
78 group, 1,101 Ps and 1,101 Sp samples. RF correlation lengths are fixed at
0.4 s and 1 s. Synthetic noise is generated from that same covariance, so
the statistical experiment is internally consistent.

Reference b5 uses dataset `dof`, its own residual normalization and absolute
sigma parameters, averaging over records (`b5_CALC_LIKELIHOOD.m:28–62`);
the example specifies Sp dof=43 and H–κ dof=15. Those are not substitutes for
our N and covariance. Sample counts alone do not prove RF domination; inspect
whitened energies and per-block likelihood/noise histories. Estimate the
appropriate covariance for real processed data instead of borrowing the
reference effective degrees of freedom without its likelihood.

The 10–30, 30–60 and 60–100 s RF band curves are diagnostic views of the same
trace, not separate independent likelihoods (`receiver_function_bands.m:25`).
Thus long-period content is included, but it is not independently upweighted
as a separate dataset. Adding overlapping bands as independent observations
would require accounting for their dependence.

## 6. Acceptance deliberately corrects both legacy second-attempt rules

**Necessary algorithm difference; do not restore the old acceptance formulas.**

Old HBI's second attempt is judged by `delta_logL + log(tau)` without explicit
forward/reverse proposal or first-rejection factors
(`Functions/delayed_rejection_damp.m:58–84`). Its cyclic driver reheats during
burn-in and runs post-burn-in probes, flagging some saved samples invalid.

The local reference THBI second attempt calls `perturb_model` from the
rejected candidate and multiplies two `p_bd` terms, then uses its ordinary
acceptor (`functions/master_subfuncs/delay_reject.m:43–54`;
`run_one_chain.m:367`). Its b6 includes the temperature multiplier but no
two-stage rejection-probability correction. These are implementation facts
about the inspected checkout, not a claim about every version of the paper.

Current code uses normalized layer-count/coordinate priors, explicit
split/merge densities and unit Jacobians. Both second-stage rejection
factors are present, and hypothetical first-stage acceptance reuses b6.
The random second kernel can select birth/death as well as scalar moves.
Invalid-support first candidates and native failures are not retried.

Fresh checks passed for reciprocal transforms, absolute proposal densities,
1,080 physical two-stage paths and an exact 36-state transition matrix.
The latter has stationarity error 1.39e-17; its deliberately naive compound
retry has error 8.20e-4. This supports the correction. It is not proof that
the physical production posterior is well mixed.

## 7. Old HBI did not always use the same Earth or data mask for all observables

**Legacy inconsistencies avoided in the current engine; important for comparisons.**

Old `Codes/Mod1_forward_models/forward_model_damp.m:75,138` extends models
shallower than its 135-km threshold for dispersion, appending material with
Vs=4.4, Vp=1.78×Vs and density=3.0. HVSR and ellipticity are called earlier
with the original model. Therefore the sampled half-space is not necessarily
the same deep Earth used for the old dispersion calculation.

Current b3 refines once and checks that every prediction used the identical
model (`THBI/b3_FORWARD_MODEL.m:51,315`). Do not restore the hidden
dispersion-only extension. Old and new dispersion curves need not coincide
even with nominally equal sampled h/Vs. Regenerate all observables together
when comparing synthetic engines.

Old c3 also allows finite-prefix/partial predictions, changing retained
observations with solver output (`c3_compute_likelihood_damp.m:80` onward).
Current a2 fixes the data vector, and b3/b4 require every prediction. This
changes failure handling and prevents likelihood comparisons over changing
data subsets. The extended launcher stops on the first native failure and
preserves evidence; this is a debug setting, not a missing retry feature.

Reference THBI has surface-wave kernels, exact MINEOS refreshes and reset/
backtracking machinery. Current structural proposals call native CPS/HVf
directly, so the kernel-reset stages are intentionally absent. Both current
truth and inversion are flat, isotropic and elastic; replacing MINEOS does
not establish numerical equivalence to it at long periods. A native benchmark
remains necessary. There is no GPU implementation in this workflow.

## 8. Posterior summaries and convergence need different interpretation

**Current storage is appropriate; successful execution is not recovery.**

Reference c1 discards per-data likelihood outliers and can select/rank a
subset (`c1_PROCESS_RESULTS.m:76–117`). Old HBI thins and flags probe/stagnation
samples invalid (`c7_run_mcmc_cyclic_damp.m:319–355`). Current collation keeps
every post-burn-in state, including repeats. Accepted-only overlays are useful
visual diagnostics but are not the probability-weighted posterior.

Current c4 chooses the retained model with minimum sum(E2) under base C0.
That selection omits sampled gamma and prior terms by design: it is a
base-error best fit, not a MAP Earth. Posterior summaries use all retained
states and remain separate. Do not judge recovery only from its thick curve.

Current diagnostics include classical split R-hat for Vs(depth), counts,
boundaries and noise, but not rank-normalized R-hat or effective sample size.
Both chains start near the same reference Earth. Add dispersed-start
comparisons and stronger mixing diagnostics before claiming stable posterior
uncertainty; the 15,000/8,000 run is still an execution/recovery experiment.

## Recommended order of subsequent changes

1. Retain the corrected split/merge and full random delayed-rejection rules,
   fixed data vectors, shared Earth and repeated-state storage.
2. Diagnose shallow recovery using per-data fits/gamma, interface and K
   histories, and dispersed starts. Test additional reversible coupled
   proposals only against a defined target; do not transplant spline factors.
3. Add a separate deep-contrast synthetic and native depth-sensitivity checks.
   Preserve the current homogeneous continuation as a control.
4. Decide whether real RF input will be traces, CCP stacks or H–κ products;
   implement the matching data adapter, processing and likelihood.
5. Test sensitivity to Vp/Vs/material relations and revise zone bounds for
   CONUS stations. Decide explicitly whether a layered near-Moho restriction
   is wanted.
6. Benchmark native long-period forward responses and assess convergence
   before production or GPU migration.

## Validation and limits of this audit

Executed under Octave, without native CPS/HVf calls:

- `test_thbi_rj`: passed.
- `test_thbi_layered_dr`: 1,080 paths passed.
- `test_thbi_dr_stationarity('layered_path')`: passed, including negative control.
- `test_thbi_extended100`: truth support, grids, RF duration, 80-s passband
  gain and signed Sp checks passed.
- Additional source-comparison calculations: allowed crustal reversal,
  gamma-coordinate equivalence, physical step sizes, refined homogeneous
  continuation, data counts and the split/merge transit-time example.

An initial audit invocation incorrectly requested an output from the void
`test_thbi_layered_dr`; it was corrected and the test passed. A temporary
calculation initially compared a floating-point sum to 70 with exact equality;
it was corrected to a 1e-10 tolerance. Neither required a production-code fix.
The [raw test transcript](THBI_HBI_REFERENCE_CONFLICT_TESTS.log) preserves the
invocation error instead of concealing it.

No new native inversion, native RF/dispersion equivalence benchmark, real-data
fit, posterior convergence assessment or causation diagnosis of the earlier
shallow-recovery failure was performed. Previously recorded validation is
additional evidence, not a replacement for those outstanding checks.
