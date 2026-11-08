# Review against HBI and Zach/Brennan THBI

Current extension: [2026-11-07 Ps/Sp, 60-s bands and delayed-rejection review](THBI_NEXT_RUN_REVIEW.md).
The audit below describes the earlier five-block profile; its exclusions and counts are historical.

Reviewed 2026-11-03. Scope: the first **thin, isotropic, layered,
five-observable synthetic RJMCMC experiment**. Source review and static
analysis are complete; MATLAB tests, native calculations and recovery
experiments were not executed during this review.

The necessary model, prior, proposal, likelihood, forward, storage and
summary stages are present for that experiment. This establishes source
coverage, not runtime correctness or convergence. It is an adapted layered
implementation, not full feature parity with either source repository.
The subsequent [dependency audit](../THBI/DEPENDENCIES.md) traced all 45
MATLAB files, checked the 18 packaged helpers and inspected native binary
and direct library presence. Startup now checks every helper's resolved path.

## References inspected

- Original HBI production engine:
  `/Users/birotimi/Downloads/hbi_all_codes/Codes/Mod5_HBI/functions_hbi_new/`:
  empirical prior, parameter bounds, likelihood, feasibility, cyclic chain
  and delayed-rejection files. The Mod1 damped forward code and thin ECM1
  synthetic-generation script were also inspected.
- Current five-observable fixed-layer `CONUS-VS-UP/HBI/`: stage equations,
  native interfaces, noise handling, chain execution and results.
- `/Users/birotimi/Downloads/THBI_MCMC`: a0/a1/a3, b1/b2/b3/b5/b6/b9,
  model/prior assembly, chain execution, posterior construction and
  `SYNTHETICS/parms/bayes_inv_parms.m`. Reference commit:
  `585ee81f55961f31cf018163eabdbf972f5d0652`.

`THBI/IMPORTS.json` records source hashes. The original repositories and
the fixed-layer HBI implementation were not edited by this review.

## Findings corrected

| Finding | Correction in the new THBI |
|---|---|
| Moho condition did not match the active reference code: 30% relative to crustal Vs, with no absolute cap | a1 now uses 25% relative to the mean across the interface and a 0.8 km/s absolute cap. Derived `fdVSmoh` follows the reference percentage convention. The synthetic prior-count check disables both caps. |
| b3 caught unexpected MATLAB exceptions and turned them into generic rejected forward models | Only explicitly identified native-output failures are collected. Other exceptions propagate. The chain distinguishes the native-failure identifier from programming/configuration errors. |
| Failed truth or starting-model files could disappear during scratch cleanup | Preserve Earth, data, configuration, RNG, failure context, exception stack/report and native files. Preserve the first failed native proposal and any fatal proposal as well. |
| Sampled observable predictions were missing, although the original HBI and reference THBI store them | Save cached responses at a fixed post-burn-in iteration stride, including repeated states. Produce pointwise response bands without extra solver calls. |
| Mixing diagnostics omitted the inferred Vs profile | Calculate classical split R-hat at every output depth, plot shallow/full-depth diagnostics and write values to the report. Counts and scalar parameters remain checked too. |
| Combined birth/death acceptance could hide an immobile sediment or mantle zone | Report attempted, prior-valid and accepted moves, plus birth/death acceptance separately by zone. |
| Changes to the prior could silently exclude the prescribed truth | a2 explicitly assembles its zoned coordinates, checks support and verifies that the physical matrix is preserved before synthetic forward calculation. |
| Complex observations/covariances/predictions were not explicitly rejected | Add real-valued checks in data preparation and misfit calculation. |
| Path guards covered only some stages | Check all main stage implementations, parameter settings and model assembly to detect mixed HBI/THBI paths. |

The checkpoint schema is now `zoned_splitmerge_v2`. Code/configuration
snapshots deliberately prevent continuation from pre-review checkpoints;
start a fresh run with these changes.

## Coverage and deliberate adaptations

| Requirement | Source relationship | Current implementation |
|---|---|---|
| Inspectable numbered stages and separate helpers | Zach/Brennan organization | Full main-body a/b/c blocks; workflow `parms/` and `tests/`, plus shared repository `functions/` |
| Thin layered physical Earth | Original HBI h/Vs and empirical material relations | Same prescribed truth; piecewise-constant layers; shared Vp/density for every observable |
| Sediment/crust/mantle zones | Reference THBI zoned structures | Variable layer counts in all three zones, movable sediment base and Moho, fixed 70 km half-space top |
| Dimension changes | Reference b2 categories, rederived layered math | Exact invertible split/merge map with explicit move/zone/layer/auxiliary densities and unit Jacobian |
| Structural prior | Reference zone bounds/count idea, rederived density | Normalized bounded count, ordered-interface and independent/ordered-velocity factors; cross-zone conditioning |
| Hierarchical errors | HBI scale multipliers, THBI data-error separation | Five sampled log standard-deviation multipliers; full normalized Gaussian likelihood |
| Prior exploration | Both sources have empirical-prior stages | Independent conditioned-prior draws, plus an optional prior-only MCMC test |
| Diffuse-field HVSR | Original HBI damped HV-DFA calculation | HVf, 10 Rayleigh/10 Love modes, 1,000 body-wave integration points, apsv=ash=0.001 |
| Ellipticity | Existing five-observable HBI interface | HVf fundamental Rayleigh H/V, no Love/body-wave contribution |
| Phase/group velocity | Original/current HBI CPS | CPS surf96, fundamental Rayleigh phase and group velocity |
| Receiver function | Current five-observable HBI | CPS hrftn96 P RF, including conversions/reverberations; one configured ray parameter/filter |
| Same Earth for all five | Current HBI forward test | Identical flat elastic model, including half-space and numerical subdivisions; input equality asserted |
| Noise-only efficiency | THBI separation/current HBI cache | Reuse predictions and whitened energies; recompute noise-dependent normalization |
| Multi-chain runs and restart | Both sources, current HBI checkpoint design | Serial or process workers, separate scratch, independent RNG, source/data/configuration checks |
| Correct retained-state multiplicity | Required by the sampled target | Keep rejected repeats; no likelihood-based posterior culling or automatic chain deletion |
| Posterior profiles and curves | Both source result workflows | Common-depth Vs bands, layer-count masses, saved response bands, scalar/depth mixing and move reports |

The original production HBI wrapper interpolates HVSR onto observation
frequencies and permits partial prediction masks. The new five-observable
workflow supplies requested HVf frequencies directly and requires complete,
fixed data grids. Thus source physics and thin truth are retained, but the
whole original production workflow is not numerically identical.

The original HBI uses bounds attached to ten fixed layer indices and
global Vs monotonicity. The new sampler uses zone bounds and sediment
monotonicity, with crust/mantle reversals allowed, following the reference
THBI organization. These are different priors. Fixed-layer and variable-layer
results therefore compare both parameterization and prior choices.

The reference's spline-knot birth/death expressions cannot be substituted
for layer split/merge probabilities. The new target retains all
dimension-dependent normalization factors and the complete proposal ratio,
following the reversible-jump formulation in
[Hastie and Green](https://groups.seas.harvard.edu/courses/cs281/papers/hastie-green-2011.pdf).
See [the derivation](TRANSDIMENSIONAL_DESIGN.md). The reference's b7/b8 kernel
and likelihood resets are unnecessary here because every structural
proposal invokes the full forward calculations; no MINEOS linearized
kernels are used.

## Features intentionally outside this first experiment

- **Vp/Vs and anisotropy:** Vp and density are deterministic HBI functions
  of Vs. Reference THBI can vary crustal Vp/Vs and radial anisotropy. The
  present synthetic truth uses the same deterministic relations, so this
  is internally consistent but does not test Vp/Vs uncertainty.
- **Scope and depths:** this is the 70 km thin example. Sediment-base
  bounds are 0.05-2 km; zero-sediment rock sites, deep basins and deeper
  mantle structure need an expanded model/prior before CONUS applications.
  The new sediment zone ends at 130 m; the original 50 m label denotes
  the soft cap comprising the first two layers, not that zone boundary.
- **Data ingestion:** a2 currently generates synthetics. Observed-station
  ingestion, station elevation conventions, RF preprocessing and multiple
  RF traces/ray parameters are not implemented. Real-data RF filtering,
  time alignment, amplitude convention and errors must match the calculator.
- **Additional data/physics:** Love-wave observations, anisotropic
  dispersion, Sp RFs, HK stacks, spherical geometry and inferred material
  attenuation are not included. HVf integration damping is not a sampled Q.
- **Noise complexity:** five scale multipliers are sampled; covariance
  correlation lengths are fixed. No cross-observable covariance is modeled
  in this controlled synthetic test.
- **Proposal acceleration:** no cyclic reheating, delayed rejection,
  adaptation, parallel tempering, surrogate kernels or GPU code. These
  need separate derivations and verification if added to this RJ target.
  Process parallelism across chains is available.
- **Production diagnostics:** classical split R-hat is included; rank-based
  diagnostics, effective sample size and repeated noise-realization
  coverage experiments remain to be added for scientific production.

None of these is a missing dependency of the current controlled thin test.
They are material differences from a full regional production system.

## Validation remaining in MATLAB

1. Run `DATA/check_thbi_implementation.m`: algebraic RJ checks, Gaussian
   likelihood/cache/continuation tests, response storage/summaries and
   deliberate failure handling. These use no native solvers by default.
2. Enable `run_prior_check`: check count sampling against analytic 1/K
   masses with cross-zone conditions disabled. Inspect mixing if it fails;
   longer runs/several seeds may be needed. Reciprocity tests alone cannot
   establish correct stationary frequencies.
3. Enable `run_native_check`: exercise the five native solvers before a
   birth, after it, and after its inverse death in reused scratch folders.
4. Run `DATA/run_ecm1_thin_thbi.m` for a fresh pilot. Inspect failures,
   acceptance by zone, layer-count transitions, noise scales, all five fits
   and depth-dependent mixing before interpreting posterior bands.
5. Test longer independent chains, multiple starts and independent noise
   realizations before claiming recovery or uncertainty coverage.

The standalone MATLAB Code Analyzer reports no syntax errors. Remaining
notices concern unused warning-suppression comments and the two deliberately
disabled optional test blocks. Source inspection is not execution of these
checks. No new synthetic results or convergence claims were produced.
