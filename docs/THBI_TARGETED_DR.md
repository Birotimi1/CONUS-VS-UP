# Targeted layered-Earth delayed rejection: derivation and numerical checks

> Configuration superseded by the user's clarification: the current production
> profiles use the paper's random second move (`layered_path`), with no zone
> preference. The full acceptance derivation below still applies with Q2=Q1
> from the rejected Earth and zero targeting correction. Targeted-mode results
> below are historical/optional comparisons. See [current setup and tests](THBI_RANDOM_SECOND_MOVE.md).


Updated 2026-11-08. This supersedes the random-second-proposal profile in
THBI_LAYERED_DR.md. The numbered scripts retain their calculation blocks.
The acceptance kernel has now been tested numerically in Octave. MATLAB
could not check out its license (error -96); MATLAB-specific summary and
failure-report checks and native CPS/HVf validation remain separate.

## What the Brennan supporting information says

Source: Brunsvik, Eilon & Lynner (2024), *Plate-scale imaging of eastern US
reveals ancient and ongoing continental deformation*,
[10.1029/2024GL109041](https://doi.org/10.1029/2024GL109041), supporting
information S1, PDF page 4, manuscript lines 60-71. The local source inspected
and rendered was
`/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Papers/Brennan 2024 supp.pdf`.

S1 motivates a second perturbation by parameter covariance. After rejection,
it describes randomly perturbing another parameter of the rejected model and
comparing the twice-perturbed result with the last accepted model. It does
not supply a targeted proposal density or a complete second-stage RJ
acceptance equation. The local `delay_reject.m` TODO at line 54 concerns
combining birth/death factors; ordinary `b6_IFACCEPT` is then used there.

Our targeted-zone proposal below is an explicit adaptation, not a method
claimed to appear in that paper. The general delayed-rejection framework is
[Green & Mira (2001)](https://people.maths.bris.ac.uk/~mapjg/abstracts/delayed.html),
DOI [10.1093/biomet/88.4.1035](https://doi.org/10.1093/biomet/88.4.1035).

## One first-stage rule, reused without consuming random numbers

Let pi(m) be our normalized hierarchical likelihood times the layered-model
prior, with its dimension-dependent count, interface and velocity factors.
For a labelled reversible move a->b, define the ordinary base kernel Q1:

```text
log R1(a,b) = log pi(b)-log pi(a)
           + log q1(b->a)-log q1(a->b) + log|J_ab|
alpha1(a,b) = min(1, exp(log R1(a,b)))
```

`THBI/b6_IFACCEPT.m` contains this calculation. Its optional fourth argument
`false` evaluates the ratio without drawing a random number. Both
alpha1(x,y) and alpha1(z,y) inside DR call this same function. The reverse
z->y descriptor swaps the SECOND leg's base forward/reverse densities and
negates its log Jacobian. No second, divergent copy of the first-stage
acceptance equation is maintained.

Using alpha1 consistently inside DR is essential. Simply using alpha1's
ordinary MH formula again on a twice-perturbed candidate omits the selection
created by the first rejection and generally changes the stationary target.
The exact finite-state negative control below measures this failure.

## Full second-stage rule for variable dimensions

Start at accepted x. Draw y with Q1 and reject it. Draw z from y using the
possibly targeted kernel Q2(z | x,y). Pair this path with the reverse path
z->y (rejected)->x. The labelled move includes its zone/index and the split
auxiliaries, so dimensions match in the augmented spaces.

```text
F = pi(x) q1(y|x) [1-alpha1(x,y)] q2(z|x,y)
B = pi(z) q1(y|z) [1-alpha1(z,y)] q2(x|z,y)

alpha2(x,y,z) = min(1, B/F * |J_path|)

log R2 = log pi(z)-log pi(x)
       + log q1(y|z)-log q1(y|x)
       + log q2(x|z,y)-log q2(z|x,y)
       + log[1-alpha1(z,y)]-log[1-alpha1(x,y)]
       + log|J_path|
```

Reverse the second RJ map, then the first. Their augmented-space dimensions
match even for two births, two deaths, or different-zone moves, and
`|J_path|=|J_xy|*|J_yz|`. Our physical split is

```text
r = (s-a)/(b-a)
v_left = v-(1-r)*delta
v_right = v+r*delta
```

Its inverse is the thickness-weighted merge. With auxiliaries (s,delta),
`|J|=1`. Both log Jacobians remain explicit. Birth density includes family,
zone, eligible-parent selection, the split-interval width and the complete
contrast-mixture density. Death uses the inverse map and its interface
selection probability. No spline-knot factors enter the calculation.

The first-stage probability in the reverse rejection is evaluated with Q1,
NOT Q2. `log(-expm1(logalpha1))` evaluates rejection probabilities stably.
If the reverse first stage would accept with certainty, alpha2 is zero.
Out-of-prior first proposals and native failures do not trigger this retry.

## Defined targeting law and its reverse density

The `extended60` and `extended100` profiles use `dr.mode='targeted_layered'` and
`dr.same_zone_weight=.5`. This means half of the second zone-selection mass
is placed on the first move's zone; the remaining half follows the usual
base zone probabilities. It does not force a particular daughter layer.
The weight may also be a three-vector for sediment/crust/mantle. Every
weight must lie in [0,1), retaining reverse support.

For base zone probabilities p_j and a first localized move in zone s:

```text
p2(j | s) = (1-w_s)*p_j + w_s*indicator(j=s)
t(s,j) = p2(j|s)/p_j
q2(z|x,y) = q1(z|y)*t(s,j)
```

Move-family probabilities, layer/index selection, scalar step scales and
birth contrasts are unchanged. For global moves (noise, zone boundaries or
half-space velocity), the targeting multiplier is one. In the reverse path,
the first move has the SECOND leg's zone j, and the second move has the
original first leg's zone s. Thus

```text
log R2 = log pi(z)-log pi(x) + H_xy + H_yz
       + log t(j,s)-log t(s,j)
       + log[1-alpha1(z,y)]-log[1-alpha1(x,y)]
H_ab = log q1(a|b)-log q1(b|a)+log|J_ab|
```

For equal targeting weights in every zone, t(j,s)=t(s,j), so this additional
ratio cancels. For different zones with unequal weights it becomes
`log(1-w_j)-log(1-w_s)`, which does NOT generally cancel. Tests exercise both
cases. Setting the weight to zero recovers the paper's random second zone;
`layered_path` remains an explicit untargeted comparison mode. `fixed_scalar`
remains the smaller same-coordinate retry from the accepted x.

b2 records **base Q1** densities in `log_q_forward/reverse` and the actual
Q2 multiplier separately in `log_targeting_forward`. The acceptance helper
constructs all four proposal densities explicitly. The chain additionally
stores `dr_log_targeting_ratio` and `dr_reverse_first_logalpha`.

## Numerical checks and implementation correction

- The original split/merge tests passed after correcting the empty interface
  array when a zone merges from two layers to one. The vector is now
  explicitly 0-by-1; this was exposed by executing the code in Octave.
- **1,080 physical paths passed**: 36 family pairs over zone combinations,
  tested with an untargeted kernel and targeting weights 0, .5, .9, and
  [.2 .5 .8]. Checks include actually generated inverse maps, absolute
  proposal densities, non-negligible rejection factors, detailed path
  balance, invalid reverse support and RNG-free first-stage evaluation.
- An independent **36-state variable-dimension transition matrix** enumerates
  648 second-stage paths. The full production acceptance rule has maximum
  stationarity error 2.0817e-17 and detailed-balance error 1.7347e-18.
- Negative controls fail: reusing ordinary compound MH without rejection
  corrections gives stationarity error 5.7019e-4; omitting the unequal-zone
  targeting correction gives 4.0247e-4. These controls demonstrate that the
  tests can detect the missing terms, not just reciprocal bookkeeping.

The finite-state experiment is an acceptance-kernel test with binary layer
parameters, not a geophysical resolution test. Continuous priors and the
actual layered proposals are checked separately by an independent prior
ensemble and the prior-count sampler.
These checks do not establish native forward accuracy or synthetic recovery.

## Run scripts and main implementation files

Run [DATA/check_thbi_targeted_dr.m](../DATA/check_thbi_targeted_dr.m) in MATLAB.
It also runs in Octave, explicitly omitting MATLAB-only posterior summaries
and exception-report formatting there. Optional 30k prior-only sampling is
controlled by `run_prior_counts`; no native solver or production inversion
is launched by this check script.

| File | Responsibility |
|---|---|
| [b2_PERTURB_MODEL.m](../THBI/b2_PERTURB_MODEL.m) | Layered transformations, base Q1 densities and targeted Q2 zone draw |
| [b6_IFACCEPT.m](../THBI/b6_IFACCEPT.m) | Single first-stage rule; deterministic evaluation option |
| [delayed_rejection_layered.m](../functions/delayed_rejection_layered.m) | Four proposal densities, two alpha1 evaluations, rejection factors and Jacobians |
| [run_one_chain.m](../THBI/run_one_chain.m) | Rejected-model second proposal, correct prediction reuse, histories and checkpoints |
| [test_thbi_layered_dr.m](../THBI/tests/test_thbi_layered_dr.m) | Continuous physical path tests |
| [test_thbi_dr_stationarity.m](../THBI/tests/test_thbi_dr_stationarity.m) | Exact finite-state invariance and deliberate negative controls |
| [test_thbi_chain.m](../THBI/tests/test_thbi_chain.m) | Analytic likelihood, every-state cache/count checks and exact restart |
| [test_thbi_prior_counts.m](../THBI/tests/test_thbi_prior_counts.m) | Actual continuous prior-only sampler and analytic count marginals |

Checkpoint schema is `zoned_splitmerge_dr_v6`; start a fresh inversion. The
requested two chains, 10k iterations and 4k burn-in remain. The later
`extended100` profile extends the depth and bands; see
[THBI_EXTENDED100_REVIEW.md](THBI_EXTENDED100_REVIEW.md). Original Downloads repositories and existing results are untouched.


## Continuous prior checks and sampling limitations

The test `test_thbi_prior_invariance` drew 6,000 independent Earths from the
exact continuous prior, then applied one production targeted RJ/DR transition
to each. Twenty-five count/velocity/boundary/noise metrics had maximum paired
drift **2.068 standard errors** (threshold 4.5); first/second stages accepted
3,997/41 moves. This avoids assuming that a random walk has mixed. It is a
stochastic check, not a proof.

A separate 30,000-step, 5,000-burn-in prior-only chain **failed** its count
frequency check with production Vs/interface steps: sediment K=1 was 0.2812
versus 0.4800 expected. A controlled diagnostic with wider prior-scale steps
improved sediment K=1 to 0.4743 but still failed the crustal count tolerance
(K=1: 0.3305 versus 0.4380 expected). Cross-zone constraints were disabled in
these diagnostics so the analytic 1/K count target applies. Neither run is
hidden or used to change production priors or proposal tuning. Long-chain
mixing and scientific posterior recovery are separate from balance tests.

The count diagnostic now has optional `niter`, `seed`, and `driver` arguments.
`driver='chain'` uses the full production histories/checkpoints. For large
prior-only experiments, `driver='kernel'` stores only counts/logP/stages and
calls the same b1/b2/a3/b6/DR functions. It first checks exact agreement of
500 states and acceptance stages with `run_one_chain` under the same seed.
This avoids large history-copy costs in Octave without changing the kernel.

The `extended100` profile also passed the 250-km truth/support checks,
0.01-Hz/100-s grid checks, 60–100-s band-isolation test and signed Sp test.
Its toy chain passed the analytic likelihood, every-state prediction cache,
layer-count accounting, repeated-state retention and exact restart checks.
Native deep RF controls have been written but not executed.


The longer controlled experiment completed **150,000 steps, 25,000 burn-in**
with `wide_prior_steps` and seed 42000. The compact driver first matched 500
production states exactly. Maximum absolute count-frequency errors were
0.004736 (sediment), 0.0326682 (crust), and 0.006840 (mantle), all below the
pre-existing 0.08 check threshold. This supports slow sampling as an
explanation of the short-run failures; it does not certify production-chain
mixing. The run used an identical copy of the relevant source on local
scratch storage to avoid shared-filesystem overhead; hashes are recorded.
Two initial full-history 150k attempts were stopped before completion because
of diagnostic storage/runtime costs; neither is counted as a passing run.

The complete `DATA/check_thbi_targeted_dr.m` run then completed with exit 0.
See [the numerical manifest](THBI_TARGETED_VALIDATION.json),
[full suite output](THBI_TARGETED_NUMERICAL_LOG.txt),
[long prior output](THBI_PRIOR_COUNTS_150K_LOG.txt), and
[static source inventory](THBI_TARGETED_STATIC.json). Native solver tests
and MATLAB-only summary/report/plot checks remain outstanding.

The final MATLAB Code Analyzer scan completed on all 122 source files
(48 helpers), with no syntax errors. Optional disabled branches, obsolete
suppressions and existing style/performance diagnostics are preserved in
[the analyzer log](THBI_TARGETED_MLINT.txt); static analysis is not execution.
