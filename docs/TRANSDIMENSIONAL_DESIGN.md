# Layered RJMCMC target and split/merge proposal

**Subsequent update (2026-11-08):** extended60 now uses two full layered
proposal stages, including birth/death in stage two. See
[THBI_LAYERED_DR.md](THBI_LAYERED_DR.md). Fixed-scalar-only descriptions
and earlier static manifests below describe the previous implementation.


This note specifies the target implemented in `THBI/`, so that proposal
and prior edits can be checked against the sampled measure. This is a
new piecewise-constant layered kernel guided by the organization of the
local Zach/Brennan THBI checkout; spline-knot proposal factors are not
reused as layer-split factors.

## State and density

For each zone z, K_z is the number of finite physical layers, f_z contains
K_z-1 ordered internal fractions, and v_z contains K_z velocities. Zone
boundaries are zsed and zmoh; maxz is fixed at 70 km for legacy profiles
and 200 km for `extended60`. There is a separate half-space velocity and
D eta = ln(gamma) variables (D=5 legacy, D=6 with Ps and Sp).
The continuous dimension is 2*sum(K_z)+D. A birth or death changes it by 2.

The target density uses counting measure on K and ordinary Lebesgue
measure on f, v, zsed, zmoh, half-space Vs and eta. Internal depths are
derived from fractions and zone thickness, so moving a zone boundary
while holding fractions fixed introduces no depth-coordinate Jacobian.
If the sampled coordinates change to absolute internal depths, the prior
and transformations must be rederived in that measure.

Before physical cross-zone conditions, the implemented priors are:

- Count mass p(K) = (1/K) / sum_{kmin}^{kmax}(1/k), or a uniform mass.
- Let epsilon be the minimum fractional layer thickness and
  S = 1-K*epsilon. The ordered interface density is
  (K-1)! / S^(K-1). It follows by transforming f_i to f_i-i*epsilon,
  which leaves K-1 ordered points in [0,S].
- For velocity range width W, independent uniform density is W^(-K).
  If within-zone ordering is imposed, the normalized ordered-uniform
  density is K!*W^(-K).
- Boundary and half-space densities are normalized uniform distributions
  on their configured ranges.
- Each eta has a normalized truncated normal density in eta coordinates.
  No additional gamma-coordinate Jacobian is included.

`a1_TEST_CONDITIONS` multiplies that joint density by the indicator for
positive basement/Moho jumps and the permitted Moho contrast. The global
conditioning constant cancels in MH. This choice can change the marginal
distribution of K; do not renormalize each count separately unless that
different prior is intentionally desired. `a3_BUILD_EMPIRICAL_PRIOR`
uses independent base-prior draws rejected by the same conditions to
estimate the actual conditioned prior.
The current Moho checks use `(v_mantle-v_crust)/mean([v_crust,v_mantle])`
at most 0.25 and `v_mantle-v_crust` at most 0.8 km/s, matching the active
reference THBI implementation. The base count diagnostic disables both
limits along with positive-jump conditions to recover an analytic marginal.

For data type d, C_d = gamma_d^2*C0_d and

```text
E2_d  = residual_d' * inv(C0_d) * residual_d
logL_d = -0.5 * [N_d*log(2*pi) + logdet(C0_d)
                 + 2*N_d*eta_d + E2_d*exp(-2*eta_d)]
```

The implementation uses Cholesky whitening, not an explicit inverse.
Data grids and C0 are fixed before MCMC. The normalization in eta remains
in the likelihood. These controlled synthetics assume independent data
types, with within-RF correlation represented in C0.

## Birth and its inverse death

Move-family probabilities are fixed. A zone has probability pz. For a
birth from K to K+1 layers, select one of E splittable layers uniformly.
For its fractional interval [a,b], choose

```text
w = b-a-2*epsilon
s ~ Uniform(a+epsilon, b-epsilon)
delta ~ g_z(delta)   # legacy Normal; extended60 fixed Gaussian mixture
r = (s-a)/(b-a)
v_left  = v - (1-r)*delta
v_right = v + r*delta
```

The split preserves the thickness-weighted mean velocity. The Jacobian
of (v,s,delta) -> (v_left,s,v_right) has determinant 1; the remaining
parameters pass through unchanged. A death selects one of K-1 internal
interfaces and recovers

```text
v = r*v_left + (1-r)*v_right
delta = v_right-v_left
```

With pb and pd the birth/death family probabilities, the birth densities
are

```text
q_forward = pb * pz / E / w * g_z(delta)
q_reverse = pd * pz / K
```

For death, q_forward = pd*pz/(K-1). Its reverse birth density uses E and
w computed in the merged model and the recovered delta. This includes
the selection of a splittable layer, not all layers. The auxiliary
contrast has an untruncated Gaussian or Gaussian-mixture density: out-of-prior velocities and negative
contrasts in ordered zones are rejected. Do not take abs(delta), clip
velocities, retry invalid draws, or omit a truncation normalizer.

The log acceptance ratio is

```text
logalpha = logL_new-logL_old + logP_new-logP_old
           + log(q_reverse)-log(q_forward) + log(abs(J))
```

This is the standard reversible-jump target/proposal/Jacobian construction;
see equations (2) and (4) in
[Hastie and Green's primary review](https://groups.seas.harvard.edu/courses/cs281/papers/hastie-green-2011.pdf).
The unit Jacobian and selection densities above are derived for this
layered transformation, not asserted for the reference spline kernel.

With DR disabled, within-dimension Gaussian walks use reflection at fixed
conditional bounds. With DR enabled they use unbounded Gaussian draws and
prior rejection, so the delayed-rejection density ratio is analytic. The bounds for a selected velocity depend only on its unchanged
neighbors; those for an interface depend on its unchanged adjacent
interfaces. These proposals are symmetric. Global physical conditions
are enforced by rejection. At count limits, or when no split is possible,
the attempted move is a self-transition. Family probabilities are never
renormalized based on which moves happen to be feasible.

For the current full two-stage layered proposal, including birth/death in
both stages, see [the layered DR derivation](THBI_LAYERED_DR.md). The earlier
fixed-scalar kernel remains selectable as a comparison mode.

## Storage and planned checks

Every iteration saves the current state, including rejected repeats.
Fraction and velocity vectors are cells because dimension changes.
Summaries evaluate each physical model on a common depth grid; layer
index is not averaged across changing counts. Interface summaries are
expected numbers of interfaces per depth bin, which can exceed one, not
probabilities of a unique interface.

A separate, fixed iteration stride saves cached responses after burn-in,
including rejected repeats. Response quantiles are computed from these
draws, not by evaluating a mean layer list. They describe uncertainty in
model responses, without an additional observational-noise realization.
Velocity mixing is checked at every fixed output depth; this remains
well-defined when K changes.

The MATLAB tests cover inverse transformations, reciprocal full acceptance
ratios, a finite-difference Jacobian, count-limit self-transitions, ordered
velocity normalization, correlated likelihood, noise-only prediction reuse,
variable-size history, sampled-response alignment and exact checkpoint
continuation. Absolute proposal-density checks also use unequal birth/death
probabilities and a zone where only two of three layers can split, so
reciprocity alone cannot conceal a shared missing factor. Failure-injection
checks verify exception propagation and diagnostic preservation. The optional
prior-only experiment disables cross-zone conditions and compares sampled
counts with analytic 1/K masses; it is a stochastic diagnostic whose
mixing must be inspected. The optional native test exercises all five
solvers through a split and inverse merge in reused scratch directories.

The user completed the earlier five-block two-chain synthetic; recovery and
mixing were not established. The new extended60/DR numerical checks remain
unexecuted by the agent. No claim of correct runtime behavior, posterior mixing
or synthetic recovery follows from source inspection alone.
