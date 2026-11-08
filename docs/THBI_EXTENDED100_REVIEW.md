# Thin synthetic model to 250 km, with periods through 100 seconds

Updated 2026-11-08. The current profile is `extended100` in
[bayes_inv_parms.m](../THBI/parms/bayes_inv_parms.m). The earlier `extended60`
profile and old result directories remain available. Changes retain explicit
calculation blocks in the numbered stages and separate DATA launch scripts.

## What was checked in the reference implementation

Sources: the local `/Users/birotimi/Downloads/THBI_MCMC` checkout, especially
`b3_FORWARD_MODEL_BW.m`, `b3_FORWARD_MODEL_RF_ccp.m`, `b4_CALC_MISFIT.m`,
`b5_CALC_LIKELIHOOD.m`, `example_github/parms/bayes_inv_parms.m` and
`functions/master_subfuncs/delay_reject.m`; and Brunsvik, Eilon & Lynner (2024),
[10.1029/2024GL109041](https://doi.org/10.1029/2024GL109041), supporting
information S1–S2, PDF pages 3–10. The local supporting PDF is
`/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Papers/Brennan 2024 supp.pdf`.
Pages describing the DR update, Sp processing and Ps/forward modelling were
rendered and inspected as well as extracted as text.

| Dataset | Brennan paper / local reference | Our synthetic implementation |
|---|---|---|
| Ps | The 2024 study uses an H–kappa–xi energy volume from reprocessed EARS RFs. The error surface is 1-E/Emax+0.002; Ps/PpPs/PpSs stacking weights are .5/.3/.2. This is not direct Ps waveform least squares. The reference code also retains older waveform branches. | CPS `hrftn96 -P` signed R/Z waveform, Gaussian and explicit period filtering. Fixed correlated Gaussian residual likelihood. No H–kappa or anisotropy energy volume. |
| Sp | The observed CCP stack used waveforms filtered over 4–100 s. Reference forward modelling uses a propagator matrix (local wrapper also permits telewavesim), rotates ZR into P–SV, aligns the parent and migrates predictions to depth. Example configuration compares depths 20–250 km and downweights the shallow part of the profile. It also substitutes a simplified parent pulse to suppress parent sediment reverberations. | CPS incident-S R and Z components; regularized signed Z/R deconvolution. Negative precursor lags, no time reversal, polarity flip or peak normalization. Direct time-domain residuals for one ray parameter, with the same processing on truth and candidate. No CCP migration, depth weighting, P–SV rotation or substituted parent pulse. |
| Phase velocity | Multiple Rayleigh and Love sources have separate hierarchical errors. SI Table S1 includes 5–40/6–40 s ambient, 25–180 s earthquake and additional regional ranges. MINEOS/reference kernels reduce forward cost. | Exact CPS `surf96` fundamental Rayleigh calculation at each structural proposal; one synthetic phase block. No Love dataset or kernel approximation. |
| Ellipticity | Rayleigh H/V is a surface-wave polarization observable, distinct from diffuse-field HVSR. SI Table S1 lists 20–90 s. Reference calculations use the Tanimoto approach and perturbation kernels. | HVf fundamental Rayleigh H/V only: one Rayleigh mode, no Love or body contribution. Exact calculation at each structural proposal. |
| Group velocity | The generic local suite has a group-velocity branch; it is not among the main 2024 study datasets in SI Table S1. | Exact CPS fundamental Rayleigh group velocity; its own hierarchical block. |
| HVSR | Not a diffuse-field HVSR constraint in the 2024 setup. | HVf diffuse-field H/V, including the configured surface-wave modes and body-wave contribution. Separate from ellipticity. |

The reference b5 uses dataset sigmas and configured degrees of freedom,
including an average over records of a datatype. Our b5 uses the normalized
Gaussian likelihood with C=gamma^2*C0, the actual retained sample count and
fixed covariance shape. Synthetic noise is drawn from this same C0. These
likelihoods are not interchangeable, especially for an H–kappa energy surface
or a CCP profile. No reference weighting or RF normalization was copied
without its corresponding data definition.

## Current coverage and model

| Quantity | extended100 setting |
|---|---|
| Variable model base / half-space top | 250 km; output Vs grid to 270 km |
| Thin true interfaces | 0.010, 0.050, 0.130 and 1.000 km preserved |
| True sediment-zone base / Moho | 0.130 km / 40 km |
| True zone layer counts | [3 4 3], plus a separate half-space |
| Deep true material | Existing 4.54 km/s material continued from 70 to 250 km and below |
| HVSR | 0.01–25 Hz (100–0.04 s); all previous frequency samples retained |
| Ellipticity / group velocity | 3–100 s; all previous samples retained |
| Phase velocity | 9–100 s; all previous samples retained |
| Ps full period band | 1.25–100 s |
| Sp full period band | 3–100 s |
| RF diagnostic bands | Previous three bands plus 60–100 s; diagnostic only |
| Native RF record | 4,096 samples at 0.2 s = 819.2 s; delay 300 s |
| Retained RF cadence | 0.2 s, unchanged from extended60's 0.1-s native / stride 2 |
| RF likelihood lag windows | Ps -100 to +120 s; Sp -120 to +100 s |
| Chains | Current 100-s launcher: two serial chains, 15,000 iterations each, 8,000 burn-in; profile defaults remain 10,000/4,000 |
| DR | Random second proposal from the rejected model (`layered_path`), without zone targeting; complete reverse-path acceptance |

CPS's existing supported sample-count limit is retained. The longer duration
comes from dt=.2, with both full RF bands still below retained-data Nyquist.
Native preflight checks compare a doubled-duration Sp record and a smaller
water level, since valid sampling settings alone do not establish accuracy.

`a2_LOAD_DATA(par,true)` now returns only the actual known Earth and its
zoned representation, with `forward.status='not_computed'` and empty data.
The default `a2_LOAD_DATA(par)` still executes the full native synthetic
workflow. This permits testing the actual 250-km truth without a forward
solver or a duplicate model definition.

For that truth, travel-time sums to the model base are +26.229 s (Ps) and
-30.626 s (Sp), so the existing data windows cover primary conversions from
250 km. The extended period bands are frequency settings, not conversion
lags. The base is not an impedance contrast in the true model, so these
numbers are coverage checks, not predictions of a distinct base arrival.

The identical deep true segment and half-space are merged by numerical
layer refinement, as they should be. Candidate models can still vary all
the way to 250 km. The truth therefore tests extension into homogeneous
mantle; it does not test recovery of a newly introduced deep low-velocity
zone or interface. Period coverage alone cannot guarantee sensitivity or
resolution at a particular depth. Actual response changes, signal/noise,
kernels and posterior recovery must establish that. The current forward
geometry remains flat, isotropic and elastic.

## What to run in MATLAB

1. [check_thbi_layered_dr.m](../DATA/check_thbi_layered_dr.m): solver-free
   physical RJ, exact stationarity, independent prior ensemble, RF filtering
   and toy-chain/restart checks.
2. [check_thbi_extended100.m](../DATA/check_thbi_extended100.m): native RF
   homogeneous/bilayer controls, independent 250-km conversion controls,
   long-record/water-level comparisons, true six-block forward calculation,
   split/merge response checks, and saved truth/observable/RF-band figures.
   It stops before a production inversion.
3. [run_ecm1_thin_thbi_100s_ps_sp.m](../DATA/run_ecm1_thin_thbi_100s_ps_sp.m):
   repeats preflight, then runs the requested two chains, summarizes all
   retained states including repeats, and recomputes the best retained Earth.

Each native script creates a new timestamped folder beneath
`DATA/transdimensional/ECM1_thin_THBI_extended100_*`, and saves its complete
settings and launcher. Neither script resumes an older checkpoint.

## Validation status

Executed in Octave: split/merge algebra, 1,080 physical DR paths, exact
36-state stationarity with negative controls, the independent 6,000-Earth
prior-transition experiment, legacy/extended profiles' kernel-only chain
checks, and the 250-km/100-s truth/grid/filter tests. The new 60–100-s filter
preserves its expected 80-s response and suppresses 10-s energy. Signed Sp
lag/amplitude also passes at dt=.2.

MATLAB could not obtain a license (-96). Native CPS/HVf calculations,
MATLAB-specific result/exception/plot checks, and the new full inversion
have not been executed by the agent. In particular, the low-frequency HVf
modal/quadrature settings and all six real solver curves still need the
native preflight. No new scientific result or recovered 250-km structure is
claimed. See [THBI_TARGETED_DR.md](THBI_TARGETED_DR.md) for the acceptance
rule, numerical evidence and the failed short prior-count mixing checks.


The longer 150k controlled prior-count diagnostic passed (maximum error
0.0327 across zones); the full solver-free launcher exited successfully.
The two failed 30k checks remain documented. Exact commands/runtime, logs
and source hashes are in [THBI_TARGETED_VALIDATION.json](THBI_TARGETED_VALIDATION.json).


## User clarification: random second move

The user selected the paper's random second perturbation. Both extended
profiles now use `layered_path`, with `same_zone_weight=0`. The 250-km truth,
100-s periods, RF processing and chain lengths are unchanged. The historical
150k prior-count experiment above used the earlier targeted kernel; its
result is not relabelled as a test of this new configuration. See
[THBI_RANDOM_SECOND_MOVE.md](THBI_RANDOM_SECOND_MOVE.md) for current tests.
