# Next thin-model test: birth/death, delayed rejection, 60-s Ps/Sp

> Historical review. Current settings and executed tests are in [THBI_TARGETED_DR.md](THBI_TARGETED_DR.md) and [THBI_EXTENDED100_REVIEW.md](THBI_EXTENDED100_REVIEW.md).


**Subsequent update (2026-11-08):** extended60 now uses two full layered
proposal stages, including birth/death in stage two. See
[THBI_LAYERED_DR.md](THBI_LAYERED_DR.md). Fixed-scalar-only descriptions
and earlier static manifests below describe the previous implementation.


Subsequent review: [2026-11-08 workflow and fixes](THBI_WORKFLOW_AUDIT.md).
The earlier static record below describes the source at that review date.

Updated 2026-11-07. These edits have source/static checks only. The agent has
not run MATLAB, Octave, CPS, HVf, numerical tests or a new inversion for this
update. The new MATLAB checks below must pass on the user's installation;
passing them does not establish posterior convergence or recovery.

## What the reference comparison found

The local reference is `/Users/birotimi/Downloads/THBI_MCMC`, commit
`585ee81f55961f31cf018163eabdbf972f5d0652`, Zach Eilon/Brennan Brunsvik's
GPL-3.0 code. `b2_PERTURB_MODEL.m` lines 552-730 insert/delete **spline
knots** in crust and mantle, refit spline coefficients, and set active
`P_bd=k/(k+1)` or `k/(k-1)` factors. More elaborate proposal-density blocks
are commented out. That implementation is not a physical layer split and
its factors cannot simply be transplanted into our different coordinates.

Our `THBI/b2_PERTURB_MODEL.m` splits a finite layer into two daughters in
any of the three zones. It samples an interface fraction and velocity
contrast; the inverse death recovers their thickness-weighted parent.
It includes move/zone selection, eligible-layer selection, split width,
contrast density, normalized dimension-dependent priors and the unit
Jacobian. Existing tests check absolute densities as well as round trips.
This is an adaptation guided by THBI, not a claim that Zach's kernel was
copied verbatim. See [the target and split equations](TRANSDIMENSIONAL_DESIGN.md).

The old birth contrast was always N(0,0.08²) km/s. In `extended60`, it is a
fixed 70% narrow / 30% broad Gaussian mixture; broad standard deviations
are [0.4,0.3,0.2] km/s for sediment/crust/mantle. **Both birth and reverse
death use the full mixture density**, not the density of the drawn
component. Negative sediment contrasts and any out-of-prior proposal are
rejected; no clipping or redraw changes the proposal distribution.

The transdimensional engine previously had **no delayed rejection**.
The original HBI file inspected was
`/Users/birotimi/Downloads/hbi_all_codes/Codes/Mod5_HBI/functions_hbi_new/Functions/delayed_rejection_damp.m`.
Its second-try likelihood/tau test is not the full delayed-rejection ratio:
it omits the reverse first-stage proposal and rejection factors. The local
THBI `functions/master_subfuncs/delay_reject.m` also contains an unfinished
second-stage adjustment. Neither is imported as an exact posterior kernel.

## New delayed rejection

`functions/delayed_rejection_damp.m` implements the second-stage acceptance;
`THBI/run_one_chain.m` contains the actual proposal/forward/cache blocks.
After a rejected x->y, propose z from x on the **same scalar coordinate**,
with standard deviation 0.25 times the original. The selected family,
zone, coordinate and dimension stay fixed. With log pi = logL + logP:

```text
log r2 = log pi(z) - log pi(x)
       + log q1(y|z) - log q1(y|x)
       + log(1-alpha1(z,y)) - log(1-alpha1(x,y))
alpha2 = min(1, exp(log r2))
```

The smaller q2(z|x,y)=Normal(z;x,s2²) is symmetric and independent of y,
so its reverse/forward terms cancel. For this implementation q1 is an
unbounded Gaussian; outside-prior draws are rejected. Reflected proposals
are retained only in legacy profiles with DR disabled. DR is never applied
to a birth/death or after a first-stage native solver failure. A known
zero prior density is different from an unknown target caused by a failed
numerical solve. The new debug profile stops at the first native failure.

The acceptance is based on [Green and Mira (2001)](https://doi.org/10.1093/biomet/88.4.1035).
`test_thbi_extended60` compares forward/reverse **two-stage path fluxes**,
including both Gaussian proposals; it also checks outside-support and
certain-first-acceptance cases. `test_thbi_chain('extended60')` checks six
likelihood blocks, DR stage accounting, noise-only caching, saved responses
and exact checkpoint continuation with an inexpensive injected forward.

Saved histories now include `accepted_stage` (0/1/2), `dr_attempted`,
`dr_proposal_valid`, `dr_logalpha`, `first_proposal_valid`, and
`first_logalpha`. `proposal_valid` means at least one attempted candidate
passed structural support. Reports include retained DR attempts/accepts.

## Period coverage and RF convention

| Quantity | `extended60` settings |
|---|---|
| HVSR | Original 0.2-25 Hz grid |
| Rayleigh ellipticity | Original 3-25 s samples plus 26:2:60 s |
| Rayleigh phase velocity | Original 9-35 s samples plus 36:2:60 s |
| Rayleigh group velocity | Original 3-25 s samples plus 26:2:60 s |
| Ps full-band RF | 1.25-60 s; diagnostic bands 1.25-10, 10-30, 30-60 s |
| Sp full-band RF | 3-60 s; diagnostic bands 3-10, 10-30, 30-60 s |
| Ps settings | p=0.06 s/km; Gaussian alpha=2.5; selected lag -100 to +120 s |
| Sp settings | p=0.10 s/km; Gaussian alpha=1; selected lag -120 to +100 s |
| Native RF records | 4096 samples at 0.1 s; Ps delay 140 s, Sp delay 160 s |
| RF likelihood sampling | Every second native sample (0.2 s) |

These are **frequency/period bands**, not a relabeling of conversion-time
windows. Zero-phase fourth-order Butterworth amplitude factors add explicit
low/high corners; they roll off rather than impose hard spectral cutoffs.
The native Gaussian exp[-(pi*f/alpha)²] remains part of the response.
The 409.6-s records and margins around the selected windows accommodate
the long-period filter. Duration/sampling checks reject incompatible edits.

The likelihood has six blocks: HVSR, ellipticity, phase, group, Ps and Sp
(five observable families). Each full-band RF enters **once** with its own
gamma. The three diagnostic bands are correlated views of that same
waveform, not three additional independent observations. Synthetic RF
noise is added after filtering using the declared fixed exponential
covariance (Ps 0.4-s, Sp 1-s correlation); it is a controlled noise model,
not a claim about the covariance of an observed/deconvolved field RF.

Ps uses CPS `hrftn96 -P`, signed R/Z. Sp uses CPS `hrftn96 -S -r` and
`-S -z` component waveforms and then regularized spectral Z/R:

```text
RF_Sp(f) = Z(f)*conj(R(f)) / max(|R(f)|², 0.01*max(|R|²))
```

The component Gaussian alpha is 10, an edge taper is 5%, and the target
Gaussian alpha is 1. Sp conversions remain at **negative lag**, with CPS's
upward-positive Z convention. There is no time reversal, polarity flip
or peak normalization. All processing is identical for truth and trial
models. R/Z or Z/R is a transfer response; neither is a measured RF here.

CPS's [official change history](https://www.eas.slu.edu/eqc/eqc_cps/cpsbug.html)
records instability of direct layered S-receiver-function division (2011).
The official [CPS source distribution](https://rbherrmann.github.io/ComputerProgramsSeismology/)
was inspected (`NP330.Oct-12-2026.tgz`, `PROGRAMS.330/VOLVI/src/hrftn96.f`)
to verify component flags, signs, Gaussian and sample conventions. The
installed binary was not replaced or executed by the agent. Direct `-S`
division is used only for the **homogeneous control** in the manual native
test; the production layered Sp path uses component deconvolution.

## Depth and shallow recovery

The thin cap and all original velocity contrasts remain unchanged. To let
the longer-period data constrain a variable deeper mantle, this profile
moves the fixed model base from 70 to 200 km. The true old 4.54-km/s half-space
is continued unchanged from 70 to 200 km and below: **no new true deep
contrast** is introduced. The zoned truth has K=[3,4,3] finite layers and
a separate identical bottom half-space. Numerical refinement removes
identical adjacent material, so this extension does not change the true
elastic Earth. The legacy profiles still use a 70-km model base.

The new crustal interface proposal scale is 0.004 of zone thickness
(previously 0.04); sediment-base steps are 0.005 km (previously 0.03 km).
These are fixed proposal settings, not bounds centered on the true answers.
No truth is supplied to the chain. The broad birth mixture helps propose
larger contrasts, but none of these changes guarantees recovery.

The prior completed run's best retained interfaces near 1 km were about
1.071 km and 1.013 km. The clear error in the supplied **upper-200-m** panel
is the true 130-m interface shifted to about 166/159 m, with excess shallow
layers. Strong between-chain differences and high split R-hat mean those
chains do not establish recovery. Longer periods primarily add deeper
constraints; they do not by themselves cure shallow resolution or mixing.

## Run in MATLAB

From the repository root, first run:

```matlab
run('DATA/check_thbi_extended60.m')
```

It runs algebraic/toy checks, independent native Ps/Sp controls, the actual
thin truth on the extended grids, and native split/merge/round-trip checks
for all six blocks. It stops before production MCMC. Inspect
`true_observables.png`, `true_rf_bands.png`, and the saved native check MAT/SAC
files under the printed `DATA/transdimensional/...extended60_preflight.../`
directory. Native failures stop the script and retain diagnostic files.

Then run:

```matlab
run('DATA/run_ecm1_thin_thbi_60s_ps_sp.m')
```

This repeats the preflight and runs two independent 10,000-iteration chains
with 4,000 burn-in. It retains all 12,000 post-burn-in states including
repeats, saves six prediction histories and creates joint-fit, RF-band,
shallow structure, layer-count, noise, DR and mixing reports/figures.
It uses a new timestamped directory. The old results remain intact;
**old checkpoints cannot be resumed with these source/configuration edits**.

Native timing/amplitude controls and sensitivity to record length and water
level are acceptance gates, not tests to bypass. A failure needs inspection
of saved records and settings before sampling. All new numerical checks
remain pending until the user runs them in MATLAB.

## Static verification record

Standalone MATLAB R2022a Code Analyzer inspected 113 project MATLAB files,
including the legacy HBI, new THBI, 47 shared helpers, DATA launchers and tests.
There were no syntax errors. Remaining messages concern historical optional
branches, formatting/performance suggestions and redundant suppressions.
`THBI_EXTENDED60_MLINT.txt` records the diagnostics;
`THBI_EXTENDED60_STATIC.json` records the file inventory, hashes, reference
source hashes and the explicit unexecuted numerical-test status. Call-inventory
inspection confirmed the new DR/RF stages refer to packaged helpers and
MATLAB functions. This does not replace MATLAB/native execution.
