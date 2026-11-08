# Current delayed rejection: random second move from the rejected Earth

Updated 2026-11-08 after the user's explicit clarification to follow the
paper's random second perturbation. This supersedes the earlier default
same-zone targeting choice.

`THBI/parms/bayes_inv_parms.m` now sets both `extended60` and `extended100` to:

```matlab
par.proposal.dr.enabled = true;
par.proposal.dr.mode = 'layered_path';
par.proposal.dr.same_zone_weight = 0;
```

After a valid first candidate y is rejected from accepted Earth x, draw z
FROM y with the same ordinary move/zone mixture as stage one. No first-zone
preference or forced daughter selection is applied. Birth/death, velocity,
interface, boundary and noise moves remain available. Invalid first proposals
and native failures are not sent into the second stage.

This follows the proposal description in Brennan's supporting information
and the local `functions/master_subfuncs/delay_reject.m` control flow. The
complete layered acceptance correction is retained, with both first-stage
rejection probabilities explicit:

```text
log R2 = log pi(z) - log pi(x) + H_xy + H_yz
       + log[1-alpha1(z,y)] - log[1-alpha1(x,y)]
H_ab = log q1(a|b) - log q1(b|a) + log|J_ab|
```

Both alpha1 values use `b6_IFACCEPT(...,false)`. The physical split/merge maps
and their proposal/Jacobian terms are unchanged. There is no targeting
correction in this mode, since Q2(z|x,y)=Q1(z|y). The more general derivation
is in [THBI_TARGETED_DR.md](THBI_TARGETED_DR.md); the targeting machinery now
serves only optional comparison tests.

## MATLAB scripts

- [check_thbi_layered_dr.m](../DATA/check_thbi_layered_dr.m): current solver-free
  entry point; physical path balance, exact random-mode stationarity,
  continuous-prior ensemble, 60/100-s profiles, RF filters, cache/counts and
  checkpoint restart. MATLAB also runs result/exception checks.
- [check_thbi_extended100.m](../DATA/check_thbi_extended100.m): native forward
  preflight and truth/observable figures, without production inversion.
- [run_ecm1_thin_thbi_100s_ps_sp.m](../DATA/run_ecm1_thin_thbi_100s_ps_sp.m): the
  same 250-km, 100-s test, now using random
  second moves. The launcher asserts this setting before running.

The known Earth, measurements and RF bands are unchanged. The latest local
launcher override is 15,000 iterations and 8,000 burn-in per chain. Profile
defaults remain 10,000/4,000. The output run label and printed counts are
derived from the actual launcher settings. Start a new result directory; checkpoint comparison detects
the changed configuration/source. Existing results are not overwritten.

The optional `check_thbi_targeted_dr.m` now explicitly requests its targeted
comparison settings instead of relying on the production default. Historical
30k/150k targeted prior-count results retain their original labels.


## Executed validation after this configuration change

`DATA/check_thbi_layered_dr.m` completed with exit code 0 in Octave 11.3.0:
1,080 physical proposal paths; exact random-mode finite-state stationarity;
60-s and 100-s chain likelihood/cache/count/restart checks; 250-km truth and
100-s RF filter checks; and 6,000 independent continuous-prior transitions.
The prior-ensemble maximum paired drift was 1.976 standard errors (threshold
4.5), with 3,969 first-stage and 17 second-stage acceptances.

Code Analyzer completed on all 11 modified MATLAB files with no syntax
errors (two optional-branch and three obsolete-suppression messages remain).
See [numerical output](THBI_RANDOM_DR_NUMERICAL_LOG.txt),
[validation manifest](THBI_RANDOM_DR_VALIDATION.json), and
[analyzer output](THBI_RANDOM_DR_MLINT.txt).

No native CPS/HVf calculation or production inversion was launched. MATLAB
summary/exception/plot checks remain user-run; prior license checkout failed.
These checks establish implementation behavior, not geological recovery or
convergence of the requested two 10k chains.
