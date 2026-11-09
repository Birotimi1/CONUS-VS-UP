# Synthetic transdimensional model plots

Run `DATA/plot_ecm1_thin_thbi_colored_chains.m` in MATLAB. The script reads the
completed run referenced by `DATA/transdimensional/LATEST_RUN.txt`; set its
`result_dir` to an absolute run folder to select another experiment. It uses
the saved `par`, `chains`, `TRUEmodel` and data from `results.mat`, so the saved
burn-in applies automatically. It does not require `postprocessing.mat`.

The layout adapts our `plot_model_unified_colored_chains.m` and
`plot_model_unified_all_2.m` from the shared forward project's `Figure_codes`
folder. Chain 1 is blue and chain 2 red. The main-body code exposes model
selection, reconstruction, statistics and plotting for editing/debugging.
The originals are unchanged and are not MATLAB dependencies.

## Selection and transdimensional representation

- Thin colored curves show **every accepted proposal after burn-in** by
  default. `include_noise_accepts=true` includes accepted noise-only moves,
  which repeat an existing Vs model. Set it to `false` to show only accepted
  structural moves.
- `max_models_per_chain=Inf` shows all selected models. A finite display limit
  uses deterministic, evenly spaced selections; it changes only the drawn
  ensemble. The selected iteration IDs are saved in `plot_summary.mat`.
- Each model is reconstructed from its own zone fractions and velocities.
  Ensemble, true, initial, best and last models use **exact step coordinates**,
  including a half-space extension. No spline or common layer-index average
  is introduced.
- Mean/median curves, 68%/95% intervals, count probabilities and interface
  summaries use **all retained states**, including rejected repeats. They
  have equal sample weights; there is no extra likelihood weighting or
  likelihood-percentile filtering. Depth summaries are evaluated on the
  saved common-depth grid augmented by the true interfaces.
- “Best retained” is the actual post-burn-in model with minimum `sum(E2)`,
  selected separately for each chain. It is not labeled MAP. “Last state”
  means the state at the final iteration; it can repeat an earlier accepted
  state. Pointwise means and medians need not describe a sampled layering.
- Interface panels report the expected number of physical interfaces **per
  km**, accounting for nonuniform bin widths. They exclude the fixed
  half-space top and include all inferred interfaces without a contrast cut.

## Outputs

Each selected run receives a `Figures_THBI_Chains` subfolder. PNG (200 dpi)
and raster PDF (300 dpi) exports are made for these figures:

| Output stem | Contents |
|---|---|
| `accepted_vs_full` | Accepted physical models by chain, best retained Earths, truth and interface occurrence over the full depth range |
| `accepted_vs_200m`, `accepted_vs_500m`, `accepted_vs_5km` | Shallow views with the same selection and chain colors |
| `accepted_vs_1km` | View to 1.2 km so recovered discontinuities slightly deeper than the true 1 km interface remain visible |
| `accepted_vs_deep` | Crust/mantle view from 20 km to the output depth limit |
| `final_vs_by_chain` | One row per chain: start, last state, best retained Earth, median and 95% interval at full depth, 200 m and 5 km |
| `retained_vs_summary` | Pooled pointwise mean/median, 68%/95% intervals, best retained Earth and truth |
| `layer_counts_by_chain` | Per-chain count probabilities and full count histories for sediment/crust/mantle |
| `saved_five_observables_by_chain` | Optional five-observable overlays from cached retained responses, with per-chain medians and synthetic truth/data |

The response figure uses the saved `prediction_iterations`, not accepted-model
row positions. The 10k test saved 600 responses per chain; not every accepted
model has a cached prediction. At most 150 curves per chain are displayed by
default, while each response median uses every available post-burn-in saved
draw. Change `max_responses_per_chain` or disable `plot_saved_responses` at the
top of the script. No native solver is called.

`plot_summary.mat` saves the exact selected iterations, physical best/last
models, depth summaries, interface/count distributions and displayed response
iteration IDs. `plot_notes.txt` records selection definitions, counts, native
failure totals and the existing split R-hat diagnostics. Enable
`save_matlab_figures` to save editable `.fig` files as well.

The completed run inspected on 2026-11-06 is
`ECM1_thin_THBI_2chains_10k_burn4k_20261105_003201`. Its saved report records
three native proposal failures in chain 1 and none in chain 2. It also shows
large split R-hat values for several quantities (sediment-base depth about
6.69 and Moho depth about 5.93). The plots therefore expose chain differences;
pooled intervals must not be taken as evidence of convergence.

Only static source analysis and read-only inspection of saved result metadata
were performed while writing this script. MATLAB execution and visual
verification of these new figures remain for the user.
