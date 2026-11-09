# THBI function-suite comparison

Current extension: [2026-11-07 Ps/Sp, 60-s bands and delayed-rejection review](THBI_NEXT_RUN_REVIEW.md).
The audit below describes the earlier five-block profile; its exclusions and counts are historical.

Reference: `/Users/birotimi/Downloads/THBI_MCMC`, commit `585ee81f55961f31cf018163eabdbf972f5d0652`. Audit: 2026-11-04.

The previous 30-helper consolidation checked the dependencies of the implemented pilot. It did **not** establish feature parity with the full reference suite. This comparison accounts for all **127 MATLAB files** under the original `functions/`, including subfolders, and adds the applicable retained-model/posterior diagnostics.

Every reference file was inventoried by file/interface and workflow role. The original implementations of adapted helpers and c2-c5 stages were read in full. Excluded specialty files have not all received full line-by-line validation. This is an applicability audit, not a claim that all reference code is debugged.

The project now has **42 helper files**, with plotting helpers in `functions/plotting/`. New main-body stages are `THBI/c2_BUILD_POSTERIOR`, `c3_BUILD_MODEL_SUITE`, `c4_FINAL_MODEL`, and `c5_FINAL_FORWARD_MODEL`. No original source checkout was modified. No MATLAB or native solver was executed.

## How to run the added suite

The existing `DATA/run_ecm1_thin_thbi.m` now builds/saves the posterior, material-profile suite and final summary after sampling, then produces the extended plots. For a completed run, use `DATA/postprocess_ecm1_thin_thbi.m`. Its `run_final_forward=false` default avoids native calculations; setting it true recomputes all five responses of the explicitly selected retained model. `DATA/plot_ecm1_thin_thbi.m` reloads saved postprocessing outputs when available.

`postprocessing.mat` uses MATLAB v7.3 for the full Vs/Vp/density suite; storage scales with retained iterations times depth points. All retained states, including rejected repeats, contribute. Display subsets are deterministic and do not alter sampler RNG or statistical weights.

The pointwise median profile is a summary. The separately identified representative Earth is the retained state with minimum base-covariance residual energy (first occurrence breaks ties); it is not a cross-dimension MAP. c5 forwards that physical model, rather than inventing interfaces in a pointwise median.

## Reference issues handled explicitly

- `count_N_layers_in_mod` layerises spline profiles and never assigns its declared output in this checkout; the new API counts physical finite layers, zone counts and the half-space explicitly.
- `calc_K_in_model` calculates fractions of sensitivity kernels, not a number of layers.
- `make_mod_from_parms_no_spline` still accesses spline weights/matrices and declares the same function name as the spline implementation.
- The original collation helper returns early for a single chain before removing burn-in; the adapted helper applies the same retained-iteration rule to every chain.
- The original heatmap uses natural `log` while labeling the bar `log10`; the adapted plot uses `log10` consistently.
- The original `q_v_given_V` has a positive squared exponent. The current birth density is explicitly written with the negative squared exponent in b2; that helper is not imported.

## Dispositions

`Adapted` means an integrated layered counterpart was added. `Covered` means the needed role is already in a numbered stage/helper, with differences stated. `Different formulation` means importing the original would change the selected model/likelihood. `Not used` covers the spline/kernel/legacy route. `Future scope` is an actual unimplemented extension, not hidden completion.

**Future scope: 35**, **Not used: 30**, **Different formulation: 17**, **Covered: 33**, **Adapted: 12**.

## File-by-file inventory

| Reference file under functions/ | Disposition | Project implementation / reason |
|---|---|---|
| `BW_2_RF.m` | Future scope | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Observed/raw BW-to-RF deconvolution and processing require a separate data pipeline. Current CPS RF is already a Gaussian-filtered signed R/Z transfer response. |
| `T_conv_max_Ps.m` | Future scope | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Observed/raw BW-to-RF deconvolution and processing require a separate data pipeline. Current CPS RF is already a Gaussian-filtered signed R/Z transfer response. |
| `back_compatible.m` | Not used | Legacy original data/misfit schema conversion; no claim that original THBI saved runs load into the new layered sampler. |
| `calc_Ce.m` | Different formulation | [THBI/a2_PREPARE_DATA.m](../THBI/a2_PREPARE_DATA.m). Fixed within-data covariance and Cholesky whitening; five independent data types. No inverse/determinant arithmetic or station-pair covariance imported. |
| `calc_K_in_model.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `calc_Pm_prior.m` | Covered | [THBI/a3_CALC_LOG_PRIOR.m](../THBI/a3_CALC_LOG_PRIOR.m). Explicit normalized count/interface/velocity/boundary/noise prior in log coordinates. |
| `calc_RF_misfit.m` | Covered | [THBI/b4_CALC_MISFIT.m](../THBI/b4_CALC_MISFIT.m). Whitened quadratic residuals for every one of the five observables; separate wrappers would duplicate this visible stage. |
| `calc_SW_misfit.m` | Covered | [THBI/b4_CALC_MISFIT.m](../THBI/b4_CALC_MISFIT.m). Whitened quadratic residuals for every one of the five observables; separate wrappers would duplicate this visible stage. |
| `calc_Vperturbation.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `check_data.m` | Covered | [THBI/a2_PREPARE_DATA.m](../THBI/a2_PREPARE_DATA.m). Check required data, dimensions, finite values, errors and covariance before likelihood evaluation. |
| `clip_main_arrival.m` | Future scope | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Observed/raw BW-to-RF deconvolution and processing require a separate data pipeline. Current CPS RF is already a Gaussian-filtered signed R/Z transfer response. |
| `collate_allmodels_perchain.m` | Adapted | [functions/collate_allmodels_perchain.m](../functions/collate_allmodels_perchain.m). Collate all retained column/cell states, preserving repeats and chain/iteration identity; burn-in also removed for one chain. |
| `compare_kernel_full_calc.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `compare_kernel_full_calc_for_testing_2.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `compare_kernel_full_calc_for_testing_3.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `count_N_layers_in_mod.m` | Adapted | [functions/count_N_layers_in_mod.m](../functions/count_N_layers_in_mod.m). Exact finite/zone/half-space counts for one layered model. Original estimates numerical layerisation counts and never assigns avNmod; signatures intentionally differ. |
| `dofTHBI.m` | Different formulation | [THBI/b5_CALC_LIKELIHOOD.m](../THBI/b5_CALC_LIKELIHOOD.m). Normalized full-covariance Gaussian uses the observed vector dimension N. Replacing N with heuristic effective DOF would change the current target. |
| `dofTHBI_SW_low.m` | Different formulation | [THBI/b5_CALC_LIKELIHOOD.m](../THBI/b5_CALC_LIKELIHOOD.m). Normalized full-covariance Gaussian uses the observed vector dimension N. Replacing N with heuristic effective DOF would change the current target. |
| `dtype_code.m` | Different formulation | [THBI/parms/bayes_inv_parms.m](../THBI/parms/bayes_inv_parms.m). Five explicit data keys hvsr/ellip/cph/ugr/rf; original parser does not implement our HVSR key and carries BW/Sp/CCP/HK variants. |
| `duplicate_data_distribute.m` | Future scope | Multiple waveform traces/processing bands and trace-specific errors. Current synthetic has one dataset per observable and five noise multipliers. |
| `eqcluster.m` | Future scope | Observed event clustering, polarity estimation or multiple-ray-parameter grouping; current synthetic RF has one specified ray parameter in s/km. |
| `expand_dathparms_to_data.m` | Future scope | Multiple waveform traces/processing bands and trace-specific errors. Current synthetic has one dataset per observable and five noise multipliers. |
| `functionsExternal/zivlab_CellReg/gaussfit.m` | Not used | [THBI/c4_FINAL_MODEL.m](../THBI/c4_FINAL_MODEL.m). Gaussian fitting dependency not required; empirical quantiles preserve non-Gaussian sampled distributions. |
| `getPaths.m` | Covered | [THBI/a0_STARTUP_BAYES.m](../THBI/a0_STARTUP_BAYES.m). Explicit returned paths and par.forc; no generated pathsAutoGen script or source-checkout path dependency. |
| `get_freqs.m` | Covered | [THBI/a2_LOAD_DATA.m](../THBI/a2_LOAD_DATA.m). Explicit synthetic frequency/period arrays. Original reads a particular gridvalue filename convention using a hard-coded default path. |
| `get_sigma_min_prior.m` | Different formulation | [THBI/parms/bayes_inv_parms.m](../THBI/parms/bayes_inv_parms.m). Explicit synthetic base uncertainties plus bounded log-noise priors; no data-derived prior estimation in this controlled experiment. |
| `hKappaFunk/HKstack_anis_wrapper.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/hKappaError.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/hk_forward_model.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/hk_maximum_possible_value.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/plot_HK_stack.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/replace_with_hk_package/HK_stack_anis_wrapper_no_grid.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/replace_with_hk_package/HKstack_anis_wrapper_fast.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/replace_with_hk_package/HKstack_no_grid.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/replace_with_hk_package/christof_radial_anis.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `hKappaFunk/replace_with_hk_package/hk_christof.m` | Future scope | HK/anisotropic stacking is a different observable/model extension; not part of the five-observable thin test. |
| `handle_errors/assess_timeout_results.m` | Covered | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Explicit command status, fixed-grid/finite output checks and cause-preserving exceptions; chain records failures and preserves native diagnostics. No timeout wrapper is enabled. |
| `ifaccept.m` | Covered | [THBI/b6_IFACCEPT.m](../THBI/b6_IFACCEPT.m). Full reversible-jump log acceptance includes likelihood, prior, forward/reverse proposals and Jacobian. |
| `ifforwardfail.m` | Covered | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Explicit command status, fixed-grid/finite output checks and cause-preserving exceptions; chain records failures and preserves native diagnostics. No timeout wrapper is enabled. |
| `initiate_Kbase.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `kchain_addcount.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `layerise_test_sediment_sp_ccp.m` | Not used | [functions/refine_layer_model.m](../functions/refine_layer_model.m). Experimental sediment-removing Sp/CCP layerisation would violate the identical five-observable Earth. Our refinement preserves material interfaces. |
| `make_SW_weight.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `make_allkernels.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `make_mod_from_parms.m` | Different formulation | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Piecewise-constant zoned layers and HBI material equations. Even the original no_spline file still references spline matrices/weights. |
| `make_mod_from_parms_no_spline.m` | Different formulation | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Piecewise-constant zoned layers and HBI material equations. Even the original no_spline file still references spline matrices/weights. |
| `make_ram_drive.m` | Covered | [THBI/run_one_chain.m](../THBI/run_one_chain.m). Per-chain temporary scratch and direct persistent checkpoints/failure copies; no RAM-disk creation or broad shell cleanup required. |
| `make_splines.m` | Not used | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Spline basis/weight operations do not apply to physical piecewise-constant layers. |
| `make_splines_fraction.m` | Not used | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Spline basis/weight operations do not apply to physical piecewise-constant layers. |
| `mantle_vs2rho.m` | Different formulation | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Preserve the user HBI Vp/density equations and floor order in both truth and inversion. Original zone-specific/AK135 conversions would change the experiment. |
| `mantle_vs2vp.m` | Different formulation | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Preserve the user HBI Vp/density equations and floor order in both truth and inversion. Original zone-specific/AK135 conversions would change the experiment. |
| `master_subfuncs/clean_ram.m` | Covered | [THBI/run_one_chain.m](../THBI/run_one_chain.m). Per-chain temporary scratch and direct persistent checkpoints/failure copies; no RAM-disk creation or broad shell cleanup required. |
| `master_subfuncs/delay_reject.m` | Not used | [THBI/b6_IFACCEPT.m](../THBI/b6_IFACCEPT.m). Single-stage RJ is the current validated-target design. A second-stage kernel-dependent rejection scheme cannot be imported without its own acceptance derivation and tests. |
| `master_subfuncs/forward_model.m` | Covered | [THBI/run_one_chain.m](../THBI/run_one_chain.m). Visible forward/cache/failure block calls five-observable b3 only for structural proposals. |
| `master_subfuncs/initiate_model.m` | Covered | [THBI/b1_INITIATE_MODEL.m](../THBI/b1_INITIATE_MODEL.m). Explicit reference/prior initialization and bounded retries; full native initial forward in run_one_chain. |
| `master_subfuncs/make_accept_info.m` | Covered | [THBI/b9_SAVE_RESULT.m](../THBI/b9_SAVE_RESULT.m). Proposal type, zone, acceptance, prior validity, forward failure and log-ratio diagnostics already stored per iteration. |
| `master_subfuncs/perturb_model.m` | Covered | [THBI/b2_PERTURB_MODEL.m](../THBI/b2_PERTURB_MODEL.m). All layered proposals, support tests and RJ ratios remain visible in b2 and the chain loop. |
| `master_subfuncs/ram_to_HD.m` | Covered | [THBI/run_one_chain.m](../THBI/run_one_chain.m). Per-chain temporary scratch and direct persistent checkpoints/failure copies; no RAM-disk creation or broad shell cleanup required. |
| `misc/for_seizmo/getwords.m` | Not used | [functions/sample_layer_model.m](../functions/sample_layer_model.m). Discontinuous spline/polynomial interpolation package not needed for exact step-layer evaluation. |
| `misc/for_seizmo/interpdc1.m` | Not used | [functions/sample_layer_model.m](../functions/sample_layer_model.m). Discontinuous spline/polynomial interpolation package not needed for exact step-layer evaluation. |
| `misc/for_seizmo/isequalsizeorscalar.m` | Not used | [functions/sample_layer_model.m](../functions/sample_layer_model.m). Discontinuous spline/polynomial interpolation package not needed for exact step-layer evaluation. |
| `misc/for_seizmo/ppdcval.m` | Not used | [functions/sample_layer_model.m](../functions/sample_layer_model.m). Discontinuous spline/polynomial interpolation package not needed for exact step-layer evaluation. |
| `mkchainstr.m` | Covered | [THBI/run_one_chain.m](../THBI/run_one_chain.m). Unique chain_XX output directory and local scratch directory; letter-based chain names not required. |
| `model_NVG_info.m` | Future scope | Smooth-profile Moho gradients/LVZ width rules require a separately defined interpretation for stepwise layers; no silent gradient approximation added. |
| `moho_adj_gradients.m` | Future scope | Smooth-profile Moho gradients/LVZ width rules require a separately defined interpretation for stepwise layers; no silent gradient approximation added. |
| `pairwise_Ce.m` | Different formulation | [THBI/a2_PREPARE_DATA.m](../THBI/a2_PREPARE_DATA.m). Fixed within-data covariance and Cholesky whitening; five independent data types. No inverse/determinant arithmetic or station-pair covariance imported. |
| `par_equiv_prior.m` | Covered | [THBI/run_one_chain.m](../THBI/run_one_chain.m). Exact configuration/data/source comparison on resume; conditioned empirical prior regenerated for each new run. |
| `parse_dtype.m` | Different formulation | [THBI/parms/bayes_inv_parms.m](../THBI/parms/bayes_inv_parms.m). Five explicit data keys hvsr/ellip/cph/ugr/rf; original parser does not implement our HVSR key and carries BW/Sp/CCP/HK variants. |
| `parse_dtype_all.m` | Different formulation | [THBI/parms/bayes_inv_parms.m](../THBI/parms/bayes_inv_parms.m). Five explicit data keys hvsr/ellip/cph/ugr/rf; original parser does not implement our HVSR key and carries BW/Sp/CCP/HK variants. |
| `plotting/dist2line_geog.m` | Future scope | Station mapping/cross-section analysis follows multi-station inversions; not required for the single thin synthetic test. |
| `plotting/plot_DATAFITS.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Five observable fits, noisy/clean truth, initial/best model and saved-response intervals. |
| `plotting/plot_FIG1_FIT_DATA.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Five observable fits, noisy/clean truth, initial/best model and saved-response intervals. |
| `plotting/plot_FIG2_FIT_MODEL.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Layered truth/start/best profiles, shallow views, retained envelopes and counts; c4 adds material summary outputs. Original optional regional/reference-model overlays not ported. |
| `plotting/plot_FINAL_MODEL.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Layered truth/start/best profiles, shallow views, retained envelopes and counts; c4 adds material summary outputs. Original optional regional/reference-model overlays not ported. |
| `plotting/plot_HEATMAP_ALLMODELS.m` | Adapted | [functions/plotting/plot_HEATMAP_ALLMODELS.m](../functions/plotting/plot_HEATMAP_ALLMODELS.m). All retained Vs/Vp histogram masses, actual log10 labels and nonuniform physical depth grid. |
| `plotting/plot_HEATMAP_ALLMODELS_shallow.m` | Adapted | [functions/plotting/plot_HEATMAP_ALLMODELS.m](../functions/plotting/plot_HEATMAP_ALLMODELS.m). Same helper with depth_limit=.3 km for the thin shallow structure. |
| `plotting/plot_HYPERPARAMETER_TRENDS.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). All-chain noise, log-likelihood and base residual energy histories with burn-in markers. |
| `plotting/plot_KNOT_TRENDS.m` | Adapted | [functions/plotting/plot_LAYER_TRENDS.m](../functions/plotting/plot_LAYER_TRENDS.m). Physical interface locations and layer counts, including rejected repeats. |
| `plotting/plot_MISFIT_TRENDS.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). All-chain noise, log-likelihood and base residual energy histories with burn-in markers. |
| `plotting/plot_MODEL_SUMMARY.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Layered truth/start/best profiles, shallow views, retained envelopes and counts; c4 adds material summary outputs. Original optional regional/reference-model overlays not ported. |
| `plotting/plot_MOD_TRUEvsTRIAL.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Layered truth/start/best profiles, shallow views, retained envelopes and counts; c4 adds material summary outputs. Original optional regional/reference-model overlays not ported. |
| `plotting/plot_P2P_recovery.m` | Future scope | Reference-specific recovery metrics need a separate calibrated synthetic recovery study; the new plots do not certify recovery or uncertainty coverage. |
| `plotting/plot_PARAMETERISATION.m` | Adapted | [functions/plotting/plot_PARAMETERISATION.m](../functions/plotting/plot_PARAMETERISATION.m). Colored physical-layer step parameterization in three zones, not spline basis weights. |
| `plotting/plot_PERTURBATIONS.m` | Adapted | [functions/plotting/plot_PERTURBATIONS.m](../functions/plotting/plot_PERTURBATIONS.m). Retained attempted/prior-valid/accepted moves and zone birth/death totals. |
| `plotting/plot_PRIORvsPOSTERIOR.m` | Adapted | [functions/plotting/plot_PRIORvsPOSTERIOR.m](../functions/plotting/plot_PRIORvsPOSTERIOR.m). Conditioned prior/posterior count masses and velocity envelopes; no unsupported anisotropy/HK/NVG panels. |
| `plotting/plot_PRIORvsPOSTERIOR_sed_inprogress.m` | Covered | [functions/plotting/plot_PRIORvsPOSTERIOR.m](../functions/plotting/plot_PRIORvsPOSTERIOR.m). Current conditioned-prior count/envelope comparison supersedes scratch/duplicate plotting variants. |
| `plotting/plot_RF.m` | Adapted | [functions/plotting/plot_RF.m](../functions/plotting/plot_RF.m). Signed RF wiggles with explicit axes, finite time-grid checks and no external maxab dependency. |
| `plotting/plot_SUITE_of_MODELS.m` | Adapted | [functions/plotting/plot_SUITE_of_MODELS.m](../functions/plotting/plot_SUITE_of_MODELS.m). All-state suite with a deterministic display subset by chain; no RNG changes. |
| `plotting/plot_TRADEOFF.m` | Adapted | [functions/plotting/plot_TRADEOFF.m](../functions/plotting/plot_TRADEOFF.m). Cross-dimension scalar tradeoffs; no layer-index or uninferred Vp/Vs/anisotropy correlations. |
| `plotting/plot_TRU_WAVEFORMS.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/plot_TRUvsPRE.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Five observable fits, noisy/clean truth, initial/best model and saved-response intervals. |
| `plotting/plot_TRUvsPRE_WAVEFORMS.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/plot_TRUvsPRE_WAVEFORMS_indiv_figs.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/plot_clust_radplot.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/plot_data_heatmap/plot_data_heatmap.m` | Covered | [THBI/c1_PROCESS_RESULTS.m](../THBI/c1_PROCESS_RESULTS.m). Saved five-response distributions are currently shown as 5/50/95 percentiles, not a full response-density heatmap. |
| `plotting/plot_h_kappa_progress2.m` | Future scope | HK-stack-specific diagnostics; HK stack is not one of the selected five observables. |
| `plotting/plot_invtime.m` | Covered | [functions/write_thbi_report.m](../functions/write_thbi_report.m). Total elapsed seconds and forward counts per chain are reported; original detailed timing breakdown figure is not reproduced. |
| `plotting/plot_prior.m` | Covered | [functions/plotting/plot_PRIORvsPOSTERIOR.m](../functions/plotting/plot_PRIORvsPOSTERIOR.m). Current conditioned-prior count/envelope comparison supersedes scratch/duplicate plotting variants. |
| `plotting/plot_quickmodel.m` | Covered | [THBI/c2_PLOT_RESULTS.m](../THBI/c2_PLOT_RESULTS.m). Layered truth/start/best profiles, shallow views, retained envelopes and counts; c4 adds material summary outputs. Original optional regional/reference-model overlays not ported. |
| `plotting/plot_rec_sec_align.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/plot_rough_sta_map.m` | Future scope | Station mapping/cross-section analysis follows multi-station inversions; not required for the single thin synthetic test. |
| `plotting/progress/plot_progress_chain.m` | Adapted | [functions/plotting/plot_progress_chain.m](../functions/plotting/plot_progress_chain.m). Complete/partial saved-chain likelihood, cumulative acceptance/failures and noise; no dependence on a kernel-reset structure. |
| `plotting/pseudoXsection_all.m` | Future scope | Station mapping/cross-section analysis follows multi-station inversions; not required for the single thin synthetic test. |
| `plotting/receiver_functions/plot_receiver_function_sp.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/receiver_functions/plot_receiver_functions.m` | Future scope | Observed parent/daughter/event/Sp waveform displays require data absent from the current direct single-P-RF synthetic experiment. |
| `plotting/sensitivity/plot_all_sensitivity_kernels.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel sensitivity plots need kernels not computed by the current exact CPS/HVf workflow. |
| `plotting/sensitivity/plot_sensitivity_kernel_hv.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel sensitivity plots need kernels not computed by the current exact CPS/HVf workflow. |
| `plotting/sensitivity/plot_sensitivity_kernel_hv_grid.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel sensitivity plots need kernels not computed by the current exact CPS/HVf workflow. |
| `plotting/sensitivity/plot_sensitivity_kernel_love.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel sensitivity plots need kernels not computed by the current exact CPS/HVf workflow. |
| `plotting/sensitivity/plot_sensitivity_kernel_ray.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel sensitivity plots need kernels not computed by the current exact CPS/HVf workflow. |
| `polarity_est.m` | Future scope | Observed event clustering, polarity estimation or multiple-ray-parameter grouping; current synthetic RF has one specified ray parameter in s/km. |
| `populate_Kbase.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `predat_process.m` | Future scope | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Observed/raw BW-to-RF deconvolution and processing require a separate data pipeline. Current CPS RF is already a Gaussian-filtered signed R/Z transfer response. |
| `predict_atten.m` | Future scope | Attenuation/time/gain waveform comparison; current forward experiment is elastic and does not infer these quantities. |
| `prior_sigma_distribute.m` | Different formulation | [THBI/parms/bayes_inv_parms.m](../THBI/parms/bayes_inv_parms.m). Explicit synthetic base uncertainties plus bounded log-noise priors; no data-derived prior estimation in this controlled experiment. |
| `q_v_given_V.m` | Different formulation | [THBI/b2_PERTURB_MODEL.m](../THBI/b2_PERTURB_MODEL.m). Birth auxiliary Gaussian is written as a normalized log density with negative squared exponent. Original helper has a positive squared exponent and is not copied. |
| `rayp_vals.m` | Future scope | Observed event clustering, polarity estimation or multiple-ray-parameter grouping; current synthetic RF has one specified ray parameter in s/km. |
| `reset_likelihood.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `reverse_inversion_to_last_kbase.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). Kernel prediction/weight/reset/rollback infrastructure. Structural proposals currently use full CPS/HVf calculations; K here often means sensitivity kernels, not inferred layer count. |
| `run_propmat_or_telewavesim.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). PropMat/telewavesim waveform engine; replaced by CPS hrftn96 for the selected synthetic RF definition. |
| `save_inv_state.m` | Covered | [functions/save_checkpoint.m](../functions/save_checkpoint.m). Checkpoint stores full sampler/RNG/prediction/configuration/source state with temporary-write/rename; fuller continuation contract than original save. |
| `sed_vs2rho.m` | Different formulation | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Preserve the user HBI Vp/density equations and floor order in both truth and inversion. Original zone-specific/AK135 conversions would change the experiment. |
| `sed_vs2vp.m` | Different formulation | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Preserve the user HBI Vp/density equations and floor order in both truth and inversion. Original zone-specific/AK135 conversions would change the experiment. |
| `source_time_func.m` | Not used | [THBI/b3_FORWARD_MODEL.m](../THBI/b3_FORWARD_MODEL.m). PropMat/telewavesim waveform engine; replaced by CPS hrftn96 for the selected synthetic RF definition. |
| `spline_weight.m` | Not used | [functions/make_mod_from_parms.m](../functions/make_mod_from_parms.m). Spline basis/weight operations do not apply to physical piecewise-constant layers. |
| `start_hk_ax.m` | Future scope | HK-stack-specific diagnostics; HK stack is not one of the selected five observables. |
| `xconv_misfit.m` | Different formulation | [THBI/b4_CALC_MISFIT.m](../THBI/b4_CALC_MISFIT.m). Cross-convolution of parent/daughter waveforms is a different likelihood from the current direct RF residual; not silently substituted. |

## New main-body stage counterparts

| Reference / project stage | Layered implementation |
|---|---|
| `c2_BUILD_POSTERIOR` | All retained scalar/zone/interface/count/noise quantities with chain/iteration identity. |
| `c3_BUILD_MODEL_SUITE` | Exact stepwise Vs/Vp/density values at common depths; no spline smoothing. |
| `c4_FINAL_MODEL` | Empirical material/boundary/noise quantiles and an actual retained representative Earth. |
| `c5_FINAL_FORWARD_MODEL` | Optional fresh five-observable calculation, misfit and likelihood for that Earth, using CPS/HVf. |

b0/c0 output setup/save roles remain visible in the DATA launcher and chain/checkpoint code. Kernel reset stages b7/b8 are not used because structural proposals are evaluated with full forward calculations.

## Remaining scope and verification

This is still a controlled single-site thin synthetic workflow. Observed-event preparation/deconvolution, multiple RF traces/ray parameters, Sp/CCP, HK stacking, anisotropy, attenuation, regional products, kernel approximations and GPU solvers are not implemented by copying these files. The table makes those boundaries explicit.

`THBI/tests/test_thbi_results.m` adds manual checks for single-chain burn-in, rejected-repeat weighting, exact interface sampling, counts, and selection of a retained physical final model. It is included in `DATA/check_thbi_implementation.m`. The checks were written but not executed. Static analyzer/call findings are recorded in `functions/THBI_SUITE_AUDIT.json`. Source hashes trace the reference and adaptations; `functions/MIGRATION.json` remains the historical record of the earlier 30-helper relocation.
