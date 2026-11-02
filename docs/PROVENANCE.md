# Imported forward-solver provenance

These helpers were copied from the existing, previously exercised local HBI forward test.
External HVf/CPS executables are used in place and are not redistributed.

Source directory: `/Users/birotimi/Downloads/hbi_all_codes/AGU_models/five_observable_test`

| File | Original SHA256 |
|---|---|
| forward_five_observables.m | `bab0ac3e0f80746a19bded30e01055cabc8653ccb408a003a7f492f70b043a44` |
| forward_receiver_function_cps.m | `4e53f30c6bf80ddcba5544b55b9f0ede65f6b1cb142e7d24a94c0d733a76d62c` |
| quote_shell_path.m | `2076888706200497168eadb5ea542f4d88ecb315e195d36fb27ad11f17967e89` |
| read_cps_sac.m | `231cf17a378c832193120e1e024425e1b1cd384154bd95b5681244289bec8c82` |
| refine_layer_model.m | `8fae0d699e16d8b9d6627633aba3d5521fbc2103fa9fd5a12b19e74cfb448a18` |
| validate_layer_model.m | `c0e4aa47606a06e0ece6a56e7bc614620b3651c5a8b0569ccc576d3753db689f` |
| write_elastic_cps_model.m | `364c3e965ca28304acfd660fa885875ef07adc5fb8ab91a9ff1a6e7daa085be2` |

The prescribed ECM1 thin layer thicknesses, Vs values, and Brocher conversion equations are retained from `generate_ecm1_five_observables.m` in that source directory.
Only the thin case is included here. No observed dataset, prior/noise estimation, or MCMC is required by this forward-only example.
The source/refined layer matrices and exact solver inputs/logs are written for each new run.
