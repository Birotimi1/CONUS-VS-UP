# Dependency audit for the thin five-observable workflow

Current extension: [2026-11-07 Ps/Sp, 60-s bands and delayed-rejection review](../docs/THBI_NEXT_RUN_REVIEW.md).
The audit below describes the earlier five-block profile; its exclusions and counts are historical.

Static audit updated: 2026-11-04. Scope includes both HBI and THBI,
the standalone thin forward demonstration, seven DATA scripts and manual tests.
The broader [127-file reference-suite comparison](../docs/THBI_FUNCTION_SUITE.md)
separately addresses feature coverage; a resolved call graph alone does not
establish parity with Zach/Brennan's entire workflow.

The standalone R2022a Code Analyzer call listing (`mlint -calls`) covered
102 MATLAB files. The `remove_scratch` anonymous cleanup callback was also
checked separately. Every
project call has a repository implementation. Remaining calls are MATLAB
functions/language constructs, plus the optional Parallel Computing Toolbox
calls below. No unresolved custom function or execution dependency on either
Downloads source checkout was identified. This is a static inventory, not a
runtime dependency trace.

The repository-level `functions/` tree now contains **42 helpers**, including
nine in `functions/plotting/`. The initial 30-file consolidation and subsequent
12 helper additions cover the selected layered synthetic workflow. No helper copies remain in
project `HBI/functions`, `THBI/functions`, or `src`. The numbered stages retain
their main calculation blocks. The complete inventory, each file's role and
source/migration hashes are in the [functions guide](../functions/README.md)
and [migration record](../functions/MIGRATION.json).

Both a0 stages check MATLAB resolution of every helper against its exact
project file, as well as the selected numbered stages and parameter function.
Both parallel chain runners select the relevant startup on process workers.
This catches shadowing by another checkout when the user runs the launcher.
These runtime path checks have been written and statically checked, not run.
The fixed-layer compatibility helpers no longer change the search path.

## Dependencies outside the functions folder

- Base MATLAB. The source was statically analyzed using the installed
  R2022a standalone analyzer. Plotting uses MATLAB graphics/export functions.
- Parallel Computing Toolbox only when `par.inv.parallel=true`:
  `gcp`, `parcluster`, `parpool`, and process-worker execution. The default
  serial workflow does not select that branch. License availability was
  not tested.
- HVf in
  `/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/AGU_models/bin_HV_DFA/`.
- CPS `surf96` and `hrftn96` in the sibling `bin_v3.30/` directory, with
  the CPS companion programs. The installed `srfpre96`, `srfdis96`,
  `srfdrl96`, `srfdrr96`, `srfinv96` and `srfker96` were also inspected;
  some support CPS operations beyond the selected Rayleigh forward mode.

The 2026-11-03 inspection found all nine native files present with executable
permission on this machine. Reading their Mach-O dependency tables with `otool -L`
identified direct non-system dependencies on `libgfortran.5.dylib`,
`libgomp.1.dylib` and `libquadmath.0.dylib` under
`/opt/homebrew/opt/gcc/lib/gcc/current/`; those files exist. This checks
installation metadata, not executable loading or numerical operation.
These binaries and libraries are external installations, not copied into
the MATLAB functions folder.

## What this audit establishes

No missing custom function was identified for the current controlled thin
experiment. Native file presence and direct library presence were checked.
The two updated startup files have no Code Analyzer findings. Across the
102 files, findings were style/performance/redundant-suppression notices and the explicitly
disabled optional test branches; no syntax errors were reported.

The earlier relocation's static source comparison checked all 30 helper files: 27 retained their
calculation bodies with only documented path/call-name adjustments; three
had explicit path or snapshot edits. Six affected main stage bodies were
compared, and all seven former HBI b3 local helpers match the common versions.
The subsequent extension adds actual layered counterparts of the retained-
model/posterior/plotting roles, including c2-c5 stages. Its reference hashes,
adaptation hashes, current source hashes and analyzer findings are in
`functions/THBI_SUITE_AUDIT.json`. Every helper has a caller. These source
checks do not demonstrate that the new numerical checks or plots execute.

THBI's source snapshot now includes the shared helpers recursively, including
the plotting subfolder, and excludes DATA.
Its changed source paths invalidate pre-consolidation checkpoints; use a new
run directory for this code layout.

This does **not** establish successful MATLAB execution, graphics export,
parallel operation, detailed balance in execution, native model evaluation,
mixing or recovery. No MATLAB test, inversion or native solver was launched.
Use `DATA/check_thbi_implementation.m`, then its optional prior/native
checks, before the synthetic pilot. The broader differences from the
source HBI/THBI systems are listed in [the review](../docs/THBI_REVIEW.md).
