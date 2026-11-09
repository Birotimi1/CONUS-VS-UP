# Project work preservation

The user has explicitly authorized ongoing daily work logging, local Git
commits, and pushes to Birotimi1/CONUS-VS-UP. Follow this procedure for work
in this repository; do not repeatedly request permission for routine daily
checkpoints and normal pushes that are already authorized.

- Append a factual, timestamped entry to `docs/daily/YYYY-MM-DD.md` whenever
  meaningful project work is completed. Use America/New_York dates and
  ISO-8601 timestamps with the UTC offset. Record changes, decisions, tests
  actually run, unresolved issues, and result locations. Never invent
  activity, execution results or past timestamps.
- Preserve useful work with a new local commit at the end of a work session
  and in the daily scheduled checkpoint. Unfinished work may be committed
  as an explicit WIP checkpoint; do not confuse a backup with validation.
- Review the staged scope and follow `.gitignore`. Preserve user edits.
  Keep native binaries, credentials, scratch files and large generated
  inversion outputs out of source commits. Record relevant output paths in
  the log; those ignored outputs need separate shared-drive preservation.
- Commit locally BEFORE attempting a push. Do not make a local checkpoint
  depend on internet access. Include the actual date/time and timezone offset
  in the daily commit message. Let Git record the true author/committer times;
  never backdate them or amend/squash/rebase a checkpoint to hide its history.
- Verify `origin` is Birotimi1/CONUS-VS-UP, then push the checked-out project
  branch normally. Preserve the current branch; do not switch/merge branches
  or force-push as part of a backup. Verify publication before claiming it.
- If pushing fails, keep the existing local commits unchanged and retry them
  during the next scheduled run or work session. Git preserves their original
  dates and hashes even when they are pushed on a later day. Report failures
  accurately; do not delete work, reset history or silently discard changes.
- If the shared repository is unavailable, do not create a replacement folder
  at the mount path. Report that no checkpoint could be made. If the tree has
  merge conflicts, preserve the files and report the blocked checkpoint.
- With no new work and no pending push, do not create empty commits or invent
  a daily activity entry. See `docs/DAILY_WORKFLOW.md` for the operational steps.
