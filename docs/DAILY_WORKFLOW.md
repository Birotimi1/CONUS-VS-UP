# Daily work log and Git preservation

Repository: `Birotimi1/CONUS-VS-UP`.
Working checkout:
`/Volumes/AS-Filer/EES/jbrussel/SharedData/US_TA/CONUS-VS-UP`.
Schedule: daily at **22:00 America/New_York** in the existing Codex chat.
The scheduled task supplements checkpoints made at the end of work sessions.
Automation ID: `conus-daily-work-log-and-git-checkpoint` (ACTIVE, attached to
the existing project chat).

## What is recorded

Use `docs/daily/YYYY-MM-DD.md`, append-only by timestamp within that day.
Record the actual work date, a timestamp with UTC offset, the current branch,
changes/decisions, checks that really ran, unresolved work and relevant result
paths. Cite known commit hashes; the commit containing a new entry is found
from Git history, so the entry does not need to contain its own future hash.
Do not reconstruct undocumented prior days or replace real Git dates with
estimates. If work is recovered later, record the recovery time separately
from any reliably known original work time.

The daily check reviews uncommitted changes plus recent commits not yet
covered by the log. Already committed work can still need a factual log
entry. Repeated checks on an unchanged tree must not create empty commits.

## Local checkpoint, then GitHub

1. Confirm that the shared drive and this checkout exist. Read `AGENTS.md`,
   inspect the current branch, `git status`, recent commits and the day's log.
2. Review actual changes, including user edits. Append the work log before
   staging relevant source, scripts, tests, documentation and the log itself.
   Do not run a long inversion just to make a preservation checkpoint.
3. Create a local commit with the actual timestamp. Example message format:
   `Daily checkpoint YYYY-MM-DDTHH:MM:SS-05:00 — <brief summary>`.
   Use America/New_York for display; its UTC offset changes with daylight
   saving time. Do not set artificial author/committer dates. Git inherently
   stores a timestamp and timezone offset in each commit.
4. Verify the expected `origin`, then attempt a normal push of the current
   branch. Do not require a successful fetch/network connection before
   making the local commit. Record whether the push actually succeeded.
5. If a push fails, the work remains committed in the shared checkout's
   `.git` directory. Retry those same commits later. Do not amend, squash,
   rebase, reset or force-push them as part of this process. If remote history
   has diverged, report it and keep both histories intact for review.

To inspect preserved dates and pending commits:

```sh
git log --format='%h | author %aI | committed %cI | %s'
git status --short --branch
git log --oneline '@{upstream}..HEAD'
```

The upstream comparison reflects the last known remote state until refreshed.
A successful push and an independent remote-head check establish publication.
Use the automation's run result or a later log entry to record a push outcome;
do not claim success in a pre-push entry or create endless status-only commits.

## Scope and availability

This preserves versioned source and work logs. `DATA/transdimensional/`,
`DATA/inversions/`, `DATA/synthetics/`, temporary files and native binaries
remain outside Git as configured. Important experiment outputs stay on the
shared drive; the log should identify their paths. This daily Git procedure
is not an additional backup of those ignored outputs or of unsaved editor
buffers.

The computer must be on, the desktop app running, and the shared drive mounted
for a local scheduled task to access the checkout. If a run is missed, perform
the checkpoint at the next opportunity using its actual timestamp; do not
fabricate a commit for a time when the job did not execute. Official behavior:
[scheduled tasks documentation](https://learn.chatgpt.com/docs/automations?surface=app).

## Resume after usage interruptions

At the user's request, the chat had a heartbeat named
`CONUS resume work after usage reset`, automation ID
`conus-resume-work-after-usage-reset`, scheduled every 5 minutes (refreshed
2026-11-09 after an overdue, never-dispatched retry was found). It was deleted
on 2026-11-10 after the native execution audit completed successfully; there
is no remaining audit work to retry. The daily 22:00 America/New_York
checkpoint remains active as a separate automation. The full production
synthetic inversion remains user-run.

Any future retry should respect pauses, cancellations, normal approval review
and required user input. It must not redeem resets or purchase credits. If work is complete
or awaiting the user's native MATLAB validation/inversion results, it takes
no project action. It must not generate new work or rerun completed tests just
because the schedule fired. The same local-app, storage and permissions
requirements described above apply.
