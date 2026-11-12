#!/usr/bin/env bash
# Run one native command, echo liveness, and bound its wall-clock duration.
# Solver stdout/stderr redirects belong to COMMAND, not to this monitor.
set -u
interval=$1; limit=$2; label=$3; command=$4
[[ $interval =~ ^[0-9]+$ && $limit =~ ^[0-9]+$ ]] || exit 64
child=''
stop_child() {
    if [[ -n $child ]]; then
        kill -TERM -- "-$child" 2>/dev/null || true
        sleep 0.2
        kill -KILL -- "-$child" 2>/dev/null || true
        wait "$child" 2>/dev/null || true
    fi
}
trap 'stop_child; exit 130' INT
trap 'stop_child; exit 143' TERM HUP
# A separate process group permits cleanup of the calculator and its children.
set -m
bash -c "$command" </dev/null &
child=$!
set +m
started=$SECONDS; next=$interval
while kill -0 "$child" 2>/dev/null; do
    elapsed=$((SECONDS-started))
    if (( limit > 0 && elapsed >= limit )); then
        printf '[native] TIMEOUT %s after %ds; stopping run, not rejecting a model.\n' "$label" "$elapsed"
        stop_child
        exit 124
    fi
    if (( interval > 0 && elapsed >= next )); then
        printf '[native] %s still running (%ds), pid=%s, scratch=%s\n' "$label" "$elapsed" "$child" "$PWD"
        next=$((elapsed+interval))
    fi
    sleep 0.2
done
wait "$child"
result=$?
child=''
exit "$result"
