#!/usr/bin/env bash
# Re-run the jobs from the repository root and compare their logs with git diff.
# Usage: ./verify_all.sh [-j JOBS] [LOG_PATH_FILTER]
# A full run costs roughly 6 CPU hours (the archive-only logs are not rerun); start with a subset.
# Magma jobs run first (JOBS at a time), then the Python jobs one at a time, in the order of
# jobs.txt, since they read the logs and data files written by earlier jobs.
# Jobs whose command starts with "archived:" are not run (see jobs.txt and the README).
# Exit 0 means all selected jobs passed their completion and result checks.
# Failed jobs give a nonzero status (their count, capped at 125).
set -uo pipefail

usage_error() {
    echo "ERROR: $*" >&2
    echo "Usage: ./verify_all.sh [-j JOBS] [LOG_PATH_FILTER]" >&2
    exit 2
}

JOBS=1
while getopts "j:" opt; do
    case $opt in j) JOBS=$OPTARG;; *) usage_error "invalid option";; esac
done
shift $((OPTIND-1))
[[ $JOBS =~ ^[1-9][0-9]*$ ]] || usage_error "-j requires a positive integer"
[ "$#" -le 1 ] || usage_error "expected at most one log-path filter"
FILTER="${1:-}"
[ -r jobs.txt ] || usage_error "jobs.txt is missing or unreadable; run from the repository root"

# Validate the whole manifest before overwriting any logs. The fourth column is optional.
JOB_RE=$'^([^\t]+)\t([^\t]+)\t([^\t]+)(\t([^\t]+))?$'
declare -A SEEN=()
ALL=(); MAGMA_JOBS=(); OTHER_JOBS=(); ARCHIVED=(); NSEL=0; NEED_MAGMA=0; LINE=0
while IFS= read -r job || [ -n "$job" ]; do
    LINE=$((LINE+1))
    [[ $job =~ ^[[:space:]]*(#|$) ]] && continue
    [[ $job =~ $JOB_RE ]] || usage_error "jobs.txt:$LINE: expected three or four tab-separated fields"
    log=${BASH_REMATCH[1]}; cmd=${BASH_REMATCH[2]}; marker=${BASH_REMATCH[3]}
    expected=${BASH_REMATCH[5]:-}
    [ -z "${SEEN[$log]+present}" ] || usage_error "jobs.txt:$LINE: duplicate log path $log"
    SEEN[$log]=1
    if [ -n "$expected" ]; then
        rc=0; grep -E -- "$expected" < /dev/null > /dev/null || rc=$?
        [ "$rc" -le 1 ] || usage_error "jobs.txt:$LINE: invalid result expression"
    fi
    ALL+=("$job")
    if [[ -z $FILTER || $log == *"$FILTER"* ]]; then
        if [[ $cmd == archived:* ]]; then
            ARCHIVED+=("$log")
        elif [[ $cmd == magma\ * ]]; then
            MAGMA_JOBS+=("$job"); NEED_MAGMA=1; NSEL=$((NSEL+1))
        else
            OTHER_JOBS+=("$job"); NSEL=$((NSEL+1))
        fi
    fi
done < jobs.txt
[ "${#ALL[@]}" -gt 0 ] || usage_error "jobs.txt contains no jobs"
[ $((NSEL+${#ARCHIVED[@]})) -gt 0 ] || usage_error "no log path matches '$FILTER'"
[ "$NSEL" -gt 0 ] || usage_error "only archive-only logs match '$FILTER'; nothing to run"

PATCH=mdmagma-magma2.30.patch
if [ "$NEED_MAGMA" -eq 1 ]; then
    command -v magma > /dev/null || usage_error "magma is not on PATH"
    for d in mdmagma/v2/mdmagma.spec mdmagma/Magma/magma.spec gl2/magma.spec fast_hecke.m "$PATCH"; do
        [ -r "$d" ] && continue
        echo "Missing $d. Initialise the pinned packages from the repository root:" >&2
        echo "    git submodule update --init" >&2
        echo "    git -C mdmagma submodule update --init Magma" >&2
        exit 1
    done
fi

STATUSDIR=$(mktemp -d) || exit 1
trap 'rm -rf -- "$STATUSDIR"' EXIT
export STATUSDIR

# The closing line of a log: the last non-blank line, ignoring carriage returns, a final
# "Script done on ..." line written by script(1), and final lines holding only a CPU time.
closing_line() {
    tr -d '\r' < "$1" | grep -av '^[[:space:]]*$' | awk '
        { line[NR] = $0 }
        END {
            n = NR
            while (n > 0 && (line[n] ~ /^Script done on / || line[n] ~ /^[[:space:]]*[0-9]+\.[0-9]+[[:space:]]*$/)) n--
            if (n > 0) print line[n]
        }'
}
fail_job() {
    local record
    echo "!!! $2"
    record=$(mktemp "$STATUSDIR/fail.XXXXXX") || return 1
    printf '%s\n' "$1" > "$record"
    return 1
}
run_one() {
    local log cmd marker expected rc last
    IFS=$'\t' read -r log cmd marker expected <<< "$1"
    mkdir -p "$(dirname "$log")" || { fail_job "$log" "cannot create log directory: $log"; return 1; }
    echo "==> $log : $cmd"
    # Several scripts end without quit; Magma then stops at the end of its standard input.
    rc=0; (eval "$cmd") < /dev/null > "$log" 2>&1 || rc=$?
    [ "$rc" -eq 0 ] || { fail_job "$log" "FAILED (exit $rc): $log"; return 1; }
    # Magma reports most errors on standard output, carries on and exits with status 0.
    if grep -aEq '(Runtime|User|System|Internal) error|Assertion failed' "$log"; then
        fail_job "$log" "MAGMA ERROR: $log reports an error"
        return 1
    fi
    last=$(closing_line "$log")
    if [[ $last != "$marker" && $last != "$marker "* ]]; then
        fail_job "$log" "INCOMPLETE: $log should end in '$marker', last line is '${last:-<empty>}'"
        return 1
    fi
    if [ -n "$expected" ] && ! grep -aEq -- "$expected" "$log"; then
        fail_job "$log" "UNEXPECTED RESULT: $log does not match '$expected'"
        return 1
    fi
    mktemp "$STATUSDIR/ok.XXXXXX" > /dev/null
}
export -f closing_line fail_job run_one

# Magma V2.30 rejects two statements of the pinned mdmagma; apply the patch if needed.
patch_mdmagma() {
    if git -C mdmagma apply --reverse --check "../$PATCH" > /dev/null 2>&1; then
        echo "==> $PATCH is applied to mdmagma"
    elif git -C mdmagma apply "../$PATCH"; then
        echo "==> applied $PATCH to mdmagma"
    else
        echo "!!! could not apply $PATCH to mdmagma; see README, Getting the code" >&2
        return 1
    fi
}

# Compile every package once, in separate processes, before several jobs attach them at the
# same time: concurrent first-time attachment can cause Magma signature-cache errors.
# The three variants of fast_hecke_11.m define the same intrinsics, hence separate processes.
init_packages() {
    local set rc last i=0
    echo "==> initialising the Magma packages"
    for set in \
        "mdmagma/v2/mdmagma.spec fast_hecke.m fast_hecke_11.m fast_hecke_cyclic_commonzero.m hecke_two_commonzero.m fast_diamond_hash_m2.m fast_diamond_hash_roots.m class_rational_reduction.m model_6_12.m" \
        "mdmagma/v2/mdmagma.spec fast_hecke.m fast_hecke_11_m2_guarded.m fast_hecke_cyclic_commonzero.m" \
        "gl2/magma.spec"; do
        i=$((i+1))
        {
            printf 'SetColumns(0);\n'
            for f in $set; do
                case $f in *.spec) printf 'AttachSpec("%s");\n' "$f";; *) printf 'Attach("%s");\n' "$f";; esac
            done
            printf 'print "PACKAGE_READY";\nquit;\n'
        } > "$STATUSDIR/init.m" || return 1
        echo "    $set"
        rc=0; magma -n -b "$STATUSDIR/init.m" < /dev/null > "$STATUSDIR/init.log" 2>&1 || rc=$?
        last=$(closing_line "$STATUSDIR/init.log")
        if [ "$rc" -ne 0 ] || [ "$last" != PACKAGE_READY ]; then
            echo "!!! could not attach $set (exit $rc):" >&2
            tail -20 "$STATUSDIR/init.log" >&2
            return 1
        fi
    done
}

echo "$NSEL of ${#ALL[@]} jobs selected (${#MAGMA_JOBS[@]} Magma, ${#OTHER_JOBS[@]} Python); $JOBS Magma jobs at a time"
for log in "${ARCHIVED[@]}"; do echo "--- skipping archive-only $log"; done
if [ "$NEED_MAGMA" -eq 1 ]; then
    patch_mdmagma || exit 1
    init_packages || exit 1
fi

RUN_STATUS=0
if [ "${#MAGMA_JOBS[@]}" -gt 0 ]; then
    if [ "$JOBS" = 1 ]; then
        for job in "${MAGMA_JOBS[@]}"; do run_one "$job" || RUN_STATUS=$?; done
    else
        printf '%s\n' "${MAGMA_JOBS[@]}" | xargs -d '\n' -P "$JOBS" -I{} bash -uo pipefail -c 'run_one "$@"' _ {} || RUN_STATUS=$?
    fi
fi
for job in "${OTHER_JOBS[@]}"; do run_one "$job" || true; done

shopt -s nullglob
FAILED=("$STATUSDIR"/fail.*); PASSED=("$STATUSDIR"/ok.*)
NFAIL=${#FAILED[@]}
MISSING=$((NSEL-${#FAILED[@]}-${#PASSED[@]}))
if [ "$MISSING" -ne 0 ]; then
    echo "!!! execution did not record a result for every selected job" >&2
    [ "$MISSING" -gt 0 ] && NFAIL=$((NFAIL+MISSING)) || NFAIL=$((NFAIL+1))
elif [ "$RUN_STATUS" -ne 0 ] && [ "$NFAIL" -eq 0 ]; then
    echo "!!! job execution failed (exit $RUN_STATUS)" >&2
    NFAIL=1
fi
if [ "$NFAIL" -gt 0 ]; then
    echo "$NFAIL of $NSEL jobs failed:"
    [ "${#FAILED[@]}" -eq 0 ] || sort "${FAILED[@]}" | sed 's/^/    /'
    [ "$NFAIL" -le 125 ] || NFAIL=125
    exit "$NFAIL"
fi
echo "all $NSEL jobs completed (${#ARCHIVED[@]} archive-only logs not rerun); now run  git diff  to compare with the committed logs"
