#!/usr/bin/env bash
# Usage: ./sweep_to_zero.sh ./symnmf_wrap norm ./tests/input_1.txt
# Env (optional):
#   VG_LOG=valgrind_all.log   # combined log file (default)
#   START_AT=1                # first FAIL_MALLOC_AT to try
#   MAX_LOOPS=100000          # safety cap
#   APPEND=0                  # 1 to append to existing log, 0 to truncate

set -uo pipefail

BIN="${1:?Usage: $0 BIN [ARGS...]}"
shift
ARGS=("$@")

VG_LOG="${VG_LOG:-valgrind_all.log}"
START_AT="${START_AT:-1}"
MAX_LOOPS="${MAX_LOOPS:-100000}"
APPEND="${APPEND:-0}"

# Prepare the log and FD 3 (append-only)
if [[ "$APPEND" != "1" ]]; then : > "$VG_LOG"; fi
exec 3>>"$VG_LOG"
trap 'exec 3>&-' EXIT

log() { printf "%s\n" "$*" >&3; }

# Ensure valgrind exists (log only, no console output)
if ! command -v valgrind >/dev/null 2>&1; then
  log "FATAL: valgrind not found in PATH."
  exit 3
fi

# Rich memcheck options; nonzero exit on any memcheck issue
VG_OPTS=(
  --tool=memcheck
  --leak-check=full
  --show-leak-kinds=all
  --errors-for-leak-kinds=all
  --track-origins=yes
  --read-var-info=yes
  --num-callers=40
  --error-exitcode=99
  --error-limit=no
  --verbose
  --log-fd=3
)

i="$START_AT"
while (( i <= MAX_LOOPS )); do
  log ""
  log "===== BEGIN VALGRIND RUN | FAIL_MALLOC_AT=${i} | $(date -Iseconds) ====="
  log "cmd: $BIN ${ARGS[*]}"

  # Run fully quiet; only Valgrind writes into FD 3 (the combined log).
  FAIL_MALLOC_AT="$i" valgrind "${VG_OPTS[@]}" "$BIN" "${ARGS[@]}" >/dev/null 2>&1
  status=$?

  log "===== END VALGRIND RUN   | FAIL_MALLOC_AT=${i} | exit=${status} ====="
  log ""

  # Stop on first clean run (exit code 0)
  if (( status == 0 )); then
    exit 0
  fi

  ((i++))
done

# If we got here, never reached a clean run within MAX_LOOPS
exit 2
