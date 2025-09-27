#!/usr/bin/env bash
# Usage: ./test_malloc_failures.sh ./symnmf_wrap goal data.txt
# Env toggles:
#   USE_VG=1          (default) run under Valgrind
#   MAX_LOOPS=10000   safety cap

set -uo pipefail

BIN="${1:?Usage: $0 BIN GOAL DATA}"
GOAL="${2:?Usage: $0 BIN GOAL DATA}"
DATA="${3:?Usage: $0 BIN GOAL DATA}"

USE_VG="${USE_VG:-1}"
MAX_LOOPS="${MAX_LOOPS:-10000}"

VG_OPTS=(
  --leak-check=full
  --show-leak-kinds=all
  --errors-for-leak-kinds=all
  --track-origins=yes
  --error-exitcode=99
)

run_case() {
  local i="$1"
  local out="out_${i}.txt" err="err_${i}.txt" vglog="valgrind_${i}.log"

  if [ "$USE_VG" = "1" ]; then
    FAIL_MALLOC_AT=$i valgrind "${VG_OPTS[@]}" \
      --log-file="$vglog" \
      "$BIN" "$GOAL" "$DATA" >"$out" 2>"$err"
  else
    FAIL_MALLOC_AT=$i "$BIN" "$GOAL" "$DATA" >"$out" 2>"$err"
  fi
  return $?
}

i=1
while [ "$i" -le "$MAX_LOOPS" ]; do
  echo ">>> Testing FAIL_MALLOC_AT=$i"
  run_case "$i"
  status=$?

  # Did our wrapper report an injected failure? (relies on your fprintf in wrap_alloc.c)
  if grep -q "wrap.*failing" "err_${i}.txt"; then
    # Failure path exercised
    if [ "$USE_VG" = "1" ] && [ $status -ne 0 ]; then
      echo "  Valgrind reported issues at i=$i (exit=$status). See valgrind_${i}.log"
    else
      echo "  Program failed as expected at i=$i (error path hit)."
    fi
  else
    echo "No failure injected at i=$i → reached normal execution path."
    # Final: check success path under Valgrind too (no injection)
    if [ "$USE_VG" = "1" ]; then
      echo ">>> Final check: success path (no injection) under Valgrind"
      valgrind "${VG_OPTS[@]}" \
        --log-file="valgrind_success.log" \
        "$BIN" "$GOAL" "$DATA" >"out_success.txt" 2>"err_success.txt" || true
      echo "  See valgrind_success.log"
    fi
    echo "Done. Maximum failure site tested: $((i-1))"
    exit 0
  fi

  i=$((i+1))
done

echo "Reached MAX_LOOPS=$MAX_LOOPS without finding a non-injected run."
echo "Check wrapper output or increase MAX_LOOPS."
exit 2
