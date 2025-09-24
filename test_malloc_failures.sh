#!/bin/bash
# Usage: ./test_malloc_failures.sh ./symnmf_wrap goal data.txt

BIN="$1"
GOAL="$2"
DATA="$3"

i=1
while true; do
    echo ">>> Testing FAIL_MALLOC_AT=$i"
    FAIL_MALLOC_AT=$i $BIN "$GOAL" "$DATA" >out.txt 2>err.txt
    status=$?

    if [ $status -eq 0 ]; then
        echo "Program succeeded at FAIL_MALLOC_AT=$i"
        echo "=> Stopping: no more failure points beyond this"
        break
    else
        echo "Program failed as expected at FAIL_MALLOC_AT=$i"
    fi

    i=$((i+1))
done

echo "Done. Maximum failure site tested: $((i-1))"
