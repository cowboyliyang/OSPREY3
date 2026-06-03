#!/bin/bash
# Build the Option A native log-sum-exp kernel into the JNA resource dir.
# NOT run as part of the gradle build (no gradle-native integration yet);
# run manually on the target node, then re-run/submit benchmarks.
#
#   bash src/main/c/osprey_logsumexp/build.sh          # safe (inf/nan-correct, ~1.5x)
#   OSPREY_LSE_FAST=1 bash .../build.sh                # -ffast-math (~2.8x, breaks -inf!)
#
# Safe flags: -fno-math-errno + -funsafe-math-optimizations let gcc vectorize
# the reduction while -fno-finite-math-only KEEPS exp(-inf)=0 working, so the
# bounds sentinels survive. The FAST variant adds -ffast-math (=> finite-only)
# and is ONLY valid if inputs are preconditioned finite.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/../../resources/linux-x86-64"
mkdir -p "$OUT"
CC=${CC:-gcc}
COMMON="-O3 -march=native -fPIC -shared -fno-math-errno"
if [ "${OSPREY_LSE_FAST:-0}" = "1" ]; then
    FLAGS="$COMMON -ffast-math"
    echo "WARNING: FAST build (-ffast-math) breaks -inf sentinels; preconditioned inputs only."
else
    FLAGS="$COMMON -funsafe-math-optimizations -fno-finite-math-only"
fi
echo "$CC $FLAGS -o $OUT/libOspreyLogSumExp.so $HERE/logsumexp.c"
$CC $FLAGS -o "$OUT/libOspreyLogSumExp.so" "$HERE/logsumexp.c"
echo "built: $OUT/libOspreyLogSumExp.so"
