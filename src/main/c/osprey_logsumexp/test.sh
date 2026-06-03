#!/bin/bash
# Build + run the native log-sum-exp kernel tests (AVX-512 vs scalar, inf/nan).
# Run on an AVX-512 node (e.g. fennario or this login node). NOT part of gradle.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CC=${CC:-gcc}
SAFE="-O3 -fno-math-errno -funsafe-math-optimizations -fno-finite-math-only"
echo "== 1. AVX-512 build + correctness test =="
$CC $SAFE -mavx512f -mavx512dq -o /tmp/lse_test_avx512 "$HERE/logsumexp.c" "$HERE/lse_test.c" -lm
/tmp/lse_test_avx512
echo "== 2. scalar-only build compiles (no AVX-512; avx512 symbol falls back) =="
$CC $SAFE -mno-avx512f -o /tmp/lse_test_scalar "$HERE/logsumexp.c" "$HERE/lse_test.c" -lm
/tmp/lse_test_scalar
echo "== 3. shared lib builds (what build.sh ships) =="
$CC $SAFE -mavx512f -mavx512dq -fPIC -shared -o /tmp/libOspreyLogSumExp.so "$HERE/logsumexp.c"
echo "OK: /tmp/libOspreyLogSumExp.so"
rm -f /tmp/lse_test_avx512 /tmp/lse_test_scalar /tmp/libOspreyLogSumExp.so
echo "NATIVE KERNEL TEST SCRIPT PASSED"
