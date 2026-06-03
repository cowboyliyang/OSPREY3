/* Test harness for the Option A native log-sum-exp kernel.
 * Compares the explicit AVX-512 path against the scalar reference (and checks
 * inf/nan semantics) across many lengths incl. non-multiples of 8. Built+run by
 * test.sh; NOT part of the gradle build. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

double osprey_logsumexp_scalar_f64(const double *, int64_t);
double osprey_logsumexp_avx512_f64(const double *, int64_t);
double osprey_logsumexp_f64(const double *, int64_t);
int osprey_logsumexp_version(void);

static int fails = 0;

static void check_close(const char *name, int n, double ref, double got, double tol) {
    int ok = 0;
    if (isnan(ref) && isnan(got)) ok = 1;
    else if (isinf(ref) && isinf(got) && ((ref > 0) == (got > 0))) ok = 1;
    else ok = (fabs(ref - got) / (1.0 + fabs(ref))) <= tol;
    if (!ok) {
        printf("FAIL %s n=%d: scalar=%.17g avx=%.17g\n", name, n, ref, got);
        fails++;
    }
}

int main(void) {
    printf("osprey_logsumexp_version=%d\n", osprey_logsumexp_version());
    srand(12345);
    int ns[] = {0, 1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 33, 100, 257, 1000};
    for (size_t t = 0; t < sizeof(ns) / sizeof(ns[0]); t++) {
        int n = ns[t];
        double *v = malloc(sizeof(double) * (n ? n : 1));
        for (int i = 0; i < n; i++) v[i] = ((double)rand() / RAND_MAX) * 100.0 - 50.0;
        check_close("finite", n, osprey_logsumexp_scalar_f64(v, n), osprey_logsumexp_avx512_f64(v, n), 1e-9);
        for (int i = 0; i < n; i += 3) v[i] = -INFINITY;   /* sprinkle -inf sentinels */
        check_close("with_neg_inf", n, osprey_logsumexp_scalar_f64(v, n), osprey_logsumexp_avx512_f64(v, n), 1e-9);
        free(v);
    }
    { double v[8]; for (int i = 0; i < 8; i++) v[i] = -INFINITY;
      double a = osprey_logsumexp_avx512_f64(v, 8);
      if (!(isinf(a) && a < 0)) { printf("FAIL all_neg_inf: %g\n", a); fails++; } }
    { double v[8] = {1, 2, NAN, 4, 5, 6, 7, 8};
      double a = osprey_logsumexp_avx512_f64(v, 8);
      if (!isnan(a)) { printf("FAIL nan_present: %g\n", a); fails++; } }
    { double a = osprey_logsumexp_avx512_f64(NULL, 0);
      if (!(isinf(a) && a < 0)) { printf("FAIL n0: %g\n", a); fails++; } }
    { double v[10] = {0, -1, -800, -1e6, -INFINITY, 5, 5.0001, -2, -3, -700};
      check_close("large_spread", 10, osprey_logsumexp_scalar_f64(v, 10), osprey_logsumexp_avx512_f64(v, 10), 1e-9); }

    if (fails == 0) printf("ALL NATIVE LSE TESTS PASSED\n");
    else printf("%d NATIVE LSE FAILURES\n", fails);
    return fails ? 1 : 0;
}
