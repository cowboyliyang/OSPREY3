/*
 * Option A (BranchMARK* §3.4): native SIMD log-sum-exp reduction kernel.
 *
 * Java SIMD is unusable on the JVM 17 line (C2 TypeVect::xmeet bug; Vector API
 * hits the same path). This kernel moves the hot log-sum-exp reduction over
 * lambda-states off the C2 vector path into native code.
 *
 * Three implementations share one inf/nan-safe contract (matches
 * RootedTreeEdge.LogSumExpAccumulator):
 *   - NaN anywhere               -> NaN
 *   - all entries -inf (or n==0) -> -inf
 *   - otherwise                  -> max + log(sum_i exp(vals[i] - max))
 *
 *   osprey_logsumexp_scalar_f64  : portable scalar reference (auto-vectorizes
 *                                  to ~1.5x with the safe build flags).
 *   osprey_logsumexp_avx512_f64  : explicit masked AVX-512 (AVX512F+DQ). -inf
 *                                  lanes are clamped so exp() underflows to ~0
 *                                  WITHOUT -ffinite-math-only, so the bounds
 *                                  sentinels survive (the ~2.8x "safe" path).
 *                                  Falls back to scalar if not compiled w/AVX512.
 *   osprey_logsumexp_f64         : dispatches to AVX-512 when available.
 *
 * Two-pass form: mathematically equivalent to the Java streaming accumulator
 * but NOT bit-identical (Option A caveat) — validate before trusting bounds.
 * The AVX-512 exp is a degree-11 range-reduced polynomial (~1e-13 rel on finite
 * inputs), also not bit-identical to libm exp; lse_test.c bounds the error.
 */
#include <stddef.h>
#include <stdint.h>
#include <math.h>

#define OSPREY_LSE_VERSION 2

int osprey_logsumexp_version(void) {
    return OSPREY_LSE_VERSION;
}

/* ------------------------- scalar reference ------------------------- */
double osprey_logsumexp_scalar_f64(const double *vals, int64_t n) {
    double maxv = -INFINITY;
    int has_nan = 0;
    for (int64_t i = 0; i < n; i++) {
        double v = vals[i];
        has_nan |= (v != v);
        maxv = (v > maxv) ? v : maxv;
    }
    if (has_nan) return NAN;
    if (maxv == -INFINITY) return -INFINITY;
    double sum = 0.0;
    for (int64_t i = 0; i < n; i++) {
        sum += exp(vals[i] - maxv);   /* exp(-inf)=0 (no -ffinite-math-only) */
    }
    return maxv + log(sum);
}

/* ------------------------- masked AVX-512 --------------------------- */
#if defined(__AVX512F__) && defined(__AVX512DQ__)
#include <immintrin.h>

/* exp(x) for x in [-700, 0] (pass-2 inputs are <= 0; -inf clamped to -700 -> ~0).
 * Range reduction x = k*ln2 + r, exp(x) = 2^k * exp(r); Horner Taylor for exp(r). */
static inline __m512d osprey_exp512_pd(__m512d x) {
    const __m512d L2E   = _mm512_set1_pd(1.4426950408889634074);
    const __m512d LN2HI = _mm512_set1_pd(6.93145751953125e-1);
    const __m512d LN2LO = _mm512_set1_pd(1.42860682030941723212e-6);
    x = _mm512_max_pd(x, _mm512_set1_pd(-700.0));   /* keep 2^k in normal range */
    x = _mm512_min_pd(x, _mm512_set1_pd(0.0));
    __m512d k = _mm512_roundscale_pd(_mm512_mul_pd(x, L2E),
                    _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    __m512d r = _mm512_fnmadd_pd(k, LN2HI, x);
    r = _mm512_fnmadd_pd(k, LN2LO, r);
    __m512d p = _mm512_set1_pd(2.5052108385441720e-08);                 /* 1/11! */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(2.7557319223985893e-07));   /* 1/10! */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(2.7557319223985888e-06));   /* 1/9!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(2.4801587301587302e-05));   /* 1/8!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(1.9841269841269841e-04));   /* 1/7!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(1.3888888888888889e-03));   /* 1/6!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(8.3333333333333329e-03));   /* 1/5!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(4.1666666666666664e-02));   /* 1/4!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(1.6666666666666666e-01));   /* 1/3!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(5.0000000000000000e-01));   /* 1/2!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(1.0));                      /* 1/1!  */
    p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(1.0));                      /* 1/0!  */
    __m512i ki = _mm512_cvtpd_epi64(k);                                 /* AVX512DQ */
    ki = _mm512_add_epi64(ki, _mm512_set1_epi64(1023));
    ki = _mm512_slli_epi64(ki, 52);
    return _mm512_mul_pd(p, _mm512_castsi512_pd(ki));                   /* p * 2^k */
}

double osprey_logsumexp_avx512_f64(const double *vals, int64_t n) {
    const __m512d NEG_INF = _mm512_set1_pd(-INFINITY);
    __m512d vmax = NEG_INF;
    __mmask8 nanmask = 0;
    int64_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512d v = _mm512_loadu_pd(vals + i);
        nanmask |= _mm512_cmp_pd_mask(v, v, _CMP_UNORD_Q);
        vmax = _mm512_max_pd(vmax, v);
    }
    if (i < n) {
        __mmask8 m = (__mmask8)((1u << (n - i)) - 1u);
        __m512d v = _mm512_mask_loadu_pd(NEG_INF, m, vals + i);  /* inactive = -inf */
        nanmask |= (_mm512_cmp_pd_mask(v, v, _CMP_UNORD_Q) & m);
        vmax = _mm512_max_pd(vmax, v);
    }
    if (nanmask) return NAN;
    double maxv = _mm512_reduce_max_pd(vmax);
    if (maxv == -INFINITY) return -INFINITY;

    __m512d vmaxb = _mm512_set1_pd(maxv);
    __m512d vsum = _mm512_setzero_pd();
    for (i = 0; i + 8 <= n; i += 8) {
        __m512d v = _mm512_loadu_pd(vals + i);
        vsum = _mm512_add_pd(vsum, osprey_exp512_pd(_mm512_sub_pd(v, vmaxb)));
    }
    if (i < n) {
        __mmask8 m = (__mmask8)((1u << (n - i)) - 1u);
        __m512d v = _mm512_maskz_loadu_pd(m, vals + i);
        __m512d e = osprey_exp512_pd(_mm512_sub_pd(v, vmaxb));
        vsum = _mm512_add_pd(vsum, _mm512_maskz_mov_pd(m, e));   /* zero inactive */
    }
    return maxv + log(_mm512_reduce_add_pd(vsum));
}
#else
double osprey_logsumexp_avx512_f64(const double *vals, int64_t n) {
    return osprey_logsumexp_scalar_f64(vals, n);   /* no AVX-512 at compile time */
}
#endif

/* ------------------------- public dispatch -------------------------- */
double osprey_logsumexp_f64(const double *vals, int64_t n) {
#if defined(__AVX512F__) && defined(__AVX512DQ__)
    return osprey_logsumexp_avx512_f64(vals, n);
#else
    return osprey_logsumexp_scalar_f64(vals, n);
#endif
}

/* Fused lower+upper: one JNA call per M-state row amortizes call overhead. */
void osprey_logsumexp2_f64(const double *lower, const double *upper,
                           int64_t n, double *out2) {
    out2[0] = osprey_logsumexp_f64(lower, n);
    out2[1] = osprey_logsumexp_f64(upper, n);
}
