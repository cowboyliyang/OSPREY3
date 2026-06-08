/*
 * BranchMARK* full-DP CUDA kernels.
 *
 * full_dp_n_children: one non-leaf edge with N (>=1) F-set children, dense
 * parent/child DP tables, full lambda enumeration. Per-child fold plans are
 * concatenated CSR-style (per-child term offsets/counts + per-child table base).
 * One CUDA block per parent M-state; mStateOffset + local output indexing lets
 * the M-state range be split across multiple GPUs.
 *
 * Compile with:
 * nvcc -fatbin
 *   -gencode=arch=compute_86,code=sm_86
 *   -gencode=arch=compute_86,code=compute_86
 *   dp.cu -o dp.bin
 */

#include <math.h>

#define MAX_EDGE_POSITIONS 32

/* Online log-sum-exp: add one term x to a running (max=*m, sum=*s) where the
 * invariant is s == sum_i exp(term_i - m). Single-pass replacement for the old
 * max-pass + exp-pass. NaN in any term poisons the row (matches old max_nan);
 * -inf terms contribute nothing. Invariant: m == -inf  <=>  s == 0 (no finite
 * term yet), so the rescale never hits exp(-inf - -inf). */
__device__ void lse_update(double *m, double *s, double x) {
	if (isnan(x)) { *m = NAN; *s = NAN; return; }
	if (x == -INFINITY) return;
	if (isnan(*m)) return;
	if (x > *m) {
		*s = (*m == -INFINITY) ? 1.0 : (*s) * exp(*m - x) + 1.0;
		*m = x;
	} else {
		*s += exp(x - *m);
	}
}

/* NaN-poisoning max for the block max-reduce: a NaN energy anywhere makes the
 * whole row NaN (matches the Java reference). */
__device__ double max_nan(double a, double b) {
	if (isnan(a) || isnan(b)) { return NAN; }
	return fmax(a, b);
}

__device__ void decode_state(long long idx, const int *counts, int count, int *out) {
	for (int i=count - 1; i>=0; i--) {
		int n = counts[i];
		out[i] = (int)(idx % (long long)n);
		idx /= (long long)n;
	}
}

__device__ double local_energy(
	int lIdx,
	const int *mLocal,
	int *lambdaLocal,
	const int *lambdaCounts,
	int lambdaCount,
	const double *lambdaOnly,
	const double *lmPairs,
	const int *lmLamSlots,
	const int *lmMSlots,
	const int *lmMCounts,
	const long long *lmOffsets,
	int lmPairCount
) {
	decode_state((long long)lIdx, lambdaCounts, lambdaCount, lambdaLocal);

	double e = lambdaOnly[lIdx];
	for (int p=0; p<lmPairCount; p++) {
		int lamSlot = lmLamSlots[p];
		int mSlot = lmMSlots[p];
		long long off = lmOffsets[p]
			+ (long long)lambdaLocal[lamSlot]*(long long)lmMCounts[p]
			+ (long long)mLocal[mSlot];
		e += lmPairs[off];
	}
	return e;
}

/* Mixed-radix child index for one child slice (M terms then lambda terms). */
__device__ long long child_index_slice(
	const int *mLocal,
	const int *lambdaLocal,
	const int *childMSrcAll,
	const long long *childMStrideAll,
	int mOff, int mCnt,
	const int *childLSrcAll,
	const long long *childLStrideAll,
	int lOff, int lCnt
) {
	long long idx = 0;
	for (int t=0; t<mCnt; t++) {
		idx += (long long)mLocal[childMSrcAll[mOff + t]]*childMStrideAll[mOff + t];
	}
	for (int t=0; t<lCnt; t++) {
		idx += (long long)lambdaLocal[childLSrcAll[lOff + t]]*childLStrideAll[lOff + t];
	}
	return idx;
}

/* Sum of one child table (childTable = lower or upper) over all N children
 * at the current (mLocal, lambdaLocal) state. */
__device__ double child_sum(
	const double *childTable,
	const long long *childTableBase,
	int numChildren,
	const int *mLocal,
	const int *lambdaLocal,
	const int *childMSrcAll,
	const long long *childMStrideAll,
	const int *childMTermOff,
	const int *childMTermCnt,
	const int *childLSrcAll,
	const long long *childLStrideAll,
	const int *childLTermOff,
	const int *childLTermCnt
) {
	double s = 0.0;
	for (int c=0; c<numChildren; c++) {
		long long fIdx = child_index_slice(mLocal, lambdaLocal,
			childMSrcAll, childMStrideAll, childMTermOff[c], childMTermCnt[c],
			childLSrcAll, childLStrideAll, childLTermOff[c], childLTermCnt[c]);
		s += childTable[childTableBase[c] + fIdx];
	}
	return s;
}

extern "C" __global__ void full_dp_n_children(
	const int *mCounts,
	const int *lambdaCounts,
	const double *lambdaOnlyRigid,
	const double *lambdaOnlyMin,
	const double *lmRigid,
	const double *lmMin,
	const int *lmLamSlots,
	const int *lmMSlots,
	const int *lmMCounts,
	const long long *lmOffsets,
	const int *childMSrcAll,
	const long long *childMStrideAll,
	const int *childMTermOff,
	const int *childMTermCnt,
	const int *childLSrcAll,
	const long long *childLStrideAll,
	const int *childLTermOff,
	const int *childLTermCnt,
	const long long *childTableBase,
	const double *childLowerAll,
	const double *childUpperAll,
	double *outLower,
	double *outUpper,
	long long mStateCount,
	long long mStateOffset,
	int totalLambdaStates,
	int mCount,
	int lambdaCount,
	int lmPairCount,
	int numChildren,
	double invRT
) {
	long long localIdx = (long long)blockIdx.x;
	long long mIdx = mStateOffset + localIdx;
	if (mIdx >= mStateCount) {
		return;
	}

	extern __shared__ double scratch[];
	double *lowerMax = scratch;
	double *upperMax = lowerMax + blockDim.x;
	double *lowerSum = upperMax + blockDim.x;
	double *upperSum = lowerSum + blockDim.x;

	int tid = threadIdx.x;
	int mLocal[MAX_EDGE_POSITIONS];
	int lambdaLocal[MAX_EDGE_POSITIONS];

	decode_state(mIdx, mCounts, mCount, mLocal);

	/* Single pass: one (max,sum) running pair per thread, child tables + decode
	 * touched once (was twice). The old max-pass then exp-pass are fused. */
	double tLowerMax = -INFINITY, tLowerSum = 0.0;
	double tUpperMax = -INFINITY, tUpperSum = 0.0;
	for (int lIdx=tid; lIdx<totalLambdaStates; lIdx += blockDim.x) {
		double eRigid = local_energy(lIdx, mLocal, lambdaLocal, lambdaCounts, lambdaCount,
			lambdaOnlyRigid, lmRigid, lmLamSlots, lmMSlots, lmMCounts, lmOffsets, lmPairCount);
		double fLower = child_sum(childLowerAll, childTableBase, numChildren, mLocal, lambdaLocal,
			childMSrcAll, childMStrideAll, childMTermOff, childMTermCnt,
			childLSrcAll, childLStrideAll, childLTermOff, childLTermCnt);
		double eMin = local_energy(lIdx, mLocal, lambdaLocal, lambdaCounts, lambdaCount,
			lambdaOnlyMin, lmMin, lmLamSlots, lmMSlots, lmMCounts, lmOffsets, lmPairCount);
		double fUpper = child_sum(childUpperAll, childTableBase, numChildren, mLocal, lambdaLocal,
			childMSrcAll, childMStrideAll, childMTermOff, childMTermCnt,
			childLSrcAll, childLStrideAll, childLTermOff, childLTermCnt);

		lse_update(&tLowerMax, &tLowerSum, -eRigid*invRT + fLower);
		lse_update(&tUpperMax, &tUpperSum, -eMin*invRT + fUpper);
	}

	/* Block reduction in two cheap halves to minimise fp64 exp (A5000 has no
	 * fp64 SFU, so exp is the costly op): max-reduce (no exp), then exactly ONE
	 * rescale per thread, then a plain additive sum-reduce. A full (max,sum)
	 * merge tree would instead pay 2 exp per combine (~8x more exp here). */
	lowerMax[tid] = tLowerMax;
	upperMax[tid] = tUpperMax;
	__syncthreads();
	for (int stride=blockDim.x/2; stride>0; stride >>= 1) {
		if (tid < stride) {
			lowerMax[tid] = max_nan(lowerMax[tid], lowerMax[tid + stride]);
			upperMax[tid] = max_nan(upperMax[tid], upperMax[tid + stride]);
		}
		__syncthreads();
	}
	double gLowerMax = lowerMax[0];
	double gUpperMax = upperMax[0];

	/* Rescale each thread's running sum from its local max to the block max.
	 * gMax == -inf means no finite term in the whole row -> sum stays 0. A finite
	 * thread paired with gMax == NaN yields NaN, which the final guard maps to NaN. */
	lowerSum[tid] = (gLowerMax == -INFINITY) ? 0.0 : tLowerSum * exp(tLowerMax - gLowerMax);
	upperSum[tid] = (gUpperMax == -INFINITY) ? 0.0 : tUpperSum * exp(tUpperMax - gUpperMax);
	__syncthreads();
	for (int stride=blockDim.x/2; stride>0; stride >>= 1) {
		if (tid < stride) {
			lowerSum[tid] += lowerSum[tid + stride];
			upperSum[tid] += upperSum[tid + stride];
		}
		__syncthreads();
	}

	if (tid == 0) {
		outLower[localIdx] = (gLowerMax == -INFINITY || isnan(gLowerMax))
			? gLowerMax
			: gLowerMax + log(lowerSum[0]);
		outUpper[localIdx] = (gUpperMax == -INFINITY || isnan(gUpperMax))
			? gUpperMax
			: gUpperMax + log(upperSum[0]);
	}
}
