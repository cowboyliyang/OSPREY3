/*
 * Batched residue CCD kernel.
 *
 * Keep this translation unit separate from residueCcd.cu.  The scalar module
 * is used by the long-standing residue-CCD path and should retain its exact
 * code generation; adding a second large entry point to that translation unit
 * can perturb floating-point reassociation/register allocation even when the
 * scalar source is unchanged.
 */

#include "residueCcd.cu"

/* Four offsets are stored per independent problem by the Java packer. */
typedef struct __align__(8) {
	unsigned long dataOffset;   // bytes into rawdata
	unsigned long coordsOffset; // doubles into coords
	unsigned long xinOffset;    // doubles into xin
	unsigned long outOffset;    // doubles into out
} BatchInstance;

/*
 * Run the same update sequence as the scalar ccd entry point, but keep it in
 * a noinline helper so the batch wrapper remains a small, launchable kernel.
 * Every block owns one helper invocation and its own dynamically allocated
 * shared-memory region.
 */
__device__ __noinline__ void runCcdBatch(
	const Data & data,
	double * const coords,
	const double * const xin,
	double * const out,
	byte * const shared
) {
	double * const threadEnergies = (double *)shared;
	double * const nextx = threadEnergies + blockDim.x;
	double * const firstSteps = nextx + data.header.numDihedrals;
	double * const lastSteps = firstSteps + data.header.numDihedrals;

	double & herefx = out[data.header.numDihedrals];
	double * const herex = out;

	for (int d = threadIdx.x; d < data.header.numDihedrals; d += blockDim.x) {
		firstSteps[d] = OneDegree;
		lastSteps[d] = OneDegree;
	}
	__syncthreads();

	herefx = calcEnergy(coords, data, NULL, threadEnergies);

	for (int d = threadIdx.x; d < data.header.numDihedrals; d += blockDim.x) {
		herex[d] = xin[d];
	}
	__syncthreads();

	for (int iter = 0; iter < MaxIterations; iter++) {
		copyx(herex, nextx, data.header.numDihedrals);

		for (int d = 0; d < data.header.numDihedrals; d++) {
			const Dihedral & dihedral = data.getDihedral(d);
			double xd = nextx[d];
			double step;
			{
				double firstStep = firstSteps[d];
				double lastStep = lastSteps[d];
				if (fabs(lastStep) > Tolerance && fabs(firstStep) > Tolerance) {
					step = InitialStep * fabs(lastStep / firstStep);
				} else {
					step = InitialStep / pow(iter + 1.0, 3.0);
				}
				while (dihedral.xdmax > dihedral.xdmin
						&& xd - step < dihedral.xdmin
						&& xd + step > dihedral.xdmax) {
					step /= 2;
				}
			}

			LinesearchOut lsout = linesearch(
					coords, data, dihedral, threadEnergies, xd, step);
			if (threadIdx.x == 0) {
				if (iter == 0) firstSteps[d] = lsout.step;
				lastSteps[d] = lsout.step;
				nextx[d] = lsout.xdstar;
			}
			__syncthreads();
		}

		double nextfx = calcEnergy(coords, data, NULL, threadEnergies);
		double improvement = herefx - nextfx;
		if (improvement > 0) {
			copyx(nextx, herex, data.header.numDihedrals);
			herefx = nextfx;
			if (improvement < ConvergenceThreshold) break;
		} else {
			break;
		}
	}
}

extern "C" __global__ void ccdBatch(
	const byte * const rawdata,
	double * const coords,
	const double * const xin,
	double * const out,
	const BatchInstance * const instances
) {
	const BatchInstance instance = instances[blockIdx.x];
	const Data data(rawdata + instance.dataOffset);
	extern __shared__ byte shared[];
	runCcdBatch(
		data,
		coords + instance.coordsOffset,
		xin + instance.xinOffset,
		out + instance.outOffset,
		shared);
}

/*
 * Pointer/offset smoke probe used only while diagnosing a new batch launch.
 * It deliberately does not construct Data or read any forcefield payload;
 * success therefore proves that the Java-side instance table and output
 * addressing are valid independently of the CCD implementation.
 */
extern "C" __global__ void ccdBatchProbe(
	const BatchInstance * const instances,
	double * const out
) {
	const BatchInstance instance = instances[blockIdx.x];
	if (threadIdx.x == 0) {
		out[instance.outOffset] = (double)instance.dataOffset;
	}
}
