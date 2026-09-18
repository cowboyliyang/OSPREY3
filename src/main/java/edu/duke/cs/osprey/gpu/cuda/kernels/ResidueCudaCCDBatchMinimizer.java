/*
** This file is part of OSPREY 3.0
*/

package edu.duke.cs.osprey.gpu.cuda.kernels;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache.ResPair;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.Gpu;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.MathTools;
import jcuda.Pointer;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.LongBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Runs several independent residue-forcefield CCD problems in one CUDA grid.
 *
 * <p>Each block owns one problem and therefore has the same numerical update
 * order as {@link ResidueCudaCCDMinimizer}.  The static forcefield/geometry
 * data for each problem is concatenated into one buffer; block-local offsets
 * make heterogeneous RC assignments safe.  This class is intentionally
 * limited to the ordinary residue CUDA forcefield path.  Unsupported DOFs or
 * clash-recovery modes are rejected by the caller and use the existing
 * per-conformation fallback.</p>
 */
public final class ResidueCudaCCDBatchMinimizer
        extends Kernel implements AutoCloseable {

    /** One independent objective in the batch. */
    public static final class Problem {
        public final MoleculeObjectiveFunction mof;
        public final ResidueForcefieldEnergy efunc;

        public Problem(MoleculeObjectiveFunction mof,
                       ResidueForcefieldEnergy efunc) {
            this.mof = mof;
            this.efunc = efunc;
        }
    }

    private static final int HeaderBytes = Integer.BYTES * 4 + Double.BYTES * 2;
    private static final int DihedralBytes = Long.BYTES + Integer.BYTES * 4
            + Short.BYTES * 4 + Double.BYTES * 2;
    private static final int ResPairBytes = Long.BYTES * 3 + Double.BYTES * 2;
    private static final int AtomPairBytes = Short.BYTES * 2 + Integer.BYTES
            + Double.BYTES * 9;

    /** Four unsigned-long offsets consumed by residueCcdBatch.cu. */
    private static final int InstanceWords = 4;

    /**
     * Resource probing launches one real CCD block.  Triple correction creates
     * many short-lived batch objects, so probing every object would add a
     * hidden CCD-like cost.  A block size found safe for a larger shared-memory
     * footprint is also safe for every smaller footprint on the same GPU.
     */
    private static final class LaunchConfig {
        final int maximumDihedrals;
        final int blockThreads;

        LaunchConfig(int maximumDihedrals, int blockThreads) {
            this.maximumDihedrals = maximumDihedrals;
            this.blockThreads = blockThreads;
        }
    }

    private static final Map<String, LaunchConfig> launchConfigs =
            new HashMap<>();

    private static final class DihedralInfo {
        final int dof;
        final Residue residue;
        final int[] angleAtoms;
        final int[] rotatedAtoms;
        final int residueIndex;
        final int[] residuePairIndices;
        final double min;
        final double max;

        DihedralInfo(int dof, FreeDihedral source, int residueIndex,
                     int[] residuePairIndices, double min, double max) {
            this.dof = dof;
            this.residue = source.getResidue();
            this.angleAtoms = source.getResidue().template
                    .getDihedralDefiningAtoms(source.getDihedralNumber());
            this.rotatedAtoms = source.getResidue().template
                    .getDihedralRotatedAtoms(source.getDihedralNumber())
                    .stream().mapToInt(Integer::intValue).toArray();
            this.residueIndex = residueIndex;
            this.residuePairIndices = residuePairIndices;
            this.min = min;
            this.max = max;
        }
    }

    private static final class ProblemInfo {
        final Problem problem;
        final List<DihedralInfo> dihedrals;
        final int[] atomOffsetsByResidue;
        final int numAtoms;
        final int maxNumAtoms;
        final int dataBytes;
        long dataOffset;
        long coordsOffset;
        long xinOffset;
        long outOffset;

        ProblemInfo(Problem problem, List<DihedralInfo> dihedrals,
                    int[] atomOffsetsByResidue, int numAtoms,
                    int maxNumAtoms, int dataBytes) {
            this.problem = problem;
            this.dihedrals = dihedrals;
            this.atomOffsetsByResidue = atomOffsetsByResidue;
            this.numAtoms = numAtoms;
            this.maxNumAtoms = maxNumAtoms;
            this.dataBytes = dataBytes;
        }
    }

    private final GpuStreamPool streams;
    private final GpuStream stream;
    private final List<ProblemInfo> infos;
    private CUBuffer<ByteBuffer> data;
    private CUBuffer<DoubleBuffer> coords;
    private CUBuffer<DoubleBuffer> xin;
    private CUBuffer<DoubleBuffer> out;
    private CUBuffer<LongBuffer> instances;
    private Kernel.Function func;
    private boolean cleaned;

    public ResidueCudaCCDBatchMinimizer(GpuStreamPool streams,
                                        List<Problem> problems) {
        super(checkStream(streams), "residueCcdBatch");
        this.streams = streams;
        this.stream = getStream();
        this.infos = new ArrayList<>();

        try {
        if (problems == null || problems.isEmpty()) {
            throw new IllegalArgumentException(
                    "CCD batch must contain at least one problem");
        }
        long totalDataBytes = 0L;
        long totalCoordValues = 0L;
        long totalXinValues = 0L;
        long totalOutValues = 0L;
        int maxDihedrals = 0;

        for (Problem problem : problems) {
            if (problem == null || problem.mof == null
                    || problem.efunc == null) {
                throw new IllegalArgumentException(
                        "CCD batch contains a null problem");
            }
            if (problem.efunc.isBroken) {
                throw new UnsupportedOperationException(
                        "broken conformations use the scalar CCD fallback");
            }

            ProblemInfo info = makeProblemInfo(problem);
            info.dataOffset = align8(totalDataBytes);
            totalDataBytes = info.dataOffset + info.dataBytes;
            info.coordsOffset = totalCoordValues;
            totalCoordValues += (long) info.numAtoms * 3L;
            info.xinOffset = totalXinValues;
            totalXinValues += info.dihedrals.size();
            info.outOffset = totalOutValues;
            totalOutValues += info.dihedrals.size() + 1L;
            maxDihedrals = Math.max(maxDihedrals, info.dihedrals.size());
            infos.add(info);
        }

        if (totalDataBytes > Integer.MAX_VALUE
                || totalCoordValues > Integer.MAX_VALUE
                || totalXinValues > Integer.MAX_VALUE
                || totalOutValues > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(
                    "CCD batch buffers exceed the direct-buffer limit");
        }

        data = stream.byteBuffers.checkout((int) totalDataBytes);
        coords = stream.doubleBuffers.checkout((int) totalCoordValues);
        xin = stream.doubleBuffers.checkout((int) totalXinValues);
        out = stream.doubleBuffers.checkout((int) totalOutValues);
        instances = stream.longBuffers.checkout(
                Math.multiplyExact(infos.size(), InstanceWords));

        writeStaticBuffers();
        data.uploadAsync();
        instances.uploadAsync();
        // Function.calcMaxBlockThreads() probes launchability by executing one
        // block.  Upload the deterministic center pose before that probe so
        // the first problem is not evaluated against uninitialized cudaMalloc
        // contents (the real minimize() upload still refreshes these buffers).
        coords.uploadAsync();
        xin.uploadAsync();

        boolean probe = Boolean.getBoolean("osprey.cuda.debugBatchProbe");
        func = makeFunction(probe ? "ccdBatchProbe" : "ccdBatch");
        func.numBlocks = infos.size();
        final int fMaxDihedrals = maxDihedrals;
        if (probe) {
            func.sharedMemCalc = new Kernel.SharedMemCalculator.None();
            func.setArgs(Pointer.to(
                    instances.getDevicePointer(),
                    out.getDevicePointer()));
        } else {
            func.sharedMemCalc = blockThreads -> blockThreads * Double.BYTES
                    + fMaxDihedrals * Double.BYTES * 3;
            func.setArgs(Pointer.to(
                    data.getDevicePointer(),
                    coords.getDevicePointer(),
                    xin.getDevicePointer(),
                    out.getDevicePointer(),
                    instances.getDevicePointer()));
        }
        // The launch-resource probe uses the same grid and shared-memory
        // footprint as the real batch, but does not need output contents.
        func.blockThreads = getCachedBlockThreads(
                getContext().getGpu(), maxDihedrals, func);
        } catch (RuntimeException | Error ex) {
            clean();
            throw ex;
        }
    }

    private static synchronized int getCachedBlockThreads(
            Gpu gpu, int maximumDihedrals, Kernel.Function func) {
        // Debug-only escape hatch for compute-sanitizer: resource probing
        // executes the kernel repeatedly and can obscure the first illegal
        // address with a poisoned-context cleanup error.  A conservative warp
        // size is sufficient for the diagnostic launch and never used by the
        // production path unless explicitly requested.
        if (Boolean.getBoolean("osprey.cuda.debugBatchFixedThreads")) {
            return Math.min(gpu.getWarpThreads(), gpu.getMaxBlockThreads());
        }
        String key = gpu.getName() + ":"
                + gpu.getComputeVersion()[0] + "."
                + gpu.getComputeVersion()[1] + ":"
                + gpu.getMaxBlockThreads();
        LaunchConfig cached = launchConfigs.get(key);
        if (cached != null && cached.maximumDihedrals >= maximumDihedrals) {
            return cached.blockThreads;
        }

        int probed = func.calcMaxBlockThreads();
        int blockThreads = cached == null
                ? probed : Math.min(cached.blockThreads, probed);
        int maximum = cached == null
                ? maximumDihedrals
                : Math.max(cached.maximumDihedrals, maximumDihedrals);
        launchConfigs.put(key, new LaunchConfig(maximum, blockThreads));
        return blockThreads;
    }

    private static GpuStream checkStream(GpuStreamPool streams) {
        if (streams == null) {
            throw new IllegalArgumentException("CCD batch stream pool is null");
        }
        return streams.checkout();
    }

    private static long align8(long value) {
        return (value + 7L) & ~7L;
    }

    private ProblemInfo makeProblemInfo(Problem problem) {
        MoleculeObjectiveFunction mof = problem.mof;
        ResidueForcefieldEnergy efunc = problem.efunc;

        int[] atomOffsetsByResidue = new int[efunc.residues.size()];
        int atomOffset = 0;
        int numAtoms = 0;
        int maxNumAtoms = 0;
        for (int i = 0; i < efunc.residues.size(); i++) {
            Residue residue = efunc.residues.get(i);
            atomOffsetsByResidue[i] = atomOffset;
            atomOffset += residue.atoms.size() * 3;
            numAtoms += residue.atoms.size();
            maxNumAtoms = Math.max(maxNumAtoms, residue.atoms.size());
        }

        ObjectiveFunction.DofBounds bounds =
                new ObjectiveFunction.DofBounds(mof.getConstraints());
        List<DihedralInfo> dihedrals = new ArrayList<>();
        for (int d = 0; d < mof.getNumDOFs(); d++) {
            DegreeOfFreedom dof = mof.pmol.dofs.get(d);
            if (!(dof instanceof FreeDihedral)) {
                throw new UnsupportedOperationException(
                        "CCD batch supports only FreeDihedral DOFs, got "
                                + dof.getClass().getSimpleName());
            }
            FreeDihedral dihedral = (FreeDihedral) dof;
            int residueIndex = efunc.residues.findIndex(dihedral.getResidue());
            if (residueIndex < 0) continue;
            dihedrals.add(new DihedralInfo(
                    d, dihedral,
                    residueIndex,
                    efunc.makeResPairIndicesSubset(dihedral.getResidue()),
                    bounds.getMin(d), bounds.getMax(d)));
        }

        int paddedRotatedAtoms = 0;
        int paddedResPairIndices = 0;
        for (DihedralInfo dihedral : dihedrals) {
            paddedRotatedAtoms += MathTools.roundUpToMultiple(
                    dihedral.rotatedAtoms.length, Long.BYTES / Short.BYTES);
            paddedResPairIndices += MathTools.roundUpToMultiple(
                    dihedral.residuePairIndices.length,
                    Long.BYTES / Integer.BYTES);
        }
        int totalAtomPairs = 0;
        for (ResPair pair : efunc.resPairs) {
            totalAtomPairs = Math.addExact(totalAtomPairs,
                    pair.info.numAtomPairs);
        }
        long dataBytes = HeaderBytes
                + Long.BYTES * (long) dihedrals.size()
                + Long.BYTES * (long) efunc.resPairs.length
                + DihedralBytes * (long) dihedrals.size()
                + Short.BYTES * (long) paddedRotatedAtoms
                + Integer.BYTES * (long) paddedResPairIndices
                + ResPairBytes * (long) efunc.resPairs.length
                + AtomPairBytes * (long) totalAtomPairs;
        if (dataBytes > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(
                    "one CCD batch problem exceeds the data-buffer limit");
        }
        return new ProblemInfo(problem, dihedrals, atomOffsetsByResidue,
                numAtoms, maxNumAtoms, (int) dataBytes);
    }

    private void writeStaticBuffers() {
        ByteBuffer dataBuf = data.getHostBuffer();
        LongBuffer instanceBuf = instances.getHostBuffer();
        DoubleBuffer coordsBuf = coords.getHostBuffer();
        DoubleBuffer xinBuf = xin.getHostBuffer();

        dataBuf.clear();
        for (ProblemInfo info : infos) {
            int start = Math.toIntExact(info.dataOffset);
            dataBuf.position(start);
            writeProblemData(dataBuf, info);
            int expectedEnd = Math.toIntExact(info.dataOffset
                    + info.dataBytes);
            if (dataBuf.position() != expectedEnd) {
                throw new IllegalStateException("CCD batch data layout mismatch: "
                        + "expected end " + expectedEnd + ", got "
                        + dataBuf.position());
            }

            DoubleMatrix1D center = DoubleFactory1D.dense.make(
                    info.problem.mof.getNumDOFs());
            info.problem.mof.pmol.dofBounds.getCenter(center);
            info.problem.mof.setDOFs(center);

            coordsBuf.position((int) info.coordsOffset);
            for (Residue residue : info.problem.efunc.residues) {
                coordsBuf.put(residue.coords);
            }

            xinBuf.position((int) info.xinOffset);
            for (DihedralInfo dihedral : info.dihedrals) {
                xinBuf.put(Math.toRadians(center.get(dihedral.dof)));
            }

            instanceBuf.put(info.dataOffset);
            instanceBuf.put(info.coordsOffset);
            instanceBuf.put(info.xinOffset);
            instanceBuf.put(info.outOffset);
        }
        dataBuf.clear();
        coordsBuf.clear();
        xinBuf.clear();
        instanceBuf.flip();
    }

    private void writeProblemData(ByteBuffer buf, ProblemInfo info) {
        ResidueForcefieldEnergy efunc = info.problem.efunc;
        ForcefieldParams ffparams = efunc.resPairCache.ffparams;
        List<DihedralInfo> dihedrals = info.dihedrals;

        int flags = ffparams.hElect ? 1 : 0;
        flags = (flags << 1) | (ffparams.hVDW ? 1 : 0);
        flags = (flags << 1) | (ffparams.distDepDielect ? 1 : 0);
        flags = (flags << 1)
                | (ffparams.solvationForcefield == SolvationForcefield.EEF1
                ? 1 : 0);
        buf.putInt(flags);
        buf.putInt(dihedrals.size());
        buf.putInt(efunc.resPairs.length);
        buf.putInt(info.maxNumAtoms);
        double coulombFactor = ForcefieldParams.coulombConstant
                / ffparams.dielectric;
        buf.putDouble(coulombFactor);
        buf.putDouble(coulombFactor * ffparams.forcefld.coulombScaling);

        int dihedralOffsetsPos = buf.position();
        for (int i = 0; i < dihedrals.size(); i++) buf.putLong(0L);
        int resPairOffsetsPos = buf.position();
        for (int i = 0; i < efunc.resPairs.length; i++) buf.putLong(0L);

        for (int d = 0; d < dihedrals.size(); d++) {
            DihedralInfo dihedral = dihedrals.get(d);
            buf.putLong(dihedralOffsetsPos + d * Long.BYTES,
                    buf.position() - info.dataOffset);
            buf.putInt(dihedral.residueIndex);
            buf.putInt(dihedral.residue.atoms.size());
            buf.putLong(info.atomOffsetsByResidue[dihedral.residueIndex]);
            for (int atom : dihedral.angleAtoms) {
                buf.putShort((short) (atom * 3));
            }
            buf.putInt(dihedral.rotatedAtoms.length);
            buf.putInt(dihedral.residuePairIndices.length);
            buf.putDouble(Math.toRadians(dihedral.min));
            buf.putDouble(Math.toRadians(dihedral.max));

            int n = MathTools.roundUpToMultiple(
                    dihedral.rotatedAtoms.length, Long.BYTES / Short.BYTES);
            for (int i = 0; i < n; i++) {
                buf.putShort((short) (i < dihedral.rotatedAtoms.length
                        ? dihedral.rotatedAtoms[i] * 3 : 0));
            }
            n = MathTools.roundUpToMultiple(
                    dihedral.residuePairIndices.length,
                    Long.BYTES / Integer.BYTES);
            for (int i = 0; i < n; i++) {
                buf.putInt(i < dihedral.residuePairIndices.length
                        ? dihedral.residuePairIndices[i] : 0);
            }
        }

        for (int i = 0; i < efunc.resPairs.length; i++) {
            ResPair pair = efunc.resPairs[i];
            buf.putLong(resPairOffsetsPos + i * Long.BYTES,
                    buf.position() - info.dataOffset);
            buf.putLong(pair.info.numAtomPairs);
            buf.putLong(info.atomOffsetsByResidue[pair.resIndex1]);
            buf.putLong(info.atomOffsetsByResidue[pair.resIndex2]);
            buf.putDouble(pair.weight);
            buf.putDouble(pair.offset + pair.solvEnergy);
            for (int j = 0; j < pair.info.numAtomPairs; j++) {
                buf.putLong(pair.info.flags[j]);
            }
            for (int k = 0; k < pair.info.numPrecomputedPerAtomPair; k++) {
                for (int j = 0; j < pair.info.numAtomPairs; j++) {
                    buf.putDouble(pair.info.precomputed[
                            j * pair.info.numPrecomputedPerAtomPair + k]);
                }
            }
        }
    }

    /** Execute the batch and return one minimized energy per problem. */
    public double[] minimize() {
        if (cleaned) {
            throw new IllegalStateException("CCD batch has been cleaned");
        }
        coords.uploadAsync();
        xin.uploadAsync();
        func.runAsync();
        DoubleBuffer outBuf = out.downloadSync();

        double[] energies = new double[infos.size()];
        for (int i = 0; i < infos.size(); i++) {
            ProblemInfo info = infos.get(i);
            int outIndex = Math.toIntExact(info.outOffset);
            // CUDA writes the dihedrals first, followed by the energy.
            energies[i] = outBuf.get(outIndex + info.dihedrals.size());

            // Keep the temporary molecule pose synchronized with the kernel
            // result.  The triple-correction caller only consumes energies,
            // but this preserves the normal Minimizer contract for future
            // callers and makes debugging output meaningful.
            DoubleMatrix1D result = DoubleFactory1D.dense.make(
                    info.problem.mof.getNumDOFs());
            info.problem.mof.pmol.dofBounds.getCenter(result);
            for (int d = 0; d < info.dihedrals.size(); d++) {
                result.set(info.dihedrals.get(d).dof,
                        Math.toDegrees(outBuf.get(outIndex + d)));
            }
            info.problem.mof.setDOFs(result);
        }
        return energies;
    }

    @Override
    public void close() {
        clean();
    }

    public void clean() {
        if (cleaned) return;
        cleaned = true;
        if (data != null) stream.byteBuffers.release(data);
        if (coords != null) stream.doubleBuffers.release(coords);
        if (xin != null) stream.doubleBuffers.release(xin);
        if (out != null) stream.doubleBuffers.release(out);
        if (instances != null) stream.longBuffers.release(instances);
        if (stream != null) streams.release(stream);
    }
}
