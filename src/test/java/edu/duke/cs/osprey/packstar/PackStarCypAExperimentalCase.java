package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import java.nio.file.Files;
import java.nio.file.Path;
import java.io.PrintWriter;
import java.util.*;

/** Frozen local-model experiment. Experimental NMR state mapping is a hypothesis. */
public class PackStarCypAExperimentalCase {
    private static final double T = 283.15;
    private static final double RT = BoltzmannCalculator.RClassic * T;
    private static final String[] VARIANTS = {"WT", "S99T", "S99T_C115S", "S99T_C115S_I97V"};
    private static final double[] EXP = {Double.NaN, .0122, .0193, .0432};
    private static final double[] EXP_ERR = {Double.NaN, .0003, .0005, .0009};
    private static String phe;

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT);
        if (System.getenv("SLURM_JOB_ID") == null) throw new IllegalStateException("Slurm required");
        if (BoltzmannCalculator.TClassic != T) throw new IllegalStateException("Need fully rebuilt 283.15 K experimental overlay");
        Path pdb = Path.of(args[0]), mapping = Path.of(args[1]), out = Path.of(args[2]);
        int requestedVariant = args.length >= 4 ? Integer.parseInt(args[3]) : -1;
        if (requestedVariant < -1 || requestedVariant >= VARIANTS.length)
            throw new IllegalArgumentException("variant index must be -1 or in [0,3]: " + requestedVariant);
        List<Integer> variantIndices = new ArrayList<>();
        if (requestedVariant < 0) {
            for (int variant = 0; variant < VARIANTS.length; variant++)
                variantIndices.add(variant);
        } else {
            variantIndices.add(requestedVariant);
        }
        Files.createDirectories(out);
        Map<Integer,String> ids = new HashMap<>();
        for (String line : Files.readAllLines(mapping)) {
            String[] p = line.trim().split("\\s+");
            if (p.length >= 5 && p[1].equals("A")) ids.put(Integer.parseInt(p[2]), "A" + p[4]);
        }
        for (int n : new int[]{55,61,97,98,99,113,115})
            if (!ids.containsKey(n)) throw new IllegalArgumentException("Missing mapped residue " + n);
        phe = ids.get(113);
        List<SimpleConfSpace> spaces = new ArrayList<>();
        try (PrintWriter preflight = new PrintWriter(out.resolve("preflight.tsv").toFile())) {
            preflight.println("variant\trc_count\tphe_residue\ttemperature_K");
            for (int variant : variantIndices) {
                ForcefieldParams ff = new ForcefieldParams();
                Strand strand = new Strand.Builder(PDBIO.readFile(pdb.toString()))
                    .setTemplateLibrary(new ResidueTemplateLibrary.Builder(ff.forcefld).build())
                    .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames).build();
                for (int n : new int[]{55,61,98})
                    strand.flexibility.get(ids.get(n)).setNoRotamers().addWildTypeRotamers().setContinuous();
                String[] types = {variant == 3 ? "VAL" : "ILE", variant == 0 ? "SER" : "THR",
                                  "PHE", variant >= 2 ? "SER" : "CYS"};
                int[] nums = {97,99,113,115};
                for (int i = 0; i < nums.length; i++)
                    strand.flexibility.get(ids.get(nums[i])).setNoRotamers().setLibraryRotamers(types[i]).setContinuous();
                SimpleConfSpace space = new SimpleConfSpace.Builder().addStrand(strand).build();
                int count = new RCs(space).getNumConformations().intValueExact();
                if (count <= 1 || count > 65536) throw new IllegalArgumentException("RC cap exceeded: " + count);
                spaces.add(space);
                preflight.printf("%s\t%d\t%s\t%.2f%n", VARIANTS[variant], count, phe, T);
            }
        }
        // Check the independent signed-dihedral convention before any reference calculation.
        double[][] fixture = {{1,0,0},{0,0,0},{0,1,0},{.5,1,-.866025403784}};
        if (Math.abs(Protractor.measureDihedral(fixture) - independentAngle(fixture)) > 1e-8)
            throw new IllegalStateException("Dihedral convention mismatch");
        try (PrintWriter results = new PrintWriter(out.resolve("results.tsv").toFile())) {
            results.println("variant\tseed\trc_count\tP_ref\tP_ref_30_90\tP_ref_m15_135\tP_C\tpc_status\tabs_error_ref\tNMR_pB\tNMR_reported_error\tabs_error_NMR\tESS\tmax_weight\tsamples\thits\tZ_status\twithin_model_check\texperimental_mapping\tdiagnostic");
            for (int i = 0; i < spaces.size(); i++) {
                int v = variantIndices.get(i);
                SimpleConfSpace space = spaces.get(i);
                Parallelism parallelism = Parallelism.makeCpu(Integer.getInteger("cypa.threads", 8));
                try (EnergyCalculator min = new EnergyCalculator.Builder(space, new ForcefieldParams())
                         .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true).setParallelism(parallelism).build();
                     EnergyCalculator rigid = new EnergyCalculator.SharedBuilder(min).setIsMinimizing(false).build()) {
                    ConfEnergyCalculator calc = new ConfEnergyCalculator.Builder(space, min).build();
                    ConfEnergyCalculator rcalc = new ConfEnergyCalculator(calc, rigid);
                    RCs rcs = new RCs(space);
                    var rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rcalc).build().calcEnergyMatrix();
                    var minEmat = new SimplerEnergyMatrixCalculator.Builder(calc).build().calcEnergyMatrix();
                    List<double[]> census;
                    try (PrintWriter audit = new PrintWriter(out.resolve(VARIANTS[v] + "_census.tsv").toFile())) {
                        audit.println("assignment\tenergy_kcal_mol\tphe113_chi1_deg");
                        audit.flush();
                        census = enumerate(calc, rcs, audit);
                    }
                    double emin = census.stream().mapToDouble(x -> x[0]).min().orElseThrow();
                    if (!Double.isFinite(emin)) throw new IllegalStateException("No finite conformations");
                    double total = 0, hit = 0, narrow = 0, wide = 0;
                    for (double[] row : census) {
                        double w = Math.exp((emin - row[0]) / RT);
                        total += w;
                        if (inside(row[1], 0, 120)) hit += w;
                        if (inside(row[1], 30, 90)) narrow += w;
                        if (inside(row[1], -15, 135)) wide += w;
                    }
                    double ref = hit / total;
                    System.out.printf("CYPA %s RC=%d P_ref=%.10g%n", VARIANTS[v], census.size(), ref);
                    for (int seed : new int[]{42,43,44}) {
                        System.setProperty("packstar.pac.randomSeed", Integer.toString(seed));
                        try (var contexts = min.tasks.contextGroup();
                             PackStarPartitionFunction pf = new PackStarPartitionFunction(space, rigidEmat, minEmat,
                                 calc, rcs, parallelism, "cypa-" + VARIANTS[v] + "-" + seed)) {
                            pf.setCorrections(new UpdatingEnergyMatrix(space, minEmat, calc));
                            pf.setReduceMinimizations(true);
                            pf.setCorrectionTighteningEnabled(true);
                            pf.setInstanceId(v * 100 + seed);
                            pf.setFunctionalEvent("Phe113_in_chi1_0_120", sample ->
                                inside(Protractor.measureDihedral(coords(sample.getMolecule())), 0, 120));
                            pf.init(.2);
                            pf.putTaskContexts(contexts);
                            pf.compute(Integer.MAX_VALUE);
                            PackStarResult result = pf.makeResult();
                            PackStarFunctionalObservableResult pc = result.getFunctionalObservableResult();
                            double error = Math.abs(pc.getProbability() - ref);
                            String check = !pc.isResolved() ? "UNRESOLVED" : error <= .02 ? "PASS_ABS_0.02" : "FAIL_ABS_0.02";
                            results.printf("%s\t%d\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%s\t%.12g\t%.12g\t%.12g\t%.12g\t%.8g\t%.8g\t%d\t%d\t%s\t%s\t%s\t%s%n",
                                VARIANTS[v], seed, census.size(), ref, narrow/total, wide/total,
                                pc.getProbability(), pc.getStatus(), error, EXP[v], EXP_ERR[v],
                                Math.abs(pc.getProbability()-EXP[v]), pc.getEffectiveSampleSize(), pc.getMaxNormalizedWeight(),
                                pc.getSampleCount(), pc.getHitCount(), result.status, check,
                                "PROVISIONAL_PHE113_IN_TO_NMR_B", pc.getDiagnostic().replace('\t',' ').replace('\n',' '));
                            results.flush();
                        }
                    }
                }
            }
        }
        Files.writeString(out.resolve("COMPUTATION_FINISHED"), "Inspect results; completion is not biological validation success.\n");
    }

    private static List<double[]> enumerate(ConfEnergyCalculator calc, RCs rcs, PrintWriter audit) {
        List<int[]> assignments = new ArrayList<>();
        enumerateAssignments(rcs, new int[rcs.getNumPos()], 0, assignments);
        System.out.printf("CYPA census assignments=%d workers=%d%n", assignments.size(), calc.tasks.getParallelism());
        System.out.flush();
        double[][] rows = new double[assignments.size()][];
        for (int i = 0; i < assignments.size(); i++) {
            int index = i;
            int[] assignment = assignments.get(i);
            calc.calcEnergyAsync(new RCTuple(assignment), energy -> {
                if (Double.isNaN(energy.energy) || energy.energy == Double.NEGATIVE_INFINITY)
                    throw new IllegalStateException("Invalid reference energy");
                if (energy.params != null)
                    for (int j=0; j<energy.params.size(); j++) energy.pmol.dofs.get(j).apply(energy.params.get(j));
                double angle = independentAngle(coords(energy.pmol.mol));
                rows[index] = new double[]{energy.energy, angle};
            });
        }
        calc.tasks.waitForFinish();
        List<double[]> census = new ArrayList<>(rows.length);
        for (int i = 0; i < rows.length; i++) {
            if (rows[i] == null) throw new IllegalStateException("Missing reference energy at assignment " + i);
            census.add(rows[i]);
            audit.printf("%s\t%.12g\t%.12g%n", Arrays.toString(assignments.get(i)), rows[i][0], rows[i][1]);
        }
        return census;
    }

    private static void enumerateAssignments(RCs rcs, int[] a, int pos, List<int[]> assignments) {
        if (pos < a.length) {
            for (int rc : rcs.get(pos)) {
                a[pos] = rc;
                enumerateAssignments(rcs, a, pos + 1, assignments);
            }
        } else {
            assignments.add(a.clone());
        }
    }

    private static boolean inside(double a, double lo, double hi) {
        if (!Double.isFinite(a)) throw new IllegalStateException("Nonfinite dihedral");
        return a >= lo && a <= hi;
    }

    private static double[][] coords(Molecule m) {
        Residue r = m.getResByPDBResNumberOrNull(phe);
        if (r == null) throw new IllegalStateException("Missing Phe113");
        double[][] out = new double[4][];
        String[] names = {"N","CA","CB","CG"};
        for (int i=0;i<4;i++) {
            var atom = r.getAtomByName(names[i]);
            if (atom == null) throw new IllegalStateException("Missing Phe113 atom " + names[i]);
            out[i] = atom.getCoords().clone();
            for (double x : out[i]) if (!Double.isFinite(x)) throw new IllegalStateException("Invalid coordinate");
        }
        return out;
    }

    private static double dot(double[] a, double[] b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
    private static double[] sub(double[] a, double[] b) { return new double[]{a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
    private static double[] cross(double[] a, double[] b) { return new double[]{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}; }
    private static double independentAngle(double[][] p) {
        double[] axis = sub(p[2],p[1]), v=sub(p[0],p[1]), w=sub(p[3],p[2]);
        double norm=Math.sqrt(dot(axis,axis));
        for(int i=0;i<3;i++) axis[i]/=norm;
        double av=dot(v,axis), aw=dot(w,axis);
        for(int i=0;i<3;i++) { v[i]-=av*axis[i]; w[i]-=aw*axis[i]; }
        if(dot(v,v)<1e-20 || dot(w,w)<1e-20) throw new IllegalStateException("Degenerate dihedral");
        return Math.toDegrees(Math.atan2(dot(cross(axis,v),w),dot(v,w)));
    }
}
