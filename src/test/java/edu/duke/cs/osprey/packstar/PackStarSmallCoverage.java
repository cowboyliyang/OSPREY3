package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import java.nio.file.*;
import java.math.BigDecimal;
import java.util.*;
import java.io.*;

/** Z-only finite CCD-target validation; no functional-event calculation. */
public final class PackStarSmallCoverage {
    static final String[] SYSTEMS = {"1gwc", "2p4a", "2rl0", "4wyu"};
    static final String[][] FLEX = {{"B314","B309","B313"}, {"B226","B225","B228"},
        {"A85","A75","A65"}, {"A94","A70","A90"}};
    static final double RT = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
    static final List<int[]> assignments = new ArrayList<>();
    static void visit(RCs rcs, int p, int[] a) {
        if (p == a.length) { assignments.add(a.clone()); return; }
        for (int rc : rcs.get(p)) { a[p] = rc; visit(rcs, p + 1, a); }
    }
    static String key(int[] a) { return Arrays.toString(a).replace(" ", ""); }
    static double log(BigDecimal x) {
        if (x.signum() == 0) return Double.NEGATIVE_INFINITY;
        if (x.signum() < 0) return Double.NaN;
        int exponent = x.precision() - x.scale() - 1;
        return Math.log(x.movePointLeft(exponent).doubleValue()) + exponent * Math.log(10);
    }
    public static void main(String[] args) throws Exception {
        String mode = args[0]; int caseId = Integer.parseInt(args[1]);
        Path out = Path.of(args[2]); Files.createDirectories(out);
        String system = SYSTEMS[caseId];
        var molecule = PDBIO.readFile(System.getProperty("coverage.pdb"));
        Strand strand = new Strand.Builder(molecule).build();
        for (String res : FLEX[caseId]) strand.flexibility.get(res)
            .setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        SimpleConfSpace space = new SimpleConfSpace.Builder().addStrand(strand).build();
        RCs rcs = new RCs(space);
        int count = rcs.getNumConformations().intValueExact();
        if (rcs.getNumPos() != FLEX[caseId].length || count < 2 || count > 65536)
            throw new IllegalStateException("RC preflight failed: positions=" + rcs.getNumPos() + " count=" + count);
        Files.writeString(out.resolve("space.tsv"), "system\tpositions\tassignments\n" + system + "\t"
            + String.join(",", FLEX[caseId]) + "\t" + count + "\n");
        Parallelism parallelism = Parallelism.makeCpu(Integer.getInteger("coverage.cpus", 8));
        try (EnergyCalculator ec = new EnergyCalculator.Builder(space, new ForcefieldParams())
                .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true).setParallelism(parallelism).build();
             EnergyCalculator rigid = new EnergyCalculator.SharedBuilder(ec).setIsMinimizing(false).build()) {
            // No reference offsets or shell-shell constant in either reference or estimator.
            ConfEnergyCalculator cc = new ConfEnergyCalculator.Builder(space, ec).build();
            if (mode.equals("reference")) {
                visit(rcs, 0, new int[rcs.getNumPos()]);
                double[] energies = new double[count]; double minimum = Double.POSITIVE_INFINITY;
                for (int i = 0; i < count; i++) {
                    energies[i] = cc.calcEnergy(new RCTuple(assignments.get(i))).energy;
                    if (Double.isNaN(energies[i]) || energies[i] == Double.NEGATIVE_INFINITY)
                        throw new IllegalStateException("invalid CCD energy at " + i);
                    minimum = Math.min(minimum, energies[i]);
                    if (i % 100 == 0) System.out.println("REFERENCE " + i + "/" + count);
                }
                if (!Double.isFinite(minimum)) throw new IllegalStateException("no finite energy");
                double[] mass = new double[count]; double sum = 0;
                for (int i = 0; i < count; i++) { mass[i] = Math.exp((minimum - energies[i])/RT); sum += mass[i]; }
                for (int i = 0; i < count; i++) mass[i] /= sum;
                try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out.resolve("census.tsv")))) {
                    w.println("assignment\tenergy\tmass");
                    for (int i=0;i<count;i++) w.println(key(assignments.get(i))+"\t"+energies[i]+"\t"+mass[i]);
                }
                // Independent repeat calls diagnose nondeterministic CCD targets.
                double maxRepeatError=0;
                for (int i=0;i<Math.min(count,10);i++) {
                    double again=cc.calcEnergy(new RCTuple(assignments.get(i))).energy;
                    if (Double.isFinite(energies[i])) maxRepeatError=Math.max(maxRepeatError,Math.abs(again-energies[i]));
                    else if (again!=energies[i]) throw new IllegalStateException("nonfinite repeat drift");
                }
                Files.writeString(out.resolve("reference.tsv"), "system\tcount\tlogZ\tmax_repeat_energy_error\n"
                    +system+"\t"+count+"\t"+(-minimum/RT+Math.log(sum))+"\t"+maxRepeatError+"\n");
                if (!(maxRepeatError<=1e-7)) throw new IllegalStateException("CCD target is not repeatable");
                Files.writeString(out.resolve("READY"), "reference complete\n"); return;
            }
            Path reference = Path.of(args[3]);
            if (!Files.exists(reference.resolve("READY"))) throw new IllegalStateException("reference not ready");
            Map<String,Double> exactEnergies=new HashMap<>();
            for (String line : Files.readAllLines(reference.resolve("census.tsv")).subList(1,count+1)) {
                String[] f=line.split("\t"); exactEnergies.put(f[0],Double.parseDouble(f[1]));
            }
            String[] truth=Files.readAllLines(reference.resolve("reference.tsv")).get(1).split("\t");
            double logZ=Double.parseDouble(truth[2]);
            EnergyMatrix min = new SimplerEnergyMatrixCalculator.Builder(cc).build().calcEnergyMatrix();
            EnergyMatrix rig = new SimplerEnergyMatrixCalculator.Builder(new ConfEnergyCalculator(cc,rigid)).build().calcEnergyMatrix();
            double[] energyError={0}; Set<String> sampled=new HashSet<>(); int[] draws={0};
            try (var contexts=ec.tasks.contextGroup();
                 PackStarPartitionFunction pf=new PackStarPartitionFunction(space,rig,min,cc,rcs,parallelism,"coverage-"+system)) {
                pf.setReduceMinimizations(true); pf.setCorrectionTighteningEnabled(true);
                pf.setSampleListener(s -> {
                    String k=key(s.getConf()); sampled.add(k); draws[0]++;
                    Double expected=exactEnergies.get(k);
                    if (expected==null) throw new IllegalStateException("sample outside census");
                    if (Double.isFinite(expected) && Double.isFinite(s.eTrue))
                        energyError[0]=Math.max(energyError[0],Math.abs(s.eTrue-expected));
                    else if (expected!=s.eTrue) throw new IllegalStateException("nonfinite sampled energy mismatch");
                });
                pf.init(Double.parseDouble(System.getProperty("packstar.pac.targetEpsilon")));
                pf.putTaskContexts(contexts); pf.compute(Integer.MAX_VALUE);
                PackStarResult r=pf.makeResult();
                double lo=log(r.values.calcLowerBound()), hi=log(r.values.calcUpperBound());
                boolean covered=lo<=logZ+1e-10 && hi>=logZ-1e-10;
                Files.writeString(out.resolve("result.tsv"),
                    "system\tseed\tepsilon\tstatus\tlogZ\tlogLower\tlogUpper\tcovered\ttraced_draws\tunique_traced\tmax_energy_error\n"
                    +system+"\t"+System.getProperty("packstar.pac.randomSeed")+"\t"+System.getProperty("packstar.pac.targetEpsilon")+"\t"+r.status+"\t"+logZ+"\t"+lo+"\t"+hi+"\t"+covered+"\t"+draws[0]+"\t"+sampled.size()+"\t"+energyError[0]+"\n");
                if (!(energyError[0]<=1e-7)) throw new IllegalStateException("sample CCD target differs from census");
            }
        }
    }
}
