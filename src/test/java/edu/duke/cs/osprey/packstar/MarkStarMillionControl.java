package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import java.nio.file.*;
import java.math.BigDecimal;
import java.util.*;

/** Deterministic MARK* Z control using the frozen census space and energy model. */
public final class MarkStarMillionControl {
    static double log(BigDecimal x) {
        if (x.signum() == 0) return Double.NEGATIVE_INFINITY;
        if (x.signum() < 0) throw new IllegalArgumentException("negative bound");
        int exponent = x.precision() - x.scale() - 1;
        return Math.log(x.movePointLeft(exponent).doubleValue()) + exponent*Math.log(10);
    }
    public static void main(String[] args) throws Exception {
        Path census=Path.of(args[0]), out=Path.of(args[1]);
        String system=args[2];
        String[] config=Files.readAllLines(census.resolve("screen/"+system+".tsv")).get(1).split("\t");
        var space=PackStarMillionCensus.space(census.resolve("input/"+system+".pdb"),config[3]);
        RCs rcs=new RCs(space);
        if (rcs.getNumConformations().longValueExact()!=Long.parseLong(config[2]))
            throw new IllegalStateException("frozen census RC count mismatch");
        List<String> layout=new ArrayList<>(); layout.add("index\tresidue\trcs");
        for (var p:space.positions) layout.add(p.index+"\t"+p.resNum+"\t"+Arrays.toString(rcs.get(p.index)));
        Files.write(out.resolve("layout.tsv"),layout);
        Files.copy(census.resolve("screen/"+system+".tsv"),out.resolve("config.tsv"));
        var parallel=Parallelism.makeCpu(Integer.getInteger("control.cpus",4));
        long start=System.nanoTime();
        try (EnergyCalculator ec=new EnergyCalculator.Builder(space,new ForcefieldParams())
                .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true).setParallelism(parallel).build();
             EnergyCalculator rigid=new EnergyCalculator.SharedBuilder(ec).setIsMinimizing(false).build()) {
            // Match census/PACK* exactly: no reference-energy subtraction.
            ConfEnergyCalculator cc=new ConfEnergyCalculator.Builder(space,ec).build();
            EnergyMatrix min=new SimplerEnergyMatrixCalculator.Builder(cc).build().calcEnergyMatrix();
            EnergyMatrix rig=new SimplerEnergyMatrixCalculator.Builder(new ConfEnergyCalculator(cc,rigid)).build().calcEnergyMatrix();
            long searchStart=System.nanoTime();
            try (var contexts=ec.tasks.contextGroup()) {
                var pf=new MARKStarBoundFastQueues(space,rig,min,cc,rcs,parallel);
                pf.reduceMinimizations=true;
                pf.setCorrectionTighteningEnabled(true);
                pf.setCorrections(new UpdatingEnergyMatrix(space,min,cc));
                pf.setReportProgress(true);
                pf.init(0.683); pf.putTaskContexts(contexts); pf.compute(Integer.MAX_VALUE);
                var r=pf.makeResult();
                Files.writeString(out.resolve("result.tsv"),
                    "system\tmethod\tepsilon\tcount\tstatus\tlogLower\tlogUpper\tsearchSeconds\ttotalSeconds\n"
                    +system+"\tMARKStarBoundFastQueues\t0.683\t"+config[2]+"\t"+r.status
                    +"\t"+log(r.values.calcLowerBound())+"\t"+log(r.values.calcUpperBound())
                    +"\t"+(System.nanoTime()-searchStart)/1e9+"\t"+(System.nanoTime()-start)/1e9+"\n");
                if (r.status!=PartitionFunction.Status.Estimated)
                    throw new IllegalStateException("MARK* did not reach target: "+r.status);
            }
        }
        Files.writeString(out.resolve("READY"),"MARK* Z control completed.\n");
    }
}
