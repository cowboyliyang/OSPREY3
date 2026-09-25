package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.parallelism.Parallelism;
import java.nio.file.*;
import java.math.BigDecimal;
import java.util.*;

/** Independent estimator runs; never reads census energies or reference Z. */
public final class PackStarMillionCoverage {
    static double log(BigDecimal x) {
        if(x.signum()==0)return Double.NEGATIVE_INFINITY;
        if(x.signum()<0)return Double.NaN;
        int exponent=x.precision()-x.scale()-1;
        return Math.log(x.movePointLeft(exponent).doubleValue())+exponent*Math.log(10);
    }
    public static void main(String[] args)throws Exception {
        Path census=Path.of(args[0]),out=Path.of(args[1]);
        String system=args.length>2?args[2]:"4wyu";
        String[] config=Files.readAllLines(census.resolve("screen/"+system+".tsv")).get(1).split("\t");
        // Reuse the exact already-compiled census space builder and frozen PDB.
        var space=PackStarMillionCensus.space(census.resolve("input/"+system+".pdb"),config[3]);
        RCs rcs=new RCs(space);
        if(rcs.getNumConformations().longValueExact()!=Long.parseLong(config[2]))
            throw new IllegalStateException("census workload mismatch");
        List<String> layout=new ArrayList<>();layout.add("index\tresidue\trcs");
        for(var p:space.positions)layout.add(p.index+"\t"+p.resNum+"\t"+Arrays.toString(rcs.get(p.index)));
        Files.write(out.resolve("layout.tsv"),layout);
        Files.copy(census.resolve("screen/"+system+".tsv"),out.resolve("config.tsv"));
        var parallel=Parallelism.makeCpu(4);
        try(EnergyCalculator ec=new EnergyCalculator.Builder(space,new ForcefieldParams())
                .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true).setParallelism(parallel).build();
            EnergyCalculator rigid=new EnergyCalculator.SharedBuilder(ec).setIsMinimizing(false).build()) {
            ConfEnergyCalculator cc=new ConfEnergyCalculator.Builder(space,ec).build();
            EnergyMatrix min=new SimplerEnergyMatrixCalculator.Builder(cc).build().calcEnergyMatrix();
            EnergyMatrix rig=new SimplerEnergyMatrixCalculator.Builder(new ConfEnergyCalculator(cc,rigid)).build().calcEnergyMatrix();
            try(var contexts=ec.tasks.contextGroup();
                PackStarPartitionFunction pf=new PackStarPartitionFunction(space,rig,min,cc,rcs,parallel,"million-"+system)) {
                pf.setReduceMinimizations(true);pf.setCorrectionTighteningEnabled(true);
                pf.init(0.683);pf.putTaskContexts(contexts);pf.compute(Integer.MAX_VALUE);
                PackStarResult r=pf.makeResult();
                Files.writeString(out.resolve("result.tsv"),"system\tseed\tepsilon\tcount\tstatus\tlogLower\tlogUpper\n"
                    +system+"\t"+System.getProperty("packstar.pac.randomSeed")+"\t0.683\t"+config[2]+"\t"+r.status
                    +"\t"+log(r.values.calcLowerBound())+"\t"+log(r.values.calcUpperBound())+"\n");
            }
        }
    }
}
