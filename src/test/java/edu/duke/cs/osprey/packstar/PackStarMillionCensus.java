package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import java.nio.file.*;
import java.io.*;
import java.util.*;

/** CPU-only, bounded-memory enumeration of the finite full-CCD target. */
public final class PackStarMillionCensus {
    static final String[] SYSTEMS={"1gwc","2p4a","2rl0","4wyu"};
    static final String[] CANDIDATES={
        "B314,B309,B313,B311,B317,C514,C510,C505,C446,C503,C511,C498,C500",
        "B226,B225,B228,D476,D475,A111,A69,A71",
        "A85,A75,A65,A63,G95,G97,A49,A48,A51,A50,A87,A77,A76",
        "A94,A70,A90,A96,A93,A92,A72,A69,A88,B212,B213,A12"};
    static final double RT=BoltzmannCalculator.RClassic*BoltzmannCalculator.TClassic;
    static SimpleConfSpace space(Path pdb,String flex) {
        Strand strand=new Strand.Builder(PDBIO.readFile(pdb.toFile())).build();
        for(String r:flex.split(",")) strand.flexibility.get(r).setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers().setContinuous();
        SimpleConfSpace s=new SimpleConfSpace.Builder().addStrand(strand).build();
        if(s.positions.size()!=flex.split(",").length) throw new IllegalStateException("missing requested position");
        return s;
    }
    static int[] decode(long id, RCs rcs) {
        int[] a=new int[rcs.getNumPos()];
        for(int p=a.length-1;p>=0;p--) { int[] options=rcs.get(p); a[p]=options[(int)(id%options.length)]; id/=options.length; }
        if(id!=0) throw new IllegalArgumentException("ID out of range");
        return a;
    }
    static double addLog(double a,double b) {
        if(a==Double.NEGATIVE_INFINITY)return b;
        if(b==Double.NEGATIVE_INFINITY)return a;
        double m=Math.max(a,b); return m+Math.log1p(Math.exp(Math.min(a,b)-m));
    }
    public static void main(String[] args)throws Exception {
        String mode=args[0]; Path build=Path.of(args[1]), out=Path.of(args[2]);
        Files.createDirectories(out);
        if(mode.equals("screen")) {
            List<String> report=new ArrayList<>(); report.add("system\tpositions\tcount\tflex");
            for(int c=0;c<SYSTEMS.length;c++) {
                String sys=SYSTEMS[c]; SimpleConfSpace all=space(build.resolve("input/"+sys+".pdb"),CANDIDATES[c]);
                Map<String,Integer> counts=new HashMap<>();
                for(var p:all.positions)counts.put(p.resNum,p.resConfs.size());
                long product=1; List<String> flex=new ArrayList<>(); String chosen=null; double best=Double.POSITIVE_INFINITY;
                for(String r:CANDIDATES[c].split(",")) {
                    product=Math.multiplyExact(product,counts.get(r)); flex.add(r);
                    String line=sys+"\t"+flex.size()+"\t"+product+"\t"+String.join(",",flex); report.add(line);
                    double distance=Math.abs(Math.log(product/1000000.0));
                    if(product>=500000 && product<=2000000 && distance<best) { chosen=line;best=distance; }
                }
                if(chosen!=null) Files.writeString(out.resolve(sys+".tsv"),report.get(0)+"\n"+chosen+"\n");
            }
            Files.write(out.resolve("screen.tsv"),report); return;
        }
        Path config=Path.of(args[3]); String[] f=Files.readAllLines(config).get(1).split("\t");
        String system=f[0]; long expected=Long.parseLong(f[2]);
        SimpleConfSpace s=space(build.resolve("input/"+system+".pdb"),f[3]); RCs rcs=new RCs(s);
        long total=rcs.getNumConformations().longValueExact();
        if(total!=expected)throw new IllegalStateException("frozen RC count drift");
        List<String> layout=new ArrayList<>();layout.add("index\tresidue\trcs");
        for(var p:s.positions) layout.add(p.index+"\t"+p.resNum+"\t"+Arrays.toString(rcs.get(p.index)));
        Files.write(out.resolve("layout.tsv"),layout);
        int shard=Integer.parseInt(args[4]), block=Integer.parseInt(args[5]);
        long start=mode.equals("pilot")?0:(long)shard*block;
        long end=mode.equals("pilot")?Math.min(block,total):Math.min(start+block,total);
        if(start>=end)throw new IllegalStateException("empty shard");
        Files.copy(config,out.resolve("config.tsv"));
        try(EnergyCalculator ec=new EnergyCalculator.Builder(s,new ForcefieldParams())
            .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true)
            .setParallelism(Parallelism.makeCpu(Integer.getInteger("census.cpus",4))).build();
            DataOutputStream binary=new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(out.resolve("energies.bin"))))) {
            ConfEnergyCalculator cc=new ConfEnergyCalculator.Builder(s,ec).build();
            double[] logZ={Double.NEGATIVE_INFINITY}; long[] finite={0}, done={0};
            Map<Long,Double> repeat=new LinkedHashMap<>(); long t0=System.nanoTime();
            Set<Long> restored=new HashSet<>();
            if(args.length>6 && mode.equals("shard")) {
                Path prior=Path.of(args[6]);
                if(Files.exists(prior.resolve("energies.bin"))) {
                    if(!Files.readString(prior.resolve("layout.tsv")).equals(Files.readString(out.resolve("layout.tsv")))
                        || !Files.readString(prior.resolve("config.tsv")).equals(Files.readString(config)))
                        throw new IllegalStateException("resume configuration/layout drift");
                    long bytes=Files.size(prior.resolve("energies.bin"));
                    try(DataInputStream input=new DataInputStream(new BufferedInputStream(Files.newInputStream(prior.resolve("energies.bin"))))) {
                        for(long record=0;record<bytes/16;record++) {
                            long id=input.readLong();double energy=input.readDouble();
                            if(id<start||id>=end||!restored.add(id)||Double.isNaN(energy)||energy==Double.NEGATIVE_INFINITY)
                                throw new IllegalStateException("invalid/duplicate checkpoint record");
                            binary.writeLong(id);binary.writeDouble(energy);
                            if(Double.isFinite(energy)){finite[0]++;logZ[0]=addLog(logZ[0],-energy/RT);}
                            if(id<start+8)repeat.put(id,energy);
                            done[0]++;
                        }
                    }
                    Files.writeString(out.resolve("resume.tsv"),"source\trestored\ttrailing_bytes_recomputed\n"+prior+"\t"+restored.size()+"\t"+(bytes%16)+"\n");
                    binary.flush();
                    System.out.println("RESUMED "+restored.size()+" complete records; ignored partial tail bytes="+(bytes%16));
                }
            }
            for(long batch=start;batch<end;batch+=64) {
                for(long k=batch;k<Math.min(batch+64,end);k++) {
                    final long ordinal=k;
                    final long id=mode.equals("pilot")?(k*total)/(end-start):k;
                    if(restored.contains(id))continue;
                    cc.calcEnergyAsync(new RCTuple(decode(id,rcs)),e->{
                        synchronized(binary) {
                            if(Double.isNaN(e.energy)||e.energy==Double.NEGATIVE_INFINITY)
                                throw new IllegalStateException("invalid energy for ID "+id);
                            try { binary.writeLong(id);binary.writeDouble(e.energy); }
                            catch(IOException ex){throw new UncheckedIOException(ex);}
                            if(Double.isFinite(e.energy)) { finite[0]++;logZ[0]=addLog(logZ[0],-e.energy/RT); }
                            if(ordinal<start+8)repeat.put(id,e.energy);
                            done[0]++;
                        }
                    });
                }
                ec.tasks.waitForFinish(); binary.flush();
                System.out.println("CCD "+done[0]+"/"+(end-start)+" seconds="+(System.nanoTime()-t0)/1e9);
            }
            double seconds=(System.nanoTime()-t0)/1e9, maxError=0;
            for(var row:repeat.entrySet()) {
                double actual=cc.calcEnergy(new RCTuple(decode(row.getKey(),rcs))).energy;
                if(Double.isFinite(row.getValue()))maxError=Math.max(maxError,Math.abs(actual-row.getValue()));
                else if(actual!=row.getValue())throw new IllegalStateException("nonfinite repeat mismatch");
            }
            if(done[0]!=end-start||!(maxError<=1e-7))throw new IllegalStateException("count or repeat check failed");
            Files.writeString(out.resolve("result.tsv"),"system\tmode\ttotal\tstart\tend\tcount\tfinite\tlogZ\tseconds\tmax_repeat_error\n"
                +system+"\t"+mode+"\t"+total+"\t"+start+"\t"+end+"\t"+done[0]+"\t"+finite[0]+"\t"+logZ[0]+"\t"+seconds+"\t"+maxError+"\n");
        }
        Files.writeString(out.resolve("READY"),"count and repeat checks passed\n");
    }
}
