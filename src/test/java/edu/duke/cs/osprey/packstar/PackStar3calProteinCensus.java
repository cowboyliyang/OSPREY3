package edu.duke.cs.osprey.packstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.branchdp.InteractionGraph;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import java.io.*;
import java.nio.file.*;
import java.util.*;

/** Exhaustive discrete assignment sum, using the frozen pair38 CPU CCD target. */
public final class PackStar3calProteinCensus {
    static final long TOTAL = 30481920L;
    static final double RT = BoltzmannCalculator.RClassic * BoltzmannCalculator.TClassic;
    static RCs rcs(SimpleConfSpace s) {
        RCs r = s.seqSpace.makeWildTypeSequence().set("C147", "TRP").makeRCs(s);
        int[] expected = {4,6,9,9,35,4,4,7,4};
        if (r.getNumPos() != expected.length || r.getNumConformations().longValueExact() != TOTAL)
            throw new IllegalStateException("Sequence/space count mismatch");
        for (int p=0; p<expected.length; p++)
            if (r.get(p).length != expected[p]) throw new IllegalStateException("RC count mismatch at " + p);
        return r;
    }
    static int[] decode(long id, RCs r) {
        int[] a = new int[r.getNumPos()];
        for (int p=a.length-1; p>=0; p--) { int[] opts=r.get(p); a[p]=opts[(int)(id%opts.length)]; id/=opts.length; }
        if (id!=0) throw new IllegalArgumentException("ID out of range");
        return a;
    }
    static double addLog(double a, double b) {
        if (a==Double.NEGATIVE_INFINITY) return b;
        if (b==Double.NEGATIVE_INFINITY) return a;
        double m=Math.max(a,b); return m+Math.log1p(Math.exp(Math.min(a,b)-m));
    }
    static void layout(SimpleConfSpace s, RCs r, Path out) throws IOException {
        List<String> lines=new ArrayList<>(); lines.add("index\tresidue\trcs");
        for (var p:s.positions) lines.add(p.index+"\t"+p.resNum+"\t"+Arrays.toString(r.get(p.index)));
        Files.write(out.resolve("layout.tsv"), lines);
    }
    static ResidueInteractions interactions(ConfEnergyCalculator cc, int[] conf, int[][] edges) {
        if(edges==null)return cc.makeFragInters(new RCTuple(conf));
        ResidueInteractions inters=new ResidueInteractions();
        for(int p=0;p<conf.length;p++)inters.addAll(cc.makeSingleInters(p,conf[p]));
        for(int[] e:edges)inters.addAll(cc.makePairInters(e[0],conf[e[0]],e[1],conf[e[1]]));
        if(cc.addShellInters)inters.addAll(cc.makeShellInters());
        return inters;
    }
    public static void prepare(TestKStar.ConfSpaces spaces, Parallelism parallelism, String outputDir) throws Exception {
        Path out=Path.of(outputDir); Files.createDirectories(out);
        SimpleConfSpace s=spaces.protein; RCs r=rcs(s); layout(s,r,out);
        Files.writeString(out.resolve("config.tsv"), "system\tpositions\tcount\tflex\tsequence\tstate\n"
            +"3cal\t9\t"+TOTAL+"\tC126,C128,C113,C137,C111,C139,C120,C147,C149\tTRP\tprotein\n");
        try (EnergyCalculator ec=new EnergyCalculator.Builder(spaces.complex, spaces.ffparams)
                .setType(EnergyCalculator.Type.Cpu).setParallelism(parallelism).build()) {
            SimpleReferenceEnergies eref=new SimplerEnergyMatrixCalculator.Builder(s,ec).build().calcReferenceEnergies();
            ConfEnergyCalculator cc=new ConfEnergyCalculator.Builder(s,ec).setReferenceEnergies(eref).build();
            EnergyMatrix minEmat=new SimplerEnergyMatrixCalculator.Builder(cc).build().calcEnergyMatrix();
            InteractionGraph graph;
            try(EnergyCalculator rigid=new EnergyCalculator.Builder(spaces.complex,spaces.ffparams)
                    .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(false).setParallelism(parallelism).build()) {
                ConfEnergyCalculator rigidCC=new ConfEnergyCalculator.Builder(s,rigid)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(s,rigid).build().calcReferenceEnergies()).build();
                EnergyMatrix rigidEmat=new SimplerEnergyMatrixCalculator.Builder(rigidCC).build().calcEnergyMatrix();
                graph=InteractionGraph.buildWithResidualBudget(s,rigidEmat,minEmat,r,1.0,true);
            }
            int[][] edges=graph.getEdgeList().stream().map(e->new int[]{Math.min(e[0],e[1]),Math.max(e[0],e[1])})
                .sorted(Comparator.<int[]>comparingInt(e->e[0]).thenComparingInt(e->e[1])).toArray(int[][]::new);
            int[] degree=new int[9];List<String> edgeLines=new ArrayList<>();
            for(int[] e:edges){degree[e[0]]++;degree[e[1]]++;edgeLines.add(e[0]+"\t"+e[1]);}
            if(!Arrays.equals(degree,new int[]{5,5,7,8,6,6,5,8,4}) || edges.length!=27)
                throw new IllegalStateException("Original rb=1 graph mismatch");
            Files.write(out.resolve("edges.tsv"),edgeLines);
            Files.writeString(out.resolve("target.txt"),"sparse=original rb=1 target; full=all pair interactions\n"
                +"retained_edges="+edges.length+" rho="+graph.getCutResidualUpperBound()+"\n");
            // Workers use the identical frozen builder and this exact energy gauge.
            try (ObjectOutputStream obj=new ObjectOutputStream(new BufferedOutputStream(Files.newOutputStream(out.resolve("eref.ser"))))) {
                obj.writeObject(eref);
            }
            Path samples=Path.of(System.getProperty("census.originalSamples"));
            List<String> rows=Files.readAllLines(samples);
            List<String> header=Arrays.asList(rows.get(0).split("\t"));
            int ci=header.indexOf("conf"), ei=header.indexOf("eTrueKcal");
            if (ci<0 || ei<0) throw new IllegalStateException("Unknown sample schema");
            Map<String,Double> unique=new LinkedHashMap<>();
            for (String row:rows.subList(1,rows.size())) {
                String[] f=row.split("\t"); double val=Double.parseDouble(f[ei]);
                Double prior=unique.putIfAbsent(f[ci],val);
                if (prior!=null && Double.compare(prior,val)!=0) throw new IllegalStateException("Inconsistent original energy");
            }
            List<String> audit=Collections.synchronizedList(new ArrayList<>());
            audit.add("conf\toriginal_sparse_energy\trecomputed_sparse_energy\tabs_difference\tfull_energy");
            double[] max={0};
            for (var row:unique.entrySet()) {
                int[] a=Arrays.stream(row.getKey().split(",")).mapToInt(Integer::parseInt).toArray();
                if (a.length!=r.getNumPos()) throw new IllegalStateException("Original assignment length mismatch");
                for (int p=0;p<a.length;p++) { boolean found=false; for(int opt:r.get(p)) if(opt==a[p])found=true;
                    if(!found)throw new IllegalStateException("Original assignment outside sequence RCs"); }
                ec.tasks.submit(() -> new double[]{cc.calcEnergy(new RCTuple(a),interactions(cc,a,edges)).energy,
                        cc.calcEnergy(new RCTuple(a)).energy}, e -> {
                    double diff=Math.abs(e[0]-row.getValue());
                    if (!Double.isFinite(diff) || diff>1e-7) throw new IllegalStateException("Original CCD mismatch: "+row.getKey()+" "+diff);
                    synchronized(max) {max[0]=Math.max(max[0],diff);}
                    audit.add(row.getKey()+"\t"+row.getValue()+"\t"+e[0]+"\t"+diff+"\t"+e[1]);
                });
            }
            ec.tasks.waitForFinish();
            if(audit.size()!=unique.size()+1)throw new IllegalStateException("Incomplete original-energy audit");
            Files.write(out.resolve("original_energy_check.tsv"),audit);
            Files.writeString(out.resolve("original_energy_check.txt"), "unique="+unique.size()+" max_abs_error="+max[0]+"\n");
        }
        Files.writeString(out.resolve("MODEL_READY"), "Frozen pair38 model and original final-sample energies verified\n");
    }
    public static void main(String[] args) throws Exception {
        String mode=args[0]; Path build=Path.of(args[1]), out=Path.of(args[2]);
        Files.createDirectories(out);
        System.setProperty("osprey.bench.pdbPath",build.resolve("input/3cal.pdb").toString());
        System.setProperty("osprey.bench.mutable","C147");
        System.setProperty("osprey.bench.flexible","D210;D212;C126;C128;C113;C137;C111;C139;C120;C149");
        System.setProperty("osprey.bench.outputDir",out.toString());
        // The census build script supplies this generated builder at runtime.
        // Keep the regular test-source compilation independent of that artifact.
        TestKStar.ConfSpaces spaces=(TestKStar.ConfSpaces) Class.forName(
                "edu.duke.cs.osprey.markstar.bench.Census3calBench")
                .getMethod("buildSpaces", String[].class).invoke(null, (Object) new String[0]);
        if(mode.equals("prepare")) {
            prepare(spaces,Parallelism.makeCpu(Integer.getInteger("census.cpus",4)),out.toString());
            return;
        }
        if(!Files.exists(build.resolve("MODEL_READY")))throw new IllegalStateException("Unverified model");
        SimpleReferenceEnergies eref;
        try(ObjectInputStream obj=new ObjectInputStream(new BufferedInputStream(Files.newInputStream(build.resolve("eref.ser"))))) {
            eref=(SimpleReferenceEnergies)obj.readObject();
        }
        SimpleConfSpace s=spaces.protein; RCs r=rcs(s); layout(s,r,out);
        if(!Files.readString(build.resolve("layout.tsv")).equals(Files.readString(out.resolve("layout.tsv"))))
            throw new IllegalStateException("Frozen RC layout mismatch");
        Files.copy(build.resolve("config.tsv"),out.resolve("config.tsv"));
        String target=System.getProperty("census.target","sparse");
        if(!target.equals("sparse")&&!target.equals("full"))throw new IllegalArgumentException("Unknown target");
        final int[][] edges=target.equals("full")?null:Files.readAllLines(build.resolve("edges.tsv")).stream()
            .map(line->Arrays.stream(line.split("\t")).mapToInt(Integer::parseInt).toArray()).toArray(int[][]::new);
        Files.writeString(out.resolve("target.txt"),target+"\n");
        int shard=Integer.parseInt(args[3]), block=Integer.parseInt(args[4]);
        long start=mode.equals("pilot")?0:(long)shard*block;
        long end=Math.min(start+block,TOTAL);
        if(start<0 || start>=end)throw new IllegalStateException("Empty/invalid shard");
        try(EnergyCalculator ec=new EnergyCalculator.Builder(spaces.complex,spaces.ffparams)
                .setType(EnergyCalculator.Type.Cpu).setParallelism(Parallelism.makeCpu(Integer.getInteger("census.cpus",4))).build();
            DataOutputStream bin=new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(out.resolve("energies.bin"))))) {
            ConfEnergyCalculator cc=new ConfEnergyCalculator.Builder(s,ec).setReferenceEnergies(eref).build();
            double[] logZ={Double.NEGATIVE_INFINITY}; long[] done={0},finite={0};
            Map<Long,Double> repeat=new LinkedHashMap<>(); Set<Long> restored=new HashSet<>();
            long t0=System.nanoTime();
            if(args.length>5 && mode.equals("shard")) {
                Path prior=Path.of(args[5]);
                if(!Files.readString(prior.resolve("layout.tsv")).equals(Files.readString(out.resolve("layout.tsv")))
                    || !Files.readString(prior.resolve("config.tsv")).equals(Files.readString(out.resolve("config.tsv")))
                    || !Files.readString(prior.resolve("target.txt")).equals(target+"\n"))
                    throw new IllegalStateException("Resume configuration mismatch");
                long bytes=Files.size(prior.resolve("energies.bin"));
                try(DataInputStream in=new DataInputStream(new BufferedInputStream(Files.newInputStream(prior.resolve("energies.bin"))))) {
                    for(long k=0;k<bytes/16;k++) {
                        long id=in.readLong();double energy=in.readDouble();
                        if(id<start||id>=end||!restored.add(id)||Double.isNaN(energy)||energy==Double.NEGATIVE_INFINITY)
                            throw new IllegalStateException("Invalid resume record");
                        bin.writeLong(id);bin.writeDouble(energy);done[0]++;
                        if(Double.isFinite(energy)){finite[0]++;logZ[0]=addLog(logZ[0],-energy/RT);}
                        if(id<start+8)repeat.put(id,energy);
                    }
                }
                bin.flush();
                Files.writeString(out.resolve("resume.txt"),"source="+prior+" restored="+restored.size()+" ignored_tail_bytes="+(bytes%16)+"\n");
            }
            for(long batch=start;batch<end;batch+=64) {
                for(long k=batch;k<Math.min(batch+64,end);k++) {
                    final long ordinal=k, id=mode.equals("pilot")?(k*TOTAL)/(end-start):k;
                    if(restored.contains(id))continue;
                    int[] assignment=decode(id,r);
                    cc.calcEnergyAsync(new RCTuple(assignment),interactions(cc,assignment,edges),e->{
                        synchronized(bin) {
                            if(Double.isNaN(e.energy)||e.energy==Double.NEGATIVE_INFINITY)throw new IllegalStateException("Invalid CCD energy "+id);
                            try {bin.writeLong(id);bin.writeDouble(e.energy);}catch(IOException ex){throw new UncheckedIOException(ex);}
                            done[0]++;
                            if(Double.isFinite(e.energy)){finite[0]++;logZ[0]=addLog(logZ[0],-e.energy/RT);}
                            if(ordinal<start+8)repeat.put(id,e.energy);
                        }
                    });
                }
                ec.tasks.waitForFinish();bin.flush();
                if(done[0]%1024==0 || done[0]==end-start)System.out.println("CCD "+done[0]+"/"+(end-start)+" seconds="+(System.nanoTime()-t0)/1e9);
            }
            double seconds=(System.nanoTime()-t0)/1e9,maxError=0;
            for(var row:repeat.entrySet()) {
                int[] assignment=decode(row.getKey(),r);
                double actual=cc.calcEnergy(new RCTuple(assignment),interactions(cc,assignment,edges)).energy;
                if(Double.isFinite(row.getValue()))maxError=Math.max(maxError,Math.abs(actual-row.getValue()));
                else if(actual!=row.getValue())throw new IllegalStateException("Nonfinite repeat mismatch");
            }
            if(done[0]!=end-start || !(maxError<=1e-7))throw new IllegalStateException("Count/repeat check failed");
            Files.writeString(out.resolve("result.tsv"),"system\tmode\ttotal\tstart\tend\tcount\tfinite\tlogZ\tseconds\tmax_repeat_error\n"
                +"3cal\t"+mode+"\t"+TOTAL+"\t"+start+"\t"+end+"\t"+done[0]+"\t"+finite[0]+"\t"+logZ[0]+"\t"+seconds+"\t"+maxError+"\n");
        }
        Files.writeString(out.resolve("READY"),"Count and repeat checks passed\n");
    }
}
