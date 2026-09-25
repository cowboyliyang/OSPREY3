package edu.duke.cs.osprey.packstar;

import java.nio.file.*;
import java.util.*;

/** Size-only selection, frozen before energies or estimator results are inspected. */
public final class PackStarSlowThreeScreen {
    public static void main(String[] args) throws Exception {
        Path root=Path.of(args[0]), out=root.resolve("screen");
        Files.createDirectories(out);
        String[] systems={"1a0r","4wwi","4kt6"};
        String[] candidates={
            "G390,G387,G392,G386,G393,B280,B281,B283,B279,B284",
            "B60,E382,E384,E540,E378,E538,B62,B64,B65,B58",
            "D712,D822,D815,D759,C512,C511,C509,C519,C515,C518"};
        String header="system\tpositions\tcount\tflex";
        List<String> audit=new ArrayList<>(); audit.add(header);
        List<String> layouts=new ArrayList<>(); layouts.add("system\tresidue\trcs");
        for(int i=0;i<systems.length;i++) {
            String sys=systems[i]; String[] residues=candidates[i].split(",");
            var space=PackStarMillionCensus.space(root.resolve("input/"+sys+".pdb"),candidates[i]);
            Map<String,Integer> counts=new HashMap<>();
            for(var p:space.positions) { counts.put(p.resNum,p.resConfs.size()); layouts.add(sys+"\t"+p.resNum+"\t"+p.resConfs.size()); }
            String chosen=null; double best=Double.POSITIVE_INFINITY;
            // Keep the original mutable site as WT; enumerate subsets in fixed baseline order.
            for(int mask=1;mask<(1<<residues.length);mask+=2) {
                long product=1; List<String> flex=new ArrayList<>();
                for(int j=0;j<residues.length;j++) if((mask&(1<<j))!=0) {
                    product=Math.multiplyExact(product,counts.get(residues[j])); flex.add(residues[j]);
                }
                String line=sys+"\t"+flex.size()+"\t"+product+"\t"+String.join(",",flex);
                audit.add(line);
                double distance=Math.abs(Math.log(product/1000000.0));
                if(product>=1000000 && product<=2000000 && distance<best) {best=distance;chosen=line;}
            }
            if(chosen==null)throw new IllegalStateException("No 1M..2M subset for "+sys);
            Files.writeString(out.resolve(sys+".tsv"),header+"\n"+chosen+"\n");
            System.out.println("SELECTED\t"+chosen);
        }
        Files.write(out.resolve("subsets.tsv"),audit);
        Files.write(out.resolve("rc_counts.tsv"),layouts);
    }
}
