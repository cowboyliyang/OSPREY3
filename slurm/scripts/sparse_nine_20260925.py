"""Stage, validate and submit the nine archived sparse CCD targets. Slurm only."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

BASE = Path('/usr/xtmp/lz280')
REPO = Path('/home/users/lz280/OSPREY3-fresh-packstar')
JAVA = '/home/users/lz280/java/jdk-17.0.2+8/bin/'
root = Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=True)
sources = root / 'source'
sources.mkdir(exist_ok=True)
for name in ['sparse_nine_20260925.py', 'sparse_nine_20260925.slurm', 'run_sparse_nine_20260925.slurm', 'merge_million_census_20260918.py']:
    shutil.copyfile(REPO / 'slurm/scripts' / name, sources / name)

million = BASE / 'packstar_million_census_20260918/build_12633933'
slow = BASE / 'packstar_slow_three_z_20260918/build_12634945'
scale = BASE / 'packstar_scale_validation_20260920/prep_12658935'
cases = [
    ('4wyu', million, BASE/'packstar_million_coverage_20260918/run_12634011'),
    ('1gwc', million, BASE/'packstar_million_coverage_20260918/1gwc_12634382'),
]
cases += [(s, slow, BASE/f'packstar_slow_three_z_20260918/{s}/seeds') for s in ['1a0r','4wwi','4kt6']]
cases += [(s, scale, BASE/f'packstar_scale_validation_20260920/{s}/packstar_100seeds') for s in ['2rfe','2xxm','4wyq','2xgy']]

# Preserve the archived implementation and energy gauge used for the 900 intervals.
src = (million/'source/PackStarMillionCensus.java').read_text()
src = src.replace('PackStarMillionCensus', 'PackStarSparseNineCensus')
helper = '''
    static ResidueInteractions sparse(ConfEnergyCalculator cc, int[] conf, int[][] edges) {
        ResidueInteractions inters=new ResidueInteractions();
        for(int p=0;p<conf.length;p++)inters.addAll(cc.makeSingleInters(p,conf[p]));
        for(int[] e:edges)inters.addAll(cc.makePairInters(e[0],conf[e[0]],e[1],conf[e[1]]));
        if(cc.addShellInters)inters.addAll(cc.makeShellInters());
        return inters;
    }
'''
src = src.replace('    public static void main(', helper + '    public static void main(')
needle = '        int shard=Integer.parseInt(args[4]), block=Integer.parseInt(args[5]);'
prepare = '''
        int[][] edges=Files.readAllLines(build.resolve("edges.tsv")).stream()
            .map(line->Arrays.stream(line.split("\\t")).mapToInt(Integer::parseInt).toArray()).toArray(int[][]::new);
        if(!Files.readString(build.resolve("original_layout.tsv")).equals(Files.readString(out.resolve("layout.tsv"))))
            throw new IllegalStateException("Archived layout mismatch");
        if(mode.equals("prepare")) {
            List<String> rows=Files.readAllLines(build.resolve("samples.tsv"));
            List<String> header=Arrays.asList(rows.get(0).split("\\t"));
            int ci=header.indexOf("conf"), ei=header.indexOf("eTrueKcal");
            if(ci<0||ei<0)throw new IllegalStateException("Sample schema mismatch");
            double max=0; int checked=0; Set<String> seen=new HashSet<>();
            List<String> audit=new ArrayList<>(); audit.add("conf\\toriginal\\trecomputed\\terror");
            try(EnergyCalculator ec=new EnergyCalculator.Builder(s,new ForcefieldParams())
                    .setType(EnergyCalculator.Type.Cpu).setIsMinimizing(true)
                    .setParallelism(Parallelism.makeCpu(4)).build()) {
                ConfEnergyCalculator cc=new ConfEnergyCalculator.Builder(s,ec).build();
                for(String row:rows.subList(1,rows.size())) {
                    String[] fields=row.split("\\t"); if(!seen.add(fields[ci]))continue;
                    int[] a=Arrays.stream(fields[ci].split(",")).mapToInt(Integer::parseInt).toArray();
                    double expectedE=Double.parseDouble(fields[ei]);
                    double actual=cc.calcEnergy(new RCTuple(a),sparse(cc,a,edges)).energy;
                    double error=Math.abs(actual-expectedE);
                    if(!Double.isFinite(error)||error>1e-7)throw new IllegalStateException("Archived sparse energy mismatch: "+error);
                    max=Math.max(max,error); checked++;
                    audit.add(fields[ci]+"\\t"+expectedE+"\\t"+actual+"\\t"+error);
                    if(checked>=64)break;
                }
            }
            if(checked==0)throw new IllegalStateException("No sample checks");
            Files.write(out.resolve("energy_check.tsv"),audit);
            Files.writeString(build.resolve("MODEL_READY"),"checked="+checked+" max_error="+max+"\\n");
            return;
        }
        if(!Files.exists(build.resolve("MODEL_READY")))throw new IllegalStateException("Unverified target");
'''
assert needle in src
src = src.replace(needle, prepare + needle)
src = src.replace('cc.calcEnergyAsync(new RCTuple(decode(id,rcs)),e->{',
                  'cc.calcEnergyAsync(new RCTuple(decode(id,rcs)),sparse(cc,decode(id,rcs),edges),e->{')
src = src.replace('cc.calcEnergy(new RCTuple(decode(row.getKey(),rcs))).energy',
                  'cc.calcEnergy(new RCTuple(decode(row.getKey(),rcs)),sparse(cc,decode(row.getKey(),rcs),edges)).energy')
assert 'sparse(cc,decode(id,rcs),edges)' in src
(sources/'PackStarSparseNineCensus.java').write_text(src)
(root/'classes').mkdir(exist_ok=True)
cp = (million/'classpath.txt').read_text().strip()
subprocess.run([JAVA+'javac','-J-Xmx2g','-cp',cp,'-d',str(root/'classes'),str(sources/'PackStarSparseNineCensus.java')], check=True)
(root/'classpath.txt').write_text(str(root/'classes')+':'+cp+'\n')
manifest=[]
for system, build, seeds in cases:
    model=root/'models'/system
    (model/'input').mkdir(parents=True)
    config=build/f'screen/{system}.tsv'
    row=next(csv.DictReader(config.open(), delimiter='\t'))
    n=int(row['positions'])
    seed0=seeds/'seed_10000'
    expected_edges=None
    for seed in range(10000,10100):
        seedpath=seeds/f'seed_{seed}'
        assert (seedpath/'config.tsv').read_bytes()==config.read_bytes()
        assert (seedpath/'layout.tsv').read_bytes()==(seed0/'layout.tsv').read_bytes()
        log=(seedpath/'run.log').read_text()
        assert 'energyMode=SPARSE' in log and 'residualBudget=1.0' in log and 'keepConnected=true' in log
        removed={tuple(sorted(map(int,m))) for m in re.findall(r'Cutting \((\d+),(\d+)\): residual',log)}
        assert removed, (system,seed,'No recorded cuts')
        edges={(i,j) for i in range(n) for j in range(i+1,n)}-removed
        if expected_edges is None: expected_edges=edges
        assert expected_edges==edges, (system,seed,'Graph changed across seeds')
    shutil.copyfile(config,model/'config.tsv')
    shutil.copyfile(build/f'input/{system}.pdb',model/f'input/{system}.pdb')
    shutil.copyfile(seed0/'layout.tsv',model/'original_layout.tsv')
    samples=list((seed0/'audit').glob('*/frequency_severity_final_samples.tsv'))
    assert len(samples)==1
    shutil.copyfile(samples[0],model/'samples.tsv')
    (model/'edges.tsv').write_text(''.join(f'{i}\t{j}\n' for i,j in sorted(expected_edges)))
    (model/'provenance.json').write_text(json.dumps(dict(build=str(build),seeds=str(seeds),sample=str(samples[0]),
        graph='All pairs minus archived RB=1 Cutting entries; identical across 100 seeds',edges=sorted(expected_edges)),indent=2))
    with (model/'prepare.log').open('w') as log:
        subprocess.run([JAVA+'java','-Xmx2g','-XX:-UseSuperWord','-XX:ActiveProcessorCount=4',
            '--add-modules','jdk.incubator.foreign','--add-opens','java.base/java.util=ALL-UNNAMED',
            '--add-opens','java.base/java.lang=ALL-UNNAMED','--add-opens','java.base/java.lang.invoke=ALL-UNNAMED',
            '-cp',(root/'classpath.txt').read_text().strip(),'edu.duke.cs.osprey.packstar.PackStarSparseNineCensus',
            'prepare',str(model),str(model/'check'),str(model/'config.tsv')],stdout=log,stderr=subprocess.STDOUT,check=True)
    manifest.append(dict(system=system,count=int(row['count']),model=str(model),seeds=str(seeds),block=2048 if system=='2rfe' else 8192))
    print(system,(model/'MODEL_READY').read_text().strip(),flush=True)
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
(root/'input.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p)+'\n'
    for p in sorted((root/'models').rglob('*')) if p.is_file()))
(root/'READY').touch()

def submit(args):
    return subprocess.check_output(['sbatch','--parsable','--account=grisman',*args],text=True).strip().split(';')[0]

jobs=[]
for case in manifest:
    system=case['system']; block=case['block']; shards=(case['count']+block-1)//block
    result=root/'results'/system
    result.mkdir(parents=True)
    env=f'ALL,BUILD_ROOT={case["model"]},RESULT_ROOT={result},RUN_ROOT={root},BLOCK={block},SYSTEM={system}'
    jid=submit(['--partition=grisman,compsci',f'--array=0-{shards-1}',f'--job-name=sp9_{system}',
                '--export='+env,str(sources/'run_sparse_nine_20260925.slurm')])
    case=dict(case,array=jid,shards=shards,result=str(result))
    jobs.append(case)
    (root/'jobs.json').write_text(json.dumps(jobs,indent=2))
    merge=root/f'merge_{system}.slurm'
    merge.write_text('#!/bin/bash\nset -euo pipefail\n'+
        f'python3 {sources}/merge_million_census_20260918.py {result} {jid} {case["model"]}/config.tsv {block}\n')
    case['merge']=submit(['--partition=grisman,compsci','--cpus-per-task=1','--mem=2G','--time=01:00:00',
        f'--job-name=sp9_merge_{system}',f'--dependency=afterok:{jid}',
        '--output=/usr/xtmp/lz280/slurm_logs/%x_%j.out','--error=/usr/xtmp/lz280/slurm_logs/%x_%j.err',str(merge)])
    (root/'jobs.json').write_text(json.dumps(jobs,indent=2))
    print('SUBMITTED',json.dumps(case),flush=True)
(root/'SUBMITTED').touch()
