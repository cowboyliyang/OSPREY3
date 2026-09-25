#!/bin/bash
#SBATCH --job-name=ps_slow3_launch
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
ROOT=/usr/xtmp/lz280/packstar_slow_three_z_20260918
BUILD_ROOT=$ROOT/build_12634945
test -f "$BUILD_ROOT/READY"
mkdir "$ROOT/submission"
cp "${BASH_SOURCE[0]}" "$ROOT/submission/launcher.slurm"
export BUILD_ROOT
python3 - "$ROOT" <<'PY'
import csv, json, os, subprocess, sys
from pathlib import Path
root=Path(sys.argv[1]); build=Path(os.environ['BUILD_ROOT']); scripts=build/'source'
manifest=root/'submission/jobs.jsonl'
def submit(system, role, args, script, env):
    command=['sbatch','--parsable',*args,'--export=ALL,'+','.join(f'{k}={v}' for k,v in env.items()),str(scripts/script)]
    result=subprocess.run(command,check=True,text=True,capture_output=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit(): raise RuntimeError(result.stdout)
    record=dict(system=system,role=role,job=job,command=command)
    with manifest.open('a') as f: f.write(json.dumps(record)+'\n'); f.flush(); os.fsync(f.fileno())
    print(system,role,job,flush=True)
    return job
for system in ('1a0r','4wwi','4kt6'):
    with (build/'screen'/f'{system}.tsv').open() as f: config=next(csv.DictReader(f,delimiter='\t'))
    count=int(config['count'])
    assert 0<count<=2000000
    census=root/system/'census'; seeds=root/system/'seeds'
    census.mkdir(parents=True); seeds.mkdir()
    block=2048; tasks=(count+block-1)//block
    common=dict(BUILD_ROOT=build,SYSTEM=system)
    ccd=submit(system,'ccd',['--array=0-'+str(tasks-1),'--time=04:00:00','--job-name=zccd_'+system],
        'run_million_census_20260918.slurm',dict(common,RESULT_ROOT=census,BLOCK=block,MODE='shard'))
    gpuenv=dict(common,RESULT_ROOT=seeds,CENSUS_ROOT=build)
    fenn=submit(system,'fennario',['--partition=grisman','--account=grisman','--nodelist=fennario-[03-06]',
        '--array=0-31','--job-name=zseed_'+system], 'run_million_coverage_20260918.slurm',gpuenv)
    comp=submit(system,'compsci_gpu',['--partition=compsci-gpu','--account=grisman','--array=32-99',
        '--job-name=zseed_'+system], 'run_million_coverage_20260918.slurm',gpuenv)
    merge=submit(system,'merge',['--dependency=afterany:'+ccd,'--job-name=zmerge_'+system],
        'merge_million_census_20260918.slurm',dict(common,RESULT_ROOT=census,BLOCK=block,ARRAY_ID=ccd))
    submit(system,'coverage',['--dependency=afterany:'+':'.join((fenn,comp,merge)),'--job-name=zcover_'+system],
        'summarize_million_coverage_20260918.slurm',dict(RESULT_ROOT=seeds,REFERENCE_ROOT=census/f'merged_A{ccd}'))
(root/'submission/READY').write_text('All three censuses, 300 seeds, merges and summaries submitted.\n')
PY
