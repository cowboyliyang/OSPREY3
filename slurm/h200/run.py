"""Portable Slurm entry point for the fixed PACK* experiment cohorts."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Submit through Slurm')
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', choices=['frontier', 'baseline38'], default='frontier')
    parser.add_argument('--frontier-manifest', type=Path,
                        help='Separate frozen add-on manifest; leaves the original build/cohort unchanged')
    parser.add_argument('--index', type=int, default=int(os.environ.get('SLURM_ARRAY_TASK_ID', '0')))
    parser.add_argument('--mode', choices=['preflight', 'full', 'pfunc'], default='full')
    parser.add_argument('--arm', choices=['budget-forward', 'decomposition-cost', 'pair-only', 'no-learning'], default='budget-forward')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--rb', type=int, choices=[1, 2, 4],
                        help='Override residual budget for the RB sweep')
    parser.add_argument('--gpus', type=int, default=2)
    parser.add_argument('--heap-gib', type=int, default=850)
    parser.add_argument('--host-gib', type=int, default=800)
    parser.add_argument('--gpu-gib', type=int, default=120)
    args = parser.parse_args()
    build = Path(os.environ['BUILD_ROOT']).resolve()
    data = Path(os.environ['INPUT_ROOT']).resolve()
    config = build / 'source/slurm/h200'
    if args.frontier_manifest and args.cohort != 'frontier':
        raise ValueError('--frontier-manifest requires the frontier cohort')
    if args.cohort == 'frontier':
        manifest_path = args.frontier_manifest or config / 'frontier.tsv'
        rows = list(csv.DictReader(manifest_path.open(), delimiter='\t'))
        package_rows = list(csv.DictReader((data / 'designs.tsv').open(), delimiter='\t'))
        if rows != package_rows:
            raise ValueError('Transferred frontier manifest differs from the frozen configuration')
    else:
        with (config / 'baseline38.csv').open() as stream:
            records = list(csv.reader(line for line in stream if not line.startswith('#')))
        rows = [dict(design_id=r[0], system=r[1], mutable=r[5], flexible=r[6],
                     expected_sequences=1 + 19 * len(r[5].split(';')),
                     pdb_relative=f'{r[1]}/{r[1]}.min.reduce.renum.pdb') for r in records]
    if not 0 <= args.index < len(rows):
        raise ValueError('Design index outside cohort')
    row = rows[args.index]
    pdb = data / row['pdb_relative']
    checksum = hashlib.sha256(pdb.read_bytes()).hexdigest()
    if row.get('pdb_sha256') and checksum != row['pdb_sha256']:
        raise ValueError('PDB checksum mismatch')
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    if min(cpus, args.gpus, args.heap_gib, args.host_gib, args.gpu_gib) <= 0:
        raise ValueError('Resource settings must be positive')
    if args.host_gib >= args.heap_gib:
        raise ValueError('Heap must leave room above the DP host budget')
    reserved_mib = int(os.environ.get('SLURM_MEM_PER_NODE', '0'))
    if reserved_mib and args.heap_gib * 1024 >= reserved_mib:
        raise ValueError('Heap must leave native/OS memory within the Slurm allocation')
    # Preserve the approved /usr/xtmp spelling required by Duke production builds.
    out = Path(os.environ['RESULT_ROOT']).absolute() / (
        f"{row['design_id']}_{args.mode}_{args.arm}_s{args.seed}_"
        f"J{os.environ['SLURM_JOB_ID']}_T{args.index}")
    # Refuse to reuse caches or overwrite a completed/partial experiment.
    out.mkdir(parents=True, exist_ok=False)
    for directory in ['tmp', 'dp_mmap', 'cuda_cache']:
        (out / directory).mkdir()
    env = dict(os.environ, TMPDIR=str(out / 'tmp'), CUDA_CACHE_PATH=str(out / 'cuda_cache'))
    properties = {}
    for line in (config / 'production.properties').read_text().splitlines():
        if line and not line.startswith('#'):
            key, value = line.split('=', 1)
            properties[key] = value
    properties.update({
        'java.io.tmpdir': str(out / 'tmp'),
        'branchdp.dp.mmap.dir': str(out / 'dp_mmap'),
        'branchdp.dp.mmap.thresholdBytes': str(min(512, args.host_gib // 2) * 2**30),
        'branchdp.rootSplit.hostBudgetBytes': str(args.host_gib * 2**30),
        'branchdp.rootSplit.gpuBudgetBytes': str(args.gpu_gib * 2**30),
        'branchdp.dp.gpu.maxGpus': str(args.gpus),
        'branchdp.dp.parallel.threads': str(cpus),
        'packstar.dp.parallel.threads': str(cpus),
        'packstar.pac.sampling.threads': str(cpus),
        'packstar.pac.randomSeed': str(args.seed),
        'packstar.pac.frequencySeverity.tripleEtaMaxHostBytes': f'{args.host_gib}GiB',
        'packstar.pac.frequencySeverity.tripleEtaMaxTotalTableBytes': f'{min(768, args.host_gib)}GiB',
        'packstar.pac.frequencySeverity.outputDir': str(out / 'adaptive_frequency_severity'),
        'osprey.bench.method': {'preflight': 'sequence_dump', 'full': 'packstar', 'pfunc': 'packstar_pfunc'}[args.mode],
        'osprey.bench.numCPUs': str(cpus),
        'osprey.bench.designId': row['design_id'],
        'osprey.bench.outputDir': str(out),
        'osprey.bench.pdbPath': str(pdb),
        'osprey.bench.mutable': row['mutable'],
        'osprey.bench.flexible': row['flexible'],
        'osprey.packstarPfunc.state': 'complex',
        'osprey.packstarPfunc.seqIndex': '0',
    })
    prefix = 'packstar.pac.frequencySeverity.'
    if args.rb is not None:
        properties['branchdp.cutoff.residualBudget'] = str(args.rb)
    if args.arm == 'decomposition-cost':
        properties[prefix + 'tripleEtaSelectionStrategy'] = args.arm
    elif args.arm in ('pair-only', 'no-learning'):
        properties[prefix + 'tripleEta'] = 'false'
        if args.arm == 'no-learning':
            properties[prefix + 'proposalLearning'] = 'false'
    java = str(Path(os.environ['JAVA_HOME']) / 'bin/java')
    command = [java, '--add-modules', 'jdk.incubator.foreign',
               '--add-opens', 'java.base/java.util=ALL-UNNAMED',
               '--add-opens', 'java.base/java.lang=ALL-UNNAMED',
               '--add-opens', 'java.base/java.lang.invoke=ALL-UNNAMED',
               f'-Xmx{args.heap_gib}g', '-Xms1g', '-XX:-UseSuperWord',
               f'-XX:ActiveProcessorCount={cpus}']
    command += [f'-D{k}={v}' for k, v in sorted(properties.items())]
    command += ['-cp', (build / 'test_classpath.txt').read_text().strip(),
                'edu.duke.cs.osprey.markstar.bench.GenericPDBBench']
    options = vars(args).copy()
    options['frontier_manifest'] = str(args.frontier_manifest) if args.frontier_manifest else None
    manifest = dict(design=row, pdb_sha256=checksum, options=options,
                    properties=properties, cpus=cpus, node=os.environ.get('SLURMD_NODENAME'),
                    slurm_job_id=os.environ['SLURM_JOB_ID'], build_root=str(build),
                    git_head=(build / 'git_head.txt').read_text().strip(),
                    status='RUNNING', start_epoch=time.time(), energy_matrices='fresh',
                    precision='FP64', ccd='CPU')
    manifest['runner_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.cohort == 'frontier':
        manifest['frontier_manifest_sha256'] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    shutil.copy2(build / 'source.sha256', out / 'build-source.sha256')
    (out / 'command.sh').write_text(shlex.join(command) + '\n')
    for name, cmd in [('slurm.txt', ['scontrol', 'show', 'job', os.environ['SLURM_JOB_ID']]),
                      ('lscpu.txt', ['lscpu'])]:
        with (out / name).open('w') as stream:
            subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, check=True)
    if args.mode != 'preflight':
        # Query only allocated devices: nvidia-smi by itself can expose other jobs' GPUs.
        visible = os.environ.get('CUDA_VISIBLE_DEVICES', '')
        devices = visible.split(',') if visible else []
        if len(devices) != args.gpus or any(not d or d == '-1' for d in devices):
            raise ValueError('CUDA_VISIBLE_DEVICES must match --gpus')
        info = subprocess.check_output(['nvidia-smi', '-i', visible,
            '--query-gpu=uuid,name,memory.total,driver_version', '--format=csv,noheader'], text=True)
        (out / 'gpus.csv').write_text(info)
        expected_gpu = os.environ.get('PACKSTAR_EXPECT_GPU', 'H200')
        for device in csv.reader(info.splitlines()):
            if expected_gpu not in device[1] or args.gpu_gib * 1024 >= int(device[2].strip().split()[0]):
                raise ValueError('Unexpected GPU model or excessive per-GPU budget')
    print('RUN_START', out, flush=True)
    with (out / 'run.log').open('w') as stream:
        result = subprocess.run(['/usr/bin/time', '-o', str(out / 'wall.time'),
            '-f', 'elapsed=%e maxRssKiB=%M', *command], cwd=out, env=env,
            stdout=stream, stderr=subprocess.STDOUT)
    manifest.update(end_epoch=time.time(), exit_code=result.returncode,
                    status='PROCESS_FAILED' if result.returncode else 'COMPLETED')
    suffix = {'preflight': '_sequences.tsv', 'full': '_packstar.csv', 'pfunc': '_packstar_pfunc.tsv'}[args.mode]
    result_file = out / (row['design_id'] + suffix)
    results = list(csv.DictReader(result_file.open(), delimiter=',' if args.mode == 'full' else '\t')) if result_file.exists() else []
    expected = 1 if args.mode == 'pfunc' else int(row['expected_sequences'])
    manifest.update(result_rows=len(results), expected_rows=expected)
    if not result.returncode and len(results) != expected:
        manifest['status'] = 'INCOMPLETE'
    elif not result.returncode and args.mode == 'preflight':
        log = (out / 'run.log').read_text()
        if 'WARNING: flexible residue' in log or 'WARNING: mutable residue' in log:
            manifest['status'] = 'INVALID_POSITIONS'
        if args.cohort == 'frontier':
            frozen = data / 'preflight' / row['design_id'] / result_file.name
            reference = list(csv.DictReader(frozen.open(), delimiter='\t'))
            if results != reference:
                manifest['status'] = 'SEQUENCE_MISMATCH'
    elif not result.returncode:
        statuses = ['prot_status', 'lig_status', 'comp_status'] if args.mode == 'full' else ['status']
        manifest['estimated_rows'] = sum(all(r[s] == 'Estimated' for s in statuses) for r in results)
        if manifest['estimated_rows'] != expected:
            manifest['status'] = 'INCOMPLETE_ESTIMATES'
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('RUN_END', manifest['status'], out, flush=True)
    return result.returncode or (0 if manifest['status'] == 'COMPLETED' else 5)


if __name__ == '__main__':
    raise SystemExit(main())
