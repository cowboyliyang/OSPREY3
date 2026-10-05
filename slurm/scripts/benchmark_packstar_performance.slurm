#!/bin/bash
#SBATCH --job-name=ps-performance-bench
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --constraint=a5000
#SBATCH --gres=gpu:a5000:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
CANDIDATE_BUILD=${1:?pass candidate validation directory}
BASELINE_BUILD=${2:-/usr/xtmp/lz280/packstar_push_20260929/check_12763717}
REPO=/home/users/lz280/OSPREY3-fresh-packstar
RUN_ROOT=/usr/xtmp/lz280/packstar_performance_20261003/bench_${SLURM_JOB_ID}
# Two tiny source harnesses; diagnostics expected <100 MiB / 1500 files. No downloads.
mkdir -p "$RUN_ROOT/harness" "$RUN_ROOT/tmp"
export TMPDIR="$RUN_ROOT/tmp"
export JAVA_HOME=/home/users/lz280/java/jdk-17.0.2+8
export PATH="$JAVA_HOME/bin:$PATH"
BASE_CP=$(<"$BASELINE_BUILD/test_classpath.txt")
if [[ -f "$CANDIDATE_BUILD/test_classpath.txt" ]]; then
    CANDIDATE_CP=$(<"$CANDIDATE_BUILD/test_classpath.txt")
else
    # Same dependency graph: shadow only Java outputs if a pre-existing test
    # failure prevented the legacy test-dependent classpath-export task.
    CANDIDATE_CP="$CANDIDATE_BUILD/build/osprey3-root/classes/java/test:$CANDIDATE_BUILD/build/osprey3-root/classes/java/main:$BASE_CP"
fi
cp "$REPO/src/test/java/edu/duke/cs/osprey/packstar/BenchPackStarPerformance.java" "$RUN_ROOT/"
cp "$REPO/src/test/java/edu/duke/cs/osprey/branchdp/BenchPackStarGpuPreparation.java" "$RUN_ROOT/"
javac -cp "$BASE_CP" -d "$RUN_ROOT/harness" "$RUN_ROOT/BenchPackStarPerformance.java" "$RUN_ROOT/BenchPackStarGpuPreparation.java"
JVM=(--add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED
     --add-opens java.base/java.lang.invoke=ALL-UNNAMED --add-modules jdk.incubator.foreign
     -Djava.io.tmpdir="$TMPDIR" -Xmx8g)
cd "$BASELINE_BUILD/source"
printf 'NODE=%s BASELINE=%s CANDIDATE=%s\n' "${SLURMD_NODENAME}" "$BASELINE_BUILD" "$CANDIDATE_BUILD"
for positions in 2 4; do
    for seed in 42 43 44; do
        modes=(baseline optimized)
        if (( seed % 2 )); then modes=(optimized baseline); fi
        for mode in "${modes[@]}"; do
            CP=$BASE_CP
            if [[ "$mode" == optimized ]]; then CP=$CANDIDATE_CP; fi
            out="$RUN_ROOT/cpu_${positions}_${seed}_${mode}"
            mkdir -p "$out"
            /usr/bin/time -v -o "$out/resources.log" java "${JVM[@]}" -cp "$RUN_ROOT/harness:$CP" \
                edu.duke.cs.osprey.packstar.BenchPackStarPerformance "$out" "$mode" "$positions" "$seed" \
                > "$out/run.log" 2>&1
            tail -n 1 "$out/result.tsv"
        done
    done
done
for mode in baseline optimized; do
    CP=$BASE_CP
    if [[ "$mode" == optimized ]]; then CP=$CANDIDATE_CP; fi
    out="$RUN_ROOT/gpu_${mode}"
    mkdir -p "$out"
    java "${JVM[@]}" -cp "$RUN_ROOT/harness:$CP" \
        edu.duke.cs.osprey.branchdp.BenchPackStarGpuPreparation "$out" "$mode" > "$out/run.log" 2>&1
done
python3 - "$RUN_ROOT" <<'PY'
import csv, json, pathlib, re, statistics, sys
root = pathlib.Path(sys.argv[1])
def rows(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter='\t'))
checks, cpu, gpu = [], [], []
for positions in (2, 4):
    speeds, calls, memory = [], [], []
    for seed in (42, 43, 44):
        old, new = [root / f'cpu_{positions}_{seed}_{mode}' for mode in ('baseline', 'optimized')]
        a, b = rows(old/'result.tsv')[0], rows(new/'result.tsv')[0]
        for key in ('status', 'lower', 'upper', 'epsilon'):
            assert a[key] == b[key], (positions, seed, key, a[key], b[key])
        old_samples = {p.relative_to(old/'measured'): p for p in (old/'measured').rglob('*samples.tsv')}
        new_samples = {p.relative_to(new/'measured'): p for p in (new/'measured').rglob('*samples.tsv')}
        assert old_samples and old_samples.keys() == new_samples.keys(), (positions, seed, 'sample files')
        for name, path in old_samples.items():
            assert path.read_bytes() == new_samples[name].read_bytes(), (positions, seed, str(name))
        for stage in ('eta_selected.tsv', 'frequency_severity_final_interval.tsv'):
            left, right = list((old/'measured').rglob(stage)), list((new/'measured').rglob(stage))
            assert len(left) == len(right)
            if left:
                def normalized(path):
                    return [line for line in path.read_text().splitlines() if not line.startswith('totalCcd\t')]
                assert normalized(left[0]) == normalized(right[0]), (positions, seed, stage)
        def ccd(path):
            found = re.findall(r'(\d+) unique CCD calls for (\d+) sample records', (path/'run.log').read_text())
            assert found, str(path)
            return tuple(map(int, found[-1]))
        ca, cb = ccd(old), ccd(new)
        assert ca[1] == cb[1]
        speeds.append(int(a['nanos']) / int(b['nanos']))
        calls.append([ca[0], cb[0]])
        memory.append([int(a['peakHeapBytes']), int(b['peakHeapBytes'])])
        checks.append({'positions': positions, 'seed': seed, 'status': a['status'],
                       'identicalSampleFiles': len(old_samples), 'logicalSamples': ca[1]})
    cpu.append({'positions': positions, 'pairedSpeedups': speeds,
                'medianPairedSpeedup': statistics.median(speeds),
                'ccdCallsOldNew': calls, 'peakHeapBytesOldNew': memory})
old, new = root/'gpu_baseline', root/'gpu_optimized'
for path in old.glob('*.txt'):
    assert path.read_bytes() == (new/path.name).read_bytes(), path.name
for method in ('gumbel', 'cdfDedup'):
    a, b = [[r for r in rows(p/'timings.tsv') if r['method'] == method] for p in (old, new)]
    assert len(a) == len(b) == 10
    timing = {'method': method, 'identicalBatches': len(a)}
    for field in ('prepareNanos', 'sampleNanos'):
        before, after = [statistics.median(int(r[field]) for r in values) for values in (a, b)]
        timing[field] = {'before': before, 'after': after, 'speedup': before/after}
    gpu.append(timing)
report = {'checks': checks, 'cpu': cpu, 'gpu': gpu}
(root/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
PY
echo "PERFORMANCE_BENCHMARK=$RUN_ROOT"
