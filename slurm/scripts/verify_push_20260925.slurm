#!/bin/bash
#SBATCH --job-name=packstar_push_check
#SBATCH --account=grisman
#SBATCH --partition=grisman,compsci
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
export PACKSTAR_REPO=/home/users/lz280/OSPREY3-fresh-packstar
export PACKSTAR_CHECK_ROOT=/usr/xtmp/lz280/packstar_push_check_20260925_$SLURM_JOB_ID
export PACKSTAR_JAVA=/home/users/lz280/java/jdk-17.0.2+8/bin/javac
export PACKSTAR_TEST_CP=/usr/xtmp/lz280/packstar_historical_triple_ablation_20260921/build_12677588/test_classpath.txt
# Source checks, copied sources, and compiled classes: <20 MiB, <300 files.
mkdir -p "$PACKSTAR_CHECK_ROOT/tmp" "$PACKSTAR_CHECK_ROOT/classes"
export TMPDIR="$PACKSTAR_CHECK_ROOT/tmp" PYTHONDONTWRITEBYTECODE=1
python3 - <<'PY'
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

repo = Path(os.environ['PACKSTAR_REPO'])
out = Path(os.environ['PACKSTAR_CHECK_ROOT'])
def git_paths(*args):
    return subprocess.check_output(['git', '-C', str(repo), *args, '-z']).decode().split('\0')
names = sorted((set(git_paths('diff', '--name-only', 'HEAD')) |
                set(git_paths('ls-files', '--others', '--exclude-standard'))) - {''})
counts = {'python': 0, 'shell': 0, 'json': 0, 'java': 0}
manifest = {}
java = []
for name in names:
    src = repo / name
    if not src.is_file():
        continue
    data = src.read_bytes()
    manifest[name] = hashlib.sha256(data).hexdigest()
    snapshot = out / 'source' / name
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, snapshot)
    if src.suffix == '.py':
        ast.parse(data.decode(), filename=name)
        counts['python'] += 1
    elif src.suffix in ('.sh', '.slurm'):
        subprocess.run(['bash', '-n', str(snapshot)], check=True)
        counts['shell'] += 1
    elif src.suffix == '.json':
        json.loads(data)
        counts['json'] += 1
    elif src.suffix == '.java':
        java.append(str(snapshot))
        counts['java'] += 1
classpath = Path(os.environ['PACKSTAR_TEST_CP']).read_text().strip()
subprocess.run([os.environ['PACKSTAR_JAVA'], '-J-Xmx4g', '-cp', classpath,
                '-sourcepath', '', '-d', str(out / 'classes'), *java], check=True)
result = {'checks': counts, 'source_sha256': manifest}
(out / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
print('PASS: syntax and Java compilation:', counts)
print('CHECK_ROOT=' + str(out))
PY
