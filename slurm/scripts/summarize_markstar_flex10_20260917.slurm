#!/bin/bash
#SBATCH --job-name=ms_flex10_summary
#SBATCH --account=grisman
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err
set -euo pipefail
: "${FRONTIER_JOB:?}"
: "${PREP_ROOT:?}"
ROOT=/usr/xtmp/lz280/packstar_flex_frontier10_20260917/A${FRONTIER_JOB}
mkdir -p "$ROOT/audits/$SLURM_JOB_ID"
export PYTHONPYCACHEPREFIX="$ROOT/audits/$SLURM_JOB_ID/pycache"
cp "${BASH_SOURCE[0]}" "$ROOT/audits/$SLURM_JOB_ID/summarize.slurm"
python3 - "$ROOT" "$PREP_ROOT" "$FRONTIER_JOB" "$SLURM_JOB_ID" <<'PY'
import csv,json,re,subprocess,sys
from pathlib import Path
root,prep=map(Path,sys.argv[1:3]); job,audit=sys.argv[3:]
target=root/'audits'/audit
command=['sacct','-j',job,'-X','-P','--format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS,Partition,NodeList']
raw=subprocess.check_output(command,text=True)
(target/'sacct.psv').write_text(raw)
states={r['JobID']:r for r in csv.DictReader(raw.splitlines(),delimiter='|')}
results=[]
for d in json.loads((prep/'package/designs.json').read_text()):
    run=root/f"T{d['task_id']}"; manifest=run/'run_manifest.json'
    metadata=json.loads(manifest.read_text()) if manifest.exists() else {}
    scheduler=states.get(f"{job}_{d['task_id']}",{})
    path=run/(d['design_id']+'_markstar.csv')
    rows=list(csv.DictReader(path.open())) if path.exists() else []
    w=run/'wall.time'; match=re.search(r'elapsed=([\d.]+)',w.read_text()) if w.exists() else None
    elapsed=float(match[1]) if match else None
    state=scheduler.get('State','UNKNOWN')
    partial=run/(d['design_id']+'_partial_scores.tsv')
    results.append(dict(design=d,run=str(run),scheduler=scheduler,manifest=metadata,
         wall_s=elapsed,complete=state=='COMPLETED' and len(rows)==d['expected_sequences'],
         right_censored=state=='TIMEOUT',rows=len(rows),
         all_estimated=sum(all(r[s+'_status']=='Estimated' for s in ('prot','lig','comp')) for r in rows),
         partial_score_lines=len(partial.read_text().splitlines()) if partial.exists() else 0))
(target/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
lines=['| Design | State | Partition | Node | Wall h | Success |', '|---|---|---|---|---:|---:|']
for r in results:
    d,s=r['design'],r['scheduler']; hours=f"{r['wall_s']/3600:.2f}" if r['wall_s'] is not None else ''
    lines.append(f"| {d['design_id']} | {s.get('State','UNKNOWN')} | {s.get('Partition','')} | {s.get('NodeList','')} | {hours} | {r['all_estimated']}/{r['rows']} |")
(target/'summary.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
PY
