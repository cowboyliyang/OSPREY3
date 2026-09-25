import csv, hashlib, json, math, re, statistics, sys
from collections import Counter
from pathlib import Path
root=Path(sys.argv[1]); out=Path(sys.argv[2]); out.mkdir(exist_ok=True)
rows=[]
for seed in range(10000,10100):
    d=root/f'seed_{seed}'
    with (d/'result.tsv').open() as f:r=next(csv.DictReader(f,delimiter='\t'))
    a=next((d/'audit').glob('*/frequency_severity_final_interval.tsv'))
    with a.open() as f:kv={x['key']:x['value'] for x in csv.DictReader(f,delimiter='\t')}
    log=(d/'run.log').read_text()
    wall=(d/'wall.time').read_text()
    r.update({k:kv[k] for k in ('selectedCandidate','certificateValid','assumptionConditional','targetReached','tailCount','finalCcd','totalCcd','totalSampleRecords','selectedTripleEtaPositionTriples','epsilon','proposalDpSweeps')})
    r['seconds']=float(re.search(r'elapsed=([\d.]+)',wall)[1])
    r['rssMiB']=int(re.search(r'maxRssKiB=(\d+)',wall)[1])/1024
    r['violations']=sum(map(int,re.findall(r'violations=(\d+)',log)))
    r['deterministic']='sampling skipped; first-round DP bounds' in log
    r['sampleSha']=hashlib.sha256((a.parent/'frequency_severity_final_samples.tsv').read_bytes()).hexdigest()
    rows.append(r)
def stats(key):
    v=[float(r[key]) for r in rows]
    return dict(min=min(v),median=statistics.median(v),max=max(v))
report={'count':len(rows),'statuses':dict(Counter(r['status'] for r in rows)),
        'candidates':dict(Counter(r['selectedCandidate'] for r in rows)),
        'certificateValid':dict(Counter(r['certificateValid'] for r in rows)),
        'targetReached':dict(Counter(r['targetReached'] for r in rows)),
        'unique_final_sample_files':len({r['sampleSha'] for r in rows}),
        'deterministic_runs':sum(r['deterministic'] for r in rows),
        'lower_bound_violations':sum(r['violations'] for r in rows),
        'stats':{k:stats(k) for k in ('logLower','logUpper','epsilon','finalCcd','totalCcd','totalSampleRecords','seconds','rssMiB','tailCount','selectedTripleEtaPositionTriples','proposalDpSweeps')},
        'common_log_interval':[max(float(r['logLower']) for r in rows),min(float(r['logUpper']) for r in rows)],
        'reference_available':False}
report['common_Z_interval']=[math.exp(v) for v in report['common_log_interval']]
report['upper_lower_factor']=dict(min=min(math.exp(float(r['logUpper'])-float(r['logLower'])) for r in rows),max=max(math.exp(float(r['logUpper'])-float(r['logLower'])) for r in rows))
(out/'inspection.json').write_text(json.dumps(report,indent=2))
(out/'rows.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(report,indent=2))
