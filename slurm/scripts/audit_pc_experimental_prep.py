#!/usr/bin/env python3
"""Slurm structural QC; never infer populations by counting PDB models."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np


def atoms(path):
    result = {}
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM  '):
            continue
        element = line[76:78].strip()
        if element in ('H', 'D') or line[12:16].strip().startswith('H'):
            continue
        key = (line[21], int(line[22:26]), line[12:16].strip())
        if key in result:
            raise ValueError('Duplicate atom ' + str(key))
        result[key] = np.array([float(line[i:i+8]) for i in (30,38,46)])
    return result


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Slurm required')
    out = Path(os.environ['PC_AUDIT_OUT'])
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for root_text in sys.argv[1:]:
        root = Path(root_text)
        with (root / 'prep_validation.tsv').open() as f:
            checks = list(csv.DictReader(f, delimiter='\t'))
        for check in checks:
            key = check['case']
            if check['status'] != 'PREPARED_PROTEIN_ONLY':
                rows.append(dict(run=root.name, case=key, status=check['status']))
                continue
            directory = root / 'pdbs_prepped' / key
            final = directory / (key + '.min.reduce.renum.pdb')
            raw_atoms = atoms(root / 'pdbs_raw' / (key + '.pdb'))
            ready = atoms(final)
            mapping = {}
            for line in (directory / (key + '_strip_reduce_prep_renum.txt')).read_text().splitlines():
                p = line.split()
                if len(p) >= 5:
                    mapping[(p[1], int(p[2]))] = (p[1], int(p[4]))
            pairs = [(raw_atoms[a], ready[mapping[a[:2]] + (a[2],)]) for a in raw_atoms
                     if a[:2] in mapping and mapping[a[:2]] + (a[2],) in ready and a[2] == 'CA']
            x = np.array([a for a,b in pairs]); y = np.array([b for a,b in pairs])
            x -= x.mean(0); y -= y.mean(0)
            u, _, vt = np.linalg.svd(x.T @ y)
            correction = np.eye(3); correction[2,2] = np.linalg.det(u @ vt)
            rmsd = float(np.sqrt(np.mean(np.sum((x @ u @ correction @ vt - y)**2, axis=1))))
            missing = [str(a) for a in raw_atoms if a[:2] not in mapping or mapping[a[:2]] + (a[2],) not in ready]
            row = dict(run=root.name, case=key, status='QC_PASS' if not missing else 'MISSING_HEAVY_ATOMS',
                       source_heavy_atoms=len(raw_atoms), prepared_heavy_atoms=len(ready),
                       matched_CA=len(pairs), CA_rmsd_A=rmsd, missing=';'.join(missing),
                       sha256=hashlib.sha256(final.read_bytes()).hexdigest())
            rows.append(row)
    fields = ('run','case','status','source_heavy_atoms','prepared_heavy_atoms','matched_CA','CA_rmsd_A','missing','sha256')
    with (out / 'structural_qc.tsv').open('w') as f:
        writer=csv.DictWriter(f, fieldnames=fields, delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    (out / 'summary.json').write_text(json.dumps({
        'models_checked':len(rows), 'qc_pass':sum(r['status']=='QC_PASS' for r in rows),
        'not_passed':[{k:r[k] for k in ('run','case','status')} for r in rows if r['status']!='QC_PASS'],
        'maximum_CA_rmsd_A':max(r.get('CA_rmsd_A',0) for r in rows),
        'population_inferred_from_model_count':False,
        '7QCY_prepared_scope':'96-residue protein domain; ligand coordinates must be checked/obtained separately before any complex calculation'
    }, indent=2)+'\n')


if __name__ == '__main__':
    main()
