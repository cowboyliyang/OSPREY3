#!/usr/bin/env python3
"""Slurm-only, isolated reuse of the 38-system prep; no population inference."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.request

PDB_IDS = tuple(os.environ.get('PC_PREP_PDBS', '1ake:4ake:7qcx:7qcy').split(':'))
LEGACY = Path('/usr/xtmp/lz280/dance_bench')
MAX_DOWNLOAD = 20 * 1024 * 1024


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, fields, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def download(url, path):
    request = urllib.request.Request(url, headers={'User-Agent': 'OSPREY-functional-validation/1.0'})
    with urllib.request.urlopen(request, timeout=90) as response:
        size = response.headers.get('Content-Length')
        if size and int(size) > MAX_DOWNLOAD:
            raise ValueError('Download exceeds 20 MiB budget: ' + url)
        data = response.read(MAX_DOWNLOAD + 1)
        if len(data) > MAX_DOWNLOAD:
            raise ValueError('Download exceeds 20 MiB budget: ' + url)
    path.write_bytes(data)


def models(path):
    result, current = [], []
    for line in path.read_text().splitlines(True):
        if line.startswith('MODEL'):
            if current:
                raise ValueError('Unexpected records before MODEL')
        elif line.startswith('ENDMDL'):
            result.append(current)
            current = []
        elif line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
            current.append(line)
    if current:
        result.append(current)
    if not result:
        raise ValueError('No PDB coordinates in ' + str(path))
    return result


def residue_list(lines):
    seen, result = set(), []
    for line in lines:
        if not line.startswith('ATOM  '):
            continue
        key = (line[21], line[22:27], line[17:20].strip())
        if key not in seen:
            result.append(key)
            seen.add(key)
    return result


def select_altlocs(lines):
    """Choose one complete available label per residue; never drop B/C-only residues."""
    grouped = {}
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            grouped.setdefault((line[21], line[22:27], line[17:20]), []).append(line)
    selected = {}
    for key, atoms in grouped.items():
        labels = sorted({l[16] for l in atoms if l[16] != ' '}) or [' ']
        def coverage(label):
            return len({l[12:16] for l in atoms if l[16] in (' ', label)
                        and l[76:78].strip() not in ('H', 'D')})
        selected[key] = min(labels, key=lambda label: (-coverage(label), label))
    output = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[21], line[22:27], line[17:20])
            if line[16] not in (' ', selected[key]):
                continue
            line = line[:16] + ' ' + line[17:]
        output.append(line)
    if residue_list(lines) != residue_list(output):
        raise ValueError('Alternate selection changed residue identity/order')
    return output, selected


def run(argv, logfile, cwd):
    with logfile.open('w') as log:
        proc = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
    return proc.returncode


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Must run through Slurm')
    root = Path(sys.argv[1]).absolute()
    allowed = Path('/usr/xtmp/lz280/packstar_functional_validation').resolve()
    if root.parent.resolve() != allowed or not root.name.startswith('experimental_'):
        raise ValueError('Output must be an isolated experimental run under xtmp')
    raw, split, prepped, scripts, logs = [root / name for name in
        ('raw', 'pdbs_raw', 'pdbs_prepped', 'source_snapshot', 'logs')]
    for directory in (raw, split, prepped, scripts, logs):
        directory.mkdir(parents=True, exist_ok=True)
    if (root / 'sources.tsv').exists():
        raise ValueError('Refusing to overwrite an existing run')
    for name in ('prep_all.sh', 'minimize_single.sh', 'PREP_PIPELINE.md'):
        shutil.copy2(LEGACY / name, scripts / name)
    shutil.copy2(__file__, scripts / 'prep_pc_experimental.py')
    sources = []
    cases = []
    for pdb in PDB_IDS:
        url = 'https://files.rcsb.org/download/' + pdb.upper() + '.pdb'
        path = raw / (pdb + '.pdb')
        download(url, path)
        sources.append(dict(pdb=pdb, url=url, bytes=path.stat().st_size, sha256=sha(path)))
        all_models = models(path)
        if pdb.startswith('7q') and len(all_models) != 40:
            raise ValueError('Expected 40 deposited models: ' + pdb)
        for i, lines in enumerate(all_models, 1):
            key = '{}_m{:03d}'.format(pdb, i)
            # Each deposited model remains a separate structure, never an occupancy sample.
            # Retain raw alternate occupancies. Choose by completeness/label, not population.
            selected, labels = select_altlocs(lines)
            table(raw / (key + '_altlocs.tsv'), ('chain', 'residue', 'resname', 'selected_altloc'),
                  [dict(chain=k[0], residue=k[1].strip(), resname=k[2].strip(), selected_altloc=v)
                   for k, v in labels.items()])
            (split / (key + '.pdb')).write_text(''.join(selected) + 'END\n')
            hetero = sorted({l[17:20].strip() for l in lines if l.startswith('HETATM')})
            cases.append(dict(case=key, pdb=pdb, model=i,
                              input_residues=len(residue_list(lines)),
                              hetero_resnames=','.join(hetero)))
    table(root / 'sources.tsv', ('pdb', 'url', 'bytes', 'sha256'), sources)
    table(root / 'models.tsv', ('case', 'pdb', 'model', 'input_residues', 'hetero_resnames'), cases)
    # Reuse legacy commands, force-field, restraint policy and 100-step CUDA minimization.
    # Path edits apply only to new working copies; original scripts remain untouched.
    prep_text = (scripts / 'prep_all.sh').read_text()
    prep_text = prep_text.replace(str(LEGACY / 'pdbs_raw'), str(split))
    prep_text = prep_text.replace(str(LEGACY / 'pdbs_prepped'), str(prepped))
    prep_text = prep_text.replace('cd ' + str(LEGACY), 'cd ' + str(root))
    (scripts / 'prep_isolated.sh').write_text(prep_text)
    min_text = (scripts / 'minimize_single.sh').read_text()
    min_text = min_text.replace(str(LEGACY / 'pdbs_prepped'), str(prepped))
    # Legacy code reads _prep instead of documented _rc; restore chains before minimization.
    min_text = min_text.replace('prepped="$dir/${pdb}_strip_reduce_prep.pdb"',
                                'prepped="$dir/${pdb}_strip_reduce_prep_rc.pdb"')
    # Preserve intermediates for atom/residue and parameterization audits.
    min_text = min_text.replace('rm -f ${pdb}.inpcrd mdinfo "$min_raw" "$ref_table" 2>/dev/null',
                                '# Preserve all intermediates for the experimental audit.')
    (scripts / 'minimize_isolated.sh').write_text(min_text)
    table(root / 'script_hashes.tsv', ('path', 'sha256'),
          [dict(path=p.name, sha256=sha(p)) for p in sorted(scripts.iterdir())])
    prep_rc = run(['bash', str(scripts / 'prep_isolated.sh')], logs / 'prep.log', root)
    if prep_rc:
        raise RuntimeError('Legacy prep failed; see logs/prep.log')
    checks = []
    for case in cases:
        key = case['case']
        work = prepped / key
        rc = run(['bash', str(scripts / 'minimize_isolated.sh'), key], logs / (key + '.log'), root)
        final = work / (key + '.min.reduce.pdb')
        reference = work / (key + '_strip_reduce_prep_rc.pdb')
        status, reason = 'PREP_FAILED', 'minimize exit {}'.format(rc)
        count = 0
        if rc == 0 and final.exists() and reference.exists():
            before = residue_list(reference.read_text().splitlines())
            after = residue_list(final.read_text().splitlines())
            count = len(after)
            norm = lambda r: (r[0], r[1], {'HID':'HIS', 'HIE':'HIS', 'HIP':'HIS', 'CYX':'CYS'}.get(r[2], r[2]))
            if [norm(r) for r in before] == [norm(r) for r in after] and count == case['input_residues']:
                status, reason = 'PREPARED_PROTEIN_ONLY', 'Residue identity/order retained; not experimental population validation'
                shutil.copy2(final, work / (key + '.min.reduce.renum.pdb'))
            else:
                status, reason = 'RESIDUE_AUDIT_FAILED', 'Input, prepared, or minimized residue identity/count differs'
        checks.append(dict(case=key, status=status, final_residues=count, reason=reason))
        table(root / 'prep_validation.tsv', ('case', 'status', 'final_residues', 'reason'), checks)
        print(key, status, flush=True)
    (root / 'eligibility.json').write_text(json.dumps({
        'schema': 'packstar-experimental-eligibility-v1',
        'population_calculation_run': False,
        'formal_biological_validation': 'BLOCKED_MODEL_AND_REFERENCE',
        'adk': ['1AKE is Ap5A-bound, 4AKE is apo; neither is the matched ATP-saturated ensemble.',
                'Legacy protein-only prep excludes nonprotein ligand and ions.',
                'Current fixed-backbone case does not integrate open/closed lid basins.'],
        'pdz2': ['7QCX/7QCY are 40 fitted conformers each, not 40 equally probable observations.',
                 'About 1:1 is CYANA-restraint-fit inference; a quantitative uncertainty remains to be established.',
                 'Apo/holo must remain separate chemical ensembles.',
                 'Independent backbone-basin weights and matching experimental state classifier are required.'],
        'craf': ['No matched full-length labelled cellular construct/structure and FRET forward model specified.'],
        'forbidden_shortcut': 'Do not average separate backbone P_C values equally or omit differing shell energies.'
    }, indent=2) + '\n')
    (root / 'PREP_FINISHED').write_text('Preparation attempted for all models; inspect prep_validation.tsv. Not formal success.\n')
    if any(c['status'] != 'PREPARED_PROTEIN_ONLY' for c in checks):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
