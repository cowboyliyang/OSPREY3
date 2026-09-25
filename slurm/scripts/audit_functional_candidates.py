"""Read-only structural candidate audit. Execute only in a Slurm allocation."""
import csv
import hashlib
import math
import os
from collections import OrderedDict
from pathlib import Path

if not os.environ.get("SLURM_JOB_ID"):
    raise SystemExit("This analysis requires Slurm")

ROOT = Path("/usr/xtmp/lz280/dance_bench/pdbs_prepped")
OUT = Path(os.environ["FUNCTIONAL_AUDIT_OUT"])
OUT.mkdir(parents=True, exist_ok=True)


def read_atoms(path):
    residues = OrderedDict()
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[16] not in " A":
            continue
        key = (line[21], line[22:27].strip())
        name = line[12:16].strip()
        if name.startswith("H") or line[76:78].strip() == "H":
            continue
        residue = residues.setdefault(key, {"type": line[17:20], "atoms": {}})
        residue["atoms"][name] = tuple(float(line[i:i+8]) for i in (30, 38, 46))
    return residues


with (OUT / "contacts.tsv").open("w") as contacts, (OUT / "chains.tsv").open("w") as chains:
    cw = csv.writer(contacts, delimiter="\t")
    sw = csv.writer(chains, delimiter="\t")
    cw.writerow(["pdb", "prepared_res1", "type1", "atom1", "prepared_res2", "type2", "atom2", "distance_A", "both_polar"])
    sw.writerow(["pdb", "chain", "count", "first", "last", "sequence_three_letter", "sha256"])
    for pdb in ("4wyu", "2xxm"):
        path = ROOT / pdb / f"{pdb}.min.reduce.renum.pdb"
        residues = read_atoms(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        for chain in dict.fromkeys(k[0] for k in residues):
            rows = [(k, v) for k, v in residues.items() if k[0] == chain]
            sw.writerow([pdb, chain, len(rows), rows[0][0][1], rows[-1][0][1], " ".join(v["type"] for k, v in rows), digest])
        items = list(residues.items())
        for idx, (k1, r1) in enumerate(items):
            for k2, r2 in items[idx+1:]:
                if k1[0] == k2[0]:
                    continue
                best = None
                polar = None
                for a1, c1 in r1["atoms"].items():
                    for a2, c2 in r2["atoms"].items():
                        d = math.dist(c1, c2)
                        row = (d, a1, a2)
                        if best is None or row < best:
                            best = row
                        if a1.startswith(("N", "O")) and a2.startswith(("N", "O")) and (polar is None or row < polar):
                            polar = row
                for row, flag in ((best, False), (polar, True)):
                    if row and row[0] <= (3.5 if flag else 4.5):
                        cw.writerow([pdb, "".join(k1), r1["type"], row[1], "".join(k2), r2["type"], row[2], f"{row[0]:.5f}", flag])
print(f"Candidate structure audit: {OUT}")
