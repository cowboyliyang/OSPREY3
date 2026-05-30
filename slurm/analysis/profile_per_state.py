#!/usr/bin/env python3
"""Per-state MARK* leaf profile analysis.

Splits leaf records by state (Protein/Complex/Ligand) and emits the same
hotspot/pair-stability/Z-concentration stats as profile_markstar_leaves.py
but separately per state, then prints a combined summary.
"""
import re, sys, math
from collections import defaultdict
import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/users/lz280/IdeaProjects/OSPREY3/slurm/outputs/markstar_profile_11739460.out"

RE = re.compile(
    r"\[LEAF_PROFILE\]\s+alg=(?P<alg>\S+),\s+state=(?P<state>\S+?),\s+min=(?P<m>\d+),\s+"
    r"eTrue=(?P<etrue>-?\d+\.\d+),\s+preE=\[(?P<elo>-?\d+\.\d+),(?P<eup>-?\d+\.\d+)\],\s+"
    r"lowerRaise=(?P<lraise>-?\d+\.\d+),\s+upperDrop=(?P<udrop>-?\d+\.\d+),\s+"
    r"oldZGap=(?P<gap>\S+?),\s+exactZ=(?P<exz>\S+?),\s+"
    r"bestAt=\d+,\s+bestE=-?\d+\.\d+,\s+epsBefore=(?P<eps>\d+\.\d+),\s*"
    r"(?:epsAfter=\S+,\s*)?conf=(?P<conf>[\d\s]+)\s*$"
)

all_rows = []
with open(PATH) as f:
    for line in f:
        if not line.startswith("[LEAF_PROFILE]"): continue
        if "alg=MARK*," not in line: continue
        m = RE.match(line.rstrip("\n"))
        if not m: continue
        d = m.groupdict()
        conf = tuple(int(x) for x in d["conf"].split())
        all_rows.append({
            "state": d["state"], "min": int(d["m"]),
            "etrue": float(d["etrue"]), "elo": float(d["elo"]), "eup": float(d["eup"]),
            "lraise": float(d["lraise"]), "udrop": float(d["udrop"]),
            "gap": float(d["gap"]), "exz": float(d["exz"]),
            "eps": float(d["eps"]), "conf": conf,
        })

if not all_rows: sys.exit(0)
RT = 1.9858775e-3 * 298.0
nflex = len(all_rows[0]["conf"])
states = sorted(set(r["state"] for r in all_rows))

print(f"# total records = {len(all_rows)}")
print(f"# states        = {states}")
print(f"# flexible pos  = {nflex}")
print()

def analyze(rows, tag):
    N = len(rows)
    if N == 0: return
    print("=" * 70)
    print(f"STATE = {tag}  (n = {N})")
    print("=" * 70)
    u = np.array([math.exp(-r["elo"]/RT) for r in rows])
    udrop = np.array([r["udrop"] for r in rows])
    lraise = np.array([r["lraise"] for r in rows])

    # ---- Rotamer space ----
    seen = defaultdict(set)
    for r in rows:
        for i, rc in enumerate(r["conf"]): seen[i].add(rc)
    print(f"\n[rot/pos] per-position rotamer counts (seen / max_idx+1):")
    for i in sorted(seen):
        rots = sorted(seen[i])
        print(f"  pos {i}: nrot_seen={len(rots):2d}, max_idx={max(rots):2d}, rots={rots}")

    # ---- Z+ concentration ----
    sortedU = np.sort(u)[::-1]; cum = np.cumsum(sortedU) / sortedU.sum()
    print("\n[Z+ concentration on minimized leaves]:")
    for frac in (0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50):
        k = max(1, int(frac * N))
        print(f"  top {frac*100:5.1f}% leaves ({k:5d}) capture {cum[k-1]*100:6.2f}% of sum(u_c)")
    print(f"  effective # leaves (1/sum(p^2)) = {1/((sortedU/sortedU.sum())**2).sum():.1f}")

    # ---- upperDrop distribution ----
    print("\n[upperDrop = E^+_rigid - E_min, kcal/mol]:")
    for q in (0.01, 0.05, 0.25, 0.50, 0.75, 0.90, 0.99):
        print(f"  p{q*100:4.1f} = {np.quantile(udrop, q):.4f}")
    print(f"  mean   = {udrop.mean():.4f}")
    print(f"  std    = {udrop.std():.4f}")
    wts = u / u.sum()
    mw_mean = (udrop * wts).sum()
    mw_std  = math.sqrt(((udrop - mw_mean)**2 * wts).sum())
    print(f"  u_c-weighted mean = {mw_mean:.4f}, std = {mw_std:.4f}")

    # ---- ICC per position ----
    print("\n[ICC per position] (fraction of udrop variance explained by 'which rotamer is at i'):")
    for i in range(nflex):
        groups = defaultdict(list)
        for r in rows: groups[r["conf"][i]].append(r["udrop"])
        big = [(rot, np.array(v)) for rot, v in groups.items() if len(v) > 30]
        if not big: continue
        means = np.array([v.mean() for _, v in big])
        sizes = np.array([len(v) for _, v in big])
        pooled_within = sum(s * v.var() for s, (_, v) in zip(sizes, big)) / sizes.sum()
        grand = (means * sizes).sum() / sizes.sum()
        between = (sizes * (means - grand)**2).sum() / sizes.sum()
        icc = between / (between + pooled_within) if (between + pooled_within) > 0 else 0
        print(f"  pos {i}: nrot_in_groups={len(big):2d}, grand_mean={grand:5.2f}, "
              f"between_std={math.sqrt(between):5.2f}, within_std={math.sqrt(pooled_within):5.2f}, ICC={icc:.3f}")

    # ---- Per-pair stable signals ----
    pair_drops = defaultdict(list); pair_u = defaultdict(list)
    for r in rows:
        c = r["conf"]; ud = r["udrop"]; uw = math.exp(-r["elo"]/RT)
        for i in range(nflex):
            for j in range(i+1, nflex):
                pair_drops[(i, c[i], j, c[j])].append(ud)
                pair_u[(i, c[i], j, c[j])].append(uw)

    pair_rows = []
    for k, vs in pair_drops.items():
        if len(vs) < 30: continue
        a = np.array(vs); us = np.array(pair_u[k]); wsum = us.sum()
        mw = (a * us).sum() / wsum
        mw_std = math.sqrt(((a - mw)**2 * us).sum() / wsum)
        pair_rows.append((*k, len(vs), float(a.mean()), float(a.std()),
                          float(mw), float(mw_std)))

    strict = [x for x in pair_rows if x[5] > 0.5 and x[6] < 0.1 * x[5]]
    stable = [x for x in pair_rows if x[5] > 0.5 and x[6] < 0.3 * x[5]]
    print(f"\n[stable pairs] n>=30, mean>0.5 kcal/mol:")
    print(f"  total evaluated     = {len(pair_rows)}")
    print(f"  STABLE (std/mean<0.3) = {len(stable)}")
    print(f"  STRICT (std/mean<0.1) = {len(strict)}")
    strict.sort(key=lambda x: -x[5])
    print("\n[top 20 STRICT pairs] (i:a j:b n mean std uwt_mean uwt_std):")
    for x in strict[:20]:
        print("  {:>2d}:{:<2d} {:>2d}:{:<2d} {:>4d}  {:6.3f}  {:5.3f}  {:6.3f}  {:5.3f}".format(*x))

    # ---- Z+ closure concentration ----
    zclose = u * (1 - np.exp(-udrop / RT))
    print(f"\n[Z+ closure] = u_c * (1 - exp(-udrop/RT))")
    print(f"  total closure achieved: {zclose.sum():.3e}")
    print(f"  total Z+ on minimized:  {u.sum():.3e}")
    print(f"  fraction closed: {zclose.sum()/u.sum()*100:.2f}%")
    sortedZ = np.sort(zclose)[::-1]; cumZ = np.cumsum(sortedZ) / sortedZ.sum()
    for frac in (0.001, 0.01, 0.05, 0.10, 0.25, 0.50):
        k = max(1, int(frac * N))
        print(f"  top {frac*100:5.1f}% leaves capture {cumZ[k-1]*100:6.2f}% of total closure")

    # ---- lowerRaise (corrections already applied; this is leftover) ----
    print(f"\n[lowerRaise = E_min - E^-_corrected (post-correction)]:")
    print(f"  mean = {lraise.mean():.4f}  std = {lraise.std():.4f}  p50 = {np.quantile(lraise,0.5):.4f}  p90 = {np.quantile(lraise,0.9):.4f}")
    print()

for st in states:
    rows_st = [r for r in all_rows if r["state"] == st]
    analyze(rows_st, st)

print("=" * 70)
print("COMBINED ALL STATES")
print("=" * 70)
analyze(all_rows, "ALL")
