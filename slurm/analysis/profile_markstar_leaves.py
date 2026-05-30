#!/usr/bin/env python3
"""Profile MARK* leaf minimization log to test region-atom hypothesis.

Reads a SLURM output with [LEAF_PROFILE] alg=MARK* lines and emits:
  (1) Z^+ concentration: fraction of Z^+ captured by top-k highest u_c leaves.
  (2) ΔE^+ (upperDrop) global distribution: percentiles + mass-weighted percentiles.
  (3) ΔE^+ stability per single-position rotamer assignment and per pair.
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

rows = []
with open(PATH) as f:
    for line in f:
        if not line.startswith("[LEAF_PROFILE]"): continue
        if "alg=MARK*," not in line: continue
        m = RE.match(line.rstrip("\n"))
        if not m: continue
        d = m.groupdict()
        conf = tuple(int(x) for x in d["conf"].split())
        rows.append({
            "alg": d["alg"], "state": d["state"], "min": int(d["m"]),
            "etrue": float(d["etrue"]), "elo": float(d["elo"]), "eup": float(d["eup"]),
            "lraise": float(d["lraise"]), "udrop": float(d["udrop"]),
            "gap": float(d["gap"]), "exz": float(d["exz"]),
            "eps": float(d["eps"]), "conf": conf,
        })

N = len(rows)
print(f"# records = {N}")
if N == 0: sys.exit(0)
state_set = sorted(set(r["state"] for r in rows))
print(f"# states  = {state_set}")
nflex = len(rows[0]["conf"])
print(f"# flexible positions = {nflex}")

# RT at 298K in kcal/mol
RT = 1.9858775e-3 * 298.0

# Per-leaf u_c = exp(-elo/RT), l_c = exp(-eup/RT), w_c = exactZ = exp(-etrue/RT)
u = np.array([math.exp(-r["elo"]/RT) for r in rows])
l = np.array([math.exp(-r["eup"]/RT) for r in rows])
w = np.array([r["exz"] for r in rows])
udrop = np.array([r["udrop"] for r in rows])
lraise = np.array([r["lraise"] for r in rows])

# ----- (1) Z+ concentration on MINIMIZED leaves only -----
# Note: this is concentration ACROSS the leaves MARK* has chosen to minimize, which
# are themselves selected to be high-u_c. So this gives the concentration AMONG the
# currently-most-uncertain leaves, not over all conformations.
sortedU = np.sort(u)[::-1]
cum = np.cumsum(sortedU) / sortedU.sum()
print("\n[1] u_c concentration on the minimized leaves (these are MARK*'s most-uncertain c):")
for frac in (0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50):
    k = max(1, int(frac * N))
    print(f"  top {frac*100:5.1f}% leaves ({k:5d}) capture {cum[k-1]*100:6.2f}% of sum(u_c)")
print(f"  effective # leaves (Gini-like, 1/sum(p^2)) = {1/((sortedU/sortedU.sum())**2).sum():.1f}")

# ----- (2) Global upperDrop distribution -----
print("\n[2] upperDrop (E^+_rigid - E_min) distribution (kcal/mol):")
for q in (0.01, 0.05, 0.25, 0.50, 0.75, 0.90, 0.99):
    print(f"  p{q*100:4.1f} = {np.quantile(udrop, q):.4f}")
print(f"  mean   = {udrop.mean():.4f}")
print(f"  std    = {udrop.std():.4f}")
print(f"  CV     = {udrop.std()/udrop.mean():.4f}")
# Mass-weighted (weight by u_c, because we care about Z^+ contribution)
wts = u / u.sum()
mw_mean = (udrop * wts).sum()
mw_std  = math.sqrt(((udrop - mw_mean)**2 * wts).sum())
print(f"  u_c-weighted mean = {mw_mean:.4f}   std = {mw_std:.4f}   CV = {mw_std/mw_mean:.4f}")

# How much would Z^- improve if we naively applied a global udrop_mean correction?
# new E^+_eff = E^+_rigid - constantDrop  =>  new l_c = exp(-(eup - drop)/RT)
print("\n  Hypothetical: if we could universally subtract constant Δ from E^+_rigid,")
print("  the resulting Z^- would scale as exp(+Δ/RT), so eps_new ~ 1 - (Z^-)/(Z^+) shrinks.")
print("  Current eps:", rows[-1]["eps"])
Zminus = l.sum(); Zplus = u.sum()
print(f"  Z^- (on these leaves) = {Zminus:.3e}")
print(f"  Z^+ (on these leaves) = {Zplus:.3e}")
print(f"  ratio Z^-/Z^+         = {Zminus/Zplus:.4e}   eps_local = {1 - Zminus/Zplus:.4f}")
for delta in (0.1, 0.3, 0.5, 1.0):
    z2 = (l * math.exp(delta/RT)).sum()
    z2 = min(z2, Zplus)  # can never exceed Z+
    print(f"  if all E^+ pushed down by {delta:.1f} kcal/mol: eps_local -> {1 - z2/Zplus:.4f}")

# ----- (3) ΔE^+ stability per single-position rotamer / per pair -----
# For each position i and rotamer r at i, gather udrop of all confs with conf[i]=r.
# A stable signal means: small std, large mean. Same for pairs.
print("\n[3] ΔE^+ stability across leaves -- per position-rotamer single:")
pos_rot_drops = defaultdict(list); pos_rot_u = defaultdict(list)
for r in rows:
    c = r["conf"]
    for i, rci in enumerate(c):
        pos_rot_drops[(i, rci)].append(r["udrop"])
        pos_rot_u[(i, rci)].append(math.exp(-r["elo"]/RT))

# Pick top single-rotamer signals by (mean udrop) * (count) where count > 30
single_rows = []
for (i, r), vs in pos_rot_drops.items():
    if len(vs) < 30: continue
    a = np.array(vs); us = np.array(pos_rot_u[(i, r)])
    wsum = us.sum()
    mw = (a * us).sum() / wsum
    mw_var = ((a - mw)**2 * us).sum() / wsum
    single_rows.append((i, r, len(vs), float(a.mean()), float(a.std()),
                        float(mw), math.sqrt(mw_var)))
single_rows.sort(key=lambda t: -t[3])
print("  pos rot  n     mean    std    u-weight-mean  u-weight-std")
for x in single_rows[:15]:
    print("  {:>3d} {:>3d} {:>4d}   {:6.3f} {:6.3f}     {:6.3f}       {:6.3f}".format(*x))
print("  ... (sorted by mean upperDrop desc, only n>=30)")

print("\n[3b] ΔE^+ stability per rotamer-pair (i:a, j:b) with n>=30:")
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
    pair_rows.append((*k, len(vs), float(a.mean()), float(a.std()), float(mw), float(mw_std)))
pair_rows.sort(key=lambda t: -t[5])  # by u-weight-mean
print("  i:a  j:b   n     mean   std   u-wt-mean  u-wt-std    (top 20 by u-weighted mean)")
for x in pair_rows[:20]:
    print("  {:>2d}:{:<2d} {:>2d}:{:<2d} {:>4d}  {:5.3f}  {:5.3f}    {:5.3f}      {:5.3f}".format(*x))

# Fixed filter — tuple is (i,a,j,b,n,mean,std,uwt_mean,uwt_std)
# mean=x[5], std=x[6], uwt_mean=x[7], uwt_std=x[8]
stable_pairs = [x for x in pair_rows if x[5] > 0.5 and x[6] < 0.3 * x[5]]
print(f"\n  # pairs with n>=30, mean>0.5 kcal/mol, std/mean<0.3 = {len(stable_pairs)}")
print(f"  # total pairs evaluated (n>=30) = {len(pair_rows)}")
stable_pairs.sort(key=lambda x: -x[5])
print("  TOP 20 STABLE PAIRS by mean ΔE^+ (i:a  j:b  n  mean  std  uwt_mean  uwt_std):")
for x in stable_pairs[:20]:
    print("    {:>2d}:{:<2d} {:>2d}:{:<2d} {:>4d}  {:6.3f}  {:5.3f}    {:6.3f}    {:5.3f}".format(*x))

strict = [x for x in pair_rows if x[5] > 0.5 and x[6] < 0.1 * x[5]]
print(f"\n  STRICT (std/mean < 0.1) pairs: {len(strict)}")
print("  TOP 10:")
strict.sort(key=lambda x: -x[5])
for x in strict[:10]:
    print("    {:>2d}:{:<2d} {:>2d}:{:<2d} {:>4d}  {:6.3f}  {:5.3f}    {:6.3f}    {:5.3f}".format(*x))

# ----- (3c) Variance decomposition per position i -----
print("\n[3c] ΔE^+ variance decomposition by position i:")
for i in range(nflex):
    rot_groups = {}
    for r in rows:
        rot_groups.setdefault(r["conf"][i], []).append(r["udrop"])
    big = [(rot, np.array(v)) for rot, v in rot_groups.items() if len(v) > 30]
    if not big: continue
    means = np.array([v.mean() for _, v in big])
    sizes = np.array([len(v) for _, v in big])
    pooled_within_var = sum(s * v.var() for s, (_, v) in zip(sizes, big)) / sizes.sum()
    grand_mean = (means * sizes).sum() / sizes.sum()
    between_var = (sizes * (means - grand_mean)**2).sum() / sizes.sum()
    print(f"  pos {i}: nrot_seen={len(big):2d}, grand_mean={grand_mean:5.2f}, "
          f"between_std={math.sqrt(between_var):5.2f}, within_pooled_std={math.sqrt(pooled_within_var):5.2f}, "
          f"ICC={between_var/(between_var+pooled_within_var):.3f}")

# ----- (5) Predictive power of (i:a) alone for predicting udrop -----
print("\n[5] Top single-position-rotamer tags by E[udrop | pos i=rot a] - global_mean, n>=100:")
global_mean = udrop.mean()
tags = []
for (i, a), vs in pos_rot_drops.items():
    if len(vs) < 100: continue
    arr = np.array(vs)
    tags.append((i, a, len(vs), float(arr.mean()), float(arr.std()),
                 float(arr.mean() - global_mean)))
tags.sort(key=lambda t: -t[5])
print("  pos rot  n      mean   std   delta-from-global")
for x in tags[:15]:
    print("  {:>3d} {:>3d} {:>5d}  {:6.3f} {:5.3f}    {:+6.3f}".format(*x))

# ----- (6) Mass-reduction concentration -----
print("\n[6] Actual Z^+ closure from each leaf min: u_c * (1 - exp(-udrop/RT))")
zclose = u * (1 - np.exp(-udrop / RT))
print(f"  total Z^+ closure achieved: {zclose.sum():.3e}")
print(f"  total Z^+ on minimized leaves: {u.sum():.3e}")
print(f"  fraction closed: {zclose.sum()/u.sum()*100:.2f}%")
sortedZ = np.sort(zclose)[::-1]
cumZ = np.cumsum(sortedZ) / sortedZ.sum()
for frac in (0.001, 0.01, 0.05, 0.10, 0.25, 0.50):
    k = max(1, int(frac * N))
    print(f"  top {frac*100:5.1f}% leaves capture {cumZ[k-1]*100:6.2f}% of total Z^+ closure")

print("\n[4] ΔE^- (lowerRaise) stats: how much triple corrections + minimization raise E^-")
print(f"  mean = {lraise.mean():.4f}  std = {lraise.std():.4f}  p50 = {np.quantile(lraise,0.5):.4f}  p90 = {np.quantile(lraise,0.9):.4f}")
