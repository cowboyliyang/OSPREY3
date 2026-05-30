#!/usr/bin/env python3
"""Profile MARK*/BranchMARK* leaf records for early hotspot scouting.

This reads [LEAF_PROFILE] records, groups them by state, and reports:
  - per-position rotamer-slice mass of oldZGap
  - ICC for preGap and actual upperDrop by position rotamer
  - top high-gap rotamer slices and pairs

preGap = E^+_rigid - E^-_corrected is available before minimization.
upperDrop = E^+_rigid - E_true requires minimization and is only a check.
"""

import argparse
import math
import re
from collections import defaultdict

import numpy as np


RE = re.compile(
    r"\[LEAF_PROFILE\]\s+alg=(?P<alg>\S+),\s+state=(?P<state>\S+?),\s+min=(?P<m>\d+),\s+"
    r"eTrue=(?P<etrue>-?\d+\.\d+),\s+preE=\[(?P<elo>-?\d+\.\d+),(?P<eup>-?\d+\.\d+)\],\s+"
    r"lowerRaise=(?P<lraise>-?\d+\.\d+),\s+upperDrop=(?P<udrop>-?\d+\.\d+),\s+"
    r"oldZGap=(?P<gap>\S+?),\s+exactZ=(?P<exz>\S+?),\s+"
    r"bestAt=\d+,\s+bestE=-?\d+\.\d+,\s+epsBefore=(?P<eps>\d+\.\d+),\s*"
    r"(?:epsAfter=\S+,\s*)?conf=(?P<conf>[-\d\s]+)\s*$"
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("--alg-prefix", default=None)
    parser.add_argument("--min-count", type=int, default=30)
    parser.add_argument("--top", type=int, default=20)
    return parser.parse_args()


def read_rows(path, alg_prefix):
    rows = []
    with open(path) as f:
        for line in f:
            if not line.startswith("[LEAF_PROFILE]"):
                continue
            m = RE.match(line.rstrip("\n"))
            if not m:
                continue
            d = m.groupdict()
            if alg_prefix and not d["alg"].startswith(alg_prefix):
                continue
            conf = tuple(int(x) for x in d["conf"].split())
            elo = float(d["elo"])
            eup = float(d["eup"])
            rows.append({
                "alg": d["alg"],
                "state": d["state"],
                "min": int(d["m"]),
                "etrue": float(d["etrue"]),
                "elo": elo,
                "eup": eup,
                "preGap": eup - elo,
                "lraise": float(d["lraise"]),
                "udrop": float(d["udrop"]),
                "gap": float(d["gap"]),
                "exz": float(d["exz"]),
                "eps": float(d["eps"]),
                "conf": conf,
            })
    return rows


def icc_by_position(rows, nflex, field, min_count):
    out = []
    for i in range(nflex):
        groups = defaultdict(list)
        for r in rows:
            groups[r["conf"][i]].append(r[field])
        big = [(rot, np.array(v)) for rot, v in groups.items() if len(v) >= min_count]
        if not big:
            continue
        sizes = np.array([len(v) for _, v in big], dtype=float)
        means = np.array([v.mean() for _, v in big])
        pooled_within = sum(s * v.var() for s, (_, v) in zip(sizes, big)) / sizes.sum()
        grand = (means * sizes).sum() / sizes.sum()
        between = (sizes * (means - grand) ** 2).sum() / sizes.sum()
        denom = between + pooled_within
        out.append((i, len(big), grand, math.sqrt(between), math.sqrt(pooled_within),
                    between / denom if denom > 0 else 0.0))
    return out


def analyze(rows, tag, min_count, top):
    if not rows:
        return

    nflex = len(rows[0]["conf"])
    gaps = np.array([max(0.0, r["gap"]) for r in rows])
    pregap = np.array([r["preGap"] for r in rows])
    udrop = np.array([r["udrop"] for r in rows])
    total_gap = gaps.sum()

    print("=" * 78)
    print(f"STATE = {tag}  n={len(rows)}  algs={sorted(set(r['alg'] for r in rows))}")
    print("=" * 78)
    print(f"eps first/last = {rows[0]['eps']:.6f} / {rows[-1]['eps']:.6f}")
    print(f"preGap mean/std/p50/p90 = {pregap.mean():.4f} / {pregap.std():.4f} / "
          f"{np.quantile(pregap, 0.50):.4f} / {np.quantile(pregap, 0.90):.4f}")
    print(f"upperDrop mean/std/p50/p90 = {udrop.mean():.4f} / {udrop.std():.4f} / "
          f"{np.quantile(udrop, 0.50):.4f} / {np.quantile(udrop, 0.90):.4f}")
    print(f"oldZGap total on profiled leaves = {total_gap:.6e}")

    print("\n[rot/pos] seen rotamers:")
    for i in range(nflex):
        rots = sorted({r["conf"][i] for r in rows})
        print(f"  pos {i}: nrot_seen={len(rots):2d}, rots={rots}")

    print("\n[ICC(preGap by position rotamer)] scout signal before minimization:")
    for i, nrot, grand, bstd, wstd, icc in icc_by_position(rows, nflex, "preGap", min_count):
        print(f"  pos {i}: nrot={nrot:2d}, grand={grand:6.3f}, "
              f"between_std={bstd:6.3f}, within_std={wstd:6.3f}, ICC={icc:.3f}")

    print("\n[ICC(upperDrop by position rotamer)] post-min check:")
    for i, nrot, grand, bstd, wstd, icc in icc_by_position(rows, nflex, "udrop", min_count):
        print(f"  pos {i}: nrot={nrot:2d}, grand={grand:6.3f}, "
              f"between_std={bstd:6.3f}, within_std={wstd:6.3f}, ICC={icc:.3f}")

    slice_rows = []
    for i in range(nflex):
        by_rot = defaultdict(list)
        for r in rows:
            by_rot[r["conf"][i]].append(r)
        for rot, vs in by_rot.items():
            if len(vs) < min_count:
                continue
            arr_gap = np.array([max(0.0, v["gap"]) for v in vs])
            arr_pre = np.array([v["preGap"] for v in vs])
            arr_u = np.array([v["udrop"] for v in vs])
            share = arr_gap.sum() / total_gap if total_gap > 0 else 0.0
            slice_rows.append((share, i, rot, len(vs), arr_gap.sum(), arr_pre.mean(),
                               arr_pre.std(), arr_u.mean(), arr_u.std()))

    slice_rows.sort(reverse=True)
    print("\n[top oldZGap position-rotamer slices] early table candidates:")
    print("  share    pos rot     n      oldZGap      preGap_mean/std   upperDrop_mean/std")
    for share, i, rot, n, gap_sum, pre_mean, pre_std, u_mean, u_std in slice_rows[:top]:
        print(f"  {share*100:6.2f}%  {i:3d} {rot:3d} {n:5d}  {gap_sum:11.4e}  "
              f"{pre_mean:7.3f}/{pre_std:6.3f}    {u_mean:7.3f}/{u_std:6.3f}")

    by_pos = defaultdict(float)
    for share, i, _, _, _, _, _, _, _ in slice_rows:
        by_pos[i] += share
    print("\n[position coverage over eligible slices]")
    for i, share in sorted(by_pos.items(), key=lambda x: -x[1]):
        print(f"  pos {i}: {share*100:6.2f}% of oldZGap in slices with n>={min_count}")

    pair_rows = []
    for i in range(nflex):
        for j in range(i + 1, nflex):
            by_pair = defaultdict(list)
            for r in rows:
                by_pair[(r["conf"][i], r["conf"][j])].append(r)
            for (ri, rj), vs in by_pair.items():
                if len(vs) < min_count:
                    continue
                arr_gap = np.array([max(0.0, v["gap"]) for v in vs])
                arr_pre = np.array([v["preGap"] for v in vs])
                share = arr_gap.sum() / total_gap if total_gap > 0 else 0.0
                pair_rows.append((share, i, ri, j, rj, len(vs), arr_pre.mean(), arr_pre.std()))

    pair_rows.sort(reverse=True)
    print("\n[top oldZGap pair slices] region/pair scouts:")
    print("  share    i:a   j:b      n      preGap_mean/std")
    for share, i, ri, j, rj, n, pre_mean, pre_std in pair_rows[:top]:
        print(f"  {share*100:6.2f}%  {i:2d}:{ri:<3d} {j:2d}:{rj:<3d} {n:5d}   "
              f"{pre_mean:7.3f}/{pre_std:6.3f}")
    print()


def main():
    args = parse_args()
    rows = read_rows(args.path, args.alg_prefix)
    print(f"# records = {len(rows)}")
    if not rows:
        return
    print(f"# states  = {sorted(set(r['state'] for r in rows))}")
    print(f"# algs    = {sorted(set(r['alg'] for r in rows))}")
    print()

    for state in sorted(set(r["state"] for r in rows)):
        analyze([r for r in rows if r["state"] == state], state, args.min_count, args.top)
    analyze(rows, "ALL", args.min_count, args.top)


if __name__ == "__main__":
    main()
