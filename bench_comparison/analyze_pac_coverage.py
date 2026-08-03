"""
PACK* paper "PAC coverage verification": aggregate the ground-truth run
(pac_coverage_ground_truth.csv, MARK* at very low epsilon) and the 100
per-seed PACK* runs (pac_coverage_seed_<seed>.csv) produced by
bench_pac_coverage.sh, and report the fraction of runs whose PAC interval
[lb_log10, ub_log10] contains the ground-truth q (point estimate =
ground truth's score_log10), per sequence and pooled overall.

Usage: python3 analyze_pac_coverage.py
"""
import glob
import os

import numpy as np
import pandas as pd

R = "/usr/xtmp/lz280/bench_comparison/results"
GT_CSV = f"{R}/pac_coverage_ground_truth.csv"
SEED_GLOB = f"{R}/pac_coverage_seed_*.csv"


def load(path):
    return pd.read_csv(path) if os.path.exists(path) else None


def main():
    gt = load(GT_CSV)
    if gt is None:
        print(f"missing {GT_CSV} -- run: bash bench_pac_coverage.sh ground_truth")
        return
    gt = gt[gt["method"] == "markstar"].set_index("sequence")

    seed_files = sorted(glob.glob(SEED_GLOB))
    if not seed_files:
        print(f"no seed files matching {SEED_GLOB} -- run: bash bench_pac_coverage.sh seeds")
        return

    rows = []  # (seed, sequence, covered)
    for f in seed_files:
        seed = os.path.basename(f)[len("pac_coverage_seed_"):-len(".csv")]
        df = load(f)
        if df is None or df.empty:
            continue
        df = df[df["method"] == "packstar"].set_index("sequence")
        for seq in df.index:
            if seq not in gt.index:
                continue
            r = df.loc[seq]
            g = gt.loc[seq]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            if isinstance(g, pd.DataFrame):
                g = g.iloc[0]
            lb, ub, q = r["lb_log10"], r["ub_log10"], g["score_log10"]
            if pd.isna(lb) or pd.isna(ub) or pd.isna(q) or not np.isfinite(q):
                continue
            covered = (lb - 1e-6) <= q <= (ub + 1e-6)
            rows.append((seed, seq, covered))

    if not rows:
        print("no comparable (seed, sequence) rows found -- check ground truth and seed CSVs share sequences")
        return

    res = pd.DataFrame(rows, columns=["seed", "sequence", "covered"])
    n_runs = res["seed"].nunique()

    print(f"PAC coverage verification: {n_runs} seed run(s) found (target 100), "
          f"{res['sequence'].nunique()} sequence(s) comparable to ground truth\n")

    print("per-sequence coverage (target >= 0.95):")
    per_seq = res.groupby("sequence")["covered"].agg(["mean", "sum", "count"])
    per_seq = per_seq.rename(columns={"mean": "coverage", "sum": "covered_runs", "count": "n_runs"})
    print(per_seq.to_string(float_format=lambda x: f"{x:.3f}"))

    overall = res["covered"].mean()
    print(f"\npooled coverage over all (seed, sequence) pairs: {overall:.3f} "
          f"({res['covered'].sum()}/{len(res)})")

    # Per-run coverage using only the top-ranked (rank==1) sequence, closest to the
    # paper's literal "fraction of runs whose interval contains q" (single target).
    rank1_seq = gt.index[0] if len(gt.index) else None
    if rank1_seq is not None:
        r1 = res[res["sequence"] == rank1_seq]
        if len(r1):
            print(f"\nrank-1 ground-truth sequence ({rank1_seq}) coverage: "
                  f"{r1['covered'].mean():.3f} ({r1['covered'].sum()}/{len(r1)} runs) "
                  f"-- this is the headline number for the paper's coverage claim")

    if n_runs < 100:
        print(f"\nNOTE: only {n_runs}/100 seed jobs have finished/written a CSV so far.")


if __name__ == "__main__":
    main()
