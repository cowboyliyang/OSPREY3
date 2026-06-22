"""
On-policy harvest for the DANCE leaf surrogate (DAgger-style).

The deterministic CCD audit (run for correctness) already produces exact
(conformation -> E_CCD) pairs drawn from EXACTLY the distribution the GNN-driven
search visits. This script folds those audited confs back into the per-confspace
training set, so the next model is trained on its own on-policy distribution.

It reconstructs the regression target from the emat tables that the data exporter
also used:  E_emat(conf) = sum_i onebody[i,rc_i] + sum_{i<j} pairwise[i,rc_i,j,rc_j]
            residual     = E_CCD - E_emat
(E_rigid is left NaN for harvested rows; the bracket-ceiling penalty skips them.)

Usage:
  python harvest_onpolicy.py --state_dir gnn_models/3gxu/complex \
      --audit audit_results/gnn_s11_leafonly/3gxu_complex_seq00007.csv \
      --out gnn_models/3gxu/complex/confs.csv   # writes augmented (backs up original)
  # add --validate to only check E_emat consistency without writing
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd


def build_emat(state_dir):
    ob = pd.read_csv(os.path.join(state_dir, "onebody_energies.csv"))
    pr = pd.read_csv(os.path.join(state_dir, "pairwise_energies.csv"))
    one = {(int(p), int(r)): float(e)
           for p, r, e in zip(ob["pos"], ob["rc"], ob["E_onebody"])}
    pair = {}
    for p1, r1, p2, r2, e in zip(pr["pos1"], pr["rc1"], pr["pos2"], pr["rc2"], pr["E_pair_min"]):
        pair[(int(p1), int(r1), int(p2), int(r2))] = float(e)
        pair[(int(p2), int(r2), int(p1), int(r1))] = float(e)
    return one, pair


def emat_energy(rcs, one, pair):
    """rcs: list of rc index per position (position order == meta order)."""
    e = 0.0
    miss = 0
    npos = len(rcs)
    for i in range(npos):
        e += one.get((i, rcs[i]), 0.0)
        if (i, rcs[i]) not in one:
            miss += 1
    for i in range(npos):
        for j in range(i + 1, npos):
            key = (i, rcs[i], j, rcs[j])
            if key in pair:
                e += pair[key]
            else:
                miss += 1
    return e, miss


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state_dir", required=True, help="gnn_models/<design>/<state>")
    ap.add_argument("--audit", required=True, nargs="+", help="audit_results CSV(s) (glob ok)")
    ap.add_argument("--out", default=None, help="output confs.csv (default: <state_dir>/confs.csv, backed up)")
    ap.add_argument("--validate", action="store_true", help="only check E_emat consistency, do not write")
    args = ap.parse_args()

    meta = pd.read_csv(os.path.join(args.state_dir, "meta.csv"))
    npos = meta.shape[0]
    rc_cols = [f"rc_{i}" for i in range(npos)]
    one, pair = build_emat(args.state_dir)

    # gather audit rows
    files = []
    for patt in args.audit:
        files.extend(sorted(glob.glob(patt)))
    if not files:
        print("no audit files matched", file=sys.stderr); sys.exit(1)
    arows = []
    for f in files:
        df = pd.read_csv(f)
        df = df[df.get("status", "ok").astype(str) == "ok"] if "status" in df.columns else df
        arows.append(df[["assignments", "ccd_energy_kcal"]])
    adf = pd.concat(arows, ignore_index=True).dropna()
    print(f"audit rows: {len(adf)} from {len(files)} file(s)")

    recs, misses = [], 0
    for assign, eccd in zip(adf["assignments"], adf["ccd_energy_kcal"]):
        rcs = [int(x) for x in str(assign).split(";")]
        if len(rcs) != npos:
            continue
        e_emat, m = emat_energy(rcs, one, pair)
        misses += m
        recs.append((tuple(rcs), float(eccd), e_emat, float(eccd) - e_emat))
    print(f"reconstructed {len(recs)} confs; pairwise/onebody table misses (total terms)={misses}")

    # ---- self-validation against existing confs.csv ----
    confs_path = os.path.join(args.state_dir, "confs.csv")
    cdf = pd.read_csv(confs_path)
    existing = {tuple(int(v) for v in row): (re, rr)
                for row, re, rr in zip(cdf[rc_cols].values, cdf["E_emat"].values, cdf["residual"].values)}
    diffs = [abs(e_emat - existing[rc][0]) for rc, _, e_emat, _ in recs if rc in existing]
    if diffs:
        print(f"E_emat self-check on {len(diffs)} overlapping confs: "
              f"max|Δ|={max(diffs):.6f}  mean|Δ|={np.mean(diffs):.6f} kcal/mol")
    else:
        print("no overlap with confs.csv to self-check (all harvested confs are new)")

    if args.validate:
        print("validate-only: not writing."); return

    # ---- build augmented confs.csv (dedup: keep existing rows, add new harvested) ----
    have = set(existing.keys())
    new_rows = []
    for rc, eccd, e_emat, resid in recs:
        if rc in have:
            continue
        have.add(rc)
        row = {f"rc_{i}": rc[i] for i in range(npos)}
        row.update({"E_CCD": eccd, "E_emat": e_emat, "E_rigid": np.nan, "residual": resid})
        new_rows.append(row)
    print(f"adding {len(new_rows)} new on-policy confs (skipped {len(recs) - len(new_rows)} already present)")

    out = args.out or confs_path
    if out == confs_path and not os.path.exists(confs_path + ".preharvest"):
        cdf.to_csv(confs_path + ".preharvest", index=False)
        print(f"backed up original -> {confs_path}.preharvest")
    aug = pd.concat([cdf, pd.DataFrame(new_rows)], ignore_index=True) if new_rows else cdf
    aug.to_csv(out, index=False)
    print(f"wrote {len(aug)} confs ({len(cdf)} original + {len(new_rows)} harvested) -> {out}")


if __name__ == "__main__":
    main()
