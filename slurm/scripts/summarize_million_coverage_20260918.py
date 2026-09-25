"""Join all predeclared seed outcomes to the independently enumerated truth; Slurm only."""
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

root, reference = Path(sys.argv[1]), Path(sys.argv[2])
if not (reference / "READY").exists():
    raise SystemExit("Exhaustive reference incomplete; coverage not declared")
truth = json.loads((reference / "reference.json").read_text())
logz = truth["logZ"]
rows, counts = [], Counter()
for seed in range(10000, 10100):
    p = root / f"seed_{seed}"
    row = {"seed": seed, "outcome": "MISSING"}
    if (p / "exit_code").exists() and (p / "exit_code").read_text().strip() != "0":
        row["outcome"] = "PROCESS_FAILED"
    elif (p / "result.tsv").exists():
        with (p / "result.tsv").open() as f:
            row.update(next(csv.DictReader(f, delimiter="\t")))
        if int(row["count"]) != truth["total"] or (p / "layout.tsv").read_bytes() != (reference / "layout.tsv").read_bytes():
            raise RuntimeError(f"RC layout/count mismatch seed {seed}")
        with (p / "config.tsv").open() as f:
            config = next(csv.DictReader(f, delimiter="\t"))
        if config != truth["config"]:
            raise RuntimeError(f"Frozen configuration mismatch seed {seed}")
        log = (p / "run.log").read_text(errors="replace")
        deterministic = "sampling skipped; first-round DP bounds" in log or "exact small-state base case" in log
        row["deterministic"] = deterministic
        row["outcome"] = row["status"]
        if row["status"] == "Estimated":
            lo, hi = float(row["logLower"]), float(row["logUpper"])
            if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
                row["outcome"] = "INVALID_INTERVAL"
            else:
                row["covered"] = lo <= logz+1e-10 and hi >= logz-1e-10
                row["relative_width"] = -math.expm1(lo-hi)
                row["precision_met"] = row["relative_width"] <= 0.683+1e-7
                group = "deterministic" if deterministic else "sampled"
                counts[group+"_issued"] += 1
                counts[group+"_covered"] += int(row["covered"])
                counts["precision_not_met"] += int(not row["precision_met"])
    counts[row["outcome"]] += 1
    rows.append(row)
n, k = counts["sampled_issued"], counts["sampled_covered"]
interval = None
if n:
    z = 1.959963984540054
    center = (k/n+z*z/(2*n))/(1+z*z/n)
    half = z*math.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/(1+z*z/n)
    interval = [center-half, center+half]
report = {"expected": 100, "epsilon": .683, "delta": .05, "logZ_reference": logz,
          "premise": "conditional-relative-gauge-S0-20-not-externally-recalibrated",
          "counts": dict(counts), "sampled_coverage_wilson95": interval, "rows": rows}
(root / "coverage.json").write_text(json.dumps(report, indent=2))
text = (f"# {truth['config']['system']} million-scale Z coverage\n\nExpected: 100 seeds; epsilon=0.683; delta=0.05.\n\n"
        f"Full CCD census: {truth['total']} assignments; log Z = {logz}.\n\n"
        f"Outcomes: {dict(counts)}\n\nSampled coverage: {k}/{n}; Wilson 95%: {interval}.\n\n"
        "Deterministic early exits remain separate. All failures and missing runs remain in the denominator.\n"
        "Conditional S0=20 premise retained; no P_C calculation.\n")
(root / "coverage.md").write_text(text)
print(text)
