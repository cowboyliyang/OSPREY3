"""Run only in Slurm. Include all 400 predeclared trials, never success-filter silently."""
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
gate = len(sys.argv) > 2 and sys.argv[2] == "gate"
systems = ("1gwc", "2p4a", "2rl0", "4wyu")

def wilson(k, n):
    if not n:
        return [None, None]
    z = 1.959963984540054
    p = k/n
    c = (p + z*z/(2*n))/(1+z*z/n)
    h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [c-h, c+h]

rows, groups = [], []
errors = []
for system in systems:
    counts = Counter()
    widths = []
    for seed in range(10000, 10001 if gate else 10100):
        path = root / "runs" / system / f"seed_{seed}"
        row = {"system": system, "seed": seed, "outcome": "MISSING"}
        if (path / "exit_code").exists() and (path / "exit_code").read_text().strip() != "0":
            row["outcome"] = "PROCESS_FAILED"
        elif (path / "result.tsv").exists():
            with (path / "result.tsv").open() as f:
                row.update(next(csv.DictReader(f, delimiter="\t")))
            log = (path / "run.log").read_text(errors="replace")
            deterministic = "sampling skipped; first-round DP bounds" in log or "exact small-state base case" in log
            row["deterministic"] = deterministic
            row["outcome"] = "DETERMINISTIC" if deterministic else row["status"]
            lo, hi, truth = (float(row[k]) for k in ("logLower", "logUpper", "logZ"))
            if row["status"] == "Estimated":
                if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
                    row["outcome"] = "INVALID_INTERVAL"
                elif 1-math.exp(lo-hi) > .683+1e-7:
                    row["outcome"] = "PRECISION_NOT_MET"
                else:
                    counts["issued"] += 1
                    covered = lo <= truth+1e-10 and hi >= truth-1e-10
                    counts["covered"] += int(covered)
                    if not deterministic:
                        counts["sampled_issued"] += 1
                        counts["sampled_covered"] += int(covered)
                    widths.append(hi-lo)
            if gate and (deterministic or row["outcome"] != "Estimated"):
                errors.append(f"{system}: pilot {row['outcome']}")
        elif (path / "STARTED").exists():
            row["outcome"] = "NO_RESULT"
        if gate and row["outcome"] in ("MISSING", "PROCESS_FAILED", "NO_RESULT"):
            errors.append(f"{system}: {row['outcome']}")
        counts[row["outcome"]] += 1
        rows.append(row)
    groups.append({"system": system, "counts": dict(counts),
                   "sampled_coverage_wilson95": wilson(counts["sampled_covered"], counts["sampled_issued"]),
                   "mean_log_width": sum(widths)/len(widths) if widths else None})

prefix = "pilot" if gate else "summary"
(root / f"{prefix}.json").write_text(json.dumps({"groups": groups, "rows": rows, "gate_errors": errors}, indent=2))
lines = ["# Small-system Z coverage", "", "epsilon=0.683; conditional S0=20; delta=0.05. No P_C.", "",
         "| System | Counts (all outcomes) | Sampled coverage Wilson 95% |", "|---|---|---|"]
for g in groups:
    lines.append(f"| {g['system']} | {g['counts']} | {g['sampled_coverage_wilson95']} |")
(root / f"{prefix}.md").write_text("\n".join(lines)+"\n")
print("\n".join(lines))
if gate:
    if errors:
        raise SystemExit("Pilot gate failed: " + "; ".join(errors))
    (root / "GATE_READY").write_text("All pilots sampled and produced usable intervals; coverage misses were not filtered.\n")
