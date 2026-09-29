#!/usr/bin/env python3
"""Reduce paired GPU runs, keeping output coverage alongside wall time."""
import csv
import json
import math
import os
from pathlib import Path
import sys

assert os.environ.get("SLURM_JOB_ID"), "Run analysis through Slurm"
out = Path(os.environ["AUTO_BENCH_SUMMARY"])
out.mkdir(parents=True, exist_ok=True)
summary = []
combined = {}
for folder in sys.argv[1:]:
    timing_file = Path(folder) / "timings.json"
    if not timing_file.exists():
        summary.append({"run": folder, "status": "no-timing-results"})
        continue
    for run in json.loads(timing_file.read_text()):
        if run["design"] not in combined:
            combined[run["design"]] = dict(run, phases={}, input_runs=[])
        merged = combined[run["design"]]
        assert (merged["pdb_sha256"], merged["gpus"], merged["cpus"]) == (run["pdb_sha256"], run["gpus"], run["cpus"])
        merged["input_runs"].append(folder)
        for name, phase in run["phases"].items():
            if name in merged["phases"]:
                assert merged["phases"][name] == phase, f"duplicate phase differs: {run['design']} {name}"
            merged["phases"][name] = phase

for run in combined.values():
        folder = ";".join(run["input_runs"])
        phases = run["phases"]
        row = {"design": run["design"], "run": folder, "node": run["node"],
               "gpus": run["gpus"], "cpus": run["cpus"]}
        row["phase_nodes"] = {k: v.get("node", "original-preflight") for k, v in phases.items()}
        valid = True
        for name in ("preflight", "baseline", "optimized"):
            phase = phases.get(name, {})
            row[name + "_seconds"] = phase.get("seconds")
            row[name + "_exit"] = phase.get("exit_code")
            valid &= phase.get("exit_code") == 0
        if not valid:
            row["status"] = "incomplete-or-failed-timing"
            summary.append(row)
            continue
        manifests = {name: json.loads((Path(phases[name]["output"]) / "manifest.json").read_text())
                     for name in ("preflight", "baseline", "optimized")}
        assert manifests["baseline"]["energy_matrix_sha256"] == manifests["optimized"]["energy_matrix_sha256"]
        assert manifests["baseline"]["energy_matrix_sha256"] == manifests["preflight"]["energy_matrix_sha256"]
        rows = {name: {r["sequence"]: r for r in phases[name].get("sequence_results", [])}
                for name in ("baseline", "optimized")}
        states = ("prot", "lig", "comp")
        signatures = {}
        for name, sequences in rows.items():
            signatures[name] = {seq: tuple(r[s + "_status"] for s in states)
                                for seq, r in sequences.items()}
            row[name + "_sequences"] = len(sequences)
            row[name + "_all_estimated_sequences"] = sum(
                all(status == "Estimated" for status in signature)
                for signature in signatures[name].values())
            row[name + "_estimated_pfunc_cells"] = sum(
                status == "Estimated" for signature in signatures[name].values() for status in signature)
        row["same_sequences"] = bool(rows["baseline"]) and rows["baseline"].keys() == rows["optimized"].keys()
        row["same_statuses"] = row["same_sequences"] and signatures["baseline"] == signatures["optimized"]
        row["status"] = "paired-matching-statuses" if row["same_statuses"] else "paired-different-coverage"
        pre, base, opt = (row[k + "_seconds"] for k in ("preflight", "baseline", "optimized"))
        row["preflight_plus_optimized_seconds"] = pre + opt
        row["observed_net_saved_seconds"] = base - opt - pre
        row["observed_net_speedup"] = base / (pre + opt)
        row["observed_break_even_uses"] = math.ceil(pre / (base - opt)) if base > opt else None
        # Different completeness may change the amount of executed work; retain
        # raw timings but do not label the resulting ratio a matched speedup.
        row["matched_net_speedup"] = row["observed_net_speedup"] if row["same_statuses"] else None
        summary.append(row)

(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
fields = list(dict.fromkeys(key for row in summary for key in row))
with (out / "summary.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)
print(json.dumps(summary, indent=2))
