#!/usr/bin/env python3
"""Summarize the new-proposal PACK* rescue array."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


def read_kv(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    if not path.is_file():
        return result
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if "\t" in line:
                key, value = line.split("\t", 1)
                result[key] = value
    return result


def finite(value: str | None) -> bool:
    if value is None or value == "":
        return False
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def parse_wall(path: Path) -> tuple[str, str]:
    if not path.is_file():
        return "", ""
    text = path.read_text()
    match = re.search(r"elapsed=([^ ]+) maxRssKiB=([^\s]+)", text)
    return (match.group(1), match.group(2)) if match else ("", "")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--result-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.candidates.open(newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))

    rows: list[dict[str, str]] = []
    for candidate in candidates:
        case_id = candidate["case_id"]
        dirs = sorted(args.result_root.glob(f"case_{case_id}_*"))
        run_dir = dirs[0] if dirs else None
        manifest = read_kv(run_dir / "run_manifest.tsv") if run_dir else {}
        result_path = run_dir / "result.tsv" if run_dir else None
        result: dict[str, str] = {}
        if result_path is not None and result_path.is_file():
            with result_path.open(newline="") as handle:
                result_rows = list(csv.DictReader(handle, delimiter="\t"))
            if result_rows:
                result = result_rows[-1]
        elapsed, rss = parse_wall(run_dir / "wall.time") if run_dir else ("", "")
        protocol = {}
        if run_dir:
            protocol_files = sorted(run_dir.glob(
                "adaptive_frequency_severity/**/frequency_severity_protocol.tsv"))
            if protocol_files:
                protocol = read_kv(protocol_files[0])
        status = result.get("status", manifest.get("newPackStarStatus", "NO_RESULT"))
        rows.append({
            "case_id": case_id,
            "design": candidate.get("design", ""),
            "state": candidate.get("state", ""),
            "sequence": candidate.get("sequence", ""),
            "sequence_index": manifest.get("sequenceIndex", ""),
            "old_packstar_status": candidate.get("packstar_status", ""),
            "markstar_status": candidate.get("markstar_status", ""),
            "new_packstar_status": status,
            "java_status": manifest.get("javaStatus", "NO_TASK"),
            "new_epsilon": result.get("epsilon", ""),
            "new_lower_log10": result.get("lower_log10", ""),
            "new_upper_log10": result.get("upper_log10", ""),
            "num_confs": result.get("num_confs", ""),
            "elapsed_s": result.get("elapsed_s", elapsed),
            "max_rss_kib": rss,
            "joint_moment_learning": protocol.get(
                "jointMomentLearning",
                manifest.get("proposalJointMomentLearning", ""),
            ),
            "selected_triples": protocol.get("tripleEtaSelectedPositionTriples", ""),
            "selected_fill_edges": protocol.get("tripleEtaSelectedFillEdges", ""),
            "run_dir": str(run_dir) if run_dir else "",
        })

    rows.sort(key=lambda row: int(row["case_id"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else ["case_id"]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    counts: dict[str, int] = {}
    for row in rows:
        status = row["new_packstar_status"]
        counts[status] = counts.get(status, 0) + 1
    rescued = sum(row["new_packstar_status"] == "Estimated" for row in rows)
    complete = sum(row["java_status"] == "0" for row in rows)
    finite_intervals = sum(
        finite(row["new_lower_log10"]) and finite(row["new_upper_log10"])
        for row in rows
    )
    summary_path = args.output.with_name("rescue_summary.tsv")
    with summary_path.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"candidate_rows\t{len(rows)}\n")
        handle.write(f"java_completed_zero\t{complete}\n")
        handle.write(f"new_estimated\t{rescued}\n")
        handle.write(f"new_estimated_finite_interval\t{finite_intervals}\n")
        handle.write(f"new_not_estimated\t{len(rows) - rescued}\n")
        for status in sorted(counts):
            handle.write(f"status_{status}\t{counts[status]}\n")
        joint_values = sorted({row["joint_moment_learning"] for row in rows})
        handle.write(f"joint_moment_values\t{','.join(joint_values)}\n")

    print(f"results_output={args.output}")
    print(f"summary_output={summary_path}")
    print(f"candidate_rows={len(rows)}")
    print(f"java_completed_zero={complete}")
    print(f"new_estimated={rescued}")
    print(f"new_estimated_finite_interval={finite_intervals}")
    print(f"new_not_estimated={len(rows) - rescued}")
    for status in sorted(counts):
        print(f"status_{status}={counts[status]}")


if __name__ == "__main__":
    main()
