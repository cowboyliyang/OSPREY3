#!/usr/bin/env python3
"""Select old PACK* failures that MARK* measured for the same pfunc.

This is deliberately a small, deterministic selector.  It is run by Slurm;
the login node is never used to analyze benchmark outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


STATES = (
    ("protein", "prot_status", "prot_qstar_lb_log10", "prot_qstar_ub_log10"),
    ("ligand", "lig_status", "lig_qstar_lb_log10", "lig_qstar_ub_log10"),
    ("complex", "comp_status", "comp_qstar_lb_log10", "comp_qstar_ub_log10"),
)


def finite_float(value: str | None) -> bool:
    if value is None or value == "":
        return False
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def read_csv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle)
        if rows.fieldnames is None or "sequence" not in rows.fieldnames:
            raise RuntimeError(f"missing sequence column in {path}")
        result: dict[str, dict[str, str]] = {}
        for row in rows:
            sequence = (row.get("sequence") or "").strip()
            if sequence:
                result[sequence] = row
        return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packstar-root", type=Path, required=True)
    parser.add_argument("--markstar-root", type=Path, required=True)
    parser.add_argument("--design-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    specs: dict[str, dict[str, str]] = {}
    with args.design_csv.open(newline="") as handle:
        raw_lines = handle.readlines()
    header = None
    data_lines: list[str] = []
    for line in raw_lines:
        stripped = line.lstrip()
        if "design_id," in stripped and stripped.startswith("#"):
            header = stripped[1:].lstrip()
        elif line.strip() and not stripped.startswith("#"):
            data_lines.append(line)
    if header is None:
        raise RuntimeError(f"could not find design CSV header in {args.design_csv}")
    for row in csv.DictReader([header, *data_lines]):
            design = (row.get("design_id") or row.get("design") or "").strip()
            if design:
                specs[design] = row

    selected: list[dict[str, str]] = []
    totals = {
        "designs": 0,
        "common_designs": 0,
        "packstar_failed_pfuncs": 0,
        "markstar_measurable_failed_pfuncs": 0,
        "excluded_markstar_failed_pfuncs": 0,
        "excluded_missing_markstar": 0,
    }

    pack_files = sorted(args.packstar_root.glob("rankings/*.csv"))
    for pack_path in pack_files:
        design = pack_path.stem
        totals["designs"] += 1
        mark_path = args.markstar_root / f"{design}_markstar.csv"
        if not mark_path.is_file():
            continue
        totals["common_designs"] += 1
        pack_rows = read_csv(pack_path)
        mark_rows = read_csv(mark_path)
        spec = specs.get(design, {})
        for sequence, pack in pack_rows.items():
            mark = mark_rows.get(sequence)
            if mark is None:
                continue
            for state, status_col, lb_col, ub_col in STATES:
                pack_status = (pack.get(status_col) or "").strip()
                if pack_status == "Estimated":
                    continue
                totals["packstar_failed_pfuncs"] += 1
                mark_status = (mark.get(status_col) or "").strip()
                if mark_status != "Estimated":
                    totals["excluded_markstar_failed_pfuncs"] += 1
                    continue
                if not finite_float(mark.get(lb_col)) or not finite_float(mark.get(ub_col)):
                    totals["excluded_markstar_failed_pfuncs"] += 1
                    continue
                totals["markstar_measurable_failed_pfuncs"] += 1
                selected.append({
                    "design": design,
                    "pdb": spec.get("pdb", ""),
                    "mutable": spec.get("mutable", ""),
                    "flexible": spec.get("flexible", ""),
                    "state": state,
                    "sequence": sequence,
                    "packstar_status": pack_status,
                    "markstar_status": mark_status,
                    "markstar_lb_log10": mark.get(lb_col, ""),
                    "markstar_ub_log10": mark.get(ub_col, ""),
                    "packstar_rank": pack.get("rank", ""),
                    "packstar_score_log10": pack.get("score_log10", ""),
                    "markstar_score_log10": mark.get("score_log10", ""),
                })

    selected.sort(key=lambda row: (
        row["design"], row["state"], int(row["packstar_rank"] or "0"), row["sequence"]
    ))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_id", "design", "pdb", "mutable", "flexible", "state", "sequence",
        "packstar_status", "markstar_status", "markstar_lb_log10",
        "markstar_ub_log10", "packstar_rank", "packstar_score_log10",
        "markstar_score_log10",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, row in enumerate(selected):
            row = dict(row)
            row["case_id"] = str(index)
            writer.writerow(row)

    summary_path = args.output.with_name("selection_summary.tsv")
    with summary_path.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key, value in totals.items():
            handle.write(f"{key}\t{value}\n")
        handle.write(f"selected_rows\t{len(selected)}\n")

    print(f"selection_output={args.output}")
    print(f"selection_summary={summary_path}")
    for key, value in totals.items():
        print(f"{key}={value}")
    print(f"selected_rows={len(selected)}")


if __name__ == "__main__":
    main()
