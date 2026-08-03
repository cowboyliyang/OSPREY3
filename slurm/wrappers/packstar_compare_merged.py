#!/usr/bin/env python3
"""Exact comparison of canonical PACK* merged CSV files."""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

from packstar_multinode_merge import EXPECTED_HEADER, normalize_sequence


SCHEMA = "packstar-merged-comparison-v1"


def load(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if list(reader.fieldnames or []) != EXPECTED_HEADER:
            raise ValueError(f"non-canonical CSV header: {path}")
        rows: Dict[str, Dict[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"extra fields at {path}:{line_number}")
            sequence = normalize_sequence(row["sequence"])
            if sequence in rows:
                raise ValueError(f"duplicate sequence {sequence!r} in {path}")
            normalized = {column: row[column].strip() for column in EXPECTED_HEADER}
            normalized["sequence"] = sequence
            rows[sequence] = normalized
    return rows


def atomic_tsv(path: Path, rows: Sequence[Sequence[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        tmp_path = Path(handle.name)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
    os.replace(tmp_path, path)


def compare(
    reference: Mapping[str, Mapping[str, str]],
    candidate: Mapping[str, Mapping[str, str]],
    ignored: Sequence[str],
) -> Tuple[List[List[str]], List[List[str]], bool]:
    ignored_set = set(ignored)
    unknown = sorted(ignored_set.difference(EXPECTED_HEADER))
    if unknown:
        raise ValueError(f"unknown ignored columns: {unknown}")
    compared = [column for column in EXPECTED_HEADER if column not in ignored_set]
    reference_sequences = set(reference)
    candidate_sequences = set(candidate)
    missing = sorted(reference_sequences.difference(candidate_sequences))
    unexpected = sorted(candidate_sequences.difference(reference_sequences))
    shared = sorted(reference_sequences.intersection(candidate_sequences))

    detail: List[List[str]] = [
        ["sequence", "status", "differingColumns", "referenceRank", "candidateRank"]
    ]
    mismatched = 0
    for sequence in missing:
        detail.append([sequence, "MISSING", "", reference[sequence]["rank"], ""])
    for sequence in unexpected:
        detail.append([sequence, "UNEXPECTED", "", "", candidate[sequence]["rank"]])
    for sequence in shared:
        differences = [
            column
            for column in compared
            if reference[sequence][column] != candidate[sequence][column]
        ]
        if differences:
            mismatched += 1
            detail.append(
                [
                    sequence,
                    "MISMATCH",
                    ",".join(differences),
                    reference[sequence]["rank"],
                    candidate[sequence]["rank"],
                ]
            )

    exact = not missing and not unexpected and mismatched == 0
    report = [
        ["metric", "value"],
        ["schema", SCHEMA],
        ["status", "PASS" if exact else "FAIL"],
        ["referenceRows", str(len(reference))],
        ["candidateRows", str(len(candidate))],
        ["sharedSequences", str(len(shared))],
        ["missingSequences", str(len(missing))],
        ["unexpectedSequences", str(len(unexpected))],
        ["mismatchedSequences", str(mismatched)],
        ["comparedColumns", str(len(compared))],
        ["ignoredColumns", ",".join(ignored)],
    ]
    return report, detail, exact


def parser() -> argparse.ArgumentParser:
    out = argparse.ArgumentParser(
        description="Compare canonical merged PACK* CSV files exactly"
    )
    out.add_argument("--reference", type=Path, required=True)
    out.add_argument("--candidate", type=Path, required=True)
    out.add_argument("--report", type=Path, required=True)
    out.add_argument("--detail", type=Path, required=True)
    out.add_argument(
        "--ignore-column", action="append", default=["total_time_s"]
    )
    return out


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        reference = load(args.reference.resolve())
        candidate = load(args.candidate.resolve())
        report, detail, exact = compare(
            reference, candidate, args.ignore_column
        )
        atomic_tsv(args.report.resolve(), report)
        atomic_tsv(args.detail.resolve(), detail)
    except (OSError, ValueError) as exc:
        print(f"PACK* merged comparison failed: {exc}", file=sys.stderr)
        return 2

    print(
        f"PACK* merged comparison status={'PASS' if exact else 'FAIL'} "
        f"report={args.report.resolve()} detail={args.detail.resolve()}"
    )
    return 0 if exact else 6


if __name__ == "__main__":
    raise SystemExit(main())
