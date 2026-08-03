#!/usr/bin/env python3
"""Validate and merge sequence-sharded PACK* K-star CSV results."""

from __future__ import annotations

import argparse
import base64
import csv
import json
import math
import os
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


SCHEMA = "packstar-multinode-validation-v1"
EXPECTED_STATES = {"Protein", "Ligand", "Complex"}
EXPECTED_HEADER = [
    "rank", "sequence", "method", "target_eps", "score_log10",
    "lb_log10", "ub_log10",
    "prot_qstar_lb_log10", "prot_qstar_ub_log10", "prot_status",
    "prot_eps", "prot_nconf", "prot_nscored", "prot_npartial",
    "prot_s9_leafGNN", "prot_s9_subtreeGNN", "prot_s9_ccdFromGNN",
    "prot_s9_onnxCalls",
    "lig_qstar_lb_log10", "lig_qstar_ub_log10", "lig_status",
    "lig_eps", "lig_nconf", "lig_nscored", "lig_npartial",
    "lig_s9_leafGNN", "lig_s9_subtreeGNN", "lig_s9_ccdFromGNN",
    "lig_s9_onnxCalls",
    "comp_qstar_lb_log10", "comp_qstar_ub_log10", "comp_status",
    "comp_eps", "comp_nconf", "comp_nscored", "comp_npartial",
    "comp_s9_leafGNN", "comp_s9_subtreeGNN", "comp_s9_ccdFromGNN",
    "comp_s9_onnxCalls", "total_time_s",
]
NON_ESTIMATOR_COLUMNS = {"rank", "total_time_s"}


class ValidationError(RuntimeError):
    def __init__(self, message: str, exit_code: int = 3):
        super().__init__(message)
        self.exit_code = exit_code


@dataclass
class ManifestBundle:
    ordinal: int
    global_sequence: str
    display_sequence: str
    owners: Set[int] = field(default_factory=set)
    state_names_by_owner: MutableMapping[int, Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )


@dataclass
class RankCsv:
    shard: int
    path: Path
    header: List[str]
    rows: Dict[str, Dict[str, str]]


def decode_b64(value: str, label: str) -> str:
    try:
        raw = base64.b64decode(value, validate=True)
        return raw.decode("utf-8")
    except (ValueError, UnicodeDecodeError) as exc:
        raise ValidationError(f"invalid {label} base64/UTF-8 value") from exc


def normalize_sequence(value: str) -> str:
    return " ".join(value.strip().upper().split())


def render_res_types(global_sequence: str) -> str:
    value = global_sequence.strip()
    if value == "(no mutable residues)":
        return normalize_sequence(value)
    residue_types: List[str] = []
    for assignment in value.split():
        if "=" not in assignment:
            raise ValidationError(
                f"cannot render preflight sequence as residue types: {value!r}"
            )
        _, residue_type = assignment.split("=", 1)
        if not residue_type:
            raise ValidationError(f"empty residue type in sequence: {value!r}")
        residue_types.append(residue_type.upper())
    if not residue_types:
        raise ValidationError(f"empty global sequence in manifest: {value!r}")
    return " ".join(residue_types)


def parse_int(value: str, label: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"invalid integer {label}={value!r}") from exc


def load_manifest(detail_path: Path, shard_count: int) -> List[ManifestBundle]:
    if not detail_path.is_file():
        raise ValidationError(f"missing preflight detail manifest: {detail_path}")

    required = {
        "shardIndex",
        "bundleOrdinal",
        "globalSequenceB64",
        "ordinal",
        "stateNameB64",
    }
    bundles: Dict[int, ManifestBundle] = {}
    seen_shard_ordinals: Set[Tuple[int, int]] = set()

    with detail_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise ValidationError(
                f"preflight detail is missing required columns: {missing}"
            )
        for line_number, row in enumerate(reader, start=2):
            shard = parse_int(row["shardIndex"], f"line {line_number} shardIndex")
            bundle_ordinal = parse_int(
                row["bundleOrdinal"], f"line {line_number} bundleOrdinal"
            )
            state_ordinal = parse_int(
                row["ordinal"], f"line {line_number} ordinal"
            )
            if shard < 0 or shard >= shard_count:
                raise ValidationError(
                    f"manifest shard {shard} is outside [0,{shard_count})"
                )
            shard_ordinal = (shard, state_ordinal)
            if shard_ordinal in seen_shard_ordinals:
                raise ValidationError(
                    f"duplicate manifest state ordinal {state_ordinal} on shard {shard}"
                )
            seen_shard_ordinals.add(shard_ordinal)

            global_sequence = decode_b64(
                row["globalSequenceB64"], "globalSequenceB64"
            )
            state_name = decode_b64(row["stateNameB64"], "stateNameB64")
            display_sequence = render_res_types(global_sequence)
            bundle = bundles.get(bundle_ordinal)
            if bundle is None:
                bundle = ManifestBundle(
                    bundle_ordinal, global_sequence, display_sequence
                )
                bundles[bundle_ordinal] = bundle
            elif (
                bundle.global_sequence != global_sequence
                or bundle.display_sequence != display_sequence
            ):
                raise ValidationError(
                    f"conflicting manifest identity for bundle {bundle_ordinal}"
                )
            bundle.owners.add(shard)
            bundle.state_names_by_owner[shard].add(state_name)

    if not bundles:
        raise ValidationError("preflight detail manifest is empty")
    ordinals = sorted(bundles)
    if ordinals != list(range(len(ordinals))):
        raise ValidationError(
            f"manifest bundle ordinals are not contiguous from zero: {ordinals}"
        )

    displays: Dict[str, int] = {}
    replicated: List[ManifestBundle] = []
    expected_all_owners = set(range(shard_count))
    for ordinal in ordinals:
        bundle = bundles[ordinal]
        previous = displays.setdefault(bundle.display_sequence, ordinal)
        if previous != ordinal:
            raise ValidationError(
                "CSV sequence rendering is not globally unique: "
                f"bundles {previous} and {ordinal} both render as "
                f"{bundle.display_sequence!r}"
            )
        if len(bundle.owners) > 1:
            replicated.append(bundle)
            if bundle.owners != expected_all_owners:
                raise ValidationError(
                    f"replicated bundle {ordinal} is not present on every shard: "
                    f"owners={sorted(bundle.owners)}"
                )
        elif len(bundle.owners) != 1:
            raise ValidationError(f"bundle {ordinal} has no owner")
        for owner in bundle.owners:
            state_names = bundle.state_names_by_owner[owner]
            if state_names != EXPECTED_STATES:
                raise ValidationError(
                    f"bundle {ordinal} shard {owner} has states "
                    f"{sorted(state_names)}, expected {sorted(EXPECTED_STATES)}"
                )

    if len(replicated) > 1:
        raise ValidationError(
            "only the wild-type bundle may be replicated; found bundles "
            + ",".join(str(bundle.ordinal) for bundle in replicated)
        )
    if replicated and replicated[0].ordinal != 0:
        raise ValidationError(
            f"replicated bundle must be wild-type ordinal 0, got {replicated[0].ordinal}"
        )
    return [bundles[ordinal] for ordinal in ordinals]


def load_rank_csv(path: Path, shard: int) -> RankCsv:
    if not path.is_file():
        raise ValidationError(f"missing shard CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        header = list(reader.fieldnames or [])
        if not header:
            raise ValidationError(f"empty CSV header: {path}")
        if header != EXPECTED_HEADER:
            raise ValidationError(
                f"CSV {path} does not match the canonical {len(EXPECTED_HEADER)}-column schema"
            )
        rows: Dict[str, Dict[str, str]] = {}
        local_ranks: List[int] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValidationError(
                    f"CSV {path}:{line_number} has extra fields"
                )
            sequence = normalize_sequence(row["sequence"])
            if not sequence:
                raise ValidationError(
                    f"CSV {path}:{line_number} has an empty sequence"
                )
            if sequence in rows:
                raise ValidationError(
                    f"CSV {path} contains duplicate sequence {sequence!r}"
                )
            local_ranks.append(
                parse_int(row["rank"], f"{path}:{line_number} rank")
            )
            normalized = {name: row[name].strip() for name in header}
            normalized["sequence"] = sequence
            rows[sequence] = normalized
        if sorted(local_ranks) != list(range(1, len(local_ranks) + 1)):
            raise ValidationError(
                f"CSV {path} local ranks are not exactly 1..{len(local_ranks)}"
            )
    return RankCsv(shard, path, header, rows)


def validate_rank_ownership(
    bundles: Sequence[ManifestBundle], rank_csvs: Sequence[RankCsv]
) -> None:
    by_sequence = {bundle.display_sequence: bundle for bundle in bundles}
    for rank_csv in rank_csvs:
        expected = {
            bundle.display_sequence
            for bundle in bundles
            if rank_csv.shard in bundle.owners
        }
        actual = set(rank_csv.rows)
        missing = sorted(expected.difference(actual))
        unexpected = sorted(actual.difference(expected))
        wrong_owner = sorted(
            sequence
            for sequence in actual.intersection(by_sequence)
            if rank_csv.shard not in by_sequence[sequence].owners
        )
        if missing or unexpected or wrong_owner:
            raise ValidationError(
                f"shard {rank_csv.shard} ownership mismatch: missing={missing}, "
                f"unexpected={unexpected}, wrongOwner={wrong_owner}"
            )


def validate_headers(rank_csvs: Sequence[RankCsv]) -> List[str]:
    expected = rank_csvs[0].header
    for rank_csv in rank_csvs[1:]:
        if rank_csv.header != expected:
            raise ValidationError(
                f"CSV header mismatch between {rank_csvs[0].path} and {rank_csv.path}"
            )
    return expected


def validate_replicas(
    bundles: Sequence[ManifestBundle],
    rank_by_shard: Mapping[int, RankCsv],
    header: Sequence[str],
) -> List[str]:
    comparison_columns = [
        column for column in header if column not in NON_ESTIMATOR_COLUMNS
    ]
    replicated_sequences: List[str] = []
    for bundle in bundles:
        if len(bundle.owners) <= 1:
            continue
        replicated_sequences.append(bundle.display_sequence)
        owner_rows = [
            rank_by_shard[owner].rows[bundle.display_sequence]
            for owner in sorted(bundle.owners)
        ]
        baseline = owner_rows[0]
        for owner, row in zip(sorted(bundle.owners)[1:], owner_rows[1:]):
            differences = [
                column
                for column in comparison_columns
                if baseline[column] != row[column]
            ]
            if differences:
                raise ValidationError(
                    f"replicated sequence {bundle.display_sequence!r} differs on "
                    f"shard {owner}; columns={differences}",
                    exit_code=4,
                )
    return replicated_sequences


def parse_finite(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValidationError(f"invalid floating-point {label}={value!r}") from exc
    if not math.isfinite(parsed):
        raise ValidationError(f"non-finite {label}={value!r}")
    return parsed


def find_incomplete_states(
    selected_rows: Sequence[Tuple[ManifestBundle, int, Dict[str, str]]],
    epsilon_tolerance: float,
) -> List[Dict[str, object]]:
    incomplete: List[Dict[str, object]] = []
    for bundle, source_shard, row in selected_rows:
        target = parse_finite(
            row["target_eps"], f"{bundle.display_sequence} target_eps"
        )
        if target < 0.0:
            raise ValidationError(
                f"negative target epsilon for {bundle.display_sequence}: {target}"
            )
        for prefix, state_name in (
            ("prot", "Protein"),
            ("lig", "Ligand"),
            ("comp", "Complex"),
        ):
            status = row[f"{prefix}_status"]
            epsilon_text = row[f"{prefix}_eps"]
            epsilon = None
            reason = None
            if status != "Estimated":
                reason = "status"
            elif epsilon_text:
                epsilon = parse_finite(
                    epsilon_text,
                    f"{bundle.display_sequence} {state_name} epsilon",
                )
                if epsilon > target + epsilon_tolerance:
                    reason = "epsilon"
            if reason is not None:
                incomplete.append(
                    {
                        "sequence": bundle.display_sequence,
                        "bundleOrdinal": bundle.ordinal,
                        "sourceShard": source_shard,
                        "state": state_name,
                        "status": status,
                        "epsilon": epsilon,
                        "targetEpsilon": target,
                        "reason": reason,
                    }
                )
    return incomplete


def lower_bound_sort_key(
    item: Tuple[ManifestBundle, int, Dict[str, str]]
) -> Tuple[int, float, int]:
    bundle, _, row = item
    value = row["lb_log10"]
    if not value:
        return (0, 0.0, bundle.ordinal)
    try:
        number = float(value)
    except ValueError as exc:
        raise ValidationError(
            f"invalid lower bound for {bundle.display_sequence}: {value!r}"
        ) from exc
    if math.isnan(number):
        return (0, 0.0, bundle.ordinal)
    return (1, -number, bundle.ordinal)


def atomic_write_csv(path: Path, header: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        tmp_path = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp_path, path)


def atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        tmp_path = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_path, path)


def atomic_write_provenance(
    path: Path,
    selected_rows: Sequence[Tuple[ManifestBundle, int, Dict[str, str]]],
) -> None:
    header = [
        "globalRank",
        "sequence",
        "bundleOrdinal",
        "sourceShard",
        "sourceLocalRank",
        "ownerShards",
        "sourceTotalTimeSeconds",
    ]
    rows = []
    for global_rank, (bundle, source_shard, row) in enumerate(
        selected_rows, start=1
    ):
        rows.append(
            {
                "globalRank": str(global_rank),
                "sequence": bundle.display_sequence,
                "bundleOrdinal": str(bundle.ordinal),
                "sourceShard": str(source_shard),
                "sourceLocalRank": row["rank"],
                "ownerShards": ",".join(str(owner) for owner in sorted(bundle.owners)),
                "sourceTotalTimeSeconds": row["total_time_s"],
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        tmp_path = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp_path, path)


def merge(args: argparse.Namespace) -> Tuple[int, Dict[str, object]]:
    result_root = args.result_root.resolve()
    bundles = load_manifest(args.detail.resolve(), args.shard_count)
    rank_csvs = [
        load_rank_csv(
            result_root / f"rank_{shard}" / f"{args.design_id}_packstar.csv",
            shard,
        )
        for shard in range(args.shard_count)
    ]
    header = validate_headers(rank_csvs)
    validate_rank_ownership(bundles, rank_csvs)
    rank_by_shard = {rank_csv.shard: rank_csv for rank_csv in rank_csvs}
    replicated_sequences = validate_replicas(
        bundles, rank_by_shard, header
    )

    selected_rows: List[Tuple[ManifestBundle, int, Dict[str, str]]] = []
    for bundle in bundles:
        source_shard = min(bundle.owners)
        selected_rows.append(
            (
                bundle,
                source_shard,
                dict(rank_by_shard[source_shard].rows[bundle.display_sequence]),
            )
        )
    incomplete = find_incomplete_states(selected_rows, args.epsilon_tolerance)
    selected_rows.sort(key=lower_bound_sort_key)

    output_path = args.output.resolve()
    validation_status = "PASS"
    exit_code = 0
    if incomplete and not args.allow_incomplete:
        output_path = args.incomplete_output.resolve()
        validation_status = "INCOMPLETE"
        exit_code = 5
    elif incomplete:
        validation_status = "INCOMPLETE_ALLOWED"

    merged_rows: List[Dict[str, str]] = []
    for global_rank, (bundle, _, row) in enumerate(selected_rows, start=1):
        row["rank"] = str(global_rank)
        row["sequence"] = bundle.display_sequence
        merged_rows.append(row)
    atomic_write_csv(output_path, header, merged_rows)
    atomic_write_provenance(args.provenance.resolve(), selected_rows)

    report: Dict[str, object] = {
        "schema": SCHEMA,
        "status": validation_status,
        "structuralValid": True,
        "complete": not incomplete,
        "allowIncomplete": bool(args.allow_incomplete),
        "shardCount": args.shard_count,
        "sequenceCount": len(bundles),
        "replicatedSequences": replicated_sequences,
        "incompleteStateCount": len(incomplete),
        "incompleteStates": incomplete,
        "mergedCsv": str(output_path),
        "provenanceTsv": str(args.provenance.resolve()),
        "ignoredReplicaComparisonColumns": sorted(NON_ESTIMATOR_COLUMNS),
    }
    return exit_code, report


def parser() -> argparse.ArgumentParser:
    out = argparse.ArgumentParser(
        description="Validate and merge PACK* rank CSV files"
    )
    out.add_argument("--result-root", type=Path, required=True)
    out.add_argument("--detail", type=Path, required=True)
    out.add_argument("--shard-count", type=int, required=True)
    out.add_argument("--design-id", required=True)
    out.add_argument("--output", type=Path, required=True)
    out.add_argument("--incomplete-output", type=Path, required=True)
    out.add_argument("--report", type=Path, required=True)
    out.add_argument("--provenance", type=Path, required=True)
    out.add_argument("--allow-incomplete", action="store_true")
    out.add_argument("--epsilon-tolerance", type=float, default=5.0e-7)
    return out


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.shard_count < 1:
        print("PACK* result validation failed: shard count must be positive", file=sys.stderr)
        return 2
    if args.epsilon_tolerance < 0.0 or not math.isfinite(args.epsilon_tolerance):
        print("PACK* result validation failed: invalid epsilon tolerance", file=sys.stderr)
        return 2

    try:
        exit_code, report = merge(args)
    except ValidationError as exc:
        report = {
            "schema": SCHEMA,
            "status": "FAIL",
            "structuralValid": False,
            "complete": False,
            "error": str(exc),
        }
        try:
            atomic_write_json(args.report.resolve(), report)
        except OSError as report_exc:
            print(
                f"PACK* result validation report write failed: {report_exc}",
                file=sys.stderr,
            )
        print(f"PACK* result validation failed: {exc}", file=sys.stderr)
        return exc.exit_code
    except OSError as exc:
        print(f"PACK* result validation I/O failure: {exc}", file=sys.stderr)
        return 3

    atomic_write_json(args.report.resolve(), report)
    print(
        "PACK* result validation "
        f"status={report['status']} sequences={report['sequenceCount']} "
        f"shards={report['shardCount']} incompleteStates={report['incompleteStateCount']} "
        f"mergedCsv={report['mergedCsv']}"
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
