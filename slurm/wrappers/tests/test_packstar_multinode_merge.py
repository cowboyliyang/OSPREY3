from __future__ import annotations

import base64
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import packstar_multinode_merge as merger
import packstar_compare_merged as comparator


class PackStarMultinodeMergeTest(unittest.TestCase):

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self.tmp.name)
        self.result_root = self.root / "result"
        self.detail = self.result_root / "detail.tsv"
        self.output = self.result_root / "merged_packstar.csv"
        self.incomplete_output = (
            self.result_root / "merged_packstar.incomplete.csv"
        )
        self.report = self.result_root / "validation.json"
        self.provenance = self.result_root / "provenance.tsv"
        self.result_root.mkdir(parents=True)
        self.write_manifest()
        self.write_rank_csvs()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    @staticmethod
    def b64(value: str) -> str:
        return base64.b64encode(value.encode("utf-8")).decode("ascii")

    def write_manifest(self) -> None:
        header = [
            "shardIndex", "bundleOrdinal", "globalSequenceB64", "ordinal",
            "stateNameB64", "sequenceB64", "stateKeyB64", "branchwidth",
            "rootSplitEdge", "predictedSeconds", "gpuWork",
            "oocTrafficBytes", "oocTrafficAvailable", "dpSweeps",
            "adaptiveAttempted", "adaptiveAccepted",
        ]
        ownership = {0: [0, 1], 1: [1], 2: [0]}
        sequences = {0: "A1=asp", 1: "A1=ALA", 2: "A1=CYS"}
        states = ["Protein", "Ligand", "Complex"]
        with self.detail.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=header, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for bundle, owners in ownership.items():
                for owner in owners:
                    for state_offset, state in enumerate(states):
                        ordinal = bundle * 3 + state_offset
                        writer.writerow(
                            {
                                "shardIndex": owner,
                                "bundleOrdinal": bundle,
                                "globalSequenceB64": self.b64(sequences[bundle]),
                                "ordinal": ordinal,
                                "stateNameB64": self.b64(state),
                                "sequenceB64": self.b64(sequences[bundle]),
                                "stateKeyB64": self.b64(f"{state}|{bundle}"),
                                "branchwidth": 2,
                                "rootSplitEdge": 1,
                                "predictedSeconds": "1.0",
                                "gpuWork": "10",
                                "oocTrafficBytes": "0",
                                "oocTrafficAvailable": "true",
                                "dpSweeps": "6",
                                "adaptiveAttempted": "false",
                                "adaptiveAccepted": "false",
                            }
                        )

    @staticmethod
    def result_row(
        local_rank: int, sequence: str, lower: float, total_time: float
    ) -> dict[str, str]:
        row = {column: "0" for column in merger.EXPECTED_HEADER}
        row.update(
            {
                "rank": str(local_rank),
                "sequence": sequence,
                "method": "packstar",
                "target_eps": "0.683000",
                "score_log10": f"{lower + 0.5:.6f}",
                "lb_log10": f"{lower:.6f}",
                "ub_log10": f"{lower + 1.0:.6f}",
                "total_time_s": f"{total_time:.1f}",
            }
        )
        for prefix in ("prot", "lig", "comp"):
            row[f"{prefix}_qstar_lb_log10"] = "1.000000"
            row[f"{prefix}_qstar_ub_log10"] = "1.100000"
            row[f"{prefix}_status"] = "Estimated"
            row[f"{prefix}_eps"] = "0.100000"
        return row

    def write_rank_csvs(self) -> None:
        rows = {
            0: [
                self.result_row(1, "CYS", 15.0, 11.0),
                self.result_row(2, "ASP", 10.0, 11.0),
            ],
            1: [
                self.result_row(1, "ALA", 20.0, 12.0),
                self.result_row(2, "ASP", 10.0, 12.0),
            ],
        }
        for shard, shard_rows in rows.items():
            rank_dir = self.result_root / f"rank_{shard}"
            rank_dir.mkdir()
            with (rank_dir / "case_packstar.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=merger.EXPECTED_HEADER,
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(shard_rows)

    def argv(self, *extra: str) -> list[str]:
        return [
            "--result-root", str(self.result_root),
            "--detail", str(self.detail),
            "--shard-count", "2",
            "--design-id", "case",
            "--output", str(self.output),
            "--incomplete-output", str(self.incomplete_output),
            "--report", str(self.report),
            "--provenance", str(self.provenance),
            *extra,
        ]

    def read_rank_rows(self, shard: int) -> list[dict[str, str]]:
        path = self.result_root / f"rank_{shard}" / "case_packstar.csv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def overwrite_rank_rows(
        self, shard: int, rows: list[dict[str, str]]
    ) -> None:
        path = self.result_root / f"rank_{shard}" / "case_packstar.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=merger.EXPECTED_HEADER, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_valid_shards_merge_and_recompute_global_rank(self) -> None:
        self.assertEqual(0, merger.main(self.argv()))
        with self.output.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(["ALA", "CYS", "ASP"], [row["sequence"] for row in rows])
        self.assertEqual(["1", "2", "3"], [row["rank"] for row in rows])
        report = json.loads(self.report.read_text(encoding="utf-8"))
        self.assertEqual("PASS", report["status"])
        self.assertEqual(["ASP"], report["replicatedSequences"])

    def test_replicated_wild_type_mismatch_fails(self) -> None:
        rows = self.read_rank_rows(1)
        rows[1]["comp_eps"] = "0.200000"
        self.overwrite_rank_rows(1, rows)

        self.assertEqual(4, merger.main(self.argv()))
        self.assertFalse(self.output.exists())
        report = json.loads(self.report.read_text(encoding="utf-8"))
        self.assertEqual("FAIL", report["status"])
        self.assertIn("comp_eps", report["error"])

    def test_missing_owned_mutant_fails(self) -> None:
        rows = self.read_rank_rows(1)
        remaining = [row for row in rows if row["sequence"] != "ALA"]
        remaining[0]["rank"] = "1"
        self.overwrite_rank_rows(1, remaining)

        self.assertEqual(3, merger.main(self.argv()))
        self.assertFalse(self.output.exists())

    def test_incomplete_state_is_not_published_as_canonical(self) -> None:
        rows = self.read_rank_rows(1)
        rows[0]["comp_status"] = "Estimating"
        rows[0]["comp_eps"] = "1.000000"
        self.overwrite_rank_rows(1, rows)

        self.assertEqual(5, merger.main(self.argv()))
        self.assertFalse(self.output.exists())
        self.assertTrue(self.incomplete_output.is_file())
        report = json.loads(self.report.read_text(encoding="utf-8"))
        self.assertEqual("INCOMPLETE", report["status"])
        self.assertEqual(1, report["incompleteStateCount"])

    def test_incomplete_state_requires_explicit_opt_in(self) -> None:
        rows = self.read_rank_rows(1)
        rows[0]["comp_status"] = "Estimating"
        rows[0]["comp_eps"] = "1.000000"
        self.overwrite_rank_rows(1, rows)

        self.assertEqual(0, merger.main(self.argv("--allow-incomplete")))
        self.assertTrue(self.output.is_file())
        report = json.loads(self.report.read_text(encoding="utf-8"))
        self.assertEqual("INCOMPLETE_ALLOWED", report["status"])

    def test_canonical_comparator_ignores_only_time(self) -> None:
        self.assertEqual(0, merger.main(self.argv()))
        candidate = self.result_root / "candidate.csv"
        with self.output.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            row["total_time_s"] = "999.9"
        with candidate.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=merger.EXPECTED_HEADER, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

        compare_report = self.result_root / "compare.tsv"
        compare_detail = self.result_root / "compare_detail.tsv"
        compare_argv = [
            "--reference", str(self.output),
            "--candidate", str(candidate),
            "--report", str(compare_report),
            "--detail", str(compare_detail),
        ]
        self.assertEqual(0, comparator.main(compare_argv))

        rows[0]["comp_s9_onnxCalls"] = "1"
        with candidate.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=merger.EXPECTED_HEADER, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        self.assertEqual(6, comparator.main(compare_argv))


if __name__ == "__main__":
    unittest.main()
