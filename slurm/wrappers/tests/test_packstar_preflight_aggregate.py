from __future__ import annotations

import base64
import csv
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


class PackStarPreflightAggregateTest(unittest.TestCase):

    HEADER = [
        "bundleOrdinal", "globalSequenceB64", "ordinal", "stateNameB64",
        "sequenceB64", "stateKeyB64", "branchwidth", "rootSplitEdge",
        "predictedSeconds", "gpuWork", "oocTrafficBytes",
        "oocTrafficAvailable", "dpSweeps", "adaptiveAttempted",
        "adaptiveAccepted",
    ]

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.root = Path(self.tmp.name)
        self.script = (
            Path(__file__).resolve().parents[2]
            / "scripts"
            / "aggregate_packstar_preflight_sequence.sh"
        )
        self.write_fixture()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    @staticmethod
    def b64(value: str) -> str:
        return base64.b64encode(value.encode("utf-8")).decode("ascii")

    def shard_rows(self, bundles: list[int]) -> list[dict[str, str]]:
        rows = []
        states = ["Protein", "Ligand", "Complex"]
        for bundle in bundles:
            global_sequence = f"A1={'asp' if bundle == 0 else 'MUT' + str(bundle)}"
            for offset, state in enumerate(states):
                rows.append(
                    {
                        "bundleOrdinal": str(bundle),
                        "globalSequenceB64": self.b64(global_sequence),
                        "ordinal": str(bundle * 3 + offset),
                        "stateNameB64": self.b64(state),
                        "sequenceB64": self.b64(global_sequence),
                        "stateKeyB64": self.b64(f"{state}|{bundle}"),
                        "branchwidth": "2",
                        "rootSplitEdge": "1",
                        "predictedSeconds": "1.0",
                        "gpuWork": "10",
                        "oocTrafficBytes": "0",
                        "oocTrafficAvailable": "true",
                        "dpSweeps": "6",
                        "adaptiveAttempted": "false",
                        "adaptiveAccepted": "false",
                    }
                )
        return rows

    def write_shard(self, shard: int, bundles: list[int]) -> None:
        rank_root = self.root / f"rank_{shard}"
        rank_root.mkdir(parents=True, exist_ok=True)
        rows = self.shard_rows(bundles)
        path = rank_root / f"shard_{shard}.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            handle.write("# packstar preflight shard v3\n")
            handle.write(f"# shardIndex={shard}\n")
            handle.write("# shardCount=2\n")
            handle.write("# totalStates=9\n")
            handle.write("# totalBundles=3\n")
            handle.write("# statesPerBundle=3\n")
            handle.write(f"# assignedBundles={len(bundles)}\n")
            handle.write("# replicatedBundleOrdinal=0\n")
            handle.write("# globalSlaHours=10.0\n")
            handle.write("# localSlaHours=10.0\n")
            writer = csv.DictWriter(
                handle,
                fieldnames=self.HEADER,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        (rank_root / "policy.tsv").write_text(
            "# packstar exact-policy merged dump v1\n", encoding="utf-8"
        )

    def write_fixture(self) -> None:
        self.write_shard(0, [0, 2])
        self.write_shard(1, [0, 1])

    def run_aggregate(self) -> subprocess.CompletedProcess[str]:
        env = os.environ.copy()
        env.update(
            {
                "RESULT_ROOT": str(self.root),
                "SHARD_COUNT": "2",
                "SLA_HOURS": "10",
                "ADMISSION_METRIC": "makespan",
            }
        )
        return subprocess.run(
            ["bash", str(self.script)],
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def summary(self) -> dict[str, str]:
        values = {}
        for line in (self.root / "summary.tsv").read_text(
            encoding="utf-8"
        ).splitlines():
            key, value = line.split("\t", 1)
            values[key] = value
        return values

    def test_replicated_wild_type_counts_actual_scheduled_work(self) -> None:
        result = self.run_aggregate()
        self.assertEqual(0, result.returncode, result.stderr)
        summary = self.summary()
        self.assertEqual("9", summary["totalStates"])
        self.assertEqual("12", summary["scheduledStates"])
        self.assertEqual("4", summary["scheduledBundles"])
        self.assertEqual("0", summary["replicatedBundleOrdinal"])

    def test_conflicting_replicated_preflight_row_fails(self) -> None:
        path = self.root / "rank_1" / "shard_1.tsv"
        text = path.read_text(encoding="utf-8")
        text = text.replace("\t1.0\t10\t", "\t2.0\t10\t", 1)
        path.write_text(text, encoding="utf-8")

        result = self.run_aggregate()
        self.assertNotEqual(0, result.returncode)
        self.assertIn("replicated WT", result.stderr)


if __name__ == "__main__":
    unittest.main()
