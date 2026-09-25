"""Slurm-only, complete ID audit and two independent log-sum calculations."""
import csv
import hashlib
import json
import math
import mmap
import struct
import sys
from pathlib import Path

root, array_id, config_path, block = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), int(sys.argv[4])
with config_path.open() as f:
    config = next(csv.DictReader(f, delimiter="\t"))
total = int(config["count"])
rt = 1.9891 / 1000 * 298.15
out = root / f"merged_A{array_id}"
out.mkdir(exist_ok=False)
seen = bytearray(total)
partitions = []
finite = 0
minimum = math.inf
layout = None
with (out / "energies.f64be").open("w+b") as f:
    f.truncate(8 * total)
    with mmap.mmap(f.fileno(), 0) as energies:
        for task in range((total + block - 1) // block):
            shard = root / f"shard_A{array_id}" / f"T{task}"
            if not (shard / "READY").exists() or (shard / "exit_code").read_text().strip() != "0":
                raise RuntimeError(f"Incomplete shard {task}; no complete reference declared")
            if (shard / "config.tsv").read_bytes() != config_path.read_bytes():
                raise RuntimeError(f"Configuration mismatch {task}")
            current_layout = (shard / "layout.tsv").read_bytes()
            if layout is None:
                layout = current_layout
                (out / "layout.tsv").write_bytes(layout)
            elif layout != current_layout:
                raise RuntimeError(f"RC mapping mismatch {task}")
            with (shard / "result.tsv").open() as sf:
                info = next(csv.DictReader(sf, delimiter="\t"))
            start, end = task * block, min((task + 1) * block, total)
            if (int(info["start"]), int(info["end"]), int(info["count"]), int(info["total"])) != (start, end, end-start, total):
                raise RuntimeError(f"Manifest mismatch {task}")
            n = 0
            digest = hashlib.sha256()
            with (shard / "energies.bin").open("rb") as sf:
                while record := sf.read(16):
                    if len(record) != 16:
                        raise RuntimeError(f"Truncated binary {task}")
                    digest.update(record)
                    ident, energy = struct.unpack(">qd", record)
                    if not start <= ident < end or seen[ident]:
                        raise RuntimeError(f"Duplicate/out-of-range ID {ident}")
                    if math.isnan(energy) or energy == -math.inf:
                        raise RuntimeError(f"Invalid energy {ident}")
                    seen[ident] = 1
                    energies[8*ident:8*ident+8] = struct.pack(">d", energy)
                    minimum = min(minimum, energy)
                    finite += math.isfinite(energy)
                    n += 1
            if n != end-start:
                raise RuntimeError(f"Missing records {task}")
            partitions.append({"task": task, "sha256": digest.hexdigest(), **info})
        if sum(seen) != total or not math.isfinite(minimum):
            raise RuntimeError("Incomplete/zero-mass census")
        mass = math.fsum(math.exp((minimum-struct.unpack_from(">d", energies, 8*i)[0])/rt) for i in range(total))
        log_z = -minimum/rt+math.log(mass)
        part_logs = [float(p["logZ"]) for p in partitions]
        maximum = max(part_logs)
        log_z_parts = maximum+math.log(math.fsum(math.exp(v-maximum) for v in part_logs))
        if abs(log_z-log_z_parts) > 1e-8:
            raise RuntimeError("Independent reductions disagree")
report = {"config": config, "total": total, "finite": finite, "positive_infinity_zero_weight": total-finite,
          "rt": rt, "logZ": log_z, "logZ_from_shards": log_z_parts,
          "reduction_difference": abs(log_z-log_z_parts), "shards": partitions}
(out / "reference.json").write_text(json.dumps(report, indent=2))
(out / "reference.md").write_text(f"# Complete CPU CCD census\n\nSystem: {config['system']}\n\n"
    f"Assignments: {total}; finite: {finite}.\n\nlog(Z): {log_z:.14g}.\n\n"
    f"Independent reduction difference: {abs(log_z-log_z_parts):.4g}.\n\n"
    "Every ID is present exactly once. No pruning; positive infinity contributes zero weight.\n")
(out / "READY").write_text("Complete finite CCD-target enumeration verified.\n")
print((out / "reference.md").read_text())
