#!/usr/bin/env python3
"""Slurm-only paired end-to-end timing, including the cost of automatic preflight."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time


def main():
    assert os.environ.get("SLURM_JOB_ID"), "Submit through Slurm"
    root = Path(os.environ["AUTO_BENCH_ROOT"])
    repo = Path(os.environ["AUTO_BENCH_SOURCE"])
    cp = Path(os.environ["AUTO_PREFLIGHT_CLASSPATH"]).read_text().strip()
    base = {}
    for line in (repo / "slurm/h200/production.properties").read_text().splitlines():
        if line and not line.startswith("#"):
            key, value = line.split("=", 1)
            base[key] = value
    cpus = int(os.environ["SLURM_CPUS_PER_TASK"])
    gpus = len(os.environ["CUDA_VISIBLE_DEVICES"].split(","))
    base.update({
        "branchdp.dp.gpu.maxGpus": str(gpus),
        "branchdp.dp.parallel.threads": str(cpus),
        "packstar.dp.parallel.threads": str(cpus),
        "packstar.pac.sampling.threads": str(cpus),
        "packstar.pac.frequencySeverity.tripleEta": "false",
        "packstar.pac.randomSeed": "20260926",
        "branchdp.cutoff.residualBudget": "0.5",
        "osprey.bench.method": "packstar",
        "osprey.bench.numCPUs": str(cpus),
    })
    specs = Path("/usr/xtmp/lz280/bench_comparison/design_specs_prepped.csv")
    with specs.open() as stream:
        rows = {r[0]: r for r in csv.reader(line for line in stream if not line.startswith("#"))}
    results = []

    def save():
        (root / "timings.json").write_text(json.dumps(results, indent=2) + "\n")

    for index, design in enumerate(os.environ.get("DESIGNS", "1a0r,2xgy").split(",")):
        row = rows[design]
        pdb = Path(f"/usr/xtmp/lz280/dance_bench/pdbs_prepped/{row[1]}/{row[1]}.min.reduce.renum.pdb")
        cache = Path("/usr/xtmp/lz280/bench_comparison/results/emat_cache") / design
        preflight_policy = root / design / "preflight/policy.tsv"
        record = {"design": design, "node": os.environ.get("SLURMD_NODENAME"),
                  "gpus": gpus, "cpus": cpus, "pdb_sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
                  "phases": {}}
        results.append(record)
        source_root = os.environ.get("PREFLIGHT_SOURCE")
        if source_root:
            source_runs = json.loads((Path(source_root) / "timings.json").read_text())
            source_run = next(r for r in source_runs if r["design"] == design)
            assert source_run["pdb_sha256"] == record["pdb_sha256"]
            previous = source_run["phases"]["preflight"]
            assert previous["exit_code"] == 0
            record["phases"]["preflight"] = previous
            record["preflight_source"] = source_root
            preflight_policy = Path(previous["output"]) / "policy.tsv"
            assert preflight_policy.is_file()
        optimized_first = index % 2 != 0 or os.environ.get("OPTIMIZED_FIRST") == "1"
        phases = (["preflight", "optimized", "baseline"] if optimized_first
                  else ["preflight", "baseline", "optimized"])
        if os.environ.get("PHASES"):
            phases = os.environ["PHASES"].split(",")
            assert all(p in ("preflight", "baseline", "optimized") for p in phases)
        timeout_seconds = int(os.environ.get("PHASE_TIMEOUT_SECONDS", "7200"))
        assert timeout_seconds > 0
        save()
        for phase in phases:
            out = root / design / phase
            out.mkdir(parents=True)
            (out / "tmp").mkdir()
            emat = out / "emat_cache" / design
            emat.mkdir(parents=True)
            hashes = {}
            for state in ("protein", "ligand", "complex"):
                for old, new in (("rigid", "rigid"), ("min", "minimizing")):
                    source = cache / f"export.{state}.{old}.dat"
                    target = emat / f"packstar.{state}.{new}.dat"
                    shutil.copyfile(source, target)
                    hashes[target.name] = hashlib.sha256(target.read_bytes()).hexdigest()
            props = dict(base)
            props.update({"java.io.tmpdir": str(out / "tmp"),
                          "branchdp.dp.mmap.dir": str(out / "dp_mmap"),
                          "osprey.bench.outputDir": str(out), "osprey.bench.designId": design,
                          "osprey.bench.pdbPath": str(pdb), "osprey.bench.mutable": row[5],
                          "osprey.bench.flexible": row[6],
                          "packstar.pac.frequencySeverity.outputDir": str(out / "adaptive_frequency_severity")})
            if phase != "baseline":
                props["packstar.admission.mode"] = "auto"
            if phase == "preflight":
                props["osprey.bench.packstarPreflightOnly"] = "true"
                props["packstar.admission.policyOut"] = str(preflight_policy)
            if phase == "optimized":
                props["packstar.admission.policyIn"] = str(preflight_policy)
            cmd = [str(Path(os.environ["JAVA_HOME"]) / "bin/java"), "--add-modules", "jdk.incubator.foreign",
                   "--add-opens", "java.base/java.util=ALL-UNNAMED",
                   "--add-opens", "java.base/java.lang=ALL-UNNAMED",
                   "--add-opens", "java.base/java.lang.invoke=ALL-UNNAMED",
                   "-Xmx700g", "-Xms1g", "-XX:-UseSuperWord", f"-XX:ActiveProcessorCount={cpus}"]
            cmd += [f"-D{k}={v}" for k, v in sorted(props.items())]
            cmd += ["-cp", cp, "edu.duke.cs.osprey.markstar.bench.GenericPDBBench"]
            (out / "manifest.json").write_text(json.dumps({"command": cmd, "properties": props,
                    "energy_matrix_sha256": hashes, "classpath_file": os.environ["AUTO_PREFLIGHT_CLASSPATH"]}, indent=2) + "\n")
            print(f"BEGIN design={design} phase={phase}", flush=True)
            start = time.monotonic()
            with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
                try:
                    run = subprocess.run(cmd, cwd=out, stdout=stdout, stderr=stderr, timeout=timeout_seconds)
                    code = run.returncode
                except subprocess.TimeoutExpired:
                    code = 124
            seconds = time.monotonic() - start
            stats = {"seconds": seconds, "exit_code": code, "output": str(out),
                     "node": os.environ.get("SLURMD_NODENAME"), "timeout_seconds": timeout_seconds}
            result_file = out / f"{design}_packstar.csv"
            if result_file.exists():
                with result_file.open() as stream:
                    stats["sequence_results"] = list(csv.DictReader(stream))
            record["phases"][phase] = stats
            save()
            print(f"END design={design} phase={phase} seconds={seconds:.3f} exit={code}", flush=True)
            if phase == "preflight" and code:
                break
        phases_done = record["phases"]
        if all(k in phases_done and phases_done[k]["exit_code"] == 0 for k in ("preflight", "baseline", "optimized")):
            pre, baseline, optimized = (phases_done[k]["seconds"] for k in ("preflight", "baseline", "optimized"))
            saving = baseline - optimized
            record.update({"preflight_plus_optimized_seconds": pre + optimized,
                           "net_saved_seconds": saving - pre,
                           "net_speedup": baseline / (pre + optimized),
                           "break_even_uses": math.ceil(pre / saving) if saving > 0 else None})
            save()
            print("RESULT " + json.dumps({k: v for k, v in record.items() if k != "phases"}), flush=True)


if __name__ == "__main__":
    main()
