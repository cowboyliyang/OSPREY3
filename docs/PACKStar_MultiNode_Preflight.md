# PACK* sequence-bundle multi-node preflight

The multi-node preflight unit is one global K* sequence, not one branch-DP
edge or one M-state range. A rank owns a deterministic subset of global
sequences:

```text
rank = globalSequenceOrdinal % shardCount
```

For each owned sequence, the rank reconstructs the three state-specific
requests locally:

```text
protein(filtered sequence) -> ligand(filtered sequence) -> complex(filtered sequence)
```

The preflight currently previews all three requests. This is conservative: a
formal K* run can short-circuit ligand or complex after a stability check, so
the runtime may be faster than the preflight certificate.

Each rank independently reads the shared PDB and immutable emat exports. It
does not exchange DP tables. It writes:

```text
rank_<r>/shard_<r>.tsv
rank_<r>/policy.tsv
```

The reducer validates that every sequence bundle appears exactly once, merges
non-conflicting exact policies, and reports both:

```text
serialPredictedCaseHours = sum(all rank sequence workloads)
multiNodeMakespanHours   = max(rank sequence workload)
```

`ADMISSION_METRIC=serial` is the safe default while the formal K* runner is
single-node. Use `ADMISSION_METRIC=makespan` only with a sequence-bundle formal
scheduler using the same node assignment model.

The branch decomposition and all M-state/lambda-state DP work remain
node-local. Existing CPU sharding, mmap/file-backed tables, multi-GPU M-state
splitting, and GPU child-table slicing continue to operate inside each node.

Example array submission for the 38-design × 3-budget matrix:

```bash
sbatch --array=0-113%8 slurm/scripts/run_packstar_preflight_multinode.slurm
```

Override `--nodes`, `SHARDS`, `THREADS`, or the hardware model properties in
the environment when calibrations change.
