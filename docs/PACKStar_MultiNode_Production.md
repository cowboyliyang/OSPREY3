# PACK* Multi-Node Production Contract

The multi-node runner scales K-star across complete sequence bundles. Each
rank keeps every partition-function calculation node-local and uses all GPUs
visible on that node. There is no cross-node DP reduction.

## Deterministic calculation identity

PACK* derives each estimator random stream from:

```text
base seed + K-star state role + exact allowed-RC signature + estimator version
```

The stream does not depend on JVM construction order, sequence shard count,
rank, node, or GPU identifier. Runs intended for comparison must use the same
`PAC_RANDOM_SEED` (the formal scripts default to `42`).

## Sequence ownership

`KStarSequenceSharding` is the single ownership implementation used by both
formal K-star execution and PACK* preflight.

- With a wild type, bundle ordinal `0` is computed on every rank because it is
  the local stability baseline.
- Every other bundle has one owner: `bundleOrdinal % shardCount`.
- Without a wild type, bundle `0` is not replicated.

Preflight v3 records the replicated bundle explicitly. Its serial and
makespan predictions include every scheduled WT copy.

## Frozen decomposition policy

An exact policy is keyed by the K-star state role and exact allowed-RC
signature.  Its value records the bounded branchwidth-improvement search
(`minDrop`, `maxDrop`, `maxMillis`) and the accepted prediction ceiling.  It is
not an estimator seed or a Slurm scheduling policy.

`run_packstar_multinode_production.slurm` enforces three ordered phases:

1. A canonical one-shard preflight sees the complete case and writes one
   shard-invariant policy.
2. Every requested rank previews that frozen policy.  The parent validates the
   complete sequence manifest and rejects the real multi-node makespan before
   any formal DP table is materialized.
3. Formal ranks load the admitted policy in `locked-policy-only` mode, run
   their sequence bundles, and enter the strict global merger.

Set `FROZEN_POLICY_IN` to the same policy file for all runs in a comparison
matrix.  This prevents each shard layout from selecting a different locally
expensive state for exact improvement.

## Result publication

After all ranks finish, `packstar_multinode_merge.py` treats `detail.tsv` as
the expected manifest and validates:

1. the canonical 41-column CSV schema on every rank;
2. every expected sequence and its owning shard;
3. no unknown, missing, or duplicate mutant sequence;
4. an identical replicated WT on every rank, excluding only local rank and
   wall-clock fields;
5. `Estimated` status and target epsilon for protein, ligand, and complex.

It then recomputes the global lower-bound ranking and writes:

```text
merged_packstar.csv
result_validation.json
result_provenance.tsv
```

If any state is incomplete, strict mode exits nonzero and writes
`merged_packstar.incomplete.csv` instead of publishing the canonical filename.
`PACKSTAR_ALLOW_INCOMPLETE_RESULTS=true` is an explicit diagnostic-only opt-in;
the JSON status is then `INCOMPLETE_ALLOWED`.

## Slurm workflow

First compile and run contract tests:

```bash
sbatch slurm/scripts/validate_packstar_multinode.slurm
```

The validation job prints a job-specific xtmp class root, so it never replaces
classes used by an already-running calculation.  Use that root for formal
jobs. A strict two-node run is:

```bash
sbatch \
  --export=ALL,CLASS_ROOT=/usr/xtmp/lz280/osprey3-gradle-build/packstar_JOB/osprey-root \
  slurm/scripts/run_packstar_multinode_production.slurm
```

Generate the canonical policy once with a one-node run, then reuse its
`frozen_policy.tsv` verbatim for the two-node, four-node, and repeat runs:

```bash
sbatch --nodes=4 \
  --export=ALL,CLASS_ROOT=/usr/xtmp/lz280/osprey3-gradle-build/packstar_JOB/osprey-root,FROZEN_POLICY_IN=/usr/xtmp/lz280/packstar_production/formal_production_gpu/CASE/frozen_policy.tsv \
  slurm/scripts/run_packstar_multinode_production.slurm
```

The legacy `run_packstar_single_gpu_node.slurm` and
`run_packstar_multigpu_multinode.slurm` remain useful for diagnosing historical
one-phase behavior; they are not the production entry point.

Compare two canonical outputs, ignoring only `total_time_s`:

```bash
sbatch \
  --export=ALL,REFERENCE=/usr/xtmp/lz280/path/a/merged_packstar.csv,CANDIDATE=/usr/xtmp/lz280/path/b/merged_packstar.csv,REPORT=/usr/xtmp/lz280/path/comparison.tsv \
  slurm/scripts/compare_packstar_merged.slurm
```

The comparator fails unless sequence membership, global rank, and every other
CSV field are exactly equal at output precision.

## Production acceptance

A configuration is production-ready only after the following Slurm checks
pass with a fixed seed and the same per-rank GPU configuration:

- single versus 2-shard;
- single versus 4-shard;
- repeated 2-shard runs;
- replicated WT agreement inside each multi-node run;
- strict convergence validation, or an explicitly documented incomplete
  scientific result that is not labeled as a production pass.
