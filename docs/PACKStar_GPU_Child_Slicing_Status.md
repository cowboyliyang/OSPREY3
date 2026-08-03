# PACK* GPU Child-Slicing Status

Last updated: 2026-07-11 (3k3q formal GPU runs started; SPACK\* small-M audit added)

## Scope and current status

This work addresses branch-DP lambda edges whose completed child DP tables do
not fit in one GPU. The acceptance target is now stronger than fixing the 3bua
bottleneck: **every structurally supported shape must have a bounded
out-of-core GPU path**. That includes shapes where no child can remain resident
and shapes where even one complete child row does not fit in VRAM.

The bounded arbitrary-shape GPU/OOC implementation is numerically covered:
full-resident, flattened-slice, resident/streamed, bilateral row tiling, and
row-internal lambda tiling all have GPU regression evidence. The two-A5000
regression (`12025816`) passed all 40 tests, and the focused OOC audit
(`12027760`) passed the extreme-value, NaN/`-Infinity`, allocation-audit,
bilateral/lambda-tiling, planner, and multi-GPU cases.

The host-lifecycle blocker exposed by the first six submissions is fixed in the
working tree: completed PACK* partition functions close their backends,
rooted-tree DP tables/planning arrays release their large memory, TiB-scale
tables can use file-backed `mmap`, and root candidates are exhaustively scored
from allocation-free shape estimates under both host and GPU budgets.

The remaining release blocker is now **production admission**, not OOC
correctness. A case must not start allocating large tables merely because its
shape is executable; it must first prove that the complete WT+mutant workload
fits a user-specified whole-case time SLA on the declared hardware. That gate is
implemented, compiled, and covered by synthetic whole-case regressions. All
three hardware calibrations have now completed successfully. The three real
3bua/3k3q allocation-free preflights have completed: 3bua is rejected against
the 336-hour SLA on all three hardware profiles, while 3k3q is admitted on all
three. The immutable artifact for the admitted 3k3q runs is frozen, and two
formal 3k3q GPU jobs are now running. The 3bua production replacement remains
blocked by admission. The old production
jobs are not accepted performance evidence and must not be resubmitted from
the old artifact.

Terminology used below is intentionally precise:

- **Double buffering**: alternate two streamed host/device buffers so gather +
  H2D for block N+1 overlaps the kernel for block N. This is a pipeline
  optimization, not an additional tiling dimension.
- **Bilateral/multi-child row tiling**: when no child can remain resident, tile
  two or more child M/row projections at the same time.
- **Lambda tiling**: split the lambda domain inside a child row and merge
  per-output stable log-sum-exp partials. This is what removes the 3k3q
  single-row VRAM floor.

### Scope boundary: lambda tiling is not multi-GPU lambda partitioning

The current implementation has bounded lambda tiling, but it must not be
interpreted as distributing one output row's lambda domain across GPUs.
`DPGpuFullDP.compute(...)` calls `chooseGpuCount(...)` before selecting the
full-resident, child-sliced, hybrid, or bounded-OOC plan. The device count is
still capped primarily from `mStateCount / minMStatesPerGpu` (default `4096`).
Normal multi-GPU execution splits M-output ranges. Sliced, hybrid, and OOC
execution can subsequently split union rows, free-M states, or M boxes, but
they inherit that device count. Each assigned worker processes the lambda
boxes required for its own output range and performs the stable partial merge.

This is a performance gap rather than an OOC-correctness gap. It is visible in
the older SPACK\* all-20 logs, which used the shared PACK\* CUDA backend:

- `mStates=114, lambdaStates=303807105` (about `3.46e10` loop iterations) and
  `mStates=45, lambdaStates=40507614` ran as single-GPU edges because M was
  below the default per-GPU threshold;
- larger-M edges such as `mStates=161595` and `mStates=1928934` did use all
  eight A5000s;
- the three-position log contains 74 outer single-GPU and 22 outer multi-GPU
  PACK\* DP completions; the partial four-position log contains 1,032 and 515,
  respectively.

The resolution should extend the unified planner rather than add an
SPACK\*-specific kernel:

1. Choose GPU count and partition dimension together from predicted makespan,
   using lambda work (`M * lambda`), shape-dependent throughput, child
   replication/gather bytes, output traffic, launch cost, and device memory.
2. Permit a small M range to use multiple GPUs when each M output has enough
   lambda work to amortize replicated inputs and launch overhead. A low
   `minMStatesPerGpu` override is useful for validation, but is not a general
   cost model.
3. Add a true lambda-box multi-GPU plan for cases where M has too little
   parallelism. Each GPU should compute lower and upper log-sum-exp partials
   for a disjoint lambda range, followed by one numerically stable cross-device
   merge per output. Reuse the current lambda-box and extreme-value merge
   primitives so memory remains bounded.
4. Compare M, union-row, free-M, M-box, and lambda-box plans before allocation;
   record the selected dimension, predicted/actual traffic, device occupancy,
   and makespan. Root selection and admission must consume this same plan.
5. Add multi-GPU regressions against CPU/full-resident references for finite,
   all-`-Infinity`, mixed-infinity, NaN-rejection, bilateral-child, and
   row-internal-lambda cases.

This work may materially accelerate small-M/high-lambda PACK\* states, but it
does not by itself solve SPACK\*'s denominator `-Infinity`, frontier requeue,
oracle-count, or PAC-accounting problems. The old completed three-position run
spent about 43 seconds in outer GPU DP timers out of 19,498.6 seconds total, so
SPACK\* scheduling and repeated estimator/CCD work remain higher-leverage
whole-run targets there.

## 2026-07-11 verified execution status

This is the current state verified from the working tree, Gradle reports, and
Slurm accounting/logs. It supersedes the earlier pending estimates in this
document.

| Gate | Actual evidence | State |
|---|---|---|
| Bounded OOC GPU path | Earlier GPU/OOC validation jobs `12025816` and `12027760`; current source/class/kernel hashes reused by all three calibrations | complete |
| Admission implementation | `BranchDpAdmission`, `PackStarAdmissionDecision`, `PackStarCasePreflight`, retained exact-policy ceiling, and pre-allocation dry-run path compile | complete |
| Targeted whole-case regressions | 19 selected tests passed: 6 decomposition, 8 admission accounting, 1 memory-release, 4 decision-core cases | complete |
| Hardware calibration | `12036821` A5000, `12036822` Titan V, `12036823` RTX Pro 6000; both leaf and OOC tests PASSED on every job | complete |
| Real 3bua/3k3q preflight | A5000 `12046063`, Pro 6000 `12046064`, Titan V `12046065`; all completed in `grisman` under `Account=grisman`: 3bua rejected, 3k3q admitted | complete with 3bua blocker |
| Immutable artifact and production replacement | Artifact `packstar_gpu_ooc_prod_20260711_3k3q_pf` frozen; formal 3k3q jobs `12048780` (4x Pro 6000) and `12048781` (4x Titan V) running under `Account=grisman`; old `12027706` canceled | 3k3q running; 3bua blocked |

The local regression command was:

```text
./gradlew test --offline --no-daemon --max-workers=1 \
  -Dorg.gradle.jvmargs=-Xmx3g -DtestMaxHeap=1g \
  --tests edu.duke.cs.osprey.packstar.TestPackStarAdmissionDecision \
  --tests edu.duke.cs.osprey.branchdp.TestBranchDpAdmission \
  --tests edu.duke.cs.osprey.branchdp.TestBranchDecompositionStrategies \
  --tests edu.duke.cs.osprey.branchdp.TestRootedTreeMemoryRelease
```

It completed successfully in 2m58s. The calibration jobs all built from
`git_head=ad102d1b8ea6732811696ac9c1dc52b6f05ee343` and recorded the same
source/class/CUDA artifact hashes:

```text
DPGpuFullDP.java       f083ffe130098f46240d130cf39b4f66b4d702aac8db72377e60619de7c3d553
DPGpuOutOfCore.java    2ec8e2f59dbc25118e854cf901668ff32663e6bddd60d86b14de97aa490aca73
DPGpuFullDP.class      262bb0625a9fa64d74020970d94b194d8c062a8ab3cef7fb76026986a4f1eb6a
DPGpuOutOfCore.class   0167c6adadf4b7f994b3639ec1f39c6dfaab5536dea3fb3e3c31969049742a55
cuda/dp.bin            9302f5b6831adbccfc54a22ce08cc459f93d2b3e7bd24b8150e0b322108c7884
```

### Completed hardware calibration

The leaf value is the aggregate warm rate over the visible GPUs; the
per-GPU value is what admission uses. OOC is the aggregate measured traffic
rate for the four-M-box calibration shape. Production preflight values were
rounded downward from the sustained warm measurements; the configured 1.5
safety factor is applied separately by admission.

| Hardware / job | GPUs | Leaf warm aggregate | Leaf rate used per GPU | OOC warm aggregate | OOC rate used |
|---|---:|---:|---:|---:|---:|
| A5000 / `12036821` | 8 | 4,917,944,601 work/s | 614,000,000 work/s | 170,310,296 B/s | 170,000,000 B/s |
| Titan V / `12036822` | 4 | 6,791,072,918 work/s | 1,697,000,000 work/s | 190,455,994 B/s | 190,000,000 B/s |
| RTX Pro 6000 / `12036823` | 3 | 10,478,206,258 work/s | 3,492,000,000 work/s | 349,035,825 B/s | 349,000,000 B/s |

All six calibration tests passed. The Slurm accounting state for the three
jobs is `COMPLETED 0:0`.

### Real preflight settings now in the queue

All three jobs use `Account=grisman`, the `grisman` partition,
`caseSlaHours=336`, `softStateHours=24`,
`safetyFactor=1.5`, initial exact timeout 120 s, final exact timeout 300 s,
`finalMaxStates=3`, `preflightOnly=true`, and no formal K* run. Their per-GPU
memory budgets are 21.5 GiB for A5000, 10 GiB for Titan V, and 85 GiB for
RTX Pro 6000. The jobs completed as CPU-only dry runs using the calibrated
rates and GPU-memory budgets; they did not allocate CUDA devices or formal DP
tables. The resulting whole-case estimates are:

| Hardware | 3bua result | 3k3q result |
|---|---:|---:|
| A5000 | rejected, `978.6652 h` after final pass | admitted, `53.5926 h` |
| RTX Pro 6000 | rejected, `474.7665 h` after final pass | admitted, `19.4884 h` |
| Titan V | rejected, `949.8192 h` after final pass | admitted, `48.6246 h` |

The wrapper jobs report Slurm `FAILED 1:0` because they return nonzero when
any case is rejected; this is the expected admission result, not a node or
CUDA failure. Each 3bua rejection occurred before DP-table materialization.
The 3k3q preflight completed with `preflight-only complete` on every profile.

### 3k3q formal GPU runs now active

The formal runs use the frozen artifact above, `preflightOnly=false`, adaptive
decomposition, predicted root selection, calibrated admission settings, and
four GPUs on each node:

| Job | Hardware | Slurm placement | State |
|---|---|---|---|
| `12048780` | 4x RTX Pro 6000 | `compsci-gpu`, `compsci-cluster-fitz-48` | running |
| `12048781` | 4x Titan V | `grisman`, `jerry2` | running |

Both jobs use `Account=grisman`. These are formal K*/PACK* runs, not
allocation-free preflights.

## 2026-07-11 production-admission checkpoint

### Required production behavior

Admission is a two-level bounded workflow:

1. Build a cheap weighted-Hicks decomposition for one unique K* state.
2. Exhaustively preview every root without materializing enumeration arrays or
   DP tables. Select the root with the lowest calibrated compute-plus-OOC time.
3. If that state exceeds the soft time budget, run one bounded exact
   decomposition-improvement attempt, exhaustively select the root again, and
   keep only a strictly faster plan.
4. Repeat the preview for every unique WT/mutant state (about 79 for the current
   cases, after de-duplicating filtered unbound sequences), then sum their
   predicted times.
5. If the whole case exceeds the hard SLA, run a final bounded exact attempt
   only for the largest configured contributors and sum the case again.
6. If the new sum still exceeds the SLA, abort before any TiB-scale `mmap`,
   sharded DP table, or CUDA workspace is materialized. Otherwise retain every
   accepted exact policy by exact RC identity; the real run must reproduce its
   preflight prediction ceiling or abort at the same pre-materialization gate.

The hard SLA is deliberately a **whole-case** limit. It is not applied
independently to each state, and decomposition search never retries without a
configured state count and time bound.

### Prediction accounting

`gpuWork` and `estimatedSlicedTrafficBytes` are the costs of one full-tree DP
sweep. Admission multiplies both by the conservative maximum number of sweeps
that the configured PACK* estimator can execute:

- default PACK*: 1 initial `p_m` DP + 1 corrected `p_eta` DP + at most 4
  distribution-shift refinement DPs = **6 sweeps per state**;
- `packstar.pac.iterate=false` or `etaEnabled=false`: 2 sweeps;
- `packstar.pac.restoreDP=true`: add one sweep;
- `packstar.admission.dpSweeps=N`: explicit audited override.

Early epsilon convergence, stable pilots, or empty sampling sets can make the
real count lower. Admission intentionally uses the configured maximum, not an
optimistic average. Compute and OOC time are added (a conservative assumption
when they overlap), then multiplied by `packstar.admission.safetyFactor`.
Missing throughput calibration, missing OOC throughput when traffic is needed,
unavailable traffic estimation, or any GPU-unsupported edge produces an
infinite prediction and therefore fails closed.

The primary production properties are:

```text
-Dpackstar.decomp.strategy=adaptive
-Dpackstar.rootSplit=predicted
-Dpackstar.admission.gpuWorkPerSecondPerGpu=<calibrated lower-bound rate>
-Dpackstar.admission.gpuCount=<allocated GPU count>
-Dpackstar.admission.oocBytesPerSecond=<calibrated lower-bound rate>
-Dpackstar.admission.safetyFactor=<conservative factor, >=1>
-Dpackstar.admission.softStateHours=<bounded exact trigger>
-Dpackstar.admission.caseSlaHours=336
-Dpackstar.admission.finalMaxStates=<bounded contributor count>
-Dpackstar.admission.finalExactMaxMillis=<per-contributor exact timeout>
```

Equivalent `branchdp.*` keys remain fallback aliases for PACK*. A positive
`caseSlaHours` is the opt-in switch; with it enabled, adaptive decomposition and
finite calibrated hardware values are mandatory.

### Evidence completed at this checkpoint

- Allocation-free exhaustive root preview, host-budget estimates, full-tree
  exact `BigInteger` GPU-work sums, and OOC traffic estimates are implemented
  before root materialization.
- The worst observed 3bua state used a 116-second bounded exact search to reduce
  branchwidth 12 to 11 and one-sweep GPU work from 573.85 T to 96.21 T, a
  5.96x reduction. This is the concrete reason for retaining one bounded retry.
- Java main/test compilation passes after adding whole-case preflight,
  prediction-based root scoring, exact-policy retention, DP-sweep accounting,
  and the safety factor.
- The current targeted run passes 19 tests with a 1 GiB test heap: four
  whole-case decision cases, eight admission-accounting cases, six
  decomposition-strategy cases, and one memory-release case. They cover
  initial admission, largest-contributor final retry, rejection before formal
  materialization, exact-plan retention/reproduction ceiling, compute+OOC
  addition, fail-closed calibration, whole-case rather than per-state SLA,
  contributor ordering, exact RC identity, DP-sweep accounting,
  safety-factor inflation, and the adaptive initialization regression.

### Work still required before release

1. Resolve the 3bua admission blocker: the closest profile is RTX Pro 6000
   at `474.7665 h`, so it needs about a 29% whole-case reduction to fit the
   336-hour SLA; A5000 and Titan V need roughly 2.8--2.9x reductions.
2. Collect the formal 3k3q results from jobs `12048780` and `12048781`, verify
   their artifact hashes and GPU/OOC diagnostics, and record wall time.
3. Do not submit 3bua as an accepted production case yet. After resolving its
   admission gap, update the production scripts and submit only the approved
   replacement jobs.

## Historical: 2026-07-10 validation and production checkpoint

This section records the state before the lifecycle/root-budget/admission work.
Its pre-restart task list is superseded by the 2026-07-11 checkpoint above.

Completed evidence:

- Two-A5000 full GPU regression job `12025816`: **40 tests passed**.
- Production-shape allocation-free root/OOC dry run job `12027758`: **passed**
  for 3bua and 3k3q with the current OOC planner.
- Focused two-A5000 OOC job `12027760`: **passed** in 2m08s, including
  extreme/NaN/all-`-Infinity` stable merges, bilateral row plus lambda tiling,
  runtime allocation auditing, production-shape planner checks, and multi-GPU
  partitioning.
- Runtime diagnostics now report bounded child workspace, independent bounded
  output workspace, output blocks, gathers, launches, moved bytes, and timing.
  The output workspace is no longer incorrectly capped by the child M-state
  chunk.

Production defaults used by the six scripts:

- `branchdp.dp.gpu.childSliceMaxBytes=2147483648` (2 GiB combined packed child
  lower+upper workspace). This is a conservative operational cap added with the
  slicer; it is not a correctness requirement and is not the per-buffer JVM
  direct-buffer limit.
- `branchdp.dp.gpu.outOfCore.outputWorkspaceMaxBytes=4294967296` (4 GiB bounded
  output partial workspace).
- `rootSplit=gpubytes`, which still evaluates every split edge of the fixed
  branch decomposition using integer/shape arithmetic before materializing the
  winning root. The host-budget addition preserves this cheap estimate: it does
  not allocate candidate DP tables.

Current pre-restart work, in order:

1. Compile and regression-test explicit PACK* large-memory release.
2. Regression-test `branchdp.rootSplit.hostBudgetBytes` and ensure the selected
   root fits both host heap and GPU execution budgets.
3. Set explicit production host budgets and disable the optional static DP
   cache in all six scripts; repeat the production-shape dry run.
4. Create one immutable class/resource artifact and record its hashes.
5. Cancel the surviving jobs from `12027705`--`12027710`, resubmit all six, and
   verify their artifact hash, root feasibility, OOC plan, and initial progress.
6. Replace this pre-fix checkpoint with final JobIDs and measured evidence.

## Historical: 2026-07-09 implementation checkpoint

This section is the preserved 2026-07-09 handoff. Its “pending” statements are
historical and are superseded by the 2026-07-10 checkpoint above.

The following changes are in the working tree and must be preserved if work is
continued in another session:

- `chooseHybridPlan(req, budget, gpuCount)` no longer reads or depends on the
  global forced-slicing execution policy. The policy gate is in `compute()`.
  This fixes the test-order/thread/alias-scope dependence where a GPU test could
  cause a later CPU planner call to return `null`.
- The legacy row-only slicer now checks its minimum representable allocation
  before loading a CUDA module or allocating any fixed device buffer. A shape
  requiring lambda tiling therefore gets the same diagnostic on A5000 and Titan
  V instead of sometimes reaching a raw `CUDA_ERROR_OUT_OF_MEMORY`.
- Hybrid streaming now has two pinned host/device slots; next-block gather and
  H2D are submitted while the current block kernel runs. All GPU workers share
  one bounded gather pool capped at 32 threads.
- Per-output-tile `mIdxList` host/device storage is grow-once/reused in both the
  hybrid and flattened-slice paths. `CUBuffer` supports prefix-sized async
  upload/download, and `DPTable.readPair()` has a `DoubleBuffer` variant for
  direct pinned gather.
- Isolated CPU tests passed with a 512 MiB test heap: five new
  planner/preflight/gather tests plus all eleven `TestDPTable` tests (16 tests
  total). The earlier aggregate test exit 137 happened before a new XML report
  was written; it was not one of these assertions failing.
- A new in-progress `DPGpuOutOfCore` implementation has been started. Its pure
  planner assigns an independent rectangular extent to every union-M and
  lambda dimension and estimates each child's projected Cartesian tile exactly.
  Host-side M-box, lambda-box, and packed-child enumerators/gather are present.
  The Java source compiles and a forced 64-byte child-workspace CPU test proves
  that both axes are tiled, every M/lambda state is covered exactly once, every
  packed value matches the original child table, and one byte below the
  planner's minimum budget is rejected.
- `full_dp_n_children_out_of_core` now accumulates lower/upper stable
  `(max,sumExp)` partials across lambda tiles. Java allocates each device
  workspace once at the planner-declared maximum, reuses it for all tiles, and
  partitions disjoint output regions across one or multiple GPUs. The root
  preflight includes the OOC minimum and propagates `branchdp.dp.gpu.maxBytes`
  into actual planning rather than using it only as a loose admission check.
- Java compilation and the isolated CPU OOC test pass. CUDA source compilation
  and one/two-A5000 numerical tests are being run by SLURM job `12024852`.
  Job `12024851` never compiled the kernel because it selected an incomplete
  conda CUDA toolkit without headers; the replacement job validates both nvcc
  and `cuda_runtime.h` before building.

The next steps are: fix any CUDA compile/numeric issues found by `12024852`, add
extreme-value/NaN/all-`-Infinity` partial-merge cases and allocation counters,
then run the complete GPU/production validation matrix below.

## Measured 3bua bottleneck

- Usable A5000 memory per device: about 21.5 GiB.
- Child tables consumed by the problematic edge: 345.7 GiB.
- Current flattened-union slice width: 20 union states.
- Current packed child data per slice: about 575 MB.
- Number of slices: 129,024.
- Aggregate child data gathered by the current schedule: about 67.45 TiB.
- Unique child data: 345.7 GiB.
- Re-read amplification: about 200x.
- Parallel gather baseline: about 1.05 seconds per slice after the existing
  `LongStream.parallel()` change.

The 1.05-second result is a limit of the current flattened-union schedule, not
a physical lower bound for the contraction. The schedule destroys reuse across
the Cartesian product of child row dimensions.

## Completed changes

### GPU-aware root scoring

`rootSplit=gpubytes` still exhaustively evaluates every possible split edge of
the already-computed unrooted branch decomposition. It now records separate
quantities instead of treating the parent output table as resident GPU data:

- `hostTableBytes`: final DP storage in CPU memory.
- `fullDeviceBytes`: device buffers from
  `DPGpuFullDP.estimateDeviceBytesForShape()`, including only the bounded parent
  output tile.
- `fitsSingleGpu`: whether every lambda edge fits the configured or detected
  single-device budget.
- `estimatedSlicedTraffic`: the lower of the feasible resident/streamed plan
  and a conservative projection-aware upper bound for the flattened slicer.
- `logDPWork`, plus an exact `BigInteger` sum of `MStates * LambdaStates` for
  comparison without floating-point tie errors.

Leaf edges now use the same footprint estimator. A leaf that writes a 327 GB
host table can therefore correctly report only the roughly 16 MB output tile
and its actual energy/metadata buffers as device memory.

The `gpubytes` ordering is:

1. Valid TESS candidates before invalid candidates.
2. Structurally GPU-supported candidates before unsupported candidates.
3. Fully resident candidates before candidates that require child slicing.
4. For sliced candidates, lower estimated sliced traffic.
5. Lower exact GPU DP work, then lower full-device bytes and existing tie-breaks.

This finds the exact optimum for that declared ordering among all root edges of
the fixed branch decomposition. It is not a global optimum over all possible
branch decompositions, and it is not a proof of minimum wall-clock time.

The hybrid-aware rerun changed the interpretation of the original two-child
shape. Under the old flattened-slice model it was scored as tens of TiB and the
scorer preferred split 81. Under the execution-consistent resident/streamed
model it is about **371.2 GB of single-GPU-equivalent input traffic**, so the
lower-work split **237** remains viable. Root scoring and the executor must use
the same plan model; “split 81 is mandatory” is an obsolete result from the old
flattened-slice cost model.

For an offline dry-run without visible GPUs, set an explicit budget, for
example:

```text
-Dbranchdp.rootSplit=gpubytes
-Dbranchdp.rootSplit.gpuBudgetBytes=23085449216
```

### Child gather indexing

- Added paired `DPTable.readPair()` access so a sharded table computes and
  checks the shard route once for lower and upper values.
- Power-of-two shard sizes use shift/mask routing.
- Precomputed each child's lambda-key to original-table offset once per GPU
  request.
- Computed each selected child row's original base once per slice.
- Gather tasks operate in 16K-state chunks, avoiding quotient/remainder and
  mixed-radix division in the per-element loop while retaining CPU parallelism.

### Validation completed

- Isolated `compileJava`: passed.
- `TestDPTable`, including paired reads across shard boundaries: passed.
- CPU-only child-slice packed values versus the original mixed-radix index
  formula for two different child projections: passed.
- Production-equivalent root scoring has been rerun with the hybrid-aware cost
  model. The original split 237 is no longer falsely charged flattened-slice
  traffic in the tens-of-TiB range; its single-GPU-equivalent input traffic is
  about 371.2 GB. A final dry-run record with the current code and exact command
  line still needs to be captured in this document.

### Reusable slice buffers

Each GPU worker now owns grow-once lower/upper host arrays and matching raw
device buffers. Subsequent slices reuse them instead of allocating hundreds of
MB of Java arrays and calling `cuMemAlloc`/`cuMemFree` for every slice. The
isolated build and packed-value regression tests pass.

The hybrid path additionally alternates two pinned host/device streamed-child
slots. Gather and H2D for the next block are issued through a separate CUDA
stream while the current block's kernel/output work uses the compute stream.
This code has CPU-side regressions but its new overlap/concurrency behavior has
not yet been revalidated on A5000.

### Resident plus streamed tiling

When a multi-child edge is too large for full residency, the planner evaluates
each child as the streamed child. It keeps all other child tables in canonical
device layout, chooses a streamed row-block size from the real device budget,
and minimizes aggregate streamed bytes plus per-GPU resident replication.

The streamed child's row range is partitioned across GPUs. Each row block is
gathered and uploaded once, then contracted against every state of the resident
children and all parent-only M dimensions before advancing. For the original
3bua two-child shape, the planner chooses the 327 GiB child as streamed and the
18.7 GiB child as resident. Estimated aggregate movement on eight GPUs is about
477 GiB instead of 67.45 TiB, a roughly 145x reduction.

The new `full_dp_n_children_hybrid` CUDA entry point reads resident children by
their canonical strides and the streamed child by packed row/lambda strides in
the same launch. Existing full and flattened-slice kernels remain separate
fallbacks.

Validation on an RTX A5000 (SLURM job 12021477):

- Hybrid planner/gather/output enumeration CPU tests: passed.
- Allocation-free 3bua production-shape planner test: passed.
- Hybrid CUDA output versus the Java reference for two differently projected
  children: passed with zero measured lower/upper difference.
- `TestDPTable`: passed.

Job 12021477 predates the current double-buffer implementation. It proves the
resident/streamed kernel's numerical mapping, but it is not sufficient evidence
for the new two-stream pipeline.

## Known 3k3q failure class

3k3q exposed the exact boundary that row-only slicing cannot cross. For at
least one edge, a complete child row's projected lambda data is larger than
usable VRAM even after the child-row/union slice is reduced to its minimum.

- On an A5000/fennario run, the preflight reported approximately
  `need=30.3 GiB, usable=21.5 GiB` and exited cleanly.
- On a 12 GB Titan V run, an older/allocation-order-dependent path reached
  `cudaMalloc` and failed with `CUDA_ERROR_OUT_OF_MEMORY`.
- The new early row-only preflight makes the diagnostic deterministic, but it
  does **not** make the edge computable. Lambda tiling plus stable partial merge
  is the actual fix.

3bua and 3k3q are different stress cases. 3bua primarily has too many child
rows/Cartesian combinations and benefits from resident + streamed reuse. 3k3q
has a row-internal lambda footprint that exceeds VRAM, so no choice of row slice
width or resident child can solve it.

## Historical remaining critical work (2026-07-09)

The items below are preserved to explain the OOC implementation history. The
current release checklist is the 2026-07-11 production-admission list above.

### 1. Finish and validate the double-buffer pipeline

- Recompile after the current workspace changes.
- Rerun hybrid CUDA output versus Java on one A5000.
- Rerun the streamed-row partition on two A5000 GPUs.
- Confirm that separate upload/compute streams do not reuse or unpin a slot
  before its transfer/kernel consumer completes.
- Measure gather, H2D, kernel, and download times separately and confirm actual
  overlap rather than relying only on enqueue order.

### 2. Validate and harden bilateral/multi-child row tiling

- When no useful child table can remain resident, tile each relevant union-M
  dimension independently instead of taking one interval in flattened union
  order.
- For every child, derive the row-key count from only the M dimensions in that
  child's projection. Gather each distinct row once per rectangular tile.
- Support two children and N children with overlapping or disjoint M
  projections; do not assume the two-child 3bua topology.
- Partition output tiles without duplicate or missing parent M states on one
  and multiple GPUs.
- Keep every host and device workspace allocation at or below the planner's
  declared maximum.

### 3. Validate and harden lambda tiling and stable partial reduction

- Give each parent lambda dimension an independent tile extent. For each child,
  gather only the unique projected lambda keys required by that tile; do not
  allocate a full `childLambdaOriginalOffsets` map for a billion-state domain.
- Stream the combinatorial `lambdaOnlyRigid`/`lambdaOnlyMin` arrays by lambda
  tile rather than treating them as fixed resident buffers.
- Have the CUDA kernel return/accumulate `(max,sumExp)` for lower and upper
  bounds. Merge tiles using max-rescaling and convert to `max + log(sumExp)`
  only after the final lambda tile.
- Preserve NaN and all-`-Infinity` behavior of the existing Java/full CUDA
  reference.
- If the remaining non-combinatorial metadata itself exceeds the usable
  budget, report that structural fixed-buffer limit before allocation; do not
  mislabel it as a child-row failure.

### 4. Production-scale validation

- Benchmark the hybrid-aware root result (currently split 237) on the production
  3bua configuration.
- Force the original two-child shape through the hybrid path and measure actual
  host reads, H2D bytes, kernel time, and end-to-end edge time.
- Run the real 3k3q failing edge far enough to show that it selects lambda
  tiling, stays below the device budget, and advances past the former 30.3 GiB
  preflight/OOM point.
- Check NUMA placement and parallel first-touch for production sharded child
  tables.

## Historical OOC completion criteria (2026-07-09)

The issue is considered resolved only when all of the following hold:

1. Existing full-resident and forced-slice numeric GPU tests match the Java
   reference. (Passed previously; rerun before production deployment.)
2. The mixed resident/streamed path has a CPU-only mapping test and a CUDA
   numeric comparison test. (Passed on A5000, job 12021477.)
3. The 3bua shape no longer performs flattened-union 129,024-slice gathering.
4. Measured aggregate child traffic is within a small factor of the planned
   roughly 477 GiB rather than tens of TiB.
5. A production-equivalent 3bua run shows the problematic edge near the GPU
   kernel-time regime rather than the current multi-hour gather regime.
6. A forced bilateral-row test uses at least two independently tiled child M
   projections, covers every output exactly once, and matches the Java reference
   on CPU mapping and A5000 numerical checks. (CPU coverage/packing passed;
   A5000 check pending in job `12024852`.)
7. A forced lambda-tile test makes one complete child row larger than the
   artificial test budget, uses multiple lambda tiles, and matches the Java
   result including extreme-value, NaN, and all-`-Infinity` cases.
8. Planner tests assert `estimatedDeviceBytes <= budget` for resident/streamed,
   bilateral-row, and lambda-tiled plans, and runtime allocation counters show
   no device allocation larger than or outside the declared plan. (OOC planner
   upper-bound and exact-minimum boundary assertions pass; runtime counters are
   still pending.)
9. Single- and two-GPU tests pass for every forced execution tier, with no
   missing/duplicate output indices at partition boundaries.
10. A production 3k3q run advances past the former single-row OOM boundary on
    both an A5000-class budget and an artificial 12 GB budget (real Titan V is
    preferred when available), without falling back to unbounded Java DP.

Passing only the existing 3bua hybrid tests is not final acceptance. The final
deliverable is the complete arbitrary-shape bounded out-of-core path described
by criteria 6--10.
