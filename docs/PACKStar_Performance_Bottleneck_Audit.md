# PACK* Performance Bottleneck Audit

Last updated: 2026-07-11 16:20 EDT

## Scope

This document records a performance audit of the formal 3k3q PACK*/K* run
`12048780` and turns the observations into hardware-independent optimization
work. It is deliberately separate from
`PACKStar_GPU_Child_Slicing_Status.md`, which tracks correctness, release
gates, and production-job status.

The evidence here comes from:

- formal job `12048780`, four RTX Pro 6000 GPUs, frozen artifact
  `packstar_gpu_ooc_prod_20260711_3k3q_pf`;
- log `/usr/xtmp/lz280/slurm_logs/3k3q_pro6000_full_12048780.out`;
- the runtime diagnostics emitted by normal multi-GPU DP and bounded OOC DP;
- the corresponding planner/executor source in `DPGpuFullDP`,
  `DPGpuOutOfCore`, `BranchDpAdmission`, and `BranchDpBackend`.

This is a checkpoint analysis, not a completed-run profile. At the timestamp
above, 21 of 39 K* sequences had complete Protein+Complex results, the 22nd
Protein state had completed, and its Complex state was in progress.

## Main conclusion

**OOC is not the largest whole-run bottleneck at this checkpoint.** It is the
dominant local cost in the currently running OOC-heavy Complex state, and its
implementation has clear general optimization opportunities, but the largest
global cost is high-lambda Branch-DP work, executed across the DP sweeps
required by the PACK* estimator. The second substantial cost is the estimator work outside Branch-DP,
including CCD/energy evaluation, sampling, and model/refinement control.

The optimization order should therefore be:

1. make decomposition, root selection, admission, and runtime use a complete
   per-edge time model;
2. reduce high-lambda DP work through better decomposition/root choices;
3. optimize estimator sampling/CCD and overlap independent state work;
4. optimize OOC partitioning, data reuse, and pipelining;
5. pursue distributed device residency/fusion for the remaining large table
   transfers.

## Whole-run time distribution at the checkpoint

The job began at 14:19:41 EDT. The analyzed log checkpoint was 16:20:58 EDT,
for about 7,277 seconds (121.3 minutes) of wall time.

The following categories are reconstructed from completed estimator timers,
normal multi-GPU DP timers, bounded-OOC multi-GPU timers, and wall-clock
residual. They should be read to about one percentage point, not as
cycle-accurate profiling.

| Category | Attributed time | Share of wall time | Interpretation |
|---|---:|---:|---|
| Normal Branch-DP GPU execution | 4,165.99 s | 57.3% | Full/resident or ordinary tiled multi-GPU edges |
| Bounded OOC Branch-DP execution | 258.12 s | 3.5% | First two OOC edges, both in the current Complex state |
| Estimator work outside Branch-DP | about 1,948.8 s | 26.8% | CCD/energy evaluation, sampling, fitting, bounds, and refinement control in completed estimators |
| Setup and not-yet-separated residual | about 904 s | 12.4% | Reference/emat setup, preflight, decomposition/materialization, initial CPU DP, state transitions, and in-flight work |
| **Total** | **about 7,277 s** | **100%** | Checkpoint wall time |

The reported reference-energy and energy-matrix setup steps took about 108
seconds in total. The remainder of the 904-second residual cannot be split
reliably because the current log has no wall-clock span around preflight,
decomposition, state initialization, and CPU initial DP. Adding those spans is
part of the proposed instrumentation work.

This distribution answers the central question: OOC accounted for only about
3.5% of elapsed wall time so far. All Branch-DP together accounted for about
60.8%, while estimator work outside DP accounted for about 26.8%.

## Branch-DP distribution

Lambda cardinality is only a proxy for kernel cost, but it separates the
observed normal-DP workload well:

| Lambda states per edge | Normal-DP time | Share of normal DP | Share of checkpoint wall time |
|---:|---:|---:|---:|
| `<= 16` | 386.46 s | 9.3% | 5.3% |
| `17..255` | 144.04 s | 3.5% | 2.0% |
| `256..2047` | 428.48 s | 10.3% | 5.9% |
| `>= 2048` | 3,207.01 s | 77.0% | 44.1% |

The largest repeated shapes were:

| Shape | Aggregate time | Why it matters |
|---|---:|---|
| `m=2,776,032,000`, `lambda=2688`, no child | 1,591.49 s over four sweeps | One proposal collapse/drift caused three runtime refinement sweeps after the initial sweep |
| `m=2,204,496,000`, `lambda=2688`, no child | 631.97 s over two sweeps | Compute-heavy lambda enumeration |
| `m=152,409,600`, `lambda=25920`, one child | 422.23 s over two sweeps | High-lambda child fold |
| `m=4,996,857,600`, `lambda=5`, one child | 85.18 s over two sweeps | Low arithmetic intensity; roughly 80 GB of output per sweep makes it table/output-bandwidth dominated |

The first shape alone consumed about 21.9% of checkpoint wall time in GPU DP.
Its estimator emitted refinement checks with proposal collapse, then drift, and
finally stability in rounds 1, 2, and 3. These repeated sweeps are required by
the current PACK* algorithm and are not classified here as implementation
waste. They do demonstrate that sweep count is a first-order input to the
planner and admission model, rather than merely a conservative constant.

The low-lambda 4.997-billion-state edge is also important. Its nominal DP work
is only about 25 billion work units, but one sweep took about 42.7 seconds and
produced about 79.95 GB of lower/upper output. A cost model based only on
`gpuWork / gpuRate` substantially underprices this edge. Normal DP output
traffic and host-table writeback must be modeled even when no OOC child gather
is required.

## Estimator work outside Branch-DP

For the 44 completed partition-function estimators at the checkpoint:

- completed estimator timers totaled 4,389.53 seconds;
- about 2,440.70 seconds of that was normal Branch-DP execution;
- the remaining approximately 1,948.8 seconds was outside Branch-DP;
- completed estimators generally used roughly 670--840 CCD calls each.

This non-DP envelope is about 32% of the attributable completed-state work and
about 27% of checkpoint wall time. It is especially visible on Protein and
small Complex states, where DP is cheap but sampling still takes roughly
10--70 seconds. OOC work cannot improve this part of the run.

Required follow-up profiling should separate at least:

- conformation generation;
- CCD/minimization and energy evaluation;
- pilot sampling;
- proposal/model fitting;
- Bernstein/bound calculation;
- refinement-decision overhead;
- waits for CPU and GPU energy-calculation pools.

## What the first production OOC state shows

The current `C440=ARG C444=tyr` Complex state is the first formal state in this
run to exercise the bounded OOC executor. Its initial DP sweep contained about
50.5 seconds of normal GPU edges and 258.1 seconds of OOC edges, so OOC was
about 84% of that state's initial DP sweep. This is a **local** bottleneck even
though it was not the largest **whole-run** category at the checkpoint.

### Large OOC edge

Shape: `mStates=8,291,082,240`, `lambdaStates=270`, one child.

- multi-GPU elapsed: 222.99 s;
- each GPU gathered about 79.95 GB of child data;
- slowest per-GPU elapsed: 217.08 s;
- slowest per-GPU components: about 25.9--26.9 s gather, 4.3 s enqueue,
  155.0--158.8 s GPU/transfer wait, and 19.7--20.4 s host output copy.

This edge is primarily kernel/transfer-wait limited, with meaningful gather
and output-copy costs. The executor currently serializes gather, upload,
kernel, download, and host copy at every tile.

### Gather-bound OOC edge

Shape: `mStates=205,632`, `lambdaStates=241,920`, one child.

- multi-GPU elapsed: 35.12 s;
- every GPU gathered the same approximately 132.66 GB child table;
- per-GPU gather took about 19 seconds;
- enqueue took about 8 seconds;
- GPU wait took only about 4--6 seconds;
- output copy was negligible.

This edge is host gather/layout/upload limited. Kernel tuning alone cannot
materially fix it.

### Multi-GPU replication and prediction drift

`runMultiGpuOutOfCore` currently selects `free-M` splitting whenever the free-M
state count is at least the GPU count. Each GPU then traverses all M boxes and
gathers the same child tiles. In the two observed edges, the aggregate child
gather was approximately:

```text
4 * (79.95 GB + 132.66 GB) = 850.4 GB per DP sweep
```

The preflight prediction for this state reported only 132.66 GB of OOC traffic
per sweep. It did not include four-GPU replication, and the first 79.95-GB
edge crossed into OOC at runtime because preflight/root scoring used an 85-GiB
per-GPU budget while runtime reported an approximately 73.1-GiB usable budget.

This is not merely a calibration issue. The general planner must evaluate the
actual multi-GPU execution plan, replication factor, output traffic, and
runtime memory budget. Admission and root selection should consume that same
plan rather than a single-GPU traffic proxy.

## General optimization plan

### P0: complete and unify the time model

Root selection and admission currently model calibrated GPU work plus modeled
OOC bytes. The model should be expanded to an allocation-free per-edge plan
that includes:

- lambda/M arithmetic work and shape-dependent throughput;
- output lower/upper bytes for every edge, including ordinary full DP;
- child H2D bytes and replication across GPUs;
- OOC host gather bytes and memory-layout cost;
- D2H and host-table writeback bytes;
- output-block, lambda-tile, and kernel-launch counts;
- the selected multi-GPU split dimension;
- actual runtime-usable memory, with the same headroom and resident GPU-pool
  allocations seen by execution;
- expected and worst-case DP sweep counts.

A small roofline-style model is sufficient: compute, host gather, H2D, D2H,
and host writeback each get a measured rate, and the model combines them
according to whether the executor pipelines them. Root/decomposition scoring,
admission, and runtime must use the same plan object.

This work is P0 because the current incomplete objective can select a root
that minimizes nominal arithmetic while creating a very large output table or
replicated child traffic.

### P1: optimize decomposition for total state time

Weighted Hicks and root selection should jointly minimize predicted end-to-end
state time, not branchwidth, `logTESS`, or GPU work in isolation. Candidate
scoring should use the complete model above across the whole rooted tree.

The bounded exact-improvement policy should also have a performance mode that
is independent of the admission SLA. In this run, `softStateHours=24` means no
state predicted below 24 hours receives exact improvement, even if a one- or
two-minute search could save many minutes of production runtime. A payback
policy should spend bounded search time on the top predicted contributors when
the expected saved execution time exceeds the search budget.

Evidence already recorded in the status document shows that a bounded exact
search reduced one 3bua state's one-sweep work by 5.96x. In this run, one
high-lambda shape consumed about 22% of checkpoint wall time. Long-tail
decomposition work therefore has more global leverage than the OOC executor
alone.

### Algorithmic baseline: refinement sweeps are not implementation waste

Most completed states performed one runtime DP sweep after initialization, but
the largest observed state performed three because its proposal first
collapsed and then drifted. Each refinement repeated an approximately
400-second full-tree DP. This is the defined behavior of the current PACK*
algorithm, so this audit does **not** treat the extra numerical sweeps as an
implementation optimization target.

The performance model must still price every required sweep correctly.
Implementation work may reuse immutable decomposition structure, index maps,
OOC boxes, allocated buffers, and table storage between sweeps, while
recomputing every numerical value required by the algorithm. Changing the
number of sweeps, proposal update, or stopping rule would be separate
algorithm research and would require a new correctness/statistical argument.

## Confirmed implementation optimization opportunities

The following opportunities do not rely on changing PACK*'s required sweep
semantics. They are confirmed targets for further analysis and implementation.

1. **Eliminate avoidable multi-GPU child replication.** The current free-M
   split can make every GPU gather the same child table. Choose the split from
   physical traffic and makespan, not from free-M cardinality alone.
2. **Avoid host round trips between producer and consumer DP edges.** Large
   intermediate tables are downloaded to host and immediately gathered and
   uploaded by their parent. Preserve distributed residency, co-partition
   parent/child edges, or fuse compatible edges where bounded memory permits.
3. **Price all table traffic in decomposition/root/admission decisions.** The
   current model misses ordinary DP output bytes, host writeback, replicated
   child upload, and the actual runtime OOC plan, so it can select a nominally
   cheap but bandwidth-heavy root.
4. **Pipeline OOC stages.** Overlap host gather, H2D, kernel execution, D2H,
   and host-table copy with bounded double or triple buffering instead of
   synchronizing every tile serially.
5. **Remove the extra heap-to-pinned child-tile copy.** Fill reusable pinned
   workspaces directly and reuse row/lambda metadata arrays to reduce memory
   bandwidth, allocation, and GC pressure.
6. **Investigate small-state over-solving and resource scheduling.** Many
   states perform roughly 670--840 CCD calls and often finish below the target
   epsilon. Profile the estimator phases, then evaluate safe adaptive error
   allocation, batching, and overlap of independent small/sampling-heavy work.

### P1: cost-based multi-GPU OOC partitioning

Evaluate at least `free-M`, `M-box`, and, where useful, lambda-box partitioning
before execution. Choose the split that minimizes predicted makespan subject
to balance and memory constraints.

For the observed gather-bound edge, splitting M boxes would assign disjoint
child tiles to GPUs instead of making every GPU read the complete 132.66-GB
child table. The planner should also hoist a packed child tile outside the
free-output-block loop and reuse it for all compatible free-M blocks.

The expected benefit is shape dependent:

- compute-heavy OOC edges retain similar kernel work per GPU, but save host
  gather/H2D replication;
- gather-heavy OOC edges can approach a GPU-count reduction in aggregate
  child traffic, subject to load balance and output-block reuse.

### P1: pipeline OOC and reuse pinned workspaces

The OOC executor currently synchronizes after every tile before reusing one
set of buffers. Implement two or three bounded slots so that:

1. CPU workers gather tile `N+1` directly into reusable pinned buffers;
2. an upload stream transfers tile `N+1`;
3. the compute stream executes tile `N`;
4. a completed output slot downloads and copies tile `N-1` to its DP table.

The large observed OOC edge has roughly 46 seconds per GPU in gather plus host
copy outside its approximately 159-second GPU wait. A steady-state pipeline
has plausible 15--25% headroom on that edge, before considering the reduction
from a better partition.

`buildPackedChildBlock` also allocates fresh `lower` and `upper` heap arrays for
every child tile, then copies them into the CUDA buffer's pinned host storage.
Fill reusable pinned storage directly and reuse row/lambda metadata arrays.
This removes large allocation/GC pressure and one full host-memory copy. It is
especially relevant to the gather-bound edge, where gather plus enqueue is
about four to six times the GPU wait.

### P1/P2: optimize estimator sampling and scheduling

Because work outside Branch-DP is about 27% of wall time, it needs its own
optimization track:

- add sequential stopping closer to the actual K* error budget;
- allocate Protein/Ligand/Complex error budgets jointly instead of routinely
  oversolving inexpensive states;
- batch CCD/energy evaluations and measure CPU/GPU pool occupancy;
- warm-start proposal models across closely related mutation states when this
  preserves the estimator's statistical guarantees;
- use a memory-bounded state scheduler to overlap sampling or a small Protein
  state with an independent large DP phase when their resources do not
  conflict.

This scheduler should be heterogeneous: large DP edges may use all GPUs, while
small or sampling-heavy states should not force the entire job into a serial
GPU-idle phase.

### P2: retain or fuse producer/consumer DP tables

The long-term data-movement optimization is to avoid writing a very large
intermediate table to host only to gather and upload it for its parent edge.
Possible implementations are:

- keep disjoint output shards resident across GPUs;
- choose the parent partition to preserve child-shard ownership;
- use peer transfer only for the projections that cross ownership;
- fuse compatible one-child producer/consumer edges;
- spill to host only when the next edge cannot consume the resident layout.

The 4.997-billion-state, lambda-5 edge demonstrates the opportunity: a sweep
spends about 42.7 seconds on a shape whose approximately 80-GB output dominates
its small arithmetic workload. This is a general table-dataflow problem, not
an OOC-only problem.

### P3: kernel specialization and logging

After the planner/dataflow work, profile specialized kernels for:

- very small lambda domains with huge output tables;
- very large lambda domains with no children;
- high-lambda one-child folds;
- small edges that underfill multiple GPUs.

Production progress logging should be sampled per edge or time interval rather
than emitted for every million-state tile. This is unlikely to be the primary
bottleneck, but it removes synchronized I/O noise and makes profiles easier to
interpret.

## Instrumentation required before claiming final speedups

Every state should emit one structured summary with wall-clock spans for:

- interaction graph and decomposition;
- root/OOC planning;
- table materialization;
- initial DP sweep;
- every refinement DP sweep;
- host gather, H2D, kernel, D2H, and host writeback;
- CCD/energy evaluation, fitting, and bound calculation;
- cleanup and mmap I/O;
- GPU count, split dimension, aggregate bytes, and peak host/device memory.

The aggregate byte counters must distinguish logical bytes from physical bytes
replicated across GPUs. Timers should report both the sum of per-device work
and the critical-path wall time; summing four concurrent GPU timers would
otherwise overstate elapsed time.

## Acceptance criteria for optimization work

An optimization should not be accepted from a microbenchmark alone. It should:

1. preserve PACK* bounds, epsilon/confidence behavior, and existing GPU/OOC
   numerical regressions;
2. remain bounded for every supported OOC shape;
3. reduce end-to-end state and whole-case wall time on at least one
   compute-heavy, one output-bandwidth-heavy, one gather-heavy, and one small
   sampling-heavy case;
4. report no regression in aggregate physical traffic, peak memory, or GPU
   load balance unless the wall-time tradeoff is explicitly favorable;
5. reproduce planner predictions within a documented tolerance using the same
   runtime plan and memory budget.

## Relevant implementation locations

- `DPGpuFullDP.runMultiGpuOutOfCore`: fixed free-M-first multi-GPU split.
- `DPGpuFullDP.runOnGpuOutOfCore`: serial gather/upload/kernel/download/copy
  loop and synchronization boundary.
- `DPGpuOutOfCore.buildPackedChildBlock`: per-tile heap allocation and gather.
- `BranchDpAdmission.Hardware` and `Prediction`: current GPU-work plus OOC-rate
  time model.
- `BranchDpBackend`: adaptive decomposition trigger and retained exact policy.
