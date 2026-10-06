# MARK* DCC assignment, updated 2026-10-06

Assign 91 design workloads to DCC and retain 48 locally. The original split
was 82/57; the October 6 update transfers nine unstarted local designs to DCC.
This partitions the same 139 normally completed historical MARK* workloads.
It does not include timeout or failed workloads.

## Files

- `dcc_cases.txt`: the original 82 design IDs, followed by nine additions.
- `dcc_designs.tsv`: exact residue selections, expected sequence counts,
  prepared PDB paths and checksums, and historical resources and timings.
  Original rows and their order are preserved. `dcc_lane` retains the legacy
  four-lane grouping; additions use `NA`. It is not the current worker limit.
- `local_cases.txt` and `local_designs.tsv`: the complementary 48 designs.
- `dcc_additions_20261006_cases.txt` and `dcc_additions_20261006_designs.tsv`:
  only the nine new designs, for submission without repeating the original 82.
- `rebalance_20261006.json`: transferred design IDs, local task IDs, and the
  handoff status recorded when this update was prepared.
- `summary.json`: updated counts, historical workload totals, and concurrency.

## October 6 handoff

The additions are `3ma2_flex_p0`, `4wyu_flex_p0`, `5dc0_flex_p7`,
`2q2a_flex_p3`, `2rfd_flex_p7`, `2rfe_flex_p2`, `5dc0_flex_p6`,
`2rfd_flex_p6`, and `5dc0_flex_p5`. All were pending locally; no running or
completed workload is transferred. There are no DCC-to-local transfers.

The user-reported snapshot for DCC array `57713739` has 24 completed, six
running, and 52 unstarted designs. The additions bring unstarted DCC work to
61 designs. Submit only the additions, preserving the running array's input
mapping and completed results. Keep the total concurrency across the original
array and additions at six; two independent arrays capped at six each would
exceed that limit. Waiting for the original array to finish before starting
the additions is also valid.

The user confirmed DCC submission of all nine additions on October 6, with
no newly completed designs. The nine held local duplicate tasks were then
cancelled after verifying that each was still pending and held. Running
workloads were preserved. External job IDs have not been supplied. Keep
newly transferred tiers of each system on one CPU model and record actual
CPU topology and binding.

The active local array `12812439` reads its original frozen manifest outside
this repository. Do not replace that manifest or renumber its task indices.
This repository update does not submit DCC jobs or restart existing tasks.

## Requested DCC run

Run MARK* with 64 CPU threads per task and at most six simultaneous tasks.
Configure both the scheduler allocation and the program's CPU setting to 64.
The `markstar_cpus` column records the historical allocation of 16; it is not
the requested CPU count for this rerun. Memory and heap columns likewise
describe the historical runs. Select sufficient memory for each design and
leave space outside the Java heap.

Preserve the original prepared structures, mutable and flexible residue
selections, sequence workload, epsilon 0.683, CPU CCD, FP64, disabled stability
filter, and benchmark implementation. Construct fresh energy matrices for
each design. Record full process wall time, per-sequence results, actual CPU
allocation and binding, hardware model, and software version.

The input structures and frozen executable are not included here. Obtain
them from the original prepared packages before running; historical absolute
paths are provenance and must not be assumed to exist on DCC. Verify prepared
PDB files against `pdb_sha256`. This manifest is not a standalone executable
bundle. No DCC jobs have been submitted by creating this assignment.

## Scheduling basis

The original 82-design assignment represented 677.4913 historical task hours
at 16 CPUs and was sized for four workers at historical performance. That
seven-day estimate is superseded by six concurrent 64-CPU tasks and the
October 6 handoff. Historical resource columns remain unchanged as provenance.
Do not infer physical core counts or fourfold speedup from the CPU allocation.

The updated assignment uses observed rerun timings to relieve the local queue.
Runtime estimates depend on hardware, current progress, and queue availability;
they are not measured benchmarks for the transferred designs.
