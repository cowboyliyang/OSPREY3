# MARK* DCC assignment, 2026-10-05

Assign the 82 shorter completed design workloads to DCC and retain the other
57 locally. This partitions the 139 normally completed MARK* workloads in
the current runtime cohort. It does not include timeout or failed workloads.

## Files

- `dcc_cases.txt`: 82 design IDs in ascending historical runtime order.
- `dcc_designs.tsv`: exact residue selections, expected sequence counts,
  prepared PDB paths and checksums, historical resources and timings, and
  four balanced groups in `dcc_lane`.
- `local_cases.txt` and `local_designs.tsv`: the complementary 57 designs.
- `summary.json`: historical runtime totals and scheduling estimate.

## Requested DCC run

Run MARK* with 64 CPU threads per task and at most four simultaneous tasks.
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

## Scheduling estimate

The selected designs consumed 677.4913 task hours at 16 CPUs. Scheduling the
longer jobs first within the selected set gives four groups of approximately
169.4 hours each: about 7 days 1 hour 23 minutes with four continuously
available workers at historical performance. The longest selected design
took 24.61 hours. No fourfold speedup is assumed for 64 CPUs.

The remaining 57 local designs consumed 5847.8700 task hours. These estimates
exclude queueing and reruns; CPU and machine changes require measurement.
