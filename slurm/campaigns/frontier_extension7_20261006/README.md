# Seven new frontier designs and MARK* extension plan

The user approved seven new designs outside the existing 161-design paper
cohort. Historical timeout designs are separate and are not counted as new.

| System | New tiers | Actual positions | Sequences per design |
|---|---|---|---|
| 5d68 | +3 | 15 | 20 |
| 4u3s | +1, +2 | 14, 15 | 39 |
| 2rl0 | +3, +4 | 16, 17 | 20 |
| 5dc4 | +5, +6 | 16, 17 | 39 |

`requested_designs.tsv` records the approved choices. `new_designs.tsv`
records the exact residue lists, structures, checksums, sequence counts, and
lower-tier anchors. Each new design extends its anchor using the existing
frozen expansion scan. No earlier design definition or active array mapping
was changed.

## Submitted PACK* group

Preparation job `12826503` passed all seven MARK* position/sequence preflights
and all seven preflights using the current optimized PACK* executable.
Sequences match each system's existing anchor exactly. Actual complex
position counts match every requested residue list.

PACK* job **12826506** is submitted to `compsci-gpu` with account `grisman`,
128 CPU threads, four RTX PRO 6000 GPUs, exclusive use of one node, and a
**six-hour limit for the entire group**. It is a single job, not an array.
All seven designs run sequentially within that allocation; the node is kept
between designs. It was pending for resources at the submission audit.

The execution order is `5d68_p3`, `4u3s_p1`, `2rl0_p3`, `5dc4_p5`,
`4u3s_p2`, `2rl0_p4`, `5dc4_p6` (with `_flex_` in the actual IDs).
The driver allows at most two hours per design while retaining 20 minutes
for each later design and five minutes for shutdown. The overall six-hour
limit takes precedence; completion of all seven within this budget is not
guaranteed. Partial results and per-design statuses are retained if a design
times out. There is no automatic long rerun.

The group reuses the validated optimized build from current job `12814879`:
FP64, CPU CCD, pair corrections, residual budget 1, seed 42, epsilon 0.683,
850 GiB Java heap, 800 GiB DP host budget, and 85 GiB per GPU. Energy matrices
are fresh for each design. Existing workloads and their runtime settings
are preserved.

`INCOMPLETE_ESTIMATES` is normal workload completion when Java exits normally,
all configured sequence rows are present, and a wall time is recorded. The
number of returned binding scores is recorded separately. Scheduler success
for the group requires normal workload completion of all seven designs.

- Frozen preparation: `/usr/xtmp/lz280/frontier_extension7_20261006/prep_12826503`
- Run output: `/usr/xtmp/lz280/frontier_extension7_20261006/run_12826506`
- Slurm log: `/usr/xtmp/lz280/slurm_logs/packstar_extension7_12826506.out`
- Small job and validation registry: `registry.json`

The preparation package contains four PDB files (3,105,986 bytes), manifests,
and preflight logs. All experiment outputs and temporary files live under
`/usr/xtmp/lz280`; this directory contains only source and small registries.

## Proposed MARK* assignment: not submitted

`markstar24_proposed.tsv` contains the seven new designs plus all 17 historical
14-day noncompletion designs. The five failed historical designs and the
withdrawn 3k3q higher tiers are excluded. The 139-workload assignment remains
in the separate `markstar_dcc_20261005` directory.

| Destination | Available campaign slots at audit | New tiers | Historical noncompletion | Total added designs |
|---|---:|---:|---:|---:|
| fennario, Xeon Gold 5320 | 4 | 4 | 1 | 5 |
| compsci, EPYC 9554 | 8 | 3 | 8 | 11 |
| DCC, CPU model fixed per system | 6 | 0 | 8 | 8 |

Each proposed MARK* task requests 64 CPU threads and a 14-day limit. CPU
thread counts are not physical core counts. At the time of the audit, the
local slots were already occupied by the current campaign; they represent
capacity as those jobs finish, not immediately free slots. DCC's six slots
are the user-reported total concurrency across all its campaign jobs.

This allocation spreads the added 14-day budgets over the available slots:
420, 462, and 448 task hours per slot respectively. These are budget totals,
not runtime predictions or calendar completion estimates. The historical
timeouts provide lower bounds only. The original backlog must also finish.

Initially allow at most 2/4/3 added long tasks on fennario/compsci/DCC,
respectively, leaving the other half of the slots for the existing backlog.
As that backlog clears, allow the long workloads to use the full 4/8/6 slots.
Do not create a second independent six-worker DCC array. Actual availability
can change as other users' allocations change.

The 5/11/8 proposal requires one additional handoff: **4wwi_flex_p1**,
local job `12812439_24`, was pending and had not started at the audit.
Transfer it together with the proposed 4wwi +2/+3 tasks so all three tiers
can use one DCC CPU model. **This extra handoff has not been applied:**
the pending local task is unchanged and the 139-workload manifests still
contain 91 DCC and 48 local designs. Before implementing a handoff, recheck
that +1 has not started, hold its local copy, and retire it only after external
submission is confirmed. If +1 starts first, keep all 4wwi tiers on compsci;
the added-design split becomes 5/13/6.

All other proposed new tiers follow the CPU model of their local same-system
anchors. DCC must select one actual CPU model per newly launched system and
use a partition that permits 14-day tasks. Its historical per-case CPU models
and new job IDs are not available locally. Existing completed measurements
on other hardware remain separate measurements; the plan does not change them.

Obtain the seven new PACK* full workload timings before launching their
MARK* counterparts. Keep the original 48-CPU jerry jobs `12825151` and
`12825152` separate: they use the original filtered setup and do not replace
the frontier timeout reruns.

## Confirmed earlier handoff

The user confirmed that DCC submitted the nine previously transferred designs,
with no new completions (24 completed, six running). Audit job `12826543`
verified that all nine local copies remained pending and held, cancelled only
those copies, and verified that the running jobs remained unchanged.
The acceptance record is in `../markstar_dcc_20261005/rebalance_20261006.json`.

`plan_markstar.py` records that one-time acceptance/cleanup operation as well
as the scheduling analysis. Its preconditions deliberately prevent repeating
the cleanup after the held jobs have been retired. It does not launch any
MARK* benchmark. Detailed scheduler snapshots and the planning audit are in
`/usr/xtmp/lz280/frontier_extension7_20261006/markstar_plan_12826543`.
