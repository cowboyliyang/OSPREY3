# Deadline focus: active MARK* execution policy

This policy supersedes the pending launch lists in `markstar_dcc_20261005`
and `frontier_extension7_20261006`. The user authorized changing all unstarted
work, including previously submitted DCC tasks. Preserve completed results,
all started work, and every task that starts before an adjustment reaches it.
The November 4 submission favors large measured speedups, with at least one
completed current-round 64-CPU MARK* case below 14 days per system. Near-14-day
completion is preferred but is not guaranteed by these tier choices.

**Latest instruction: keep useful CPU work ready without waiting for GPUs.**
MARK* may start as soon as inputs pass preflight and a CPU slot is available.
PACK* measurement, completion status, and the 30-minute screening threshold
are no longer launch prerequisites. Existing running jobs remain untouched.
High-potential continuation tiers supplement the primary targets; the low-value
intermediate backlog stays deferred. Separate Slurm allocations cannot
guarantee uninterrupted ownership of a node between tasks.

## Latest append-only additions

The subsequent [addition registry](../near14_20261006/README.md) preserves this
phase's selected smaller queue and appends nine designs. All nine new PACK*
measurements use GPU; 2p4a uses conservative p12 instead of the withdrawn p18.
Four new local MARK* jobs follow existing lane tails. Five DCC additions extend
the combined launch list to nineteen future entries and the new-definition
handoff to nine entries; remote submission has not been performed here.
The updated TSV files are the active lists. The phase-specific counts below
describe the earlier fourteen-entry plan and its four new definitions.

The user also confirmed `2xxm_flex_p1` completed in 6338 seconds with 20 rows,
19 fully Estimated. The external snapshot now records 30 completed designs.
This supersedes the unconfirmed p1 status described in the earlier audit below.

## DCC: adjust the existing queue

**The DCC changes have not been applied by this repository update.** Read the
current queue and task mapping before executing the handoff. The report used
here contains 29 completed designs, three running designs (`2rl0_flex_m1`,
`5em2_flex_p4`, `5em2_flex_p5`), and the previously started `2xxm_flex_p1` whose
current status is unconfirmed. Preserve that last task until its status is
known. Preserve every additional task found to have started since this report.

Use these exact files, preserving the existing array's frozen row order:

- `dcc_cancel_if_pending_cases.txt`: defer these designs only if the matching
  task is still PENDING and has never started. Retain definitions and results.
- `dcc_keep_existing_pending_cases.txt`: retain these ten existing tasks;
  do not submit duplicate copies.
- `dcc_protected_cases.txt`: completed, reported running, or previous-start
  records that are protected regardless of this plan's priorities.
- `dcc_new_designs.tsv` and `dcc_new_cases.txt`: four additional definitions,
  not previously submitted in this campaign. These are not another full cohort.
- `dcc_actions.tsv`: complete 104-row DCC policy, including the four additions.
- `dcc_launch_priority.tsv`: the current long-task-first order for the 14
  retained or new future designs. It does not change an existing array's indices.

The six retained coverage/primary entries from the initial focus plan are
listed below for inventory only. Use `dcc_launch_priority.tsv` for launch order:

1. `2xgy_flex_m1`
2. `2p4a_flex_p9`
3. `4wwi_flex_p1`
4. `4wyu_flex_p0`
5. `3ma2_flex_p0`
6. `2xgy_flex_p0`

Also retain four continuation tasks: `4wyu_flex_p1`, `3bua_flex_m2`,
`3ma2_flex_p1`, and `4wwi_flex_p2`. These are long-task candidates and can start
without waiting for a new coverage measurement. The updated cancellation
list contains 57 candidate designs, not the previous 61. These four cases
were already submitted remotely: preserve their pending copies rather than
resubmitting. If an earlier policy was already applied, check accounting
before restoring a cancelled copy.

The four additions are `3gxu_flex_p12`, `2q1e_flex_p6`, `5dc0_flex_p10`, and
`3bua_flex_m4`. The first three are long speedup targets. All four have passed
MARK* and PACK* input preflight and may be launched without waiting for GPU
timings. Retain sequence statuses for subsequent comparison. `3bua_flex_m4` is a lower-tier
coverage candidate: retain it behind the long-task queue and review early progress; it has no
guaranteed completion forecast. No new design is a measured speedup yet.

Integrate retained tasks and additions into the **existing combined six-worker
limit**, at **64 CPU threads per task**, with a maximum of 14 days. Use
`--account=grisman`, sufficient memory for the recorded Java heap, and 192 GiB
allocation for these new tasks. **The eight long candidates take the first
available slots; coverage fallbacks follow them.** Count all already-running
tasks toward the six-worker limit. If the last reported three running tasks
are still running, start the first three new long targets to bring the total
to six. As those existing short tasks finish, use `2xgy_flex_p0`,
`4wyu_flex_p1`, and `3bua_flex_m2` next. Keep `4wwi_flex_p2` and
`3ma2_flex_p1` queued as further long replacements. None waits for GPU results.
Use this priority order in the existing worker dispatcher or equivalent
pending-only scheduling; never reorder a running array's input manifest.
The lower tiers are retained to protect per-system sub-14-day coverage and
can be skipped only once that system has a qualifying completed result.
Do not start a second independent six-job array.
A full 14-day observation must begin by October 16 for an October 30
data freeze; reconsider unstarted long tasks at that cutoff.

Use one actual CPU model for each system. For new `3gxu`, `2q1e`, and `5dc0`
tiers, inspect the current-round anchor hardware (`3gxu_p9`, `2q1e_p3`, and
`5dc0_p4`) and bind the new tier to that model. For `2xgy` and `3bua`, use any
already-started same-system task as the hardware anchor; otherwise choose one
model and keep it for all retained tiers. Preserve actual hardware provenance
for existing heterogeneous results. Do not guess a CPU model from the partition.

The new TSV includes exact residue lists, prepared PDB checksums, sequence
counts, historical command anchors and heaps. Preserve epsilon 0.683, FP64,
CPU CCD, disabled stability filtering, sequence order, and the existing MARK*
implementation. Build fresh energy matrices. Duke paths are provenance; use
the corresponding already-transferred structures on DCC after checking hashes.
The new inputs use the same prepared PDBs as their anchors, so a new structure
download is unnecessary. This is a manifest handoff, not a portable executable.

`dcc_apply_pending.py` accepts a reviewed TSV with `design_id` and `job_id`,
where each job ID identifies one task. Run it in a small Slurm job under
account grisman, initially without `--apply`, with an explicit output directory
on the cluster's scratch storage. Review the resulting exact task list, then
run with `--apply` if it matches. It checks live state again, rejects array
parent/range IDs and tasks with prior runtime or restarts, and uses
`scancel --state=PENDING` to protect tasks that start during the check. It
does not submit jobs or rewrite existing array indices. Missing mappings are
listed for inspection, never inferred. All remote actions need an audit back
in this registry before remote application is claimed.

## Local changes applied

Audit job `12826891` removed 38 of 50 pending local tasks, retained 11 selected
tasks, and held one `2xxm_flex_p3` fallback while `2xxm_p1` remains unconfirmed.
All 12 running 64-CPU MARK* tasks were preserved. The two separate 48-CPU jerry
experiments and unrelated workloads were also preserved. See `local_actions.tsv`.

The three retained coverage tasks are `1b6c_flex_p2`, `2hnu_flex_p5`, and
`4znc_flex_p6`. Retained primary targets are `5d68_p3`, `2rl0_p3`, `5a6y_p2`,
`3cal_p4`, `4pxf_p6`, `4u3s_p1`, `5dc4_p6`, and `5it3_p10` (all `_flex_`).
Dependencies in `local_lane_dependencies.tsv` maintain eight compsci lanes
on EPYC 9554 and four fennario lanes on Xeon Gold 5320. Started tasks are not
rebound or restarted. Dependencies can leave a lane idle while its predecessor
finishes; the counts are concurrency caps, not a utilization guarantee.

The two new local tasks are `4wem_flex_p12` (`12827029_0`, fennario) and
`3eb6_flex_p7` (`12827030_0`, compsci). Management job `12827027` replaced their
still-pending GPU-gated predecessors `12826957_0` and `12826958_0` with newly
frozen CPU-only runners. Their sole dependencies are CPU predecessor tasks;
there is no GPU dependency or GPU-result check inside the runner. Existing
frozen packages were preserved. October 16 remains the planning cutoff for
a full 14-day observation by October 30, not a runtime assertion.

Five preflighted high-potential continuations are also submitted:
`1a0r_flex_p0` (`12827031_1`), `4wyq_flex_p3` (`12827032_2`),
`2rl0_flex_p4` (`12827033_3`), `3cal_flex_p5` (`12827034_4`), and
`4u3s_flex_p2` (`12827035_1`). These follow CPU predecessors and keep the total
caps at eight compsci and four fennario tasks. They restore selected higher
tiers, not the intermediate backlog. There are now 18 scheduled future local
designs, plus the one held `2xxm` fallback. Exact jobs and predecessor links are
in `cpu_continuation_submissions.tsv` and `local_lane_dependencies.tsv`.

## One grouped PACK* screening job

The active GPU job is **12827169**, with **12 designs in one allocation** on
four RTX PRO 6000 GPUs, 128 CPU threads, and a **six-hour total limit**. It
retains the allocation between designs, allows at most one hour per candidate,
and reserves 20 minutes for every later candidate. The optimized build and
numerical protocol match the running GPU101 job `12814879`.

Management job `12827168` added `2rl0_flex_p4` and `4u3s_flex_p2` to pair the
restored CPU continuation tasks. Their residue lists, structure checksums and
sequence counts match the active CPU commands. Fresh PACK* input preflights
matched all 20 and 39 sequences, respectively, to the existing MARK* preflights.
Both are two tiers above GPU101's highest same-system designs (`2rl0_p2` and
`4u3s_p0`). They follow the original ten entries without changing their order
or definitions. All 12 entries are in `pack_designs.tsv`; the six originally
new definitions remain in `new_designs.tsv`.

Only never-started jobs `12826949` and its screening job `12826950` were
replaced. Original preparation `12826915` and all earlier frozen packages were
preserved. CPU tasks and their GPU-independent dependencies were not changed.
See `pack_pair_completion.json` for the replacement and input verification audit.
The unused intermediate `5dc4_p5` remains outside this batch.

Screening job `12827170` writes measured times and a 30-minute screening flag
after the GPU group ends. This flag is advisory for subsequent analysis and
does not control any active CPU job. Current output is
`/usr/xtmp/lz280/markstar_deadline_focus_20261006/pack_12827169/screening.json`.
Do not infer speedup from submission or preflight success. A normally
finished workload with incomplete binding estimates is distinguished from
fully estimated sequences; only matched sequence outcomes support a reported
numerical speedup.

## Registries and provenance

`design_actions.tsv` contains all 169 existing and proposed definitions.
`NA` denotes unavailable audit values or the absence of a primary target.
`system_coverage.tsv` covers all 38 systems; 22 had a current-round completion
at the source snapshot. Coverage candidates are provisions, not promises.
`policy.json` records the deadline, priorities, source snapshot, and local audit.
`registry.json` records actual submissions and preflight checks.
Old array manifests remain intact and must not be submitted as active lists.
Full scheduler snapshots and generated artifacts remain under
`/usr/xtmp/lz280/markstar_deadline_focus_20261006`; repository files contain
only source, definitions, small registries and documentation.
