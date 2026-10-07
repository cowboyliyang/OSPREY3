# Current design inventory

The versioned [unified plan table](CURRENT_PLAN.md) contains all 38 systems
in the requested format, including the 5dc0 PACK* records. The compact
[plan TSV](current_plan.tsv) also includes actual job IDs and dependencies.
These are dated snapshots; DCC plans without receipts remain labelled.

`registry.json` identifies the latest Slurm-generated inventory snapshot.
`systems.tsv` gives all 38 systems' completed, running, local pending, remote
previous-submission and remote new-plan tiers. The full overview, per-design
timings and provenance, and the inventory of deferred definitions are in the
scratch output directory recorded by the registry.

Refresh with `refresh.py` through Slurm under account grisman. The collector
merges the verified DCC completed38 transfer, local completion artifacts,
current local scheduler status and the current DCC handoff manifest. It uses
the prior timing report as its initial design and PACK* reference inventory.
Completed output counts do not mean every sequence is Estimated; the full
design table retains row quality and jointly Estimated counts separately from
raw wall-time ratios.

DCC original submissions have no current scheduler snapshot here. Nine new
definitions originally had no remote submission receipts; six remained after
the October 7 removal and the subsequent addition of `5dc0_flex_p9` brings
the current count to seven. Both categories are kept
explicitly separate. The two existing 48-CPU jerry comparison runs are recorded
separately in the registry and excluded from the 64-CPU tier counts.

On October 7, the user removed `3k3q_flex_p3`, `4znc_flex_p14`,
`2xxm_flex_p3`, `4kt6_flex_p4`, `1gwc_flex_p4`, and `2xxm_flex_p7`.
Management job `12832018` cancelled the three never-started local tasks and
withdrew the three DCC definitions from the current handoff. Remote scheduler
state remains unverified. `excluded_designs.json` is the persistent exclusion
list used by every current inventory refresh. The current per-design and
all-definition tables omit these six; original definitions, frozen packages,
results, and earlier snapshots remain available for provenance. Full backups
and the removal audit are under `/usr/xtmp/lz280/markstar_remove_six_12832018`.

The collector reloads GPU101's latest status and wall-time files, including
results that arrived after the initial timing snapshot.

For the user's plan-review table, follow `reporting_preferences.json`.
Use italics only in the current-round MARK* cell when that result exists.
Keep additions in the main table, use only `保留` or `条件保留` in its
recommendation column, and show the 5dc0 PACK* records below it. The latest
request adds 5d68 p4 to the review plan and removes the three previously
suggested holds from that view. Existing running and completed work remains
protected. This review configuration does not itself change Slurm jobs.

October 7 management job `12832108` applied the later user decisions:
`5dc0_flex_p9` is now in the DCC handoff, without a remote submission receipt.
`3eb6_flex_p7`, `4wem_flex_p12`, and `5dc4_flex_p6` remain in the table as
`条件保留`, with MARK* state `暂不排队`; their never-started submissions were
cancelled. `conditional_deferred.json` prevents reporting them as queued.
Do not automatically resubmit them. The three newly defined missing PACK*
cases (`5dc0_p9`, `5d68_p4`, `3k3q_m1`) share PRO 6000 job `12832109`.
See [the supplement registry](../plan_supplements_20261007/registry.json).

Job `12832118` subsequently cancelled the three never-started omitted tasks
`2rl0_p4`, `3cal_p5`, and `4u3s_p1`. The retained `4u3s_p2` now follows
running `4u3s_p0` directly, preserving its serial lane. All actual local
64-CPU campaign tasks were checked against the unified table; none remains
outside it. At that snapshot, five planned additions were explicitly `待提交`; this
alignment did not submit them or apply any remote DCC scheduler changes.
See the [queue alignment audit](../queue_alignment_20261007/registry.json).

Job `12832174` raised the local concurrency limits from 8 + 4 to
9 compsci + 5 fennario tasks. `1a0r_p0` and `5it3_p10` were released from
their predecessor dependencies and both started, giving 14 running and eight
pending main 64-CPU tasks in that snapshot. No tiers were added, and the
conditional deferrals remain in effect. See the
[concurrency adjustment](../local_concurrency14_20261007/README.md).

Job `12832726` subsequently submitted all five retained additions:
`1gwc_p2` and `2xxm_p2` on fennario; `3k3q_m1`, `4znc_p5`, and `5d68_p4`
on compsci. The refreshed table has 14 local tasks running and 13 pending,
with no remaining local `待提交` rows. It also includes the completed PACK*
measurements for `5dc0_p9`, `5d68_p4`, and `3k3q_m1`. See the
[five-task submission audit](../queued_five_20261007/README.md).

The `MARK*位置` column identifies the recorded destination, with DCC in bold.
Completed DCC results, prior submissions awaiting verification, and plans
without submission receipts remain distinct in the MARK* status column.

Generate this view through Slurm with `--account=grisman` using
[`report.slurm`](report.slurm), which also exports the versioned plan snapshot.
The script writes `unified_report.html`, `unified_plan.tsv`, and
`5dc0_packstar_records.tsv` under its scratch snapshot directory.
The stable [unified table link](latest_report/unified_report.html) and
[TSV link](latest_report/unified_plan.tsv) point to the most recent successful
render. Plan counts include conditional tiers; completed MARK* cells alone
are italicized.

Export the latest small versioned snapshot by running `export_plan.py`
through Slurm after a successful render. Runtime artifacts remain in scratch.
