# Local MARK* concurrency: 9 + 5

Management job `12832174` completed successfully on October 7, 2026 under
account `grisman`. The current local limits are nine compsci tasks and five
fennario tasks. Each task retains 64 CPUs, 192 GiB memory, and its existing
14-day time limit. The DCC planning limit remains six; remote state was not
queried or changed by this operation.

The operation removed the predecessor dependency from two retained tasks:

| Design | Job | Previous predecessor | Observed running node |
| --- | --- | --- | --- |
| `1a0r_flex_p0` | `12827031_1` | `12812439_26` | `compsci-cluster-fitz-44` |
| `5it3_flex_p10` | `12826673_4` | `12812439_0` | `fennario-03` |

Before release, fitz-44 had 106 unallocated CPUs and 380948 MiB free memory;
fennario-03 had 104 unallocated CPUs and 300340 MiB free memory. Slurm's
unallocated memory also exceeded each task's request. Both tasks subsequently
entered RUNNING, bringing the main 64-CPU campaign to 14 running and eight
pending tasks in the refreshed snapshot. No design was added to the plan.

The script verified the dependency graph, account and resources, preserved
all 12 previously running tasks, and left conditional tiers deferred. Other
dependencies remain in place. The root counts enforce the current plan's
concurrency; Slurm resource availability still determines actual start times.

[`registry.json`](registry.json) records the application and its immediate
scheduler states. Full before/after snapshots, command output, and configuration
backups are in `/usr/xtmp/lz280/markstar_local_concurrency14_20261007/12832174`.
The [unified plan](../current_inventory_20261006/CURRENT_PLAN.md) records the
subsequent running states. Do not rerun the one-time application script.
