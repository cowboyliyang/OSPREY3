# Current design inventory

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
definitions have no remote submission receipts. Both categories are kept
explicitly separate. The two existing 48-CPU jerry comparison runs are recorded
separately in the registry and excluded from the 64-CPU tier counts.
