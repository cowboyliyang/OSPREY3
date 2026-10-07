# October 6 additions

Preserve the selected smaller designs and append higher tiers. All nine new
PACK* workloads use four RTX PRO 6000 GPUs in one six-hour allocation after
the existing twelve-design batch, job 12827169. Each design has at most one
hour, with time reserved for later cases. No CPU PACK* batch is submitted.
The optimized build, FP64 protocol, seed 42, residual budget 1 and pair-only
settings match GPU101. GPU DP is enabled with failure on an unavailable GPU
path; CPU CCD and sampling are retained as in GPU101.

| System | Retained smaller MARK* case | Addition | MARK* destination |
| --- | --- | --- | --- |
| 1b6c | p2, queued | p4 | compsci |
| 1gwc | p1, completed | p4 | DCC handoff |
| 2hnu | p5, queued | p11 | compsci |
| 2p4a | p9, retained | **p12** | DCC handoff |
| 2xxm | p1, completed | p7 | DCC handoff |
| 3bu8 | p9, completed | p12 | DCC handoff |
| 3k3q | p0, completed | p3 | compsci |
| 4kt6 | p1, completed | p4 | DCC handoff |
| 4znc | p6, queued | p14 | compsci |

All names use the `_flex_` separator. The already-running `4z80_p1` remains
unchanged. Existing local `2xxm_p3` remains a held fallback; confirming p1
completion does not silently release or cancel that job.

The user requested a conservative 2p4a extension. Historical p9 MARK* took
4.561761 hours; this does not establish a 64-CPU tier-growth law. The new p12
adds three flexible residues, preserves 39 sequences, and has scan branchwidth
5 versus p9's 4. The withdrawn p18 added nine residues and had branchwidth 9.
Its assumed fourfold CPU speedup and 1.8-fold growth per tier were not validated.
p18 was preflighted but never submitted for production. p12 is an exploratory
case with no near-14-day runtime forecast. Branchwidth is not a MARK* timing
model. Other near-14-day targets are also uncertain selections, not measured
runtime predictions; see `selection.json`.

The lower PACK* tiers for seven systems belong to GPU101. The lower tiers of
2p4a and 3bu8 belong to the CPU cohort, but the user explicitly requested GPU
for their additions. Preserve this hardware distinction when comparing tiers.
No speedup is established until matching MARK* and PACK* outputs exist.

## Scheduling and external handoff

`registry.json` records actual local submissions. Four new 64-CPU MARK* jobs
use EPYC 9554 and follow existing compsci lane tails, preserving its eight-job
cap and all smaller queued work. They do not wait for GPU results. No existing
Slurm job is canceled, replaced, released, or rebound by this addition.

Five DCC definitions are appended to the active sibling campaign's handoff,
after its original fourteen future entries. The combined DCC concurrency
remains six, including already-running work, at 64 CPU threads per task and a
14-day limit. These additions have **not been submitted remotely**. Use
`deadline_focus_20261006/dcc_launch_priority.tsv` for the combined nineteen-entry
order and its `dcc_new_designs.tsv` for the nine definitions added across both
phases. `handoff_index` is a manifest index, not an existing Slurm task ID;
`task_id` and `group_index` retain their source-package provenance.

For the five new DCC cases, verify the actual CPU model of the current-round
same-system anchor before launch: 1gwc p1, 2p4a p9, 2xxm p1, 3bu8 p9 and
4kt6 p1. If p9 has not started for 2p4a, choose one CPU model for both p9 and
p12. Retain all existing smaller tasks before these additions. The TSVs carry
residue lists, PDB checksums, sequence counts, heaps and historical command
references. Local paths are provenance, not a portable DCC execution package.

The corrected external result for `2xxm_flex_p1` is 6338 seconds, completed
October 6 at 15:27:31 Eastern, with 20 output rows and 19 fully Estimated rows.
The active external snapshot records 30 completed designs after this correction.

Prepared inputs and run outputs stay under
`/usr/xtmp/lz280/markstar_near14_20261006`. Earlier p18 preparations are preserved
for provenance and are not active launch inputs.
# GPU memory retry

The original PACK* job `12827321` failed its preflight memory guard before
producing any full measurements. Active retry `12827800` explicitly requests
1000 GiB with the existing four RTX PRO 6000 GPUs, 128 CPUs, and six-hour
limit; its audit is `12827801`. The first design entered full GPU computation.
The twelve-design group follows as `12827802`. Frozen inputs and previous
outputs are preserved, and MARK* tasks retain their independent dependencies.
See `../gpu_memory_retry_20261006/registry.json` for the retry history.
