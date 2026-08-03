# SPACK\* — PACK\*-Specialized Sequence-Level Pruning — Design & Implementation Plan

> Name: **SPACK\*** — a sequence-level pruning algorithm specialized for PACK\*.
> It is to PACK\* what BBK\* is to MARK\*.
> Status: implemented research prototype.  Small exhaustive tests pass and a real 8,000-sequence
> three-position run completed in 5.42 h; four-position scaling and end-to-end PAC accounting remain
> open.  The active 2026-07-16 scaling route is the deterministic binding-interaction cancellation
> bound described in the recovery checkpoint below.  Built on the
> `edu.duke.cs.osprey.spackstar` framework (renamed from `varkstar`).

## 0. 2026-07-16 recovery checkpoint — current source of truth

This section supersedes earlier performance-roadmap text when the two disagree.  It records the
state recovered from the interrupted Codex session and the changes completed after recovery.

### Active mathematical route

The original AA-level envelope is fast and uses zero PACK\*/CCD calls at query time, but it is too
loose on the corrected 3BUA N=6 state split.  The rigid-A\* denominator leaf screen is also too
loose because its unminimized denominator witnesses clash badly.  Neither is the primary scaling
route now; both remain deterministic fallbacks only.

The active bound cancels protein and ligand internal energies before optimizing.  For compatible
protein and ligand RC assignments `p,l`, if

```text
C(q) >= P_min(p) + L_min(l) + I_lower(p,l)
```

for every compatible complex conformation `q`, then

```text
K* = Z_complex / (Z_protein Z_ligand)
   <= exp(-min_(p,l) I_lower(p,l) / RT).
```

`BindingInteractionEnergyMatrixCalculator` builds a factorwise lower matrix containing only the
modeled bound-minus-unbound interaction: the exact matrix constant difference, complex-minus-state
reference-energy offsets on one-body factors, flexible-to-opposite-shell interaction minima, and
cross-partner flexible pair factors.  `BindingInteractionKStarBounder` converts those RC factors to
an AA objective and maximizes it over every completion of a partial sequence with max-sum
mini-bucket.  Query-time work is deterministic and uses zero PACK\* calls and zero CCD calls.

This route is fail-closed.  Production construction requires the traditional energy partition,
a minimizing complex calculator, matching state matrices, matching residue-entropy settings, no
unsupported higher-order complex terms, an exact protein/ligand partition of complex positions,
and aligned state/complex RC counts and residue types.  A failed check disables the interaction
bound and enters the configured fallback path; a non-finite query factor returns `+inf`, never an
optimistic finite upper.

### Measurements that remain valid, and one result that does not

- Corrected 3BUA N=6 rigid AA-envelope probes were approximately `logKUpper=62,921` at i-bound 2
  and `62,882` at i-bound 3.  They are far too loose to prune usefully.
- The deterministic rigid leaf screen was approximately `logKUpper=59,510`, dominated by a ligand
  rigid-witness lower bound near `-66,638`; increasing the witness target from 1 to 8 or 16 did not
  fix the clash-driven gap.
- Slurm job `12097451` is **invalid as a quantitative calibration**.  It compiled live SPACK\*
  sources while those sources were still changing: its output contains the old
  `crossShellPairs=...` marker although the final calculator reports `modeledConstant=...`.
  Independently, the old probe supplied rigid protein/ligand matrices and rigid reference energies
  to the interaction calculator, while production assembly supplies minimizing versions.  Its
  reported interaction values (including root `logKUpper≈442`) are qualitative evidence that
  cancellation helps, but must not be quoted as final numbers.
- The first valid frozen-version 3BUA N=6 interaction calibration completed as job `12097791` on
  2026-07-16 from snapshot
  `/home/users/lz280/tmp/spackstar_source_snapshots/3bua_probe.UdCWRTbF`, source digest
  `82786b8a63f7fc06bc3f0b6564f8d5619c01f133e820b9394c8f655d9d8c04c7`.  The verified run used
  production-equivalent minimizing inputs and did not build legacy rigid matrices.  Setup plus
  interaction-matrix construction took `444.202 s`; the interaction matrix itself reported
  `modeledConstant=0.0`, `minimizedSingles=1134`, `copiedFlexiblePairs=285768`, and `12.929 s`.
  The root was `logKUpper=127.9580775` (`log10KUpper=55.571487`) for both i-bound 2 and 3.  Both
  i-bounds returned the same upper for every reported root/child/seed query, so i-bound 3 added no
  tightness.  I-bound 2 used `374 ms` over `945` bound calls; i-bound 3 used `2315 ms`; both made
  zero PACK\* oracle calls and zero CCD calls.  This is dramatically tighter than the old
  `~59,500--62,900` rigid bounds and fast after setup, but its actual prune fraction remains
  unknown because there is still no verified incumbent lower cutoff for this exact experiment.
  The recorded
  `topNLowerCutoff=35.102578` belongs to an older 2RL0 run and must not be reused for 3BUA.

### Code state after recovery

- `SpackStarPackStar` is now interaction-first.  It attempts the interaction matrix/bound before
  constructing either rigid fallback.  When interaction succeeds, the same bound serves partial
  nodes and leaf preflight and the rigid AA envelope/rigid-A\* leaf objects are not built.  Only an
  unavailable interaction bound triggers the explicitly configured old fallbacks.
- `SpackStarDriver` already schedules a full-sequence preflight before the primary PACK\* leaf
  oracle.  A preflight prune starts no PACK\*/CCD work; a non-pruning upper falls through to the
  primary oracle.  Partial `+inf` nodes expand immediately instead of being repeatedly requeued.
- `RunSpackStarEnvelopeProbe` now matches production interaction inputs: minimizing matrices and
  minimizing reference-energy conventions for complex, protein, and ligand.  Old rigid envelope
  and leaf diagnostics are opt-in via `spackstar.probe.legacyRigidDiagnostics=true`.
- The probe now reports root and every evaluated child upper along a greedy depth path, per-node
  time, per-depth prune fraction/effective-child count when a real
  `spackstar.probe.incumbentLogKLower` is supplied, suggested leaf bounds, and aggregate
  PACK\*/CCD counters.
- `submit_spackstar_3bua_envelope_probe.sh` freezes and compiles SPACK\* Java sources at submission
  time, writes a SHA-256 source manifest, and submits the Slurm job with the immutable class
  snapshot.  The Slurm runner verifies that digest and never compiles the live repository.
- `RunSpackStarIncumbentProbe` is a leaf-only follow-up: it asks the production interaction bound
  for a small ranked candidate set, calls the existing `PackStarVerifier` on only those leaves,
  takes the kth finite PACK\* lower as an incumbent cutoff, and immediately reruns the deterministic
  depth calibration with that cutoff.  It never starts the full sequence-tree driver.  Frozen job
  `12101659` started on 2026-07-16 with three candidates, top-2, i-bound 2, a 4096-conformation
  per-state cap, configured PAC delta `0.001`, snapshot
  `/home/users/lz280/tmp/spackstar_source_snapshots/3bua_incumbent.6JtWSz7R`, and source digest
  `700370aec24faed5ce3c6b77db3c389564721743f02b340a89c7fc14858ff272`.  That first attempt
  completed matrix setup but stopped before any PACK* candidate call because its CPU-vs-GPU root
  reproducibility tolerance was incorrectly set to `1e-8`: the reference root
  `127.9580775217` and GPU root `127.9580303224` differ by only `4.72e-5`.  This is a guard
  calibration failure, not an interaction-bound failure.  The corrected rerun uses tolerance
  `1e-3` and reduces GPU streams from 832 to 32.  It started as job `12116542` from the same
  immutable snapshot/digest; the job header records both corrected settings.  Here `gpuStreams=4`
  means four energy/minimization streams per GPU, hence 32 across eight GPUs.  It does **not**
  reduce PACK* Branch-DP parallelism: the runner separately enables multi-GPU DP with
  `branchdp.dp.gpu.maxGpus=8`, whose normal execution has one main compute stream per participating
  GPU (and hybrid execution may add one upload stream per GPU).  PAC sampling is CPU-only in this
  probe (`packstar.pac.sampling.gpu=false`) with 104 configured CPU threads.
- `compileJava compileTestJava` succeeds.  Direct focused entries for assembly selection,
  driver/leaf preflight, interaction/envelope exhaustive safety, AA mini-bucket, and small-state
  bounders all pass.  A combined Gradle test worker was killed with exit 137 by the environment;
  it emitted no assertion failure, and the same compiled focused tests then passed individually.

The prototype source/test directories and `run_spackstar_*.slurm` remain hidden by local
`.git/info/exclude`, so ordinary `git status` is not a checkpoint.  Per-job immutable snapshots now
make calibration runs reproducible, but a deliberate Git checkpoint is still needed before broad
algorithm changes.

### End-of-day handoff — 2026-07-16 23:59 ET

Leave corrected incumbent job `12116542` running; do not resubmit it merely because this Codex
session ended.  It began on `fennario-05` at `2026-07-16T23:50:39-04:00` with the immutable
snapshot/digest recorded above.  At this checkpoint stderr was empty.  The first minimizing
complex energy matrix (`536,949` entries) finished in `7.1 min`, versus `10.6 min` in the invalid
832-stream attempt, and the second matrix had reached `31.6%`.  Reference-energy timings also
improved from `15.5 s` to `2.5 s` for the first 1,134-RC calculation.  Thus lowering the pool from
104 to four streams per GPU has not slowed initialization so far; the earlier 832-stream pool was
GPU oversubscription, not useful PACK* DP parallelism.  The invalid attempt also ended with
`maxRssKiB=327196140` and many leaked-stream cleanup warnings, which are additional reasons not to
restore that setting.

The live logs are:

```text
/usr/xtmp/lz280/slurm_logs/spack_inc_3bua_12116542.out
/usr/xtmp/lz280/slurm_logs/spack_inc_3bua_12116542.err
```

The first action in the next session is to inspect those two files and the Slurm terminal state.
Search the output for `SPACK_INCUMBENT_WALL`, `[INCUMBENT-PROBE]`, `[INCUMBENT-START]`,
`[INCUMBENT-CANDIDATE]`, `[INCUMBENT-FAIL]`, `[INCUMBENT-CUTOFF]`, and the post-cutoff depth
calibration.  Record:

1. total matrix/setup time and peak RSS;
2. reproduced interaction root (expected GPU value about `127.9580303224`, accepted within
   `1e-3` of the frozen CPU reference);
3. each of the three candidate sequences, per-candidate wall time, PACK* lower/upper, and state
   oracle/CCD counts;
4. the top-2 finite lower cutoff;
5. per-depth child-prune fraction and effective-child count against that real cutoff.

This is the first corrected run that can reach PACK* candidate verification, so there was no valid
per-candidate timing at handoff.  The rough remaining-time estimate was 30--90 minutes, but the
first completed candidate must replace that estimate.  If the job fails before producing two
finite candidate lowers, diagnose the recorded failure and do not treat the absence of a cutoff as
evidence that the interaction bound cannot prune.  In particular, do not reuse job `12101659` or
its lack of candidate results as a bound-quality result.

### Exact continuation order and decision gate

1. Smoke-check the snapshot submit path and record the emitted source digest and Slurm job ID here.
   **Done:** job `12097791`, digest
   `82786b8a63f7fc06bc3f0b6564f8d5619c01f133e820b9394c8f655d9d8c04c7`.
2. Run corrected 3BUA N=6 once with the same frozen snapshot for i-bounds 2 and 3, legacy rigid
   diagnostics disabled.  Record root, every depth summary, per-node cost, and zero PACK\*/CCD
   query counters.  **Done as job `12097791`; results are recorded above.**
3. Do not invent a cutoff.  If no compatible verified N=6 incumbent exists, first use the
   interaction-ranked leaf proposals to obtain one controlled PACK\*-verified lower bound, then
   evaluate the frozen child uppers against that cutoff.  Leaf-only job `12101659` was a
   no-candidate guard-calibration failure as described above; the corrected rerun retains the
   requirement that top-2 needs at least two distinct candidates with finite lower bounds.
   **Corrected rerun: job `12116542`, running at the 2026-07-16 handoff; collect it before doing
   anything else.**
4. Consider N=16 only if the measured effective child count is roughly 2--3 per expanded parent
   and the projected PACK\*/CCD leaf-oracle count falls by at least one order of magnitude.
5. If that gate fails, strengthen correlation in the cancellation bound (for example, retain
   unbound-weight/interaction coupling or partitioned witnesses).  Do not spend a larger cluster
   allocation merely to brute-force the exponential frontier, and do not blindly raise i-bound.

Scope remains strict: modify only SPACK\* main sources, SPACK\* tests, SPACK\* runners, and this
SPACK\* design document.  Do not modify `packstar` or `branchdp` for this workstream.

## 1. What this is

A sequence-level branch-and-bound / pruning algorithm that sits **on top of PACK\***
(the randomized, branch-decomposition **PAC** partition-function oracle) to find the
top-k K\* sequences **without running PACK\* to convergence on every sequence**.

It is to PACK\* what BBK\* is to MARK\*:

| sequence driver | Z-oracle | oracle bound type |
|---|---|---|
| BBK\* | MARK\* | deterministic, anytime |
| **SPACK\*** | **PACK\*** | **randomized (PAC), anytime** |

We reuse the already-built `spackstar` (formerly `varkstar`) sequence-tree B&B + pluggable
bounders. The new work is (a) a PACK\*-backed bounder and (b) the parts of the driver that a
*randomized* oracle forces: e-value pruning, demand-driven budgets, and (later) a ratio estimand.

## 2. Why PACK\* needs its own driver (not just BBK\*)

- **PAC bounds, not exact.** PACK\* returns high-probability brackets `[L,U]`. BBK\*'s prune
  `U(node) < L(incumbent)` is one deterministic test; with PAC bounds **and** adaptive
  re-bounding it is not sound → we need **anytime-valid e-value** tests.
- **Fixed per-call floor.** Each PACK\* call has a sampling+CCD floor (hundreds of seconds).
  Naively bounding every sequence pays that floor on losers — the 2q1e/ARG 8.57 h pathology.
  The driver must kill losers with the **cheapest possible loose bound** and never refine them.
- **Immune to wide-Z blow-up.** PACK\* does not explode where MARK\*/MICA time out
  (4u3s, 4wyu, 3ma2; 1a0r/2xgy 14-day timeouts). The driver's payoff is concentrated there.
- **BBK\* is forced to use an equal fixed budget per node** ("the number of confs must be the
  same for every node ... not sound to do epsilon-based iterative approximations"). The
  **e-value layer removes that constraint** → unequal, adaptive, demand-driven budgets become
  sound. *This is the core enabling change of the whole project.*

## 3. The algorithm

The core mechanism: **the pruning certificate is an e-value (test martingale / e-process).**
Everything else (estimand, information geometry) is secondary and optional.

**State.** Sequence tree; a node = partial sequence = superset over its completions. Each node
carries a K\* bracket `[logKLower, logKUpper]` produced by PACK\*.

**Node bound (subtree-safe — already implemented by `SubtreeSafeKStarBounder`).**
- complex: PACK\* **upper** over the superset RCs (`sequence.makeRCs(complexConfSpace)` puts
  *all* res-types on unassigned positions).
- protein, ligand: PACK\* **uniform lower** over the subtree.
- `logKUpper(node) = logZUpper(complex) − logZLower(protein) − logZLower(ligand)`.

  This is a valid (PAC) upper bound on `K*(s)` for *every* completion `s` (admissibility),
  because the superset Z ≥ each completion's complex Z.

**Driver loop (best-first, anytime).**
1. Pop the node with the highest `logKUpper`.
2. **e-value prune test** vs `incumbent_k`: maintain an e-process for
   `H0: K*(node) ≥ K*(incumbent_k)`; draw a few more PACK\* samples on the deciding side;
   if `e > 1/δ_node` → prune the whole subtree. (The certificate *is* this e-process.)
3. If leaf and it survives → refine via the verifier, update incumbents.
4. If internal → choose: **(a) tighten** its bound with more PACK\* samples and re-queue
   (catch-and-release), or **(b) branch** one free position and push children with a cheap
   loose bound. Pick (a)/(b) by a **branchwidth-aware** cost-vs-gap rule.

**Granularity.** Edge-granular lazy expansion: push children with an inherited optimistic
bound, compute the real (expensive) PACK\* bound only when a child reaches the queue top.
Fits PACK\*'s per-call cost — never bound all ~20 siblings eagerly.

**Sharper estimand (Stage 3, optional).** For near-top leaves vs incumbent, estimate `ΔlogK*`
*directly* by reweighting the incumbent's samples onto the competitor (MBAR / cross-sequence
bridge through the **shared** branch decomposition). Information geometry `KL/χ²` is only the
overlap diagnostic / sample scheduler here — it affects the *rate* of evidence accumulation, not
the validity of a prune. Not a core mechanism; future work.

## 4. Mapping onto the `spackstar` framework

**Reused as-is**
- `SpackStar.run()` — best-first B&B, pop by `logKUpper`, incumbent cutoff (driver skeleton).
- `SubtreeSafeKStarBounder` — assembles subtree-safe K\* upper from state bounds (superset combiner).
- `KStarBounds` / `PartitionBounds` — bound types + the K\* = complex/(protein·ligand) combination.
- `OspreyWmbMfStateBounder` — WMB upper + MF lower deterministic state bounder (kept as a cheap tightener).
- `SequenceVerifier` — leaf-refinement hook.

**New classes** (in `spackstar` / `spackstar.bound`)

| class | role |
|---|---|
| `PackStarStateBounder` (`StatePartitionBounder` + `StateSubtreeLowerBounder`) | run `PackStarBound` on `sequence.makeRCs(confSpace)`; return `[logZLower, logZUpper]`; budget-parameterized / anytime |
| `CombiningStateBounder` (`StatePartitionBounder`) | intersect two brackets `[max(L), min(U)]` — wrap PACK\* with WMB/MF to tighten for free |
| `WmbImportanceSamplingStateBounder` (`StatePartitionBounder`) | WMB-IS state sequence-subtree estimator; samples from `WeightedMiniBucket.Proposal`, scores unique conformations with CCD, reports `cv`/ESS/unique/CCD diagnostics |
| `BindingInteractionEnergyMatrixCalculator` | fail-closed construction of the factorwise bound-minus-unbound interaction lower matrix |
| `BindingInteractionKStarBounder` (`SequenceSubtreeBounder`) | deterministic internal/leaf K\* upper using interaction cancellation plus AA max-sum mini-bucket; zero query-time PACK\*/CCD |
| `SequenceKStarEnvelopeBounder` (`SequenceSubtreeBounder`) | older rigid-denominator AA envelope, retained only as an interaction-unavailable fallback |
| `DeterministicLeafKStarBounder` (`SequenceSubtreeBounder`) | older full-sequence complex-WMB/rigid-A\* preflight, retained only as an interaction-unavailable fallback |
| `EValuePruner` | per-node e-process for `node ≥ incumbent`; δ-spending schedule; `shouldPrune(node, incumbent, budget)` — **the certificate engine** |
| `PackStarVerifier` (`SequenceVerifier`) | refine a chosen full sequence with PACK\* to target ε |
| `SpackStarDriver` (extends the `SpackStar` base) | catch-and-release re-bounding, demand-driven budgets, e-value prune, edge-granular lazy expansion |
| `BranchwidthPositionOrder` | replace `findFirst()` with branchwidth-reduction position scoring (ordering-only) |

**Driver extension needed.** The current `SpackStar.run()` bounds each child once, never re-bounds
an internal node, and prunes with the deterministic `logKUpper <= cutoff`. PACK\* needs
catch-and-release (re-pop → add samples → re-queue) and the e-value prune. Implement these in
`SpackStarDriver` so the deterministic `SpackStar` base stays intact.

## 5. Implementation stages (each independently testable)

**Stage 0 — PACK\* as a SpackStar bounder ("randomized BBK\*").**
`PackStarStateBounder` wired through `SubtreeSafeKStarBounder` into stock `SpackStar`. Fixed
per-node budget, deterministic prune (treat PAC `[L,U]` as exact — *not yet sound*, but gives an
end-to-end baseline). Test on 2q1e (easy), 4u3s/4wyu (hard, MARK\* timeout), 3gxu (easy small):
ranking + wall vs enumerate-all PACK\*; log expanded/pruned/verified.

**Comparison protocol.** Keep oracle pairs separate:
- **BBK\* → MARK\*.** Run BBK\* with `MARKStarBoundFastQueues` as the leaf pfunc and compare against
  an enumerate-all `KStar` baseline using the same MARK\* pfunc. This measures BBK\*'s sequence-layer
  savings over MARK\*.
- **SPACK\* → PACK\*.** Run SPACK\* with `PackStarStateBounder` / `PackStarVerifier` and compare against
  enumerate-all PACK\* leaf verification. This measures SPACK\*'s sequence-layer savings over PACK\*.
- Cross-oracle comparison is only an overlap/ranking sanity check between MARK\* and PACK\* outputs;
  it is not the correctness baseline for either sequence driver.

**Stage 1 — e-value pruning (the core; soundness).**
`EValuePruner` + `SpackStarDriver` replacing the deterministic prune; δ-spending schedule.
Test: top-k matches brute-force PACK\* ground truth on small designs; empirical wrong-prune rate ≤ δ.

**Stage 2 — demand-driven + edge-granular lazy + catch-and-release.**
Cheap loose initial bound; refine only competitive nodes; branchwidth-aware budget. Test:
wall + #PACK\*-calls drop vs Stage 1 on hard designs; ARG-type losers must get very few samples.

**Stage 3 — ratio estimand (MBAR / bridge) + IG overlap (optional).**
Cross-sequence reweighting for near-top comparisons; `KL/χ²` overlap diagnostic; bridge insertion.
Test: variance of `ΔlogK*` vs absolute-difference baseline; overlap stats on 2q1e/2rfe.

**Stage 4 — dynamic position ordering.**
`BranchwidthPositionOrder` (and/or GNN/Fisher). Ordering-only; expect fewer expanded nodes, no
correctness change. Fills BBK\*'s `// TODO: dynamic A*?` at the sequence layer.

### Implementation status (current)

- **Assembly:** `SpackStarPackStar` (`.Builder` + `StateInputs`/`Assembly`) wires the full stack —
  complex/protein/ligand `PackStarStateBounder`, protein/ligand `OspreyWitnessStateLowerBounder`,
  `BudgetedSequenceSubtreeBounder`, `PackStarVerifier`, and the clip-bias-corrected
  `PackStarIncumbentEvidenceProvider` — into a `SpackStarDriver` in one call.  The state bounders
  can now be wrapped with WMB-IS fallback and deterministic WMB/MF tightening.
- **WMB-IS sequence-subtree fallback:** `WeightedMiniBucket.Proposal` and
  `WmbImportanceSamplingStateBounder` are implemented.  If PACK\* fails or returns a full-range
  state bound, `SpackStarPackStar` can estimate the state sequence subtree with WMB-IS.  PACK\*
  evidence batches remain fail-open because WMB-IS does not yet produce PACK\* residual
  diagnostics for e-value evidence.
- **Real WMB-IS validation:** `RunWmbExactSubtreeValidation` +
  `slurm/scripts/run_wmb_exact_subtree_validation.slurm` exactly enumerated a small real 2RL0
  `complex` sequence subtree (`A156,A164` all-20 mutable, no WT-flex): `35,721` conformations,
  every conformation scored with minimized-energy/CCD, exact `logZ=58.038234`. WMB-IS covered
  this exact value across all tested seeds and sample counts; representative 2048-sample errors
  were `+0.000322`, `-0.000130`, and `+0.000127 logZ`.
- **Large sequence-subtree WMB-IS validation:** `RunWmbSequenceSubtreeValidation` +
  `slurm/scripts/run_wmb_sequence_subtree_validation.slurm` validated two larger 2RL0 `complex`
  root/unassigned sequence subtrees without exact ground truth:
  `2 full + 6 WT-flex` (`A156,A164` mutable plus `A172,A192,A193,G649,G650,G651` flex) has
  `41,047,715,520` conformations, RC domains `[189,189,28,8,19,5,6,9]`; iBound=2 WMB table max
  `1,000,188` cells; 1024 WMB draws gave 225 unique CCD calls, `cv=0.3564`, `ESS/n=0.8873`,
  empirical logZ interval `[155.4625,155.5706]`.  `3 full + 5 WT-flex` (`A156,A164,A172`
  mutable plus `A192,A193,G649,G650,G651` flex) has `277,072,079,760` conformations, RC domains
  `[189,189,189,8,19,5,6,9]`; iBound=2 WMB table max `6,751,269` cells; 4096 WMB draws gave
  483 unique CCD calls, `cv=0.7441`, `ESS/n=0.6437`, empirical logZ interval
  `[156.8708,156.9842]`.
- **WMB-IS finite-sample status:** the estimator route now has a fail-closed fixed-range PAC
  certificate with a local weight cap.  `WeightedMiniBucket.Proposal` first tries an exact
  `max_c(theta_min(c) - log q(c))` cap when the local support is enumerable; for large supports it
  uses a max-sum mini-bucket cap over the proposal's conditional penalty factors.  This fixed a
  real sentinel/range failure where the old fallback cap was dominated by `logQ_lower` and was
  `~1.24e9` log units above observed weights.  On Slurm job `12007350` (`A156,A164` mutable plus
  five WT-flex positions, `1,465,989,840` RC conformations), the cap became finite and local:
  `logWeightUpper=156.6608`, `observedMaxLogWeight=128.6639`, gap `27.9969`.  The returned safe
  upper was still the deterministic WMB upper `139.3654`; the WMB-IS diagnostic interval was
  `[128.2623,128.3396]`, and the fixed-range certificate upper was `152.3940`.
- **WMB-IS table-size hardening:** `WeightedMiniBucket` now checks factor table sizes in `long`
  before allocation and can lower the effective WMB `iBound` to respect a configured
  `maxTableCells` cap.  `WmbImportanceSamplingStateBounder` and the 2RL0 SPACK\* runners expose
  this via `spackstar.run.wmbMaxTableCells`; the default remains unlimited/old behavior.
- **e-value prune (Stage 1):** incumbent→`m` conversion + clip-bias shift (`m' = m − B`) implemented
  and unit-tested; sample pruning is disabled unless a residual-tail bound certifies `B`.
  The clip-excess certificate now tightens the clip-probability term with the minimum of
  Hoeffding and empirical-Bernstein upper bounds, so rare observed clipping gives a smaller
  certified bias shift without weakening soundness.
- **Denominator-side e-value evidence (Stage 1 partial):** `PackStarOracleDiagnostics` can now turn
  protein/ligand PACK\* samples into `boundedMeanGreaterThan` evidence against the incumbent, and
  `PackStarIncumbentEvidenceProvider` can combine complex-side and denominator-side fresh evidence
  batches into one e-process update.  `SpackStarPackStar` wires optional protein/ligand evidence
  sources fail-open; unit tests cover the threshold conversion and provider combination.
- **Global alpha/delta accounting:** `SpackStarErrorBudget` provides named geometric spending
  streams.  `EValuePruner` now allocates one alpha lease per sequence-node e-process instead of
  reusing a depth-level delta across many nodes; clip-excess correction and WMB-IS finite-sample
  certificates consume per-batch/per-state delta leases.  The SPACK\* assembly wires separate
  streams for e-process pruning, PACK\* clip-excess, and WMB-IS certificates.
- **Small ground-truth correctness validation:** `TestSpackStarDriver` now enumerates a small
  synthetic sequence space exhaustively, computes the true top-k, runs `SpackStarDriver`, and
  asserts that the returned top-k matches the exhaustive result while verifying fewer leaves than
  enumerate-all.  `RunSpackStarDesign` now performs the same top-k comparison with real
  PACK\*/OSPREY state oracles and can fail the run on mismatch.  Slurm job `12004527`
  (`slurm/scripts/run_spackstar_ground_truth_validation.slurm`) validated the 2RL0
  one-mutable-position case (`A156`, 20 sequences): SPACK\* top-2 matched exhaustive PACK\*
  ground truth (`A156=LYS`, `A156=ARG`).  Slurm job `12004530` then validated the
  two-mutable-position, single-mutation case (`A156,A164`, 39 sequences): SPACK\* expanded
  2 nodes, pruned 33, verified 6 leaves, and matched exhaustive PACK\* top-2
  (`A156=LYS A164=tyr`, `A156=ARG A164=tyr`).
- **Budgeted PACK\* bounds (Stage 2 partial):** `BudgetSchedule` is now wired into
  `SpackStarPackStar`.  The per-node budget is passed through `BudgetedSequenceSubtreeBounder`
  and caps the PACK\* estimator's per-call sample budget (`samples`, `train`, `pilot`, and
  `maxEst` are capped by the current budget).  The default library builder remains equivalent to
  the old full-budget behavior unless a schedule is supplied; the 2RL0 runners expose explicit
  `spackstar.run.budgetedBounds`, `initialComputeMaxNumConfs`, `budgetGrowthFactor`, and
  `maxComputeMaxNumConfs` settings.
- **Safe bound reuse (Stage 2 partial):** ordinary state `bound()` calls are cached by state
  sequence and budget.  A cached result is reused only when it was computed with at least the
  requested budget; higher-budget results are intersected with earlier valid brackets.  Fresh
  evidence/verifier calls through `boundWithResult()` deliberately bypass this cache so sample
  traces and e-value evidence are not replayed.
- **State-level parallelism and instrumentation (Stage 2 partial):** fully assigned
  complex/protein/ligand state bounds can be evaluated independently when
  `stateBoundParallelism > 1`; CPU runs default to auto-splitting state tasks, while GPU runs
  default to one outer state task to avoid oversubscribing the same GPU.  Driver logs can now
  report `max_remaining_upper`, `topN_lower_cutoff`, gap, node-bound wall time, oracle-call
  deltas, CCD-call deltas, cache-hit deltas, and oracle wall-time deltas.
- **PACK\* GPU implementation is shared, not SPACK\*-specific:** `PackStarStateBounder` constructs
  `PackStarBound`, which enters `PackStarBranchDpBackend` and the shared
  `BranchDpBackend`/`RootedTreeEdge`/`DPGpuFullDP` implementation.  SPACK\* has no separate DP-table
  CUDA kernel.  Consequently the current child-slicing, hybrid, and bounded-OOC code is available
  to a newly compiled SPACK\* run, while SPACK\* frontier scheduling, denominator handling, and
  oracle batching remain separate responsibilities.  `WeightedMiniBucket` and WMB-IS are still
  Java/CPU paths and do not use this branch-DP CUDA implementation.
- **Runner:** `RunSpackStarOracleComparison` (test tree) builds one 2RL0 mutable conf space and runs
  the two proper oracle pairs: BBK\*→MARK\* vs exhaustive MARK\*, and SPACK\*→PACK\* vs exhaustive PACK\*.
  `RunSpackStarDesign` remains the smaller SPACK\*→PACK\* runner and now has strict top-k
  ground-truth validation via `spackstar.run.failOnGroundTruthMismatch`.
- **Active sample-evidence loop (Stage 2 partial):** the driver can now draw multiple fresh
  sample-evidence batches for the same popped node before re-queueing, verifying, or expanding it.
  `maxEvidenceBatchesPerNode` defaults to `1` for the old one-batch behavior; the 2RL0 runners
  expose `spackstar.run.maxEvidenceBatchesPerNode` and report `evidenceBatches`.
- **Betting-mixture and certificate diagnostics (Stage 2 partial):** sample evidence can now use a
  fixed portfolio of bounded-mean betting fractions instead of a single fraction
  (`spackstar.run.bettingFractions`, current Slurm default
  `0.1,0.25,0.5,0.75,0.9`).  Certificate logs include clip/tail diagnostics
  (`clipExcess`, `shiftedM`, `clipped`, `meanX`, `logZCorr`, `essPQ`) and e-process progress
  (`eUpdates`, `eDelta`, `eLogThreshold`, `eLogValue`, `eBatchFactor`) so a run can distinguish
  loose truncated/clipped bounds from loose tail-mass corrections.
- **Still not yet:** production-scale validation/scheduling for denominator evidence, branchwidth
  ordering (Stage 4), ratio estimand (Stage 3), and multi-position PACK\* ground-truth validation
  with nontrivial branch decompositions.
- **2026-07-11 local verification:** the focused SPACK\* suite (`49` tests covering the driver,
  budgets, e-values, evidence conversion, and state bounders) and the core
  `TestWeightedMiniBucket` suite (`18` tests) pass with `-DtestMaxHeap=1g`.  These establish
  component behavior, not the unresolved global-coverage and scaling claims below.
- **Reproducibility:** the local `.git/info/exclude` deliberately excludes the SPACK\* main/test
  source trees, Slurm scripts, and this document from ordinary Git status.  The 2026-07-16 3BUA
  submission path now compiles an immutable SHA-256-manifested source snapshot before `sbatch`, so
  a queued job cannot mix live source revisions.  A deliberate Git checkpoint is still needed
  before broad algorithm changes.

### Current diagnosis and optimization roadmap

SPACK\* should not be judged by whether it can stop as early as BBK\*/MARK\* under a loose
SS-ε certificate.  The intended SPACK\* certificate is stronger:

```text
remaining true K* <= remaining upper <= incumbent lower <= returned true K*
```

This is stricter than popping a sequence from an upper-bound heap after that sequence has an
ε-approximation.  The price is that SPACK\* cannot soundly stop just because the current popped
node has the largest upper bound.  It must either prune or exhaust the remaining search frontier
using a lower-vs-upper certificate (deterministic or e-value/PAC).

**Current pruning diagnosis.** The soundness machinery is more developed than the demonstrated
sample-evidence pruning effect, but the end-to-end global PAC claim is still under audit (see the
2026-07-11 correctness finding below).  Small real 2RL0 validation shows that SPACK\* can prune
correctly, but the observed pruning is mostly ordinary upper-bound pruning.  Sample/e-process and
WMB-IS certificate pruning have not yet contributed much in real runs.  On the 39-sequence
validation, SPACK\* pruned 33 nodes and matched exhaustive PACK\* top-2, but took about `30.5 s`
versus `15.2 s` for exhaustive PACK\*.
Older 3/4-position Slurm logs also show the main blocker more directly: visible
`[SPACK*-CERT]` lines in `spackstar_3p4p_11963832_{3,4}.out` had
`topNLowerCutoff=-inf`, so the PACK\* sample/e-process bound had no finite incumbent threshold
to test against.  This motivated the cheap-seed/early-stop work below.  That finite-cutoff
milestone is now met; the newer bottlenecks are denominator `-inf`, same-depth `+inf` scheduling,
and certificate provenance.

**2026-07-06 Slurm status.** Job `12009291` is the current 3/4-position probe with WMB fallback,
budgeted bounds, and the betting portfolio enabled.  It was submitted after cancelling
`12009288`, because `sbatch --export` had truncated the comma-separated betting-fraction list to
`0.1`; the resubmission used environment export so the logs now confirm
`bettingFractions=[0.1, 0.25, 0.5, 0.75, 0.9]`.  As of the first status check both array tasks
were still running (`12009291_3` on `fennario-01`, `12009291_4` on `fennario-02`) and had not yet
entered the `==== SPACK* -> PACK* chain` phase.  The logs showed repeated BBK\*/MARK\* seed
proposal refinement (`Refining sequence ...`) but no `[SPACK*-CERT]`, `seed-verify`, or
`sampleEvidence` lines.  That means this run is currently measuring seed-proposal latency, not
yet probability-bound tightness or SPACK\* sequence pruning.

**2026-07-06 pivot.** Because `12009291` was still in BBK\*/MARK\* seed proposal after about half
an hour, that path is too slow to diagnose SPACK\* pruning.  It was cancelled and replaced with
job `12009337`, which disables BBK seed by default and uses cheap deterministic initial seeds.
The runner now supports:

```text
spackstar.seed.cheap=true
spackstar.seed.cheapWt=true
spackstar.seed.cheapSingles=true
spackstar.seed.cheapMaxMut=1
spackstar.seed.cheapCount=80
spackstar.seed.stopAfterFiniteCutoff=true
spackstar.seed.extraAfterFiniteCutoff=1
spackstar.seed.bbk=false
```

This produces WT plus single-mutant seeds before SPACK\* starts (`58` seeds for the 3-position
array task and `77` seeds for the 4-position task under the current 2RL0 spaces).  The new job's
logs confirm it enters `==== SPACK* -> PACK* chain` immediately after listing these seeds, so the
next timing question is PACK\* seed verification cost and whether it raises a finite
`topNLowerCutoff`, not BBK\* proposal latency.  Early `12009337_3` output confirms this worked:
PACK\* energy-matrix setup took `312.644 s`; the first seed-verify was still
`topNLowerCutoff=-inf`, the second raised it to `28.283412`, and the third raised it again to
`35.102578`.  The observed seed verification walls were roughly `24-61 s` after state-bound cache
reuse began.  `12009337_4` completed PACK\* energy-matrix setup in `437.645 s` and was still in
its first seed verification at the time of this note.

**2026-07-06 incumbent status.** The incumbent is now good enough to remove the previous
`topNLowerCutoff=-inf` blocker: `12009337_3` reached a finite kth cutoff after the second cheap
seed and improved it after the third, while `12009337_4` reached a finite cutoff after its second
cheap seed (`32.109238`).  That means probability-bound tightness can finally be tested once the
main frontier runs.  The new bottleneck is that verifying all cheap seeds is itself too slow
(`58`/`77` PACK\* leaf verifications).  The driver therefore now supports early stopping initial
seed verification once the kth incumbent cutoff is finite, with a small configurable number of
extra verifications after that point.  The Slurm default is `stopAfterFiniteCutoff=true` and
`extraAfterFiniteCutoff=1`, which would have stopped the 3-position seed phase after the third
seed in the observed log and moved directly into the main search.  This is sound because initial
seeds are only incumbent proposals; unverified seed sequences remain in the normal SPACK\* search
space.

The same log also exposed the next blocker.  Because `12009337` was running the pre-fix driver,
it still verified the full cheap-seed pool and then called the sample-evidence hook on the root
node whose inherited upper bound was `+inf`.  That forced PACK\* state bounds for
`complex:unassigned`: the 3-position task predicted a `24.1 GiB` root DP table and hit GPU memory
limits after entering the main search; the 4-position task predicted `570.4 GiB`.  This is not a
probability-bound tightness failure; it is a scheduling bug.  Evidence is now skipped unless both
the incumbent cutoff and the node upper bound are finite, so the root should expand lazily instead
of invoking WMB/PACK\* on the full unassigned complex.

**2026-07-11 medium-scale outcome (`spackstar_3p4p_12009341_{3,4}`).**  The successor probe used
cheap seeds with early stopping, entered SPACK\* immediately, and produced the first useful
medium-scale throughput result:

- The three-position all-20 space contains `8,000` sequences.  SPACK\* completed in
  `19,498.557 s` (`5.42 h`) with `processed=7,641`, `generated=7,241`, `expanded=362`, and only
  `3` seed verifications.  It pruned `6,806` nodes (`39` internal, `6,767` leaves), verified
  `76` leaves, and returned a top-2.  This is a **positive operational milestone**: a real
  PACK\*-backed search over 8,000 sequences completed in hours rather than failing at setup or a
  DP table cliff.  It should not be described as a speedup yet because this run deliberately
  skipped the same-oracle exhaustive PACK\* baseline and ground-truth comparison.
- The pruning mix also shows what remains inactive: `6,805/6,806` prunes were reported as ordinary
  bound comparisons, only `1` came from sample evidence, and `rebounded=0`.  Thus the completed
  throughput is encouraging, but it does not yet demonstrate that e-process accumulation,
  catch-and-release, or WMB-IS certificate pruning is providing the savings.
- The four-position space contains `160,000` sequences.  It reached the 12 h Slurm limit at
  iteration `5,827`, with `queue=6,879`, `topNLowerCutoff=35.102578`, and
  `maxRemainingUpper=+inf`.  Its recorded actions were `362` expansions and `5,465` lazy bounds,
  with no post-seed lazy prune/verify before timeout.  This run therefore identified a scaling
  boundary rather than invalidating the already-completed 8,000-sequence regime.

**2026-07-11 PACK\*/GPU integration audit.**  The medium-scale jobs predate the current PACK\*
child-slicing/OOC and production-admission work.  Their logs report `decompStrategy=HICKS` and
`rootSplit=work`, not the current production PACK\* combination of adaptive decomposition,
predicted root selection, and calibrated admission.  They nevertheless did enable the older
shared CUDA and multi-GPU DP path, so they reveal an important workload-shape limitation that is
still present in the latest `DPGpuFullDP` planner:

- Latest PACK\* **does** implement row-internal lambda tiling, but that tiling is currently a
  memory/OOC mechanism.  It splits a lambda row into boxes and stably merges partial log-sum-exp
  values on the GPU handling that output region.  It is not a multi-GPU lambda partition.
- `chooseGpuCount(...)` still chooses the number of devices before the full-resident/hybrid/OOC
  plan, primarily from `mStateCount` and `minMStatesPerGpu` (default `4096`).  Normal multi-GPU DP
  assigns disjoint M ranges.  Sliced/hybrid/OOC executors can use union rows, free-M, or M-box
  ranges, but inherit the already chosen device count; lambda boxes are not assigned to distinct
  GPUs.
- This rule is poorly matched to SPACK\* states with small M and enormous lambda domains.  The old
  four-position log repeatedly contains `mStates=114, lambdaStates=303807105` (about `3.46e10`
  lambda-loop iterations) and `mStates=45, lambdaStates=40507614`, both as single-GPU
  `PACK*: GPU DP done` calls.  Larger-M shapes such as `mStates=161595` and `mStates=1928934` did
  use all eight A5000s.  Thus PACK\* was multi-GPU capable, but the card-selection heuristic
  classified some of the heaviest SPACK\* edge shapes as single-GPU.
- Anchored outer-call counts make the distinction explicit: the completed three-position log has
  `74` single-GPU versus `22` multi-GPU PACK\* DP calls; the partial four-position log has `1,032`
  versus `515`.  Their recorded outer GPU times total only about `43.0 s` of the `19,498.6 s`
  three-position run and `1,084.0 s` of the first `43,200 s` four-position run.  In contrast,
  `12,938` completed `[PACK*] Total` estimator spans consume `18,225.0 s` in the three-position
  log, and its initial `DP tables computed` spans total another `1,119.1 s`.  Therefore improving
  small-M GPU occupancy is worthwhile, especially for the four-position tail, but child-kernel
  speed alone cannot explain or remove the overall SPACK\* cost.

**Work predictability boundary.**  Unlike exhaustive PACK\*/K\*, SPACK\* cannot know its exact
whole-run work before starting: which frontier nodes are bounded, pruned, rebound, or verified
depends on the observed bounds and incumbent cutoff.  It can still preflight each individual
protein/ligand/complex oracle after a node is known.  The appropriate integration is therefore
per-node, dynamic admission: preview each candidate oracle without materializing DP tables,
compare its predicted cost with immediately branching the sequence node, and update the remaining
frontier forecast online.  A whole-run worst case can be reported by assuming exhaustive leaves,
but it is generally too pessimistic to be an admission rule.

**Refined `+inf` diagnosis.**  The four-position timeout was not primarily a WMB-IS or complex
branch-DP failure.  Candidate root-split edges with `mStates > Integer.MAX_VALUE` were skipped, but
PACK\* selected other decompositions and returned finite complex-state `logZ` brackets.  SPACK\*
then still reported `nodeUpper=+inf`.  In the current combiner this can happen when the uniform
protein or ligand subtree lower bound is `-inf`; `OspreyWitnessStateLowerBounder` returns `-inf`
when higher-order terms, missing witnesses, or an infinite witness energy prevent a finite
certificate.  A lazily bounded node that remains `+inf` is re-queued, so thousands of same-depth
`+inf` nodes are expensively bounded before any of them branch.  No `[SPACK*-WMBIS]` fallback was
recorded in this probe, so more WMB-IS samples would not address this particular bottleneck.

**2026-07-11 correctness audit.**  The current counters call `upper <= cutoff` a
`deterministic` prune even when the upper came from an ordinary per-call PACK\* PAC bracket.  The
medium run used `PAC_CONFIDENCE=0.05` for thousands of calls, while the global geometric streams
currently cover node e-processes, clip-excess batches, and WMB-IS certificates—not these ordinary
PACK\* bracket calls.  In addition, the three-position stderr contains `342` finite cases where a
PACK\* state bracket and the deterministic WMB/MF bracket are disjoint; the current combiner logs
the crossing and keeps the PACK\* bracket.  This does not show that the returned top-2 is wrong,
but it means the advertised global error guarantee has not yet been established.  Bound
provenance must distinguish deterministic bounds from PAC bounds, crossed certificates must fail
closed until their model/coverage mismatch is explained, and every randomized direct-prune path
must consume a global lease or be converted into valid e-process evidence.

The highest-priority near-term work is:

1. **Checkpoint and reproduce the current prototype.** Preserve the exact SPACK\* source, runner,
   configuration, and log-summary command that produced the 8,000-sequence result before changing
   certificate or scheduling semantics.
2. **Close global PAC accounting before claiming strict top-k coverage.** Add explicit bound
   provenance (`deterministic`, `PAC(delta lease)`, `heuristic/order-only`); allow direct
   `upper <= cutoff` pruning only for deterministic or globally leased bounds; make crossed
   PACK\*/WMB intervals fail closed; and obtain `residualBoundKcal` from a certified backend source
   rather than a configured assumption.
3. **Make every `+inf` node progress.** Compute the cheap denominator witness bounds before an
   expensive complex PACK\* call.  If either denominator is `-inf`, branch immediately and prefer
   a position implicated in the witness failure.  After lazy bounding, expand a node whose real
   upper remains `+inf` instead of re-queueing it behind thousands of equally infinite nodes.
4. **Strengthen the uniform denominator certificate.** Instrument whether `-inf` came from
   higher-order terms, a missing residue-type witness, or an infinite single/pair energy.  Then add
   compatible multi-witness selection or another mathematically uniform lower bound.  Branching is
   the immediate scheduling fix; a finite shallow denominator bound is the route to real internal
   subtree pruning at four positions and beyond.
5. **Integrate current PACK\* planning as a dynamic SPACK\* oracle cost model.** Enable and retain
   the current allocation-free root/decomposition preview for each requested state, but do not
   apply exhaustive PACK\*'s whole-case SLA literally.  Use the preview to choose `bound now`
   versus `branch now`, and cache compatible decomposition/planning structure across related
   parent/child state restrictions.  The SPACK\* runner must explicitly select the desired
   adaptive/predicted PACK\* policy; merely recompiling the shared kernels does not enable every
   production policy.
6. **Fix small-M/large-lambda GPU planning in the shared PACK\* backend.** Replace the
   M-cardinality-only `chooseGpuCount` rule with a plan-level makespan comparison using
   `M*lambda` work, child replication/gather traffic, output traffic, launch cost, and available
   memory.  First allow multiple GPUs to split even a small M range when each output has enormous
   lambda work.  When M parallelism is insufficient, add a true lambda-box multi-GPU plan: each
   device computes lower/upper log-sum-exp partials for a disjoint lambda range and the executor
   performs a numerically stable final merge.  Keep the existing lambda tiling as the bounded-memory
   primitive; do not describe it as lambda-parallel until this cross-device reduction exists.
7. **Benchmark the encouraging 8,000-sequence throughput.** Run same-oracle exhaustive PACK\* and
   controlled SPACK\* ablations (`deterministic-only`, `+sample evidence`, `+WMB-IS`) with matched
   seeds and budgets.  Report wall time, PACK\*/CCD calls, internal subtrees pruned before leaves,
   and certificate provenance.  Treat `5.42 h` as a useful baseline, not yet a speedup claim.
8. **Push WMB-IS certificates below deterministic WMB only after the frontier can reach them.**
   The local cap is now O(10-100) rather than O(1e9), but on `12007350` the certified WMB-IS upper
   (`152.3940`) was still looser than deterministic WMB (`139.3654`).  Next targets are
   residual/tail-source bounds and better WMB position ordering so the certificate can add pruning
   beyond the deterministic tightener.
9. **Add position ordering and retain cost instrumentation.** The default next-unassigned-position
   branching is sound but weak.  Prefer positions that raise the incumbent early, reduce subtree
   upper gaps, or reduce downstream PACK\*/DP cost.  Keep cache hit/miss, duplicate state requests,
   decomposition setup, CCD time, and oracle wall time in every scaling comparison.

#### BBK\* vs SPACK\* certificates

The key distinction is what each sequence-level driver can certify at return time.

**BBK\*/MARK\* heap certificate.** BBK\* uses admissible upper bounds as heap keys.  When a node
`p` is popped from the max heap, the heap order certifies only:

```text
for every remaining node r: upper(r) <= current upper(p)
```

Together with admissibility, this implies:

```text
true K*(r) <= upper(r) <= current upper(p)
```

This is a valid upper-bound search invariant, but it is not by itself a certificate that:

```text
true K*(r) <= true K*(p)
```

or that the popped sequence is already gap-free by true K\*.  That stronger claim would require
an additional lower-vs-remaining-upper check:

```text
max_remaining_upper <= lower(p)
```

or an equivalent proof that `current upper(p)` has collapsed to the true value.  An SS-ε
approximation for `p` gives a bracket for `p`; it does not by itself compare `p`'s lower bound
against all remaining nodes' upper bounds.  Thus BBK\*/MARK\* can be very fast when used as an
upper-bound-ordered proposal method, but its heap order should not be read as a strict true-K\*
ordering certificate unless this extra check is present.

**SPACK\* top-k certificate.** SPACK\* is designed to stop only when every remaining subtree is
certified unable to beat the kth incumbent:

```text
for every remaining node r:
    true K*(r) <= upper(r) <= kth_incumbent_lower <= true K*(kth incumbent)
```

For top-k, the cutoff is the kth best verified lower bound, not the upper bound of the most
recently popped node.  A **deterministically certified** `upper(r) <= cutoff` is the direct-prune
special case; a randomized PACK\* bracket needs a globally accounted PAC lease or sample-level
e-evidence for the same conclusion.  This stronger certificate is what SPACK\* is intended to use
to justify "all remaining sequences are worse than the returned top-k", but the current ordinary
PACK\* direct-prune path still needs the provenance/accounting fix described above.

The current implementation therefore has two separate improvement tracks:

**Algorithm scheduling**

- **Raise the cutoff early — implemented for the current runner.** Cheap WT/single-mutant seeds
  plus early stopping produced a finite top-2 cutoff after two seeds in the medium-scale probe.
  BBK\*/MARK\* remains an optional proposal source rather than a prerequisite.
- **Do not make BBK\* the only incumbent source — implemented by default.** The earlier logs showed
  BBK\*/MARK\* proposal refinement could dominate wall time before SPACK\* started.  The current
  runner instead uses a reproducible cap of cheap local candidates; soundness must come from the
  later certified SPACK\* decisions, not from the proposal source.
- **Make lazy and budgeted bounds real.** Children should first receive cheap optimistic bounds,
  and expensive PACK\* refinement should be spent only on frontier nodes that can still beat the
  incumbent cutoff.  The basic budget plumbing is now present: `BudgetSchedule` controls the
  per-call PACK\* estimator sample cap.  The medium run shows the next policy requirement: if the
  computed upper is still `+inf`, branch immediately rather than re-queueing the node; also avoid
  the complex PACK\* call entirely when a cheap denominator witness is already `-inf`.
- **Improve branching policy.** Splitting the next unassigned position is sound but not always
  efficient.  A better policy should prefer positions that are expected to reduce subtree upper
  bounds the most, shrink the largest upper-gap, or reduce branchwidth/DP cost.
- **Use BBK\* as a proposal generator, not a correctness competitor.** BBK\*/MARK\* is useful
  because it can quickly suggest promising sequences, but SPACK\* should remain responsible for
  the final strict certificate.

**Implementation work**

- **Wire `BudgetedSequenceSubtreeBounder` / `BudgetSchedule` into `SpackStarPackStar`.** Done for
  the current PACK\* estimator path: the node budget caps the estimator sample budget.  Benchmark
  work remains to choose a good default initial budget and growth factor.
- **Cache safe subtree bounds.** Ordinary `bound()` results can be memoized for repeated subtree
  requests.  Evidence-producing or verifier calls such as `boundWithResult()` should stay fresh
  unless the cache entry preserves the required sample trace and certificate semantics.  Current
  code caches only ordinary state bounds and keeps `boundWithResult()` fresh.
- **Measure and tighten WMB-IS certificates.** The fixed-range PAC path and local weight cap are
  implemented.  Real-state validation reduced the cap gap from `~1.24e9` to `~28` log units, but
  the certified WMB-IS upper still has not beaten the deterministic WMB upper on the tested
  subtree.  The 8,000-sequence run now supplies overall pruning counts but did not invoke WMB
  fallback; matched ablations are still needed to measure WMB-specific sequence pruning before
  tightening the residual/tail certificate further.
- **Reuse state-level work.** Protein, ligand, and other denominator-side subproblems repeat across
  many related sequence nodes.  Memoizing compatible state bounds can avoid repeated CCD/PACK\*
  work.  Current code reuses exact compatible state-sequence bounds with budget-aware cache keys.
- **Parallelize independent bounds.** Complex/protein/ligand state bounds, sibling-node bounds,
  and independent pfunc estimators are better parallelization targets than small Java-level
  micro-optimizations.  Current code parallelizes independent fully assigned state bounds when
  configured; sibling-node parallelism remains future work.
- **Log certificate progress directly.** Each iteration should report the quantities that decide
  pruning:

```text
max_remaining_upper
topN_lower_cutoff
gap = max_remaining_upper - topN_lower_cutoff
node bound wall time
PACK*/CCD call counts
```

This separates "we do not have a strong incumbent yet" from "the remaining upper bounds are too
loose" and from "the oracle implementation is spending too much per bound."

**Sample-evidence multi-round status.** The e-process bookkeeping is implemented and supports
multiple `update()` calls for the same sequence node.  The driver has an evidence hook and PACK\*
can provide fresh complex-state sample evidence.  The driver now supports an active sampling
policy:

```text
while node is not pruned and sample-evidence budget remains:
    draw a fresh PACK* evidence batch
    update the node e-process
```

The default remains one evidence batch per popped node, preserving the previous behavior.  Setting
`spackstar.run.maxEvidenceBatchesPerNode > 1` lets sample evidence accumulate immediately while the
node is still competitive.  In the observed pre-loop 2RL0 run, `sampleEvidence=0` and `rebounded=0`,
so the actual action was an ordinary `upper <= cutoff` bound comparison, not multi-round
sample-evidence pruning.  Whether that upper was deterministic or PAC must now be retained as
explicit provenance rather than inferred from the comparison operation.

In the current 2RL0 comparison runs, the root PACK\* bound is not the dominant cost: the observed
root complex call was about 92.6 seconds out of roughly 2400 seconds total.  Most time was spent
across many PACK\*/CCD bound calls.  That points to budgeted refinement, caching/reuse, and
parallel state evaluation as the highest-value implementation fixes.

## 6. WMB, WMB-IS, and mean-field

There are now two separate WMB uses:

- **Deterministic WMB/MF tightener:** `OspreyWmbMfStateBounder` supplies a deterministic bracket
  that can be intersected with a PACK\* state bracket (`[max(lower), min(upper)]`).  The intersection
  retains the PACK\* bracket's PAC provenance; it is not automatically zero-PAC-cost.  The finite
  crossings observed in the medium run must be explained and handled fail-closed before this
  combination is used for strict pruning.
- **WMB-IS state fallback:** `WmbImportanceSamplingStateBounder` estimates a sequence subtree by
  sampling from `q_WMB(c)` and scoring unique conformations with CCD:

```text
Z(subtree) = E_{q_WMB}[ exp(-E_CCD(c)/RT) / q_WMB(c) ]
```

WMB-IS now reports two different intervals.  The empirical interval over observed IS weights is
still diagnostic only.  The pruning-safe path is the fixed-range certificate based on a
deterministic local weight cap, intersected with deterministic WMB/MF when available.  Real 2RL0
root-subtree validation remains positive: WMB-IS avoided the branch-DP table cliff, produced
stable `cv≈0.35`, `ESS/n≈0.89` estimates on a 41-billion-conformation sequence subtree, and on
the `1,465,989,840`-conformation `2 full + 5 WT-flex` subtree tightened the certified weight-cap
gap from `~1.24e9` to `27.9969` log units.

## 7. Correctness summary

- Per-node subtree K\* upper-bound algebra is **admissible conditional on valid state bounds**
  (`SubtreeSafeKStarBounder`).  The conditioning matters: bound provenance and coverage accounting
  cannot be discarded when the bracket reaches the sequence driver.
- Ordinary PACK\* state brackets are per-call **PAC** objects, not deterministic bounds.  The current
  direct `upper <= cutoff` path does not yet lease global delta for those calls, and the 8,000-sequence
  run exposed `342` finite PACK\*/deterministic-WMB interval crossings.  Until the cause is resolved
  and direct PAC pruning is globally accounted, the completed medium run is a throughput result—not
  a proof of the advertised family-wise top-k error bound.
- PACK\*'s upper certificate also requires the coverage/residual tail `τ` to be a *true* bound;
  otherwise a rare good amino acid in a superset could be missed and a good subtree wrongly pruned.
- WMB-IS state estimates now have a **finite-sample fixed-range PAC certificate path**.  The
  implementation is fail-closed and consumes global delta leases; whether the deterministic weight
  cap is tight enough for useful pruning remains an empirical question.
- For pruning decisions based solely on valid fresh e-evidence, **e-value + global alpha spending**
  gives `P(ever prune a true top-k sequence) <= alpha_total` over adaptive node e-processes: each
  sequence-node process receives a unique geometric alpha lease and Ville's inequality applies under
  optional stopping.  This statement does not retroactively cover ordinary unleased PACK\* brackets.
- A truly deterministic MICA/WMB certificate is a degenerate e-value (`e = inf` when violated).
  Randomized PACK\* brackets need their PAC lease or sample-level e-evidence; they must not be
  relabeled deterministic merely because the final comparison is `upper <= cutoff`.

## 8. De-risk first

- Three-position end-to-end feasibility is now established: the 8,000-sequence 2RL0 run completed
  in 5.42 h.  Next obtain its matched exhaustive-PACK\* baseline and rerun certificate-provenance
  ablations before changing the performance claim.
- Before another 12 h four-position job, require a short probe to show that finite complex bounds no
  longer become `K* upper=+inf` through a `-inf` denominator, and that a bound which remains infinite
  branches immediately instead of cycling through the same-depth frontier.
- Only after those gates pass, test wide-Z designs (4u3s/4wyu) that motivate PACK\* over MARK\*.
- The overlap probe (`KL/chi^2`) on existing 2q1e/2rfe samples still decides whether Stage 3 is worth
  doing, but Stage 3 is lower priority than coverage and infinite-frontier fixes.

## 9. Open questions

- True-tail `τ` for the PACK\* upper certificate (soundness of the upper bound).
  *Partially addressed in code:* `PackStarOracleDiagnostics.worstCaseClipExcessUpperBound`
  (worst-case Hoeffding branch) + the clip-bias-shifted `evidenceVersusIncumbent`
  (`m' = m − B`) make the sample prune sound **given** a true residual-energy tail bound
  `residualBoundKcal`.  The low-level provider can disable sample pruning when no finite bound is
  supplied, but the current assembly/runners pass a configured finite assumption (`0.5` or `1.0`
  kcal/mol by default).  The empirical-Bernstein clip-probability branch and
  `withResidualTailBound(...)` assembly are implemented.  Remaining: default fail-closed behavior
  plus a real `residualBoundKcal` source from the PACK\* backend.
- Global spending for ordinary PACK\* state brackets: currently missing from the direct bound-prune
  path.  Decide whether to allocate per-call leases, build node/state confidence sequences, or use
  PACK\* brackets only for ordering while sample-level e-processes carry the pruning decision.
- Explain finite PACK\*/deterministic-WMB interval crossings and make the combiner fail closed until
  both bounders are confirmed to target the same energy model with compatible certificate semantics.
- Existing geometric-stream sample-cost/tightness may still need tuning for large sequence spaces.
- Whether the ratio estimand (Stage 3) pays off given the real overlap between adjacent sequences.
- Reusing the branch decomposition across parent→child (condition one position) inside the PACK\* backend.
