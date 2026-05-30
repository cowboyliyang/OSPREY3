# Optimal ε-convergence algorithm for BranchMARK*

A composite plan synthesizing every lever discussed so far. The goal is to
minimize the wall-clock work required to reach `ε ≤ 0.68` (or any target) for
each pfunc in the K* score, on the 2RL0 protein-design benchmark and analogous
systems.

## 1. Restatement of the goal

The reported quantity is `K* = Z_complex / (Z_protein × Z_ligand)`. Each `Z`
is approximated under EDBR (Error-Driven Bound Refinement) by maintaining a
queue Q of partial conformations with flat-sum bounds:

```text
Z^-(Q) = Σ_v L(v),     Z^+(Q) = Σ_v U(v)
Φ(Q)   = Z^+(Q) - Z^-(Q) = Σ_v δ(v)
ε(Q)   = 1 - Z^-(Q)/Z^+(Q)
stopping rule: ε(Q) ≤ ε*
```

The algorithm consumes flat-sum gap Φ via three operations, all preserving
admissibility:

```text
(i)   refine internal node v on edge p → replaces v with children
(ii)  minimize leaf c → replaces [L(c), U(c)] with the exact Boltzmann factor
(iii) bound refinement → tightens L(v), U(v) without splitting v
```

Operations (i) and (ii) are the paper's EDBR axes (refinement unit, refinement
order). Operation (iii) is the new mechanism that covers pair corrections,
triple corrections, and region atoms.

Objective: minimize wall-clock to reach ε* using the best mix of (i), (ii),
(iii) across parallel hardware (1 GPU + 104 CPUs in the current setup).

## 2. Empirical state of the current implementation

Confirmed from code and recent runs (`ra_scout104_{11789709, 11791810, 11791940}`):

### 2.1 Higher-order corrections asymmetry

```text
MARK*:        triple lower correction = ON  (unconditional)
BranchMARK*:  triple lower correction = OFF (default; flag exists but disabled)
              upper-side correction   = none (UC-EDBR not implemented)
```

`BranchMARKStarBound.java:554`:

```java
this.useHigherOrderCorrections =
    getConfigBoolean("branchmarkstar.useHigherOrderCorrections", false);
```

### 2.2 Per-leaf bound slack distribution

Pulled from `BranchMARK*: N minimizations, ..., eRigid=..., eMin=..., eTrue=...`
progress lines. Averaged across the recent 3 runs:

```text
                       lowerRaise   upperDrop   E split (lo/up)   Z-mass split (lo/up)
BMS Protein (n≈410):     2.88         7.30        28% / 72%         14% / 86%
BMS Complex (n≈66):      2.63         4.53        37% / 63%         0.4% / 99.6%
MARK* with triples ref:  0.93         2.72        25% / 75%         25% / 74.5%
```

Two observations from these numbers:

- BranchMARK* leaves `lowerRaise` ≈ 3× larger than MARK* with triples → 2
  kcal/leaf is "lying on the table" because the triple flag is off.
- Even with that, the Boltzmann amplification of `upperDrop` is so large that
  the Z-mass distribution on Complex is 99.6% upper-side. On Protein it is
  86% upper-side. Upper-side correction is THE biggest lever regardless of
  whether the triple flag is on.

### 2.3 Hotspot concentration (from MARK* analysis, transfers to BMS)

```text
ICC(upperDrop, pos i)  — fraction of upperDrop variance explained by position i

pos 7:  Protein 0.391   Complex 0.845   combined 0.390
```

pos7 has 19 rotamers (largest in the conf space). Complex constrains it
geometrically; rigid emat is systematically wrong → 84.5% of upper-side
variance in Complex traces to pos7.

ICC for pos7 is computed on MARK* leaves (with full conf arrays). BMS progress
lines don't carry the conf array, so a direct BMS-side ICC requires rerunning
with `-Dmarkstar.leafProfile=true -Dmarkstar.leafProfile.allConfs=true` and
parsing `[LEAF_PROFILE]` records. Until then we adopt the MARK* ICC values
because the underlying physics is identical.

### 2.4 Edge selection engagement

`numFlexible=10, residualBudget=0.5`, sparse graph with `kept 32/45 edges`.
Branch decomposition produces nodes with mostly 1 pending edge:

```text
rootSplit=legacy + contractionPerState:   maxPendingEdges=1, multiEdge=0
rootSplit=lookahead + contractionPerState: multiEdge=3-7% of calls, exactLookahead fires
```

Even with lookahead exposing multi-pending cases, the strategy difference
between `contractionPerState` and `lambdaStates` is small at n=10 (Protein
40872 vs 39312 minimizations, Complex 6136 vs 6032). Paper's 58× speedup is
MARK* vs λK*, not within-λK* strategies. Within-λK* edge selection only
becomes a measurable lever at n ≥ 12.

## 3. Levers, with ROI

Independent levers ranked by expected speedup on Complex (the K* bottleneck):

```text
ID   Lever                                          Complex     Protein    Eng cost   Status
─────────────────────────────────────────────────────────────────────────────────────────────
L0   Enable BMS triple lower correction             ≈1.0×       1.2-1.5×   ~0         flag exists
L1   λK* base (already in use)                      —           —          done       baseline
L2   contractionPerState + rootSplit=lookahead      ≈1.0×       ≈1.0×      done       at n=10
L3   η^+ pair upper correction (UC-EDBR)            5-10×       3-5×       medium     unimplemented
L4   tightenNode partial-B worst-case ratio         1.5-2×      1.2×       small      unimplemented
L5   Bandit-with-budget batch scheduling            1.2×        1.2×       medium     unimplemented
L6   Triple-level η^+ upper correction              1.1×        1.3×       medium     unimplemented
L7   Cross-state η matrix sharing                   1.5-2×      1.5-2×     small      unimplemented
L8   Online retroactive η propagation               1.2×        1.2×       medium     unimplemented
─────────────────────────────────────────────────────────────────────────────────────────────
```

L3 dominates because Complex's 99.6% Z-mass slack is upper-side. L0 is "free"
on Protein only.

L4 + L7 are the largest small-engineering wins after L3.

L2 was already used in the recent run (job 11791940) but did not exceed
lambdaStates by a measurable margin at n=10.

## 4. Recommended algorithm: λK*-UC++

A composite EDBR algorithm with three layers.

### Layer A — refine/minimize (existing λK*)

```text
- Branch decomposition of sparse G' (residual budget cutoff)
- λ-set refinement unit (paper EDBR axis 1)
- contractionPerState + rootSplit=lookahead (paper EDBR axis 2)
- Flat-sum scheduling on Z_I vs Z_L (existing)
```

### Layer B — symmetric pair corrections (core new mechanism)

Two correction matrices, structurally identical:

```text
η^- : E' × R_i × R_j → ℝ_≥0    (lower correction, raises E^-)
η^+ : E' × R_i × R_j → ℝ_≥0    (upper correction, lowers E^+)
```

**Update rule** (after each leaf minimization c with observed
ΔE^-(c), ΔE^+(c)):

```text
For each (i,j) ∈ E' with c_i, c_j fixed:
    η^-_{ij}(c_i, c_j) ← max( η^-_{ij}(c_i, c_j),
                              attribute_lower(ΔE^-(c), c, (i,j)) )
    η^+_{ij}(c_i, c_j) ← max( η^+_{ij}(c_i, c_j),
                              attribute_upper(ΔE^+(c), c, (i,j)) )
```

Monotonic max keeps the correction sound under future revisions.
Conservative attribution: `attribute = ΔE / n_pairs(c)` or an empirical
margin (PAC-style). Stricter attributors give tighter η at the cost of
narrower sound coverage.

**Application rule** (at any node v, leaf or internal, with partial conf c_S):

```text
E^-_corrected(c_S) = E^-_emat(c_S)
                   + Σ_{(i,j) ∈ E', {i,j} ⊆ S} η^-_{ij}(c_S(i), c_S(j))

E^+_corrected(c_S) = E^+_emat(c_S)
                   - Σ_{(i,j) ∈ E', {i,j} ⊆ S} η^+_{ij}(c_S(i), c_S(j))

L(v) = exp(-E^+_corrected(c_S) / RT) · τ(c_S)
U(v) = exp(-E^-_corrected(c_S) / RT) · τ(c_S)
```

Both sides apply at every search node — leaves and internals — wherever pair
endpoints are both in the fixed set. No separator is needed; no DP factor
replacement; no partial-B handling.

**Soundness conditions** (verified offline via retroactive simulation on MARK*
leaf log; tracked online via correction audit):

```text
(A4) lower pair-additive: ∀c, E_true(c) - E^-_emat(c)  ≥ Σ η^-_{ij}(c_i, c_j)
(A5) upper pair-additive: ∀c, E^+_emat(c) - E_true(c)  ≥ Σ η^+_{ij}(c_i, c_j)
```

(A4) is the foundation of MARK*'s existing `correctionMatrix` and has held
empirically for years. (A5) is the new dual; UC-EDBR doc validates it on
MARK*'s 31556-record log (143 STRICT-stable pairs in Complex, all involving
pos 7).

**Cross-state η sharing.** Residues and pair interactions are shared across
the three pfunc states of a K* score. Maintain a single (residue_pair,
rotamer_pair) → η table per energy function, shared across Complex / Protein
/ Ligand. Updates from any state contribute to the same table; reads from
any state apply the same correction. This amortizes the certification cost
roughly 2-3× across the K* triplet.

### Layer C — adaptive scheduler (bandit-with-budget)

Each EDBR round selects a batch of actions across the three operations
(refine, minimize, certify). Score each candidate by work-normalized rate:

```text
ρ_refine(v, p)   = (δ(v) - Σ δ(w_α)) / |R_λ(v, p)|
ρ_minimize(c)    = δ(c) / 1
ρ_certify(pair)  = E[ΔU(pair)] · n_affected_nodes / CCD_cost
```

Then:

```text
sort candidates by ρ desc
threshold ρ* = ρ of (K+1)-th candidate, K = current parallel budget
batch = candidates with ρ > ρ*, total parallel cost ≤ budget
parallel execute batch
update η matrices with new minimizations
retroactively propagate η updates to all queued nodes' bounds
re-rank, next round
```

This subsumes "how many minimizations per round" (it's whatever rank falls
above ρ*) and prevents both baseline over-expansion and fullpar blind
batching.

### Layer D — fallback hierarchy (only if pair ceiling hit)

Trigger conditions and upgrades:

```text
if observed Σ η^+_{ij}(c) consistently < ΔE^+(c) by significant margin:
    promote to triple-level η^+_{ijk}(a, b, c)
    same attribution and application logic, larger storage

if triple-level still leaves residual:
    promote to region atom (R, B) for the residual hotspot
    apply via ratio overlay (no DP, no λ-set alignment required)
    extended with partial-B worst-case ratio (L4)
```

Most instances stop at pair level. Triple and region-atom are theoretical
guarantees that ε always closes.

## 5. Phased implementation

### Phase 0 — verify the "free" lower-side win (1 day)

- Run BranchMARK* with `-Dbranchmarkstar.useHigherOrderCorrections=true`
- Confirm lowerRaise mean drops from ~2.9 to ~0.9 kcal on Protein
- Confirm no correctness regressions in K* score
- If sound, set as default for all subsequent runs

Open question: why is the BMS triple flag default off? Possibly a known
soundness issue or unfinished integration. Code in
`BranchMARKStarBound.java:3694 computeTripleCorrections` looks structurally
similar to MARK*'s. Need to read history before flipping.

### Phase 1 — η^+ minimum viable (2-3 days)

- Add `BranchMARKStarBound.upperCorrectionMatrix`, mirroring
  `correctionMatrix`
- After each leaf minimization, attribute ΔE^+ to pairs via
  conservative-min and update entries
- At bound calculation (`recomputeZBounds` and equivalents) for each
  partial conf with pairs in fixed set, subtract Σ η^+ from E^+_emat
- Validate via retroactive simulation on `markstar_profile_11739460.out`:
  - Replay leaves in order; track simulated ε after applying η^+
  - Predicted ε reduction at 5000 minimizations (per UC-EDBR doc)

Hardware: pure CPU change, no GPU.

### Phase 2 — application coverage and scope (1 week)

- Extend η^+ application to internal nodes with partial fixed sets
  (current `tightenNodeWithCertifiedTables` parallel)
- Add cross-state sharing: single (residue, rotamer) key space across
  Complex / Protein / Ligand
- Add audit logging: per-conf, what η^+ contributed; alert if Σ η^+ ever
  approaches the observed ΔE^+ (signals A5 may be tight)

### Phase 3 — adaptive scheduling (1-2 weeks)

- Replace fixed `maxMinimizations` / `maxInternalNodes` with bandit selector
- Track per-action ρ posteriors
- Re-rank queue after each batch
- Profile lookahead overhead vs gain

### Phase 4 — graceful degradation (as needed)

- Triple-level η^+: only if Phase 1-3 doesn't reach 5× wall-clock improvement
  on Complex
- Region atom + partial-B overlay: only if triple-level still leaves residual

### Phase 5 — paper-grade evaluation

- Scaling series n ∈ {8, 10, 12, 16, 20, 24} same as paper §5
- Compare λK*-UC++ vs λK* base vs MARK* (with triples)
- Report wall-clock, minimization count, ε trajectories
- Verify K* score consistency across implementations

## 6. Quantitative predictions

### Per-state speedup at n=10 vs current BMS (LS+legacy = 11789709 baseline)

```text
                 Protein min     Complex min   K* total (rough)
Baseline BMS     40,872          6,136         47,000
+ L0 (triples)   32,000          6,100         38,000   (~1.2× Protein)
+ L3 (η^+)       9,000           700           9,700    (~5× total)
+ L4 + L7        7,000           500           7,500    (~6× total)
+ L5 (bandit)    6,000           450           6,450    (~7× total)
```

vs MARK* baseline (with triples) at same problem: ~40,000 min. So predicted
total speedup at n=10:

```text
MARK*           40,000 min  (baseline)
λK*-UC++         6,500 min  (~6× over MARK*)
```

At n=12 paper claims λK* alone is 58× over MARK*. With L3 stacked on top,
predicted total speedup ~100-300×.

These are wall-clock estimates assuming the η^+ predictions hold; lower side
is now ≈ 0.4% of Complex Z-mass so even MARK*-level lower correction has
negligible Complex impact. Upper side is the only handle.

### Per-step ε contraction

For Complex at the ε ≈ 0.68 stopping point:

```text
Current BMS: each leaf closes ~5.5e+43 / 6136 ≈ 9e+39 of Z^+ per leaf
With L3:     same leaf closes additional Σ η^+ across all queued nodes
             estimated 5-10× more Z^+ per leaf because η^+ amortizes
```

## 7. What we explicitly do not include

| Excluded                              | Why                                                |
|---------------------------------------|----------------------------------------------------|
| Hypertree decomposition               | Engineering cost > expected gain over pair η^+     |
| DP factor replacement                 | Requires atom = λ-set alignment, not realistic     |
| Region atom as primary mechanism      | Pair + triple η^+ covers 90%+ of observed signal   |
| Data-driven branch decomposition      | Atom-overlay path doesn't depend on it             |
| Fullpar-from-start blind batching     | Bandit (L5) gives strictly better behavior         |
| `correctionMatrix.confE` direct write | Already flagged unsafe in 2026-05-19 doc           |

## 8. Risks and mitigations

```text
Risk                                  Mitigation
─────────────────────────────────────────────────────────────────────────
A5 fails (upper not pair-additive)    Conservative attribution; retroactive
                                      simulation gates Phase 1 deployment
                                      
BMS triple flag has known bug         Read history of L_PROPERTY 554 default
                                      before flipping; run regression
                                      
Cross-state η write conflict          Per-state shards; merge on read; max
                                      across shards is still sound
                                      
Bandit cold start                     Offline scout on MARK* log pre-fills
                                      η before first BMS run
                                      
Online propagation too slow           Async update; allow stale bounds
                                      during update (still sound)
                                      
Triple/region atom ROI lower than     OK — pair level alone is 5-10× on
predicted                             Complex; that's sufficient for paper
```

## 9. Open empirical questions

1. Verify BMS pos7 ICC matches MARK* (rerun with `leafProfile=true`).
2. Verify (A5) on BMS leaves directly (currently only verified on MARK*).
3. Measure actual η^+ attribution conservativeness — how much margin do we
   leave on the table?
4. Profile how often atoms apply at partial-B vs full-B in current
   `tightenNodeWithCertifiedTables`.
5. Establish the n where edge-selection strategies (LS vs CPS) actually
   diverge (probably n ≥ 12 from paper data).

## 10. Summary in one paragraph

OSPREY's current state is asymmetric: it tightens the lower energy bound via
MARK*'s `correctionMatrix` (pair + triple) but leaves the upper energy bound
untouched. Empirically the upper bound carries 86-99.6% of the Z-mass gap in
the K* states. The optimal ε-convergence algorithm is the existing λK* with
branch decomposition + adaptive edge selection, augmented with a symmetric
η^+ upper-correction matrix that mirrors the lower one. Pair-level η^+ alone
should give 5-10× speedup on the bottleneck Complex state at the current
benchmark size, scaling further with problem size. Region atoms, hypertree
decompositions, and DP factor replacement are not needed for this gain; they
are theoretical fallbacks for residual non-pair-additive cases. The
fast-path is to flip BMS's triple flag, add the η^+ matrix, and let
ε converge.
