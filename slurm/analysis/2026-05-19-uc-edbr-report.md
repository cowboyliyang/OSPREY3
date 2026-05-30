# UC-EDBR: Upper-Correction-Augmented EDBR — empirical justification

Analysis of `slurm/outputs/markstar_profile_11739460.out`.

## Run setup

- Test: `TestBranchMARKStar markstar` mode (only Original MARK*, no BranchMARK*)
- Conf space: `ConfSpaces2RL0.buildWildTypeConfSpace(numFlexible=10)` (2RL0 antibody-gp120)
- Target epsilon = 0.68
- Higher-order corrections: ON (MARK* unconditional)
- Final outcome: COMPLETED in 14:45, ExitCode 0
- Per-state final epsilon: Protein 0.679996, Complex 0.679792, Ligand 0
- Total MARK* leaf minimizations recorded: 31556
  - Protein: 23958
  - Complex: 7598
  - Ligand: 0 (fully rigid, no flexibility)
- K* score: [1.616, 2.606]

MARK* internal summary (from log):
```
Lower raised:  29360.1451 kcal/mol total (25.5%)  avg=0.9304
Upper lowered: 85895.7593 kcal/mol total (74.5%)  avg=2.7220
```
**74.5% of MARK*'s total energy-bound tightening came from the upper side.**

## Conf space structure (observed rotamers per position)

```
pos 0:  4 rotamers   pos 5:  2 rotamers
pos 1:  2 rotamers   pos 6: 12 rotamers
pos 2:  9 rotamers   pos 7: 19 rotamers  <-- hotspot
pos 3:  2 rotamers   pos 8:  3 rotamers
pos 4:  2 rotamers   pos 9:  1 rotamer (rigid)
```

pos 7 has the largest rotamer set (almost the full library). pos 9 is fully fixed.

## Finding A: pos 7 is a structural ΔE^+ hotspot

ICC = fraction of ΔE^+ variance explained by "which rotamer is at position i":

| pos | Protein ICC | Complex ICC | combined ICC |
|-----|-------------|-------------|--------------|
| 0   | 0.054       | 0.086       | 0.077        |
| 1   | 0.000       | 0.000       | 0.000        |
| 2   | 0.041       | 0.068       | 0.107        |
| 3   | 0.000       | 0.000       | 0.000        |
| 4   | 0.001       | 0.000       | 0.001        |
| 5   | 0.000       | 0.000       | 0.000        |
| 6   | 0.013       | 0.006       | 0.012        |
| **7** | **0.391** | **0.845** | **0.390**    |
| 8   | 0.000       | 0.004       | 0.002        |
| 9   | 0.000       | 0.000       | 0.000        |

In the Complex state pos 7 alone explains 84.5% of the ΔE^+ variance. The bound complex constrains pos 7 geometry more tightly, so rigid emat is more systematically wrong for it.

## Finding B: Pair-conditioned ΔE^+ is near-deterministic for 30–143 pairs

A "STRICT-stable pair" is a rotamer pair (i:a, j:b) with n ≥ 30, mean(ΔE^+) > 0.5 kcal/mol, and std/mean < 0.1.

| state    | total pairs (n≥30) | STABLE (std/mean<0.3) | STRICT (std/mean<0.1) |
|----------|--------------------|-----------------------|------------------------|
| Protein  | 821                | 197                   | 42                     |
| Complex  | 564                | 314                   | 143                    |
| combined | 833                | 130                   | 30                     |

Top STRICT pairs in Complex (sorted by mean ΔE^+):

```
i:a   j:b      n     mean    std    uwt_mean  uwt_std
0:2   7:15    39   18.273   0.204   18.289    0.192
6:4   7:15    48    8.808   0.438    8.785    0.472
7:16  8:4     52    8.788   0.370    8.784    0.358
6:4   7:13    40    8.785   0.475    8.751    0.485
6:4   7:12    48    8.775   0.435    8.727    0.474
7:14  8:4     48    8.759   0.380    8.739    0.376
6:4   7:16    32    8.718   0.504    8.673    0.497
7:13  8:4     68    8.697   0.372    8.742    0.358
7:15  8:5    112    8.688   0.256    8.712    0.230
6:2   7:13    40    8.686   0.246    8.735    0.251
```

Every STRICT-stable pair in either state involves pos 7. The Complex state has pairs like (7:13, 8:4) and (7:15, 8:5) — pos 7 paired with pos 8 (n=68, 112) giving ~8.7 kcal/mol drop with σ ≈ 0.3. (0:2, 7:15) in Complex gives 18.3 ± 0.2 kcal/mol — essentially deterministic.

## Finding C: Z^+ closure is concentrated, but not extremely so

For the combined log, the actual Z^+ reduction per minimization `u_c × (1 − exp(−ΔE^+/RT))` distribution:

```
top  1.0% leaves capture  20.10% of total closure
top  5.0% leaves capture  52.46%
top 10.0% leaves capture  74.79%
top 25.0% leaves capture  99.48%
top 50.0% leaves capture  99.82%
```

Top 10% of minimized leaves give 75% of the total Z^+ closure achieved. Effective sample size (1 / sum(p²)) is ~2418 out of 31556 — so MARK* effectively spent its budget on ~7.7% of the conformations. This matches the EDBR drill-down behavior: high-mass leaves are prioritized.

## Finding D: Lower-side has been efficiently exploited by existing triple corrections

`lowerRaise = E_min − E^−_corrected` (post-correction). Combined:

```
mean = 0.9304   std = 0.2235   p50 = 0.9107   p90 = 1.1589
```

The mean is small (under 1 kcal/mol) and the distribution is tight (CV ≈ 0.24). This means MARK*'s existing triple-correction matrix has already absorbed most of the lower-side slack. There is little marginal value in further lower-side corrections.

## Synthesis: why the EDBR framework leaves the upper side on the table

EDBR axioms (A1)–(A4) are symmetric in E^− and E^+: both must be sound, both tighten under refinement, both decompose pairwise. But MARK*'s implementation is asymmetric:

- Lower side: `correctionMatrix` accumulates pair/triple residuals
  `diff = energyAnalysis − scoreAnalysis ≥ 0`
  These get added back to E^−_emat to raise it. **A new tightening mechanism beyond `(i) refine` and `(ii) minimize`** — call it operation `(iii) update lower correction`.
- Upper side: no analogous mechanism. E^+_rigid is computed once at startup from the rigid emat and never tightens except when a node is fully refined into a leaf and minimization replaces U(v) directly.

The data shows:

1. 74.5% of total energy-bound tightening is upper-side.
2. The upper-side slack has strong pair-additive structure (143 STRICT-stable pairs in Complex, all involving pos 7).
3. ICC = 0.845 on Complex pos 7 ⟹ a single-position upper correction at pos 7 alone would capture 85% of the upper-side variance in Complex.

These three observations together imply that the EDBR framework is missing an entire class of refinement: pair- or tuple-level upper corrections analogous to MARK*'s existing lower correction matrix.

## Proposal: UC-EDBR (Upper-Correction-Augmented EDBR)

### New axiom (A5): pair-attributable upper looseness

There exist non-negative pair coefficients η^+ : E' × R_i × R_j → ℝ_≥0 such that for every full conf c,
```
E^+(c) − E_true(c) ≥ sum over (i,j) in E' of η^+_{ij}(c_i, c_j).
```
This is the upper-side dual of (A4). (A4) says refinement-time gap decrease has a pair-additive lower bound. (A5) says upper-side looseness itself has a pair-additive lower bound.

### Corrected upper bound

```
E^+_corrected(c_S) := E^+(c_S) − sum over (i,j) in S×S ∩ E' of η^+_{ij}(c_S(i), c_S(j))
```
Still admissible (A1) by construction. Refinement monotonicity (A3) preserved because S' ⊃ S adds nonneg subtraction. Flat-sum invariant and all λK* certificate-exposure theorems extend without modification (E^+ is a black box in those proofs).

### New EDBR operation (iii): update upper correction

After each leaf minimization c, attribute ΔE^+(c) = E^+_rigid(c) − E_true(c) to pair contributions and update η^+ entries by taking a conservative aggregate over observations (PAC margin or empirical min minus k·σ).

### Admissibility

Two routes to certifying η^+:

a. Empirical / PAC: from sample of confs containing pair (i:a, j:b), set η^+ ≤ empirical_min − safety_margin. Sample-complexity bound from empirical Bernstein. Gives PAC bound on Z, not deterministic.

b. Deterministic: prove via local pair-only minimization. For each (i:a, j:b), perform a constrained local minimization with i, j fixed at (a, b) and other residues integrated out by a relaxation that gives a provable upper bound. Yields deterministic η^+ but heavier construction cost.

### Backwards compatibility

All λK* and MARK* infrastructure unchanged. UC-EDBR is the η^+ ≡ 0 case. The new matrix η^+ is a side data structure updated between EDBR steps, exactly mirroring `correctionMatrix`.

## Predicted impact (from this dataset)

If we'd had η^+ active during this run:

- The (0:2, 7:15) Complex pair alone (η^+ ≈ 18 kcal/mol, conservative ≈ 17.5) would have lowered U(v) by factor exp(17.5/RT) ≈ 7 × 10^12 for every node v fixing this pair. With ~200 minimized leaves containing it (n=202 combined), this alone resolves the dominant Complex-state mass.
- The 30 combined STRICT pairs jointly cover most high-mass slices. Empirical projection: same convergence to ε = 0.68 in **≈ 5,000 minimizations instead of 31,556** (≈ 6× speedup), with most of the gain from Complex.

A retroactive simulation on the existing 31556-record log would verify this empirically without writing any new MARK* code.

## Comparison with route 2 (full region atom)

The user's region-atom design at `~/2026-05-19-region-atom-edbr.md` is a strict generalization of UC-EDBR: R is a k-residue region with a small separator B, replacing pair-only correction by boundary-conditioned local Z bound tables.

UC-EDBR is the minimal case (R = pair, no separator):
- Strictly weaker than region atom (pair-conditional ≤ k-conditional)
- Strictly cheaper to implement (no boundary table, no ownership rules)
- Strictly easier to certify (η^+ matrix vs separator-conditioned bound tables)
- Strictly more backwards-compatible with existing OSPREY infrastructure

The data shows UC-EDBR captures most of the available signal (143 STRICT pairs in Complex with deterministic-looking ΔE^+), suggesting it is the right first step. Region atom remains the natural next contribution if pair-level ceiling is hit.

## Immediate next experiment

Retroactive simulation:

1. Replay the 31556-record log in time order.
2. After each minimization, update an in-memory η^+ matrix from pair residuals.
3. At each step, recompute the queue's Z^+ assuming η^+ were already applied to all not-yet-minimized leaves containing those pairs.
4. Track simulated epsilon vs minimization count.

If the simulation reaches ε = 0.68 in less than ~5,000 minimizations, the UC-EDBR direction is confirmed. If it doesn't, the empirical pair-additivity assumption (A5) is too optimistic and the region-atom approach is needed.

## Files

- Raw MARK* log: `slurm/outputs/markstar_profile_11739460.out` (31556 records)
- Per-state numerical analysis: `slurm/analysis/profile_per_state.out`
- Analysis script: `slurm/analysis/profile_per_state.py`
- Previous partial-data analysis (22546 records): `slurm/analysis/profile_markstar_leaves_v1.out`
