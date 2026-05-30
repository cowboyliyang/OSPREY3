# PAC Partition Function Estimation via Rao-Blackwellized Importance Sampling

## Problem

Compute K* = Z_complex / (Z_protein × Z_ligand), where each
Z = Σ_c exp(-E_true(c) / kT).

- E_true(c): expensive CCD oracle (~seconds per call)
- E_min(c), E_rigid(c): cheap pairwise bounds from energy matrix (microseconds)
- E_min(c) ≤ E_true(c) ≤ E_rigid(c)
- Conformation space: product of discrete rotamer choices, ~394K for n=10

Current approach (BranchMARK*): A* tree search + leaf-by-leaf CCD minimization.
31,556 CCD calls for ε = 0.68 (deterministic). Bottleneck: 99.6% of the ε gap
comes from Z⁻ = Z_rigid being too small.

## Core insight

The problem is not summation (DP handles that exactly for pairwise models).
The problem is that E_min ≠ E_true — the pairwise energy model is approximate.

The correction g(c) = E_true(c) - E_min(c) is **nearly pair-decomposable**:
143 STRICT-stable pairs in Complex have std/mean < 0.1 (UC-EDBR report).
This means g(c) can be well-approximated by a pair-decomposable function f_pair(c),
with residual std ≈ 0.3 kcal/mol.

## Mathematical framework

### Step 1: Reformulation as mean estimation

Z = Σ_c exp(-E_true(c)/kT)
  = Σ_c exp(-(E_min(c) + g(c))/kT)
  = Σ_c exp(-E_min(c)/kT) × exp(-g(c)/kT)
  = Z_min × 𝔼_p[φ(c)]

where:
- Z_min = Σ_c exp(-E_min(c)/kT)  — exact via DP, free
- p(c) = exp(-E_min(c)/kT) / Z_min — known distribution, DP-sampleable
- φ(c) = exp(-g(c)/kT) ∈ (0, 1]   — requires CCD to evaluate

The problem reduces to estimating 𝔼_p[φ(c)].

NOTE: E_rigid does not appear anywhere. The entire estimation is from the E_min side.
This is why the "99.6% gap from E_rigid" problem dissolves — we bypass E_rigid entirely.

### Step 2: Rao-Blackwellization (variance reduction)

After CCD minimization of conformation c, each energy term at the optimized DOFs*
is directly evaluable:

    g(c) = Σ_i [oneBody(i, DOFs*_i) - oneBody_min(i, rᵢ)]
         + Σ_{(i,j)} [pair(i,j, DOFs*) - pair_min(i,j, rᵢ,rⱼ)]

This is an exact decomposition (not regression, not attribution). Each term is a
direct energy evaluation at DOFs*.

Define η per term from observed CCD results:

    η_onebody(i, rᵢ) = mean over observed c of [oneBody(i, DOFs*_i) - oneBody_min(i, rᵢ)]
    η_pair(i,j, rᵢ,rⱼ) = mean over observed c of [pair(i,j, DOFs*) - pair_min(i,j, rᵢ,rⱼ)]

Define f_pair(c) = Σ_i η_onebody(i, cᵢ) + Σ_{(i,j)} η_pair(i,j, cᵢ, cⱼ).

f_pair is pair-decomposable. Then:

    φ(c) = exp(-f_pair(c)/kT) × exp(-residual(c)/kT)
          = φ_pair(c) × ψ(c)

- φ_pair: pair-decomposable → modify emat_min entries (add η per term) → rerun DP
  → Z_corrected = Σ exp(-(E_min(c) + f_pair(c))/kT) computed EXACTLY
- ψ(c) = exp(-residual(c)/kT): residual std ≈ 0.3 kcal/mol → CV[ψ] ≈ 0.24

Without Rao-Blackwellization: CV[φ] ≈ 5-10 → need ~10K-50K samples → no advantage.
With Rao-Blackwellization: CV[ψ] ≈ 0.24 → need ~300-500 samples → 60× speedup.

### Step 3: PAC bound

Z = Z_corrected × 𝔼_p'[ψ(c)]

Estimate 𝔼[ψ] by sample mean ψ̄ from n CCD evaluations.
Bernstein inequality (ψ bounded, known variance):

    P(|ψ̄ - 𝔼[ψ]| > Δ) ≤ 2 exp(-nΔ² / (2Var[ψ] + 2Δ·range/3))

→ Z ∈ [Z_corrected × (ψ̄ - Δ), Z_corrected × (ψ̄ + Δ)] with confidence 1-δ.
→ ε_pac = 2Δ̃ / (ψ̄ + Δ̃) where Δ̃ accounts for per-pfunc union bound.

### Why Z⁻ improves dramatically

Current:    Z⁻ = Z_rigid ≈ Z_min × 5.7e-6
PAC:        Z⁻ = Z_corrected × (ψ̄ - Δ) ≈ Z_min × 0.007
Improvement: ~1200×

## Algorithm

### Input
- branchRigidEmat, branchMinimizingEmat (existing)
- rootedRoot, rootedRootEdge (existing DP tree)
- minimizingEcalc (existing CCD calculator)
- rcs (rotamer conformations)
- interactionGraph (existing)
- Target: ε_target, confidence δ (e.g., 0.05)

### Phase 0: DP baseline (existing, free)
- DP tables already built by initSearch()
- Z_min = exp(rootedRootEdge.logZUpper[0]) [already computed]
- Z_rigid = exp(rootedRootEdge.logZLower[0]) [already computed]

### Phase 1: Sample conformations from p(c) ∝ exp(-E_min/kT)
Ancestral sampling on the tree decomposition:
1. Start at root edge
2. For each edge (top-down): sample M-state proportional to exp(-fullEnergyMin)
3. Within each edge: sample lambda-state proportional to its contribution to logZUpper[mIdx]
4. Decode (lambda-state, M-state) → full rotamer assignment

Produce N conformations. N = 500 per pfunc is a good starting point.

Implementation: new method `sampleConformationsFromDP(int n)` in BranchMARKStarBound.
Returns List<int[]> of full conformation arrays.

### Phase 2: Parallel CCD minimization
For each sampled conformation c:
1. Build RCTuple + interaction list (same as current leaf minimization)
2. Submit to minimizingEcalc.calcEnergy() or calcEnergyAsync()
3. Collect E_true(c) AND the optimized DOFs (EnergyParametricMolecule)

All N conformations can run in parallel on 104 CPUs (~5 batches of 104).

Implementation: similar to the existing leaf minimization at line 3561,
but without the search tree overhead.

### Phase 3: Extract per-term corrections (η learning)
For each minimized conformation c with DOFs*:
- For each position i:
    actual_oneBody = evaluate oneBody energy of residue i at DOFs*_i
    η_onebody_samples[i][rᵢ].add(actual_oneBody - emat_min.getOneBody(i, rᵢ))
- For each pair (i,j) in interactionGraph:
    actual_pair = evaluate pair energy of (i,j) at (DOFs*_i, DOFs*_j)
    η_pair_samples[i][j][rᵢ][rⱼ].add(actual_pair - emat_min.getPairwise(i, rᵢ, j, rⱼ))

Then: η(term) = mean of samples for that term.

Implementation: need access to per-term energy evaluation at specific DOFs.
The EnergyParametricMolecule from CCD has the optimized coordinates.
Use the same forcefield evaluator to compute individual terms.

Key code entry points:
- minimizingEcalc.calcEnergy() returns EnergyParametricMolecule (epmol)
- epmol contains ParametricMolecule with optimized DOFs
- Need to evaluate individual energy terms at those coordinates
  (similar to how ConfAnalyzer.analyze() breaks down energy)

### Phase 4: Modify emat and recompute DP
Create a corrected energy matrix:
    emat_corrected.oneBody(i, rᵢ) = emat_min.oneBody(i, rᵢ) + η_onebody(i, rᵢ)
    emat_corrected.pair(i,j, rᵢ,rⱼ) = emat_min.pair(i,j, rᵢ,rⱼ) + η_pair(i,j, rᵢ,rⱼ)

Rebuild DP tables with emat_corrected:
    RootedTreeEdge.postOrderInitIncremental(rootedRoot,
        branchRigidEmat, emat_corrected, interactionGraph, RT);
    RootedTreeEdge.postOrderComputeFullDP(rootedRoot);

Read Z_corrected = exp(rootedRootEdge.logZUpper[0]) with corrected emat.

Implementation: create a copy of branchMinimizingEmat, add η to each entry,
pass to postOrderInitIncremental. ~seconds.

### Phase 5: Estimate residual and compute PAC bound
For each of the N minimized conformations:
    g(c) = E_true(c) - E_min(c)                 [known from Phase 2]
    f_pair(c) = computeFullConfPairwiseEnergy(c, emat_corrected)
              - computeFullConfPairwiseEnergy(c, branchMinimizingEmat)   [cheap]
    residual(c) = g(c) - f_pair(c)
    ψ(c) = exp(-residual(c) / kT)

Compute:
    ψ̄ = (1/N) Σ ψ(cᵢ)
    s² = (1/(N-1)) Σ (ψᵢ - ψ̄)²
    range_ψ = max(ψᵢ) - min(ψᵢ)

Bernstein bound for confidence δ_per_pfunc = δ/3:
    Solve for Δ: 2 exp(-N Δ² / (2s² + 2Δ·range_ψ/3)) = δ_per_pfunc

Result:
    Z_lower = Z_corrected × (ψ̄ - Δ)
    Z_upper = Z_corrected × (ψ̄ + Δ)
    ε_pac = 1 - Z_lower / Z_upper = 2Δ / (ψ̄ + Δ)

### Phase 6: K* score
    K*_lower = Z_lower_complex / (Z_upper_protein × Z_upper_ligand)
    K*_upper = Z_upper_complex / (Z_lower_protein × Z_lower_ligand)
    Report K* ∈ [K*_lower, K*_upper] with confidence 1-δ.

## Cost analysis

Per pfunc:
- Phase 1 (sampling): milliseconds (DP table lookups)
- Phase 2 (CCD): 500 calls × ~0.1s = ~50s on 104 CPUs (~5 batches)
- Phase 3 (η extraction): milliseconds (energy term evaluation)
- Phase 4 (DP rebuild): ~1-3 seconds
- Phase 5 (statistics): milliseconds

Total per pfunc: ~1 minute wall-clock.
Total for K* (2 non-trivial pfuncs): ~2 minutes.

Current: BranchMARK* takes ~10 minutes for ε = 0.68 with 31,556 CCD calls.

Predicted:
| N per pfunc |  ε_pac  | total CCD | wall-clock |
|-------------|---------|-----------|------------|
| 200         | ~0.38   | ~800      | ~1 min     |
| 500         | ~0.19   | ~1,000    | ~2 min     |
| 1000        | ~0.14   | ~2,400    | ~3 min     |

## What is provable vs empirical

| Component                              | Guarantee           |
|----------------------------------------|---------------------|
| Z_min from DP                          | exact               |
| Sampling from p(c) via DP              | exact               |
| Per-term η extraction from CCD         | exact (direct eval) |
| f_pair generalizes to unseen confs     | empirical (STRICT pairs validate) |
| Z_corrected from modified DP           | exact (given f_pair)|
| PAC bound on 𝔼[ψ]                     | probabilistic (Bernstein) |
| Final K* interval                      | 95% confidence      |

The only empirical assumption: f_pair(c) learned from observed conformations
generalizes to unobserved ones. STRICT-stable pair data (std/mean < 0.1 across
143 pairs) validates this. Online monitoring: each new CCD result checks whether
the per-term corrections are consistent with η.

## Implementation plan (code changes)

### New class: PACPartitionFunction.java
Location: edu.duke.cs.osprey.markstar.framework

Responsibilities:
1. Accept DP tree (rootedRoot), emats, ecalc, interactionGraph
2. Implement ancestral sampling from DP tables
3. Run parallel CCD minimizations
4. Extract per-term η
5. Build corrected emat, recompute DP
6. Compute PAC bound
7. Return (Z_lower, Z_upper, confidence)

### Modifications to BranchMARKStarBound.java
- Add flag: `-Dbranchmarkstar.usePAC=true`
- In computePartitionFunction(): if usePAC, call PACPartitionFunction instead
  of the main search loop (tightenFlatSum)
- Expose: rootedRoot, rootedRootEdge, branchMinimizingEmat, branchRigidEmat,
  interactionGraph, minimizingEcalc, rcs as package-private for PAC access

### Modifications to RootedTreeEdge.java
- Add method: sampleLambdaState(int mIdx, Random rng) → int lIdx
  Samples lambda state proportional to exp(-fullEnergyMin[mIdx][lIdx]/RT)
- Expose: decodeLambdaState, decodeMState as public (already done: decodeMStatePublic, decodeLambdaStatePublic)

### New utility: PerTermEnergyEvaluator.java
Evaluates individual energy terms (one-body, pairwise) at specific DOF values.
Wraps the forcefield to compute E_oneBody(i, DOFs_i) and E_pair(i,j, DOFs_i, DOFs_j)
for a given ParametricMolecule with optimized coordinates.

### Test
Add test in TestBranchMARKStar:
- Run PAC on the existing 2RL0 benchmark (numFlexible=10)
- Verify K* interval contains the MARK* deterministic result
- Compare CCD count and wall-clock time

## Relationship to prior work in this project

| Prior approach                    | Status | Why PAC is different                    |
|-----------------------------------|--------|-----------------------------------------|
| Region-atom ratio overlay         | failed | Needed R∪B in one edge; accounting bugs |
| Region-atom DP integration        | failed | M-set empty; R-owned energy spans edges |
| η^+ deterministic correction      | blocked| A5 not provable; per-pair context dep.  |
| Per-edge joint CCD (provable)     | viable | Captures only within-edge effects (~30%)|
| **PAC + Rao-Blackwellization**    | **new**| **Bypasses E_rigid, search, A5 entirely**|

## Risks

1. **f_pair generalization failure**: if η learned from 500 samples doesn't predict
   well for other conformations, residual variance is larger → ε_pac worse.
   Mitigation: monitor per-sample residual; increase N if variance is high.

2. **Ancestral sampling quality**: if the DP tree structure makes sampling inefficient
   (e.g., strong correlations between edges), samples may not cover the important
   conformations. Mitigation: check effective sample size; consider alternative
   sampling (MCMC on product space).

3. **Per-term energy evaluation**: need to evaluate individual energy terms at
   CCD-optimized coordinates. If the forcefield API doesn't expose this cleanly,
   may need to work around it. Mitigation: ConfAnalyzer already does energy
   breakdown; adapt its code.

4. **K* semantics**: switching from deterministic to probabilistic bound changes
   the guarantee. For drug design ranking, 95% confidence is sufficient. For
   rigorous proofs, it's not. Mitigation: make PAC optional (flag-controlled),
   keep deterministic as default.

## References

- UC-EDBR report: slurm/analysis/2026-05-19-uc-edbr-report.md
- Optimal ε plan: slurm/analysis/2026-05-21-optimal-eps-plan.md
- Region-atom EDBR: slurm/analysis/2026-05-19-region-atom-edbr.md
- Valiant 1984: PAC learning framework
- Wainwright et al 2005: Tree-reweighted belief propagation
- Nilsson 1998: Top-K enumeration from junction trees
