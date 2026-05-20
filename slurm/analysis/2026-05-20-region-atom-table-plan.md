# Region Atom Table / EDBR implementation plan

This note turns the current region-atom idea into an implementation and
diagnostic plan for MARK*/BranchMARK*. The short version is:

- a region atom is a separator-conditioned local partition-function table,
  not a scalar correction;
- the first useful scout can be computed before leaf minimization, using only
  RCs, the current energy matrices, and the graph/factorization;
- the deterministic certificate should be integrated as a DP/factor-graph
  replacement, not as `correctionMatrix.confE(partialConf)`;
- BranchMARK* sparse mode is the best first target because it already has the
  graph/separator structure needed to make `R/B` meaningful.

## 1. Definitions: what are R and B?

Let the pfunc variables be residue positions:

```text
V = {0, 1, ..., n-1}
c_i = rotamer/RC assignment at position i
```

Let `G` be the interaction graph used by the certificate. This detail matters:

- in BranchMARK* sparse mode, `G` is the sparse interaction graph;
- in full pairwise mode, `G` is effectively the full pairwise graph unless we
  introduce an explicit residual-bound treatment for cut edges;
- using the wrong graph makes the separator proof invalid.

A region atom is:

```text
a = (R, B)
```

where:

- `R` is the internal region whose assignments are locally summed out;
- `B` is the boundary/separator whose assignment indexes the table;
- in a pairwise graph, the safe default is:

```text
B = N_G(R) \ R
```

`B` is not summed out by the atom. It is the key:

```text
b -> [L_R(b), U_R(b)]
```

For each boundary assignment `b`, the atom certifies:

```text
L_R(b) <= Z_R(b) <= U_R(b)

Z_R(b) = sum_{x_R} exp(-E_R(x_R; b) / RT)
```

The energy ownership rule is:

```text
E_R(x_R; b)
  = one-body terms in R
  + R-R pair/factor terms
  + R-B crossing pair/factor terms
```

The atom must not own:

- one-body terms in `B`;
- `B-B` terms;
- outside-only terms;
- terms crossing `B` to outside that do not touch `R`.

This ownership rule is what prevents double counting when the table is inserted
back into the global DP/factor graph.

## 2. Can we infer useful tables from computed energy matrices alone?

Yes, but with an important boundary:

From RCs + `rigidEmat` + `minimizingEmat` + the current graph, we can infer:

1. candidate `R/B` atoms;
2. table size and expected cost;
3. a pre-minimization local Z sandwich table;
4. which boundary/local slices are high-gap;
5. a ranking of which atoms should be worth certification.

We cannot infer the final minimized local oracle improvement from emats alone.
The following require local certification/minimization:

```text
E_true
actual upperDrop = E^+_rigid - E_true
deterministic eta^+ or tightened local table values
```

So the scout should be viewed as a deterministic pre-min prioritizer, not as the
final certificate.

### Pre-min table from emats

For a candidate atom `(R, B)`, and a boundary assignment `b`, enumerate local
assignments `x_R` and compute owned local energies from the two existing emats:

```text
E^+_R(x_R; b) = owned local energy from rigidEmat
E^-_R(x_R; b) = owned local energy from minimizing/corrected emat
```

Then:

```text
L_R^pre(b) = sum_{x_R} exp(-E^+_R(x_R; b) / RT)
U_R^pre(b) = sum_{x_R} exp(-E^-_R(x_R; b) / RT)
gap_R^pre(b) = U_R^pre(b) - L_R^pre(b)
```

This table is already a valid emat-level local bound table, but it usually has
the same looseness as the current global bound. Its main purpose is ranking:

```text
score(atom) = expected epsilon drop / table cost
```

A simple first score is:

```text
score(R,B) =
  sum_b max(0, U_R^pre(b) - L_R^pre(b))
  / (prod_{i in B} |RC_i| * prod_{j in R} |RC_j|)
```

For leaf-queue scout data, use the same idea on selected leaves:

```text
score_i(a)    = selected ZGap mass with c_i = a
score_ij(a,b) = selected ZGap mass with c_i = a, c_j = b
```

where:

```text
ZGap = Z^+ - Z^-
preGap = E^+_rigid - E^-_corrected
```

These quantities are known before leaf minimization.

## 3. Why the sparse BranchMARK* run does not refute the pos7 analysis

The recent BranchMARK* sparse mode run doing more leaf minimizations than MARK*
without reaching `epsilon = 0.68` does not by itself invalidate the hotspot
analysis.

It says the current BranchMARK* sparse execution path is not yet exploiting the
hotspot as a table/factor replacement. It is still mostly doing queue expansion
and leaf minimization under the existing bound mechanics. The pos7 conclusion
was about where the pre-min and post-min gap mass is concentrated:

- pos7 has the largest rotamer entropy;
- pos7 rotamer slices already control much of `preGap`;
- high-mass selected leaves are enriched for pos7 slices;
- many stable pair-conditioned drops involve pos7.

Those are statements about where table amortization should pay off. They are
not statements that the existing BranchMARK* sparse leaf-minimization loop must
already converge faster.

In fact, if sparse BranchMARK* spends many minimizations without closing the
bound, that strengthens the case for table/factor replacement: the algorithm is
paying leaf-by-leaf for a repeated local uncertainty that should be certified
once per boundary slice.

## 4. MARK* vs BranchMARK* target

### MARK*

MARK* is useful for an offline prototype:

- existing leaf logs expose strong pos7 signals;
- `calcEnergyAsync(tuple)` can prototype CPU local oracles;
- it is easy to replay selected leaves and run retroactive what-if analysis.

But MARK* is not the clean final integration point for region atoms because its
current correction mechanism is scalar/tuple-like and assigned-conf oriented.
Putting a table into `correctionMatrix.confE(partialConf)` risks repeating the
same unsafe abstraction: treating a partial g-score correction as a subtree
partition-function bound.

### BranchMARK* sparse mode

This is the best first target.

Reasons:

- it already has an interaction graph;
- separators are meaningful;
- DP/factorization is already the natural language of the bound;
- sparse `B = N_G(R) \ R` is likely small enough for tables;
- region atoms can be inserted as separator-conditioned factor replacements.

The certificate is for `Z_sparse`. If we also need full-pfunc certification,
the sparse-to-full residual inflation must be applied after the sparse result.

### BranchMARK* full mode

Full mode should be inspected, but it is not the best first target.

In a complete pairwise graph, for a small `R`:

```text
B = V \ R
```

which usually makes the table enormous. To use full mode productively, we need
one of:

1. complete graph but very tiny `R`;
2. explicit factor ownership with bounded residual cut terms;
3. a hybrid graph where omitted full-mode edges are handled by a certified
   residual envelope.

Without that, a full-mode region table will either explode or silently miss
outside interactions.

## 5. Multiple-table route

The first implementation should support multiple candidate atoms, but it should
be conservative about how they are combined.

### Phase 0: emat-only atom scout

Inputs:

```text
RCs
branchRigidEmat
branchMinimizingEmat
interactionGraph
optional current DP tables
```

Outputs:

```text
[REGION_ATOM_SPEC]
[REGION_ATOM_TABLE]
[REGION_ATOM_BOUNDARY]
[REGION_ATOM_LOCAL]
[REGION_ATOM_WHATIF]
```

No leaf minimization is required.

What this phase answers:

- what is `R`;
- what is `B`;
- how many boundary cells;
- how many local states per boundary;
- where the pre-min local Z gap is concentrated;
- whether the table is small enough to certify.

### Phase 1: queue scout before minimization

Run BranchMARK* until it exposes high-priority leaves, but intercept before
calling leaf minimization. For each exposed leaf, record:

```text
conf
oldLowerZ
oldUpperZ
oldZGap
preGap = E^+_rigid(conf) - E^-_corrected(conf)
boundary key for each atom
local key for each atom
```

This is still pre-min data. It answers whether the emat-only atom ranking
matches the actual search pressure.

### Phase 2: table-level what-if

For each candidate atom, simulate improvements such as:

```text
lower local energy bound increases by delta
upper local energy witness decreases by delta
local U_R(b) shrinks by a factor
local L_R(b) grows by a factor
```

Report:

```text
epsilon_before
epsilon_after
epsilon_drop
gap share captured
jobs required
drop per job
```

This is the go/no-go point before expensive local oracle work.

### Phase 3: deterministic local table oracle

For a selected atom `(R,B)`, build jobs:

```text
TableJob(atomId, boundaryAssignment b, localAssignment x_R, ownedInteractions)
```

CPU path:

```text
chunk jobs -> TaskExecutor -> fixed table indices
```

GPU path:

```text
TableJob -> compiled ConfEnergyCalculator.MinimizationJob
batch by CudaConfEnergyCalculator.maxBatchSize()
feed batches like COFFEE MinimizationQueue
```

The important design point is deterministic indexing:

```text
tableIndex = atomId + boundaryIndex + localIndex
```

Each job writes exactly one fixed slot. Final merge/validation is single-threaded
and stable.

### Phase 4: DP/factor replacement

The certified table should be inserted into the branch DP/factor graph as a
separator-conditioned local factor:

```text
old local factor over R with boundary B
        replaced by
certified [L_R(b), U_R(b)]
```

This is not a correction added to every partial conformation. It is a local
partition-function replacement under the separator proof:

```text
Z(b) = Z_outside(b) * exp(-E_B(b)/RT) * Z_R(b)
```

If:

```text
L_R(b) <= Z_R(b) <= U_R(b)
```

then multiplying by the nonnegative outside factor preserves the global bound.

### Phase 5: multiple atom composition

Default safe rule:

- disjoint owned factor sets can be multiplied/combined directly;
- overlapping atoms cannot be additively combined by default;
- if two atoms overlap, use factor replacement ownership or take conservative
  max/min style certificates.

For energy certificates:

```text
E_lower(c) = max_a E_lower_a(c)
E_upper(c) = min_a E_upper_a(c)
```

For partition-function factors, prefer explicit ownership and replacement over
ad hoc additive correction.

## 6. Candidate selection algorithm from emats

A practical automatic scout can do:

1. Build the graph `G` used by the current run.
2. Compute per-position metadata:

```text
rotamer count
graph degree
sum of incident pair gap magnitudes
local one-body gap
selected leaf ZGap share if queue scout is enabled
```

3. Seed regions from high-score positions:

```text
R0 = {i}
```

4. Expand by neighbors while table cost is under budget:

```text
R <- R union {best neighbor}
B <- N_G(R) \ R
cost <- prod_{i in R} |RC_i| * prod_{j in B} |RC_j|
```

5. Keep candidates with good benefit/cost:

```text
benefit = pre-min local ZGap share
cost = total local jobs
score = benefit / cost
```

6. Deduplicate candidates with same owned factor set.

For pos7, the hand-picked first candidates are:

```text
R = {7}
R = {6,7,8}
R = {0,2,6,7,8}
```

The boundary is computed from the graph:

```text
B = N_G(R) \ R
```

The single-position atom `{7}` is useful as a sanity check, but the better
scientific hypothesis is a pos7-centered interaction region, because the signal
comes from pos7 together with positions like `0,2,6,8`.

## 7. Implementation touchpoints

### BranchMARKStarBound

Add a diagnostic mode controlled by properties:

```text
branchmarkstar.regionAtom.enabled=true
branchmarkstar.regionAtom.scoutOnly=true
branchmarkstar.regionAtom.regions=7;6,7,8;0,2,6,7,8
branchmarkstar.regionAtom.scout.maxLeaves=2000
branchmarkstar.regionAtom.scout.topCells=12
branchmarkstar.regionAtom.table.maxJobs=50000000
branchmarkstar.regionAtom.whatIfDeltas=0.1,0.3,0.5
```

Initial integration:

- initialize atom specs after `interactionGraph`, `branchRigidEmat`, and
  `branchMinimizingEmat` are ready;
- compute emat-only local pre-min tables;
- in scout-only mode, when leaves are pulled, record them and put internals
  back, but skip actual minimization;
- print summary after max scout leaves or queue exhaustion.

This mode is diagnostic only and should not set a final certified pfunc value.

### DP/factor code

The final certified implementation should live closer to:

```text
RootedTreeEdge DP tables
BranchMARK* sparse DP factorization
```

not inside:

```text
correctionMatrix.confE(partialConf)
```

The DP layer is where a separator-conditioned table is mathematically natural.

### Batch oracle

CPU prototype can use the existing task executor. The production GPU route
should follow the compiled/COFFEE path:

```text
edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator.MinimizationJob
CudaConfEnergyCalculator.minimizeEnergies(List<...>)
coffee.NodeProcessor.MinimizationQueue
```

The older MARK* pattern:

```text
calcEnergyAsync(tuple)
```

is fine for a CPU prototype, but it is not the right final GPU batching model.

## 8. SLURM diagnostics to add

Create a short-running scout job that does not attempt to finish the pfunc:

```bash
sbatch slurm/scripts/run_region_atom_table_scout.slurm
```

Suggested defaults:

```text
NUM_FLEXIBLE=10
EPSILON=0.68
BUDGET=0.5
ENERGY_MODE=sparse
REGIONS="7;6,7,8;0,2,6,7,8"
SCOUT_MAX_LEAVES=2000
TOP_CELLS=12
TABLE_MAX_JOBS=50000000
```

The output should be parseable with stable prefixes:

```text
[REGION_ATOM_SPEC]
[REGION_ATOM_TABLE]
[REGION_ATOM_SCOUT]
[REGION_ATOM_BOUNDARY]
[REGION_ATOM_LOCAL]
[REGION_ATOM_WHATIF]
```

## 9. Success criteria

The route is promising if the diagnostic shows:

1. table size is feasible:

```text
prod |RC_R| * prod |RC_B| <= tableMaxJobs
```

2. a small number of boundary/local cells capture a large fraction of selected
   `ZGap`;
3. what-if deltas produce epsilon drops comparable to thousands of leaf
   minimizations;
4. the same atom family appears in both emat-only scout and queue scout;
5. for pos7, `{6,7,8}` or `{0,2,6,7,8}` beats `{7}` on benefit/cost.

## 10. Main risks and mitigations

### Boundary explosion

Risk:

```text
prod_{i in B} |RC_i|
```

gets too large.

Mitigation:

- start in sparse mode;
- cap `table.maxJobs`;
- shrink `R` or use branch separators instead of raw graph neighbors;
- choose atoms by benefit/cost, not just raw gap.

### Wrong ownership / double counting

Risk: table includes `B` one-body or `B-B` terms and gets counted twice.

Mitigation:

- table owns only `R`, `R-R`, and `R-B`;
- DP/global factor owns `B` and outside terms;
- emit an ownership audit for each atom.

### Overlapping atoms

Risk: additive combination becomes unsound.

Mitigation:

- begin with one atom at a time;
- then allow only disjoint owned factor sets;
- for overlapping atoms, use explicit factor replacement or conservative
  max/min certificates.

### Upper-side local witness

Risk: a local minimization witness may depend on continuous outside geometry not
represented by boundary RCs.

Mitigation:

- first certify emat-level local tables;
- for tightened upper-side tables, include all necessary boundary variables or
  use a worst-case envelope;
- keep the local oracle deterministic and auditable.

## 11. Immediate next steps

1. Implement BranchMARK* sparse `regionAtom.scoutOnly` diagnostics.
2. Run scout on the pos7 candidates:

```text
7
6,7,8
0,2,6,7,8
```

3. Compare:

```text
emat-only local ZGap
queue-selected ZGap
boundary cell concentration
what-if epsilon drop
jobs per epsilon drop
```

4. If pos7-centered tables pass the benefit/cost test, implement the first
   deterministic local table oracle using CPU batching.
5. Move the oracle to compiled/CUDA batching once the certificate table format
   is stable.
6. Integrate certified tables into BranchMARK* sparse DP as factor replacements.

The key experimental question is not whether pos7 is a hotspot. The current
evidence already says it is. The key question is which separator-conditioned
atom gives the best amortized certificate:

```text
{7}
{6,7,8}
{0,2,6,7,8}
or an automatically selected sparse-graph neighborhood
```

