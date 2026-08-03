> **Separate historical branch.** This plan documents SPACK/PACK state-bound
> WMB-IS work. COHERE-IDP reuses parts of its known-logQ and certificate
> machinery, but COPPER is specified and tracked in
> [`COHERE_IDP_CURRENT.md`](COHERE_IDP_CURRENT.md). Statements labelled
> "current" below apply to this WMB-IS branch at their recorded date.

# SPACK\*/PACK\* low-cost proposal for DP-blown-up sequence subtrees — WMB-IS Plan

> Status: WMB-IS sampler + SPACK\* state-bounder fallback implemented and validated on real
> SPACK\* sequence subtrees, including exact real-CCD enumeration and large root-subtree stress
> tests.  The fixed-range certificate now has a local weight cap: exact support enumeration when
> possible, otherwise a max-sum mini-bucket cap over proposal conditional penalties.  A real
> 8,000-sequence SPACK\* run completed in 5.42 h, but did not invoke WMB-IS fallback; its current
> certificate still does not beat deterministic WMB.  Answers:
> when a sequence subtree's branch-DP table is too large to solve exactly (`TooLargeForDenseDP`),
> what should SPACK\*/PACK\* sample from instead?
> Verdict: **Weighted-Mini-Bucket importance sampling (WMB-IS)**, not mean-field. Mean-field is
> demoted to a diagnostic/init role only. See `docs/SPACKStar_Sequence_Pruning_Plan.md` §6 for
> the (different) use of WMB/MF as a deterministic bound-combiner.
>
> **2026-07-16 scope note:** WMB-IS remains a PACK\* state-oracle fallback, but it is not the
> current primary sequence-pruning route.  The active SPACK\* scaling experiment uses deterministic
> bound/unbound interaction cancellation; its recovered state, invalidated job `12097451`, frozen
> runner protocol, and N=6→N=16 decision gate are recorded in
> `docs/SPACKStar_Sequence_Pruning_Plan.md` §0.

## 1. The problem

`BranchDpBackend`'s DP table over a branch-decomposition edge's separator costs `D^branchwidth`.
Real example (`spackstar_3p4p_11963832_4.err`): 6 separator positions `[1,2,3,4,5,7]`,
`mStates=5,772,334,995` (implied `D≈42` per position) `> Integer.MAX_VALUE` →
`Status.Aborted`, `qprime=+inf`, subtree becomes unprunable. This is a hard combinatorial wall,
not a tuning problem — the fix is a proposal distribution whose cost doesn't depend on
branchwidth at all.

PACK\*'s existing estimator (`PackStarEstimator`) needs, at the single method boundary
`sampleConformationsFromDP(int n, Random rng) -> List<int[]>`, a proposal `p` that is (a) exactly
normalizable and (b) exactly sample-able, so `w(c) = exp(-E_true(c)/RT) / p(c)` gives an unbiased
`Z` estimate (`Z = E_p[w]`, standard IS identity — no PAC approximation in this step). The DP does
this exactly today; the question is what replaces it when the DP can't be built.

## 2. Verdict

**Use WMB-IS.** Do not invest further in mean-field as a primary proposal — it is a real,
measured dead end for this use case (data below). Mean-field stays useful only as (a) a coordinate
-ascent initializer for something else, and (b) the deterministic lower-bound tightener in
`OspreyWmbMfStateBounder`.

Current implementation status:

- `WeightedMiniBucket.Proposal` now provides top-down conditional sampling and exact per-sample
  `logQ`.
- `WmbImportanceSamplingStateBounder` estimates a SPACK\* state sequence subtree with WMB-IS,
  scores samples using `ConfEnergyCalculator`/CCD, and reports `cv`, `ESS/n`, unique samples, and
  CCD calls.
- `SpackStarPackStar` can wrap PACK\* state bounders with WMB-IS fallback when PACK\* fails or
  returns a full-range bound, while keeping deterministic WMB/MF tightening separate.
- The WMB-IS route now has a fail-closed fixed-range finite-sample certificate.  The empirical
  interval remains a diagnostic; the returned pruning-safe state bound is the intersection of the
  fixed-range certificate, the probability-tail diagnostic cap when applicable, and deterministic
  WMB/MF tightening.

## 3. Evidence (measured today, not just argued)

Three data sources, all grounded in this repo:

**(a) Real production baseline — PACK\*'s current DP+eta-correction proposal.**
`grep cvPsi= /usr/xtmp/lz280/slurm_logs/pac_11920687_*.out` → n=3871 real calls:
`median=0.105, p95=0.367, mean=0.139`. In IS terms, `ESS/n = 1/(1+cv²)` → **88–99% effective
sample efficiency**. This is the bar any replacement proposal is judged against.

**(b) Earlier — synthetic-but-shape-matched benchmark run against the real code** (not a
reimplementation): `WmbModel` / `WeightedMiniBucket` / `MeanFieldBound` exercised directly via a
throwaway harness (`ExploreProposalPotential`, methodology in §9, deleted after use —
reproducible from this doc). Synthetic emats use OSPREY-realistic steric-clash structure (most
RC-pairs mild, ~10–15% hard clashes with large positive energy), at three scales, the largest
matching the real abort case exactly (`numPos=6`, `D=42`):

| tier | numPos, D | exact logZ available? | WMB upper gap to exact @ iBound=2/3/4 | mean-field-IS cv | ESS/n |
|---|---|---|---|---|---|
| A | 5, 6 | yes (brute force) | 1.74 / 0.81 / 0.00 (exact) | 8.2–13.6 | **1.4% → 0.5%** |
| B | 6, 15 | yes (11.4M, brute force) | 4.01 / 2.12 / 1.07 | 15.3–19.7 | **0.4% → 0.3%** |
| C | 6, 42 (**real shape**) | no (5.5B, matches the abort) | n/a / n/a / n/a (bracket only, see §6) | 11.2–31.7 | **0.1% → 0.8%** |

Mean-field-IS's `cv` is **20–300x worse** than PACK\*'s real production `cv`, and it does not
improve with more samples (Tier C: cv=11→32→12 across n=2000/8000/20000) — this is weight
*collapse* (a few lucky non-clashing draws carry all the mass), not ordinary variance that
averages out. Required-sample count for a fixed target epsilon scales roughly with `cv²`
(Bernstein/empirical-Bernstein sizing, same code PACK\* already uses) → **~100–1000x more samples
needed**. Since each PACK\* sample costs a real CCD minimization (the dominant real cost, not the
proposal draw), mean-field is not actually a *cheap* fallback once this is accounted for — it is
frequently a net loss. This confirms, with numbers instead of intuition, the physical argument
from the design discussion: full factorization throws away exactly the pairwise steric-clash
structure OSPREY's energy model (and all of DEE/A\*/branch-DP) is built around.

**Caveat on (b):** the WMB column above is the deterministic `[Z-,Z+]` *bracket*, not a measured
sampler variance. The bracket's lower half is, by the code's own docstring, a deliberately crude
"worst-case min" bound, not a proposal-quality proxy — so the *full* gap (upper − lower, tens of
nats even at iBound=3–4 in the table) overstates how bad a real WMB sampler would be. The upper
(Hölder) half alone — the more relevant half — tightens much faster: within **1–2 nats (3–8x in
Z-space) by iBound=3–4**, even on this adversarial, fully-dense synthetic graph.

**(c) New — exact real-CCD sequence-subtree validation against ground truth.** Runner:
`src/test/java/edu/duke/cs/osprey/spackstar/RunWmbExactSubtreeValidation.java`, submitted via
`slurm/scripts/run_wmb_exact_subtree_validation.slurm`.

Case: 2RL0 `complex` state, root/unassigned sequence subtree, mutable residues `A156,A164`
with all-20 residue types and no WT-flex positions. This is small enough to enumerate exactly:

```text
state positions = 2
RC domains      = [189,189]
total confs     = 35,721
exact logZ      = 58.038234
exact CCD time  = 393.592 s
```

The exact value was computed by enumerating every real RC conformation, running minimized
energy/CCD for each, and accumulating:

```text
logZ_exact = log sum_conf exp(-E_CCD(conf)/RT)
```

WMB-IS on the same subtree covered the exact logZ across all tested seeds and sample counts.
Representative results:

| job | iBound | seed | samples | unique CCD confs | logZ interval | center error | covers exact? |
|---|---:|---:|---:|---:|---|---:|---|
| `11991963` | 1 | 1 | 128 | 34 | `[58.023862, 58.054496]` | `+0.000945` | yes |
| `11991963` | 1 | 1 | 512 | 61 | `[58.031877, 58.043200]` | `-0.000696` | yes |
| `11991963` | 1 | 1 | 2048 | 100 | `[58.036477, 58.040635]` | `+0.000322` | yes |
| `11991963` | 1 | 2 | 2048 | 97 | `[58.035947, 58.040261]` | `-0.000130` | yes |
| `11991963` | 1 | 3 | 2048 | 113 | `[58.035676, 58.041047]` | `+0.000127` | yes |

`iBound=2` produced identical numbers in this two-position case because the elimination width is
already one. This validation confirms the implemented chain `WMB proposal -> exact logQ ->
real CCD energy -> IS logZ` on a real OSPREY sequence subtree. It does **not** prove absolute
logZ correctness for the much larger subtrees below, where exact enumeration is infeasible.

**(d) New — real SPACK\* sequence-subtree WMB-IS validation.** Runner:
`src/test/java/edu/duke/cs/osprey/spackstar/RunWmbSequenceSubtreeValidation.java`, submitted via
`slurm/scripts/run_wmb_sequence_subtree_validation.slurm`.

Case: 2RL0 `complex` state, root/unassigned sequence subtree, mutable residues `A156,A164`
with all-20 residue types plus WT-flex positions `A172,A192,A193,G649,G650,G651`.

```text
state positions = 8
RC domains      = [189,189,28,8,19,5,6,9]
total confs     = 41,047,715,520
```

This is a **sequence subtree** estimate:

```text
Z(subtree) = sum_{seq in subtree} sum_{conf compatible with seq} exp(-E_CCD(conf)/RT)
```

It is not a single fully assigned sequence and not a synthetic 72-conf toy. WMB prep and WMB-IS
results on Slurm:

| job | iBound | samples | unique CCD confs | maxTableCells | logZ interval | cv | ESS/n |
|---|---:|---:|---:|---:|---|---:|---:|
| `11991904` | 1 | 64 | 37 | 35,721 | `[153.8727, 156.0547]` | 0.5328 | 0.7789 |
| `11991904` | 2 | 64 | 36 | 1,000,188 | `[154.9635, 155.8662]` | 0.3416 | 0.8955 |
| `11991905` | 2 | 256 | 94 | 1,000,188 | `[155.3660, 155.6632]` | 0.3571 | 0.8869 |
| `11991905` | 2 | 1024 | 225 | 1,000,188 | `[155.4625, 155.5706]` | 0.3564 | 0.8873 |

Interpretation:

- WMB avoids the dense branch-DP table cliff: this root sequence subtree has `4.10e10`
  conformations, but iBound=2 WMB builds million-cell mini-bucket tables in about `0.12 s` after
  the energy matrix is available.
- The real CCD-scored IS weights do **not** collapse: `cv≈0.35`, `ESS/n≈0.89`.
- Sample de-duplication matters: 1024 WMB draws required 225 unique CCD minimizations, not 1024.
- The interval shown is the current empirical Bernstein diagnostic interval over observed IS
  weights. It is useful evidence that the proposal is stable, but it is not yet a rigorous
  SPACK\* pruning certificate.

**(e) New — larger `3 full + 5 flex` root subtree stress test.** Same runner/script as (d).

Case: 2RL0 `complex` state, root/unassigned sequence subtree, all-20 mutable residues
`A156,A164,A172` plus WT-flex positions `A192,A193,G649,G650,G651`.

```text
state positions = 8
RC domains      = [189,189,189,8,19,5,6,9]
total confs     = 277,072,079,760
energy entries  = 135,247
```

This run did not perform sequence pruning and did not have exact ground truth. It measured whether
WMB-IS remains computationally usable and whether the sampled real-CCD weights are stable on a
much larger sequence subtree:

| job | iBound | samples | unique CCD confs | maxTableCells | logZ interval | cv | ESS/n | time |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| `11991964` | 1 | 256 | 86 | 35,721 | `[156.5024, 157.2394]` | 0.8464 | 0.5826 | 3.596 s |
| `11991964` | 1 | 1024 | 188 | 35,721 | `[156.6904, 157.0945]` | 0.8896 | 0.5582 | 7.311 s |
| `11991964` | 1 | 4096 | 454 | 35,721 | `[156.8069, 157.0114]` | 1.0200 | 0.4901 | 18.007 s |
| `11991964` | 2 | 256 | 90 | 6,751,269 | `[156.2794, 157.4023]` | 0.8835 | 0.5616 | 3.831 s |
| `11991964` | 2 | 1024 | 203 | 6,751,269 | `[156.7730, 157.0982]` | 0.7791 | 0.6223 | 7.955 s |
| `11991964` | 2 | 4096 | 483 | 6,751,269 | `[156.8708, 156.9842]` | 0.7441 | 0.6437 | 19.301 s |

Interpretation:

- The root sequence subtree has `2.77e11` conformations, but iBound=2 WMB prep used about
  `6.75M` table cells and took `0.349 s` after the energy matrix was available. This is the
  main practical difference from PACK\* branch-DP, whose dense separator table can hit an
  integer/table-size wall.
- The main cost is not WMB sampling. At 4096 samples, the run scored only 483 unique CCD
  conformations because duplicate samples reused the energy cache.
- The best observed estimate is centered near `logZ ~= 156.93`. From `cv=0.7441` and `n=4096`,
  the ordinary Monte Carlo standard-error scale is roughly `cv/sqrt(n) ~= 0.0116` in logZ units,
  assuming no unseen rare high-weight tail. A conservative practical uncertainty for this
  diagnostic run is closer to `0.05-0.1 logZ`, because the interval is still empirical, not a
  certified finite-sample pruning bound.

## 4. Why WMB is structurally the right unit (not a new algorithm family)

WMB-IS is not a departure from what PACK\* already does — it's the same pattern:

```text
today (DP):  Z = q_eta * E_{p_eta}[ exp(-xi/RT) ]     p_eta = exact top-down sample from branch-DP
proposed:    Z =         E_{q_WMB}[ exp(-E_true/RT) / q_WMB(c) ]   q_WMB = top-down sample from WMB messages
```

Both are ancestral samplers with an *exactly known* per-sample proposal probability, built by
construction (the sampler defines a proper distribution regardless of how tight the elimination
bound behind it is). The DP's `q_eta` happens to be an exact normalizer; WMB's `Z_WMB+` is only an
upper bound of true `Z` — but that is irrelevant to unbiasedness, only to variance (how well `q_WMB`
matches the target shape). This also means the two-stage eta-correction trick (`extractEtaCorrections`
/ `buildCorrectedEmat`) is not DP-specific — `WeightedMiniBucket` takes the same
`(EnergyMatrix, RCs, assignments, iBound, rt)` inputs the corrected-DP path already builds, so it can
run on `correctedEmat` exactly like `recomputeDP` does today (§6, Phase 3).

WMB's cost knob is decoupled from branchwidth (validated numbers, this machine, real abort-case
shape `numPos=6, D=42`):

| iBound | maxTableCells | wall time |
|---|---|---|
| 2 | 74,088 | 0.18 s |
| 3 | 3,111,696 | 2.72 s |
| 4 (not run; login node) | ~1.3×10^8 | predicted ~GB-scale, minutes |

vs. the DP table at the same shape: `5,772,334,995` cells → `Integer.MAX_VALUE` abort. The formula
`D^(iBound+1)` (validated: `42.35^6 ≈ 5.77e9` matches the real logged `mStates` exactly) degrades
*gracefully* with iBound instead of hitting a cliff — this is the actual fix for "can't legally
construct a smaller subtree."

## 5. Real limitations found while reading the code (must be designed around, not ignored)

- **`WmbModel.edges()` is hardcoded dense** ("Energy-matrix graphs are dense, so every pair of
  variables interacts"). It ignores the sparsity the rest of the codebase already computes and
  trusts — `branchdp.InteractionGraph` (`hasEdge`, `getNeighbors`, `buildWithResidualBudget`,
  already used by `PackStarEstimator.extractEtaCorrections`). Every gap number in §3 is therefore
  a **worst case** (dense graph, no locality exploited). Threading `InteractionGraph` through
  `WmbModel`/mini-bucket partitioning is free extra tightness, not new theory — do this before
  concluding WMB is "only" as good as §3 shows.
- **Mini-bucket partitioning is scope-size-greedy, not interaction-aware**
  (`WeightedMiniBucket.partition`: sorts by scope length only). For a clash-dominated energy
  model, what matters is keeping strongly-coupled (large `|pairwise energy|`) position pairs in
  the same mini-bucket, not just balancing bucket sizes. This is a likely quick win worth an
  A/B before the full sampler.
- **Fixed global `iBound` risks `int` overflow**, not just cost: a mini-bucket table is `double[]`
  of size `D^(iBound+1)`; at `iBound=3`, `D>215` overflows a 32-bit array index. Superset
  positions near the sequence-tree root can plausibly have `D` in the hundreds. This is now
  partially hardened: `WeightedMiniBucket` checks table sizes in `long` before allocation, and
  WMB-IS can lower the effective `iBound` to satisfy a configured `maxTableCells` cap. Remaining
  work is to choose a production default cap and make the policy interaction-aware rather than
  only domain-size-aware.
- **The WMB `[Z-,Z+]` bracket's lower half is not a variance proxy** (§3 caveat) — don't use the
  full bracket gap to gate "is this node's WMB proposal good enough"; use the upper-only gap, or
  better, the WMB-IS sampler's empirical `cv`/ESS diagnostics (exactly how PACK\* already
  self-monitors via `cvPsi`).

## 6. Design: WMB-IS estimator

**Implemented capability: top-down conditional sampling from WMB messages.**
`WeightedMiniBucket` now has `proposalForModel(...)`, returning a `Proposal` that samples
variables in reverse elimination order and returns `(domainValues, logQ)`. The proposal keeps the
conditional mini-bucket factors produced by the upper WMB pass; for each sampled variable it
normalizes the conditional scores under already-fixed later variables. `logQ` is the exact product
of the conditional draw probabilities, so the IS estimator is unbiased when `q_WMB` has support.

**Estimator integration.** Single-stage IS (no eta-correction initially):

```text
Z = E_{q_WMB}[ exp(-E_true(c)/RT) / q_WMB(c) ]
```

`q_WMB(c)` is known exactly per drawn sample (product of the conditional draw probabilities used
to construct it), independent of how tight `Z_WMB+` is.

Current SPACK\* integration is a state-bounder fallback, not yet a replacement inside
`PackStarEstimator.sampleConformationsFromDP(...)`: `WmbImportanceSamplingStateBounder` builds a
proposal from `(minimizingEmat, stateSequence.makeRCs(confSpace), iBound, rt)`, samples the
sequence subtree, scores each unique conformation with CCD through `ConfEnergyCalculator`, and
returns a diagnostic empirical interval.

### Finite-sample guarantee status

The WMB-IS route now has a fail-closed finite-sample state-bound certificate path.  The remaining
issue is usefulness, not basic soundness: on the measured real subtree the certified WMB-IS upper
is still looser than the deterministic WMB upper, so WMB-IS has not yet contributed additional
sequence pruning beyond the deterministic tightener.

For ordinary importance sampling:

```text
W(conf) = exp(-E_true(conf)/RT) / q_WMB(conf)
Zhat    = mean_i W(conf_i)
E[Zhat] = Z
```

If `q_WMB(conf) > 0` for every conformation in the sequence subtree and `logQ` is exact, then the
estimator is unbiased. A high-probability finite-sample interval then needs a deterministic weight
range:

```text
0 <= W(conf) <= B
```

Given such a certified `B`, Hoeffding, Bennett, or empirical-Bernstein-with-fixed-range intervals
can give:

```text
Pr[ |Zhat - Z| <= eps ] >= 1 - delta
```

and the result can be converted to a logZ interval. This is different from many Monte Carlo/IS
uses that rely only on asymptotic CLT intervals; those become unsafe when weights have heavy
tails, incomplete support, infinite/huge variance, or only an observed sample max.

What is still missing for rigorous SPACK\* pruning:

1. **Certified support. Done.** WMB proposal construction guarantees `q_WMB(conf) > 0` over the
   full subtree, including uniform fallback when a conditional has all `-inf` scores.
   `WeightedMiniBucket.Proposal.logProbabilityLowerBound()` gives a conservative global lower
   bound on proposal probability, and regression tests enumerate a small space to verify the bound
   covers every assignment.
2. **Certified weight upper bound. Local-cap version done.** The interval no longer uses observed
   range/sample max for certification.  `WmbImportanceSamplingStateBounder` first asks
   `WeightedMiniBucket.Proposal` for a local cap:

   ```text
   logW(conf) <= max_conf [ theta_min(conf) - log q_WMB(conf) ]
   ```

   For enumerable local supports, this max is exact.  For larger supports, the proposal builds a
   max-sum mini-bucket upper bound on the same objective using conditional-penalty factors that
   preserve the sampler's zero-support/uniform-fallback semantics.  If that local cap is
   unavailable, the fail-closed fallback remains:

   ```text
   E_true(conf) >= E_LB(conf)
   logW(conf) <= logZ_WMB_upper(E_LB subtree) - logQ_lower
   ```

   This fallback is still intentionally conservative; the local cap is the path expected to be
   useful on large real subtrees.
3. **Fixed-range concentration. Done.** The bounder scales weights by `logWeightUpper`, applies an
   empirical-Bernstein interval on `[0,1]`, and returns full range if the certificate cannot be
   formed or a sample exceeds the deterministic cap.
4. **Delta accounting. Done for the implementation path.** `SpackStarErrorBudget` provides named
   geometric streams; WMB-IS state certificates consume per-call delta leases.  SPACK\* assembly
   wires a shared `spackstar.wmbis` stream across complex/protein/ligand fallback state bounders.

Current conclusion: WMB-IS now has a sound, fail-closed finite-sample certificate path, and the
largest observed cap pathology is fixed.  On Slurm job `12007350`, the real 2RL0 `complex`
subtree with `A156,A164` mutable and `A192,A193,G649,G650,G651` WT-flex had
`1,465,989,840` RC conformations.  The old fallback cap was `1242117566.8368`; the local
mini-bucket cap was `156.6608`, with observed max sampled log weight `128.6639` and gap
`27.9969`.  The sampler diagnostic remained tight (`[128.2623,128.3396]`), the fixed-range
certificate upper became finite (`152.3940`), and the final returned safe upper was the
deterministic WMB upper (`139.3654`).

The standalone WMB validation does not measure sequence pruning counts.  The later 8,000-sequence
SPACK\* run measures overall pruning, but did not invoke WMB-IS fallback and therefore still does
not measure WMB-specific pruning.  The current evidence says WMB-IS is no longer blocked by a
giant sentinel/range cap, but on this state it still does not beat the deterministic WMB upper
used for sound pruning.

**2026-07-06 SPACK\* integration note.** The then-current 3/4-position SPACK\* probe (`12009291`) was
launched with WMB fallback enabled (`wmbIBound=3`, `wmbMaxTableCells=10000000`,
`wmbISFallback=true`, `wmbISSamples=1024`) and the e-value betting portfolio enabled.  Early log
checks show both array tasks are still in the BBK\*/MARK\* seed-proposal phase, with repeated
`Refining sequence ...` lines and no `==== SPACK* -> PACK* chain`, `[SPACK*-CERT]`,
`seed-verify`, or `sampleEvidence` output yet.  Therefore this run has not yet tested WMB-IS
certificate tightness in the SPACK\* search loop.  The immediate SPACK\* bottleneck is getting a
finite incumbent cutoff cheaply; WMB-IS remains the fallback for PACK\* state bounds that hit
branch-DP/table-size limits.

That probe was superseded the same day by job `12009337`: BBK seed is disabled by default and
SPACK\* starts with cheap PACK\*-verified candidate seeds (WT plus single mutants, capped at
`80`).  This deliberately isolates WMB-IS from the incumbent-seeding problem.  If `12009337`
reaches WMB fallback calls, the relevant WMB-IS question is still whether its certified state
upper is tighter than the deterministic WMB/MF upper under a finite incumbent cutoff; it should
not be judged by the earlier BBK seed latency.  Early `12009337_3` output already produced a
finite SPACK\* cutoff during seed verification (`topNLowerCutoff=28.283412` after the second seed,
then `35.102578` after the third), so subsequent WMB-IS observations from this job are finally in
the regime where state-bound tightness can matter for sequence pruning.

Follow-up from the same day: the incumbent itself is no longer the hard blocker, but verifying
all cheap seeds is too slow to use as a pre-search phase.  SPACK\* now has a driver-level
early-stop option for initial seed verification: stop once the kth incumbent cutoff is finite,
then verify only a small extra number of seeds (`spackstar.seed.extraAfterFiniteCutoff`, default
`1`) before entering the main frontier.  This keeps the WMB-IS question clean: WMB-IS should be
evaluated during subtree bounding and pruning under a finite cutoff, not during a long incumbent
proposal sweep.

The old `12009337` driver then exposed a second scheduling issue: after the seed phase it invoked
sample evidence on the root node while the root upper bound was still `+inf`, which triggered
PACK\*/WMB-style bounds for the full unassigned complex (`24.1 GiB` predicted DP table for the
3-position task and `570.4 GiB` for the 4-position task).  This should not be counted as WMB-IS
certificate tightness.  The SPACK\* driver now skips evidence calls unless both the incumbent
cutoff and the node upper bound are finite, so WMB-IS is only exercised on bounded subtrees where
it can actually contribute to pruning.

**2026-07-11 medium-scale SPACK\* follow-up.**  The successor
`spackstar_3p4p_12009341_{3,4}` probe clarifies both the promise of the current implementation and
WMB-IS's actual place in the next work:

- The three-position all-20 space (`8,000` sequences) completed in `19,498.557 s` (`5.42 h`).
  SPACK\* processed `7,641` nodes, pruned `6,806`, and verified `76` leaves.  This is an encouraging
  real throughput milestone: a PACK\*-backed 8,000-sequence search completed in hours.  Because the
  run skipped the matched exhaustive PACK\* baseline, it establishes feasibility rather than a
  speedup factor.
- Only `1/6,806` prunes was attributed to sample evidence; `6,805` were ordinary bound comparisons,
  and `rebounded=0`.  The result therefore does not yet show a WMB-IS or e-process performance
  contribution.
- The four-position space (`160,000` sequences) reached its 12 h limit at iteration `5,827`, with
  `queue=6,879` and `maxRemainingUpper=+inf`.  It had performed `5,465` lazy bounds but no
  post-seed lazy prune/verify.  PACK\* skipped several oversized candidate root edges, selected
  other decompositions, and returned finite complex-state bounds.  The remaining `K* upper=+inf`
  is therefore attributable to a `-inf` uniform protein/ligand subtree lower bound, followed by
  re-queueing of a node whose real bound is still infinite—not to WMB proposal variance.
- No `[SPACK*-WMBIS]` fallback line appears in this probe.  WMB-IS remains valuable for an actual
  PACK\*/dense-DP state-bound failure, but it was not on the critical path of this four-position
  timeout.  Increasing WMB sample count would not fix the current denominator/frontier problem.

The local 2026-07-11 verification also reran all `18` `TestWeightedMiniBucket` tests successfully,
including proposal probability, support, exact/local weight caps, iBound behavior, and the small
importance-sampling validation.  This supports the sampler implementation; it does not change the
finding that the real fixed-range certificate is still looser than deterministic WMB.

### 6.1 Next work: make the certificate useful

The current WMB-IS certificate is designed to be safe first.  It should be treated as a pruning
certificate only when the logged tightness diagnostics show the fixed range is not dominating the
interval.  The 2026-07-11 medium-scale run also establishes an ordering constraint: first make
SPACK\* progress when a uniform denominator bound is `-inf` and stop re-queueing real `+inf` nodes;
only then can an end-to-end run fairly measure whether WMB-IS helps states that actually reach its
fallback boundary.  The next concrete WMB tasks are:

1. **Keep certificate tightness instrumentation.** For every WMB-IS call, log:

   ```text
   logWeightUpper
   logWeightCapExact / logWeightCapAssignments / logWeightCapReason
   logWeightFallbackUpper / logWeightMiniBucketUpper / logWeightMiniBucketReason
   observedMaxLogWeight
   logWeightUpper - observedMaxLogWeight
   certified logZ width
   empirical logZ width
   ESS/n, CV, unique samples, CCD calls
   delta stream/index/delta
   ```

   This tells whether poor pruning comes from proposal variance, a loose deterministic cap, too
   little sampling, or overly small delta leases.
2. **Push the certified upper below deterministic WMB.** The local cap reduced the real-state gap
   to `~28` log units, but the certificate upper is still above the deterministic WMB upper on the
   tested subtree.  Next candidates are residual/tail-source bounds, better WMB position ordering,
   and branch-decomposition-aware local caps.
3. **Keep fail-closed semantics.** If the deterministic cap is unavailable or a sample exceeds it,
   return full range and report the reason.  Do not silently fall back to observed range for
   pruning.
4. **Compare diagnostic vs certified intervals.** Keep the old empirical interval as a diagnostic
   only.  Large gaps between diagnostic and certified intervals indicate the range certificate, not
   the sampler, is the limiting factor.
5. **Baseline and attribute full sequence-pruning impact.** The 8,000-sequence run now supplies a
   first complete counter set, but not a matched baseline: `6,805` ordinary bound prunes, `1`
   sample-evidence prune, no WMB fallback, and no rebound.  Run same-oracle exhaustive PACK\* plus
   SPACK\* ablations (`deterministic-only`, `+e-process`, `+WMB-IS`) on the `~400`- and
   `8,000`-sequence spaces before assigning speedup to WMB-IS.
6. **Keep WMB-IS separate from incumbent and denominator scheduling.** Cheap PACK\*-verified seeds
   now raise a finite cutoff after two candidates, so initial seeding is no longer the hard blocker.
   WMB-IS tightens or replaces a state partition bound; it does not repair a uniform denominator
   lower bound of `-inf`.  Branch such a node immediately (preferably on the witness-conflict
   position) instead of spending a complex PACK\*/WMB call and re-queueing it at `+inf`.  Likewise,
   do not call WMB-IS/sample evidence on a lazy node whose inherited upper is still `+inf`.

## 7. Implementation stages

**Stage 0 — real sequence-subtree de-risking. Done for 2RL0 root complex subtree.** See §3(c-e):
exact `35,721`-conformation real-CCD enumeration matches WMB-IS within `~1e-3-4e-3 logZ`; the
larger `41,047,715,520`-conformation subtree has iBound=2 WMB table max `1,000,188` cells,
real CCD-scored `cv≈0.35`, `ESS/n≈0.89`; the `3 full + 5 flex` `277,072,079,760`-conformation
subtree has iBound=2 WMB table max `6,751,269` cells and 4096-sample `cv≈0.74`, `ESS/n≈0.64`.
Still open: pull the exact logged `spackstar_3p4p`
`TooLargeForDenseDP` node (`positions [1,2,3,4,5,7]`) and repeat the same measurement on that
specific node.

**Stage 1 — top-down WMB sampler + single-stage IS estimator. Partially done.** The sampler and
single-stage state estimator exist:
`WeightedMiniBucket.Proposal` and `WmbImportanceSamplingStateBounder`. They are wired into
SPACK\* assembly as a fallback state bounder when PACK\* fails or returns a full-range bound.
Still open: wire WMB-IS behind `PackStarEstimator.sampleConformationsFromDP(...)` so the existing
PACK\* estimator can fall back at the exact `DPTableTooLargeException` boundary.

**Stage 2 — interaction-aware mini-bucket partitioning + adaptive iBound. Partially done.** The
current WMB graph is still dense and mini-bucket partitioning is still scope-size-greedy. The
table-size overflow path is now guarded by `long` size checks plus an optional `maxTableCells`
cap that lowers effective `iBound` for WMB-IS and deterministic WMB/MF tighteners. Remaining
work: pick a production default cap and make partitioning interaction-aware.

**Stage 3 (stretch, do only if Stage 1–2 look good) — port two-stage eta-correction onto WMB.**
Mirror `recomputeDP`: learn eta on a train split from the base WMB proposal, rebuild
`correctedEmat`, rerun WMB on it for a refined `q_WMB`, then estimate the (hopefully much smaller)
residual exactly as `computePACBoundResidual` does today. This is the Rao-Blackwellization step
that gets PACK\*'s real `cv≈0.1` — worth attempting once Stage 1 shows single-stage WMB-IS is in
the right neighborhood (say `cv < 2–3`, not `cv > 10`).

**Assembly wiring. Done for SPACK\*.** `SpackStarPackStar` now wraps each state in a fallback
bounder: try PACK\* first; if PACK\* fails or returns `[-inf,+inf]`, run WMB-IS; then optionally
intersect with deterministic WMB/MF when usable. WMB-IS fallback now reports a fixed-range PAC
certificate or full range; PACK\* evidence batches remain fail-open because WMB-IS does not yet
produce PACK\* residual diagnostics for e-value evidence.

## 8. Compute

The WMB proposal and state estimator are CPU Java code. Real validation that scores samples with
CCD should run through Slurm rather than an interactive node. Current checked-in runners:

- `slurm/scripts/run_wmb_sequence_subtree_validation.slurm` — real SPACK\* sequence-subtree
  validation; defaults to CPU-only `grisman`, `complex` root subtree, configurable with
  `STATE`, `IBOUNDS`, `SAMPLES`, `ASSIGNMENTS`, `LOCAL_CAP_MAX_ASSIGNMENTS`,
  `PROB_TAIL_QUANTILE`, and residue-list environment variables.
- `slurm/scripts/run_wmb_exact_subtree_validation.slurm` — exact real-CCD enumeration versus
  WMB-IS for small real sequence subtrees; default `A156,A164` all-20 mutable with no WT-flex
  enumerates `35,721` conformations.  It accepts the same WMB-IS local-cap and probability-tail
  environment variables.
- `slurm/scripts/run_wmb_is_small_validation.slurm` — small synthetic WMB/logQ/IS math smoke test;
  useful for sampler regressions but **not** evidence for sequence-subtree performance.

GPU is not required for the current WMB implementation. `ConfEnergyCalculator` may use whatever
energy-calculator backend the caller supplied, but the validation runs above used CPU Slurm.
The measured bottleneck is real minimized-energy/CCD scoring, not WMB table construction or WMB
sampling. For example, in the `3 full + 5 flex` run, iBound=2 WMB prep took `0.349 s`, while
4096 samples plus 483 unique CCD scores took `19.301 s`. If GPU work is added, batching the
unique sampled conformations for energy/minimization should be prioritized over rewriting WMB.

## 9. Harness for reproducing §3 (recreate on demand; not checked in)

`src/test/java/edu/duke/cs/osprey/wmb/ExploreProposalPotential.java` (throwaway, deleted after
use): builds `clashEmat` (mild background + ~10–15% hard-clash pairwise energies) at the three
scales in §3 via the same fixture pattern as `TestWeightedMiniBucket.randomEmat`, computes exact
`logZ` by brute force where feasible, sweeps `WeightedMiniBucket.boundsForModel` over iBound,
copies `MeanFieldBound`'s coordinate-ascent loop to expose `q[][]` (cross-checked against
`MeanFieldBound.logZLower` to confirm the copy is faithful), and does plain self-normalized
importance sampling from `prod_i q_i` to measure `cv`/ESS/bias against brute-force truth. Rebuild:
`./gradlew compileTestJava exportTestClasspath`, then
`java -cp "build/classes/java/main:build/classes/java/test:build/resources/main:build/resources/test:$(cat bench_logs/test_classpath.txt)" edu.duke.cs.osprey.wmb.ExploreProposalPotential`.
