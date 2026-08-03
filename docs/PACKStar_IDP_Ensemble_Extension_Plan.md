> **Historical research ledger.** This file preserves the chronological
> COHERE-IDP audit trail through §10.59. Earlier statements labelled
> "current", "next", or "only handoff" are historical at their point in
> time. The sole current status and execution handoff is
> [`COHERE_IDP_CURRENT.md`](COHERE_IDP_CURRENT.md). Do not delete or rewrite
> frozen results, job records, hashes, or failed branches in this ledger.

# COHERE-IDP: PACKStar/WMB-based IDP/IDR Observable Estimation Plan

Status: active research prototype + feasibility audit, not yet scoped/approved
by Bruce and not yet paper-ready.
Updated 2026-07-19 after the independent SASDXC6 ensemble-average holdout,
known-logQ proposal audit, MultiFoXS baseline, first generic IDP/WMB adapter,
failed locked steric-surrogate replication gate, its importance-weighted
held-out FoXS follow-up, and the corrected but negative MaxEnt/MAP blocked-CV
result. Unrestricted MaxEnt has now also failed the independent preregistered
SASDNV6/E1A FoXS endpoint and every registered candidate-pool perturbation;
the E1A secondary forward-model sensitivity is complete. The current evidence
and claim boundaries are recorded through §10.52.
Origin: 2026-07-02/03 discussion, triggered by reading PACK* (pac_draft) and
asking "what's actually still open in structural biology, and what's the best
next problem given (a) MD is banned, (b) 48x A5000 on grisman, (c) existing
PACK*/WMB* machinery."

Working name: **COHERE-IDP** — Certified Observable-guided Heterogeneous
Ensemble Refinement for IDPs. Avoid new PACK*/PackStar-derived names; the
new capability is certified observable estimation, not another K*/PACK*
variant.

## 1. Motivation / positioning (why this, why now)

Landscape read (2026):
- Computational protein design (DEE/A*/K*) is being reshaped by deep learning
  (AlphaFold3, ProteinMPNN, RFdiffusion, AlphaProteo) for backbone/sequence
  generation, but ddG / binding-specificity prediction for point mutants
  (SKEMPI-hard-case territory) is still unsolved and still the active
  frontier — this is what multi-backbone K* already targets.
- NMR structure/assignment lines ([8][7][9][16][11][13][15][6][14]) are mostly
  superseded: assignment by CYANA/FLYA/ARTINA, homology detection by
  AlphaFold, backbone-from-minimal-RDC motivation weakened by AlphaFold
  giving a starting model for free. The interesting exception: AlphaFold-
  informed NOE assignment (FAAST, RASP-based protocols, 2025) is literally
  reviving the NVR "structure-informed assignment" idea 15-20 years later,
  now that AlphaFold is a good enough prior — vindicates the idea, not the
  old implementation.
- IDP/IDR conformational ensemble determination ([10], Wang & Donald 2006) is
  the one problem in the old NMR list that is *still* open and *still* being
  actively attacked at top venues in 2025-2026 (bAIes/Nat Comm, CALVADOS/
  AF-CALVADOS, IDPFold/IDPFold2, X-EISD-style Bayesian reweighting, a Nature
  Methods paper literally titled "Toward a unified framework..." in 2026).
- Blue-ocean read relative to two other candidates (plasma-proteomics
  diagnostics, single-cell circadian rhythm detection): IDP ensembles has
  fewer groups than plasma proteomics (which is a red ocean: UKB-PPP/deCODE
  pQTL-MR mining is nearly industrialized), a bigger prize than single-cell
  circadian methods (small niche, low funding, ~zero skill transfer), and by
  far the best skill-transfer match to what we already do (K*/partition-
  function/ensemble machinery).
- Industry pull is real, not just academic: Microsoft Research AI4Science's
  BioEmu (Science 2025) is exactly "generative model for protein equilibrium
  ensembles," trained on ~200ms aggregate MD + AlphaFold DB + experimental
  stability data; Isomorphic Labs' IsoDDE (announced 2026) explicitly frames
  "beyond AlphaFold" as including conformational/dynamic modeling; condensate
  biotech (Dewpoint $287M w/ Bayer/Merck/Pfizer partnerships, Nereid $50M;
  Faze $81M but has since folded — real money, but a real failure rate too)
  needs exactly this kind of IDR/IDP ensemble tooling for "undruggable"
  targets (MYC, p53 TAD, AR-NTD, tau, alpha-synuclein, TDP-43/FUS).
- Given Bruce forbids MD (any flavor, including coarse-grained MD like
  CALVADOS — TBD how strict, see Open Questions), the CALVADOS/bAIes/IDPFold2
  lineage is partly or fully off-limits as *our own* method, but this is not
  a dead end: it pushes us toward exactly the algorithmic family we already
  own (branch-decomposition DP, provable/certified estimation), which is also
  the one differentiator none of CALVADOS/bAIes/X-EISD/IDPFold2 currently
  have (a rigorous finite-sample certificate, not a Bayesian point estimate
  or a black-box generative sample).

## 2. Problem statement

An IDP/IDR does not have one dominant conformation. Any solution NMR/SAXS/
smFRET observable is a population average over an entire conformational
ensemble:

```text
<O> = E_p[O(c)] = [ Σ_c exp(-E_t(c)/RT) O(c) ] / [ Σ_c exp(-E_t(c)/RT) ]
```

where `O(c)` is a back-calculation function (predicted chemical shift / SAXS
intensity at a given q / RDC / etc. for a frozen conformation `c`), and
`E_t` is the explicitly declared target energy/effective potential. It must
not be called the unknowable "true energy"; changing this target changes the
ensemble whose expectation is being estimated.

This is a severely underdetermined inverse problem: tens to a few hundred
independent scalar observables vs. an astronomically large conformational
space. No amount of clever reweighting recovers information that was never in
the candidate pool (garbage-in-garbage-out on the proposal distribution), and
multiple distinct ensembles can fit the same averaged data equally well
(non-identifiability). Better priors (physics-based CG force fields, AI
structural priors) narrow the search but do not remove the underdetermination
— this is why the problem is still explicitly called "unsolved" in 2026
(Nat Methods "Toward a unified framework...", Nat Comm "...remains a major
challenge").

**Operational identifiability rule.** COHERE has three distinct uses of
observables, and reports must say which one is active:

- a proposal can be trained only for sampling efficiency (for example, the
  CA-validity classifier in §10.41), without fitting any experimental
  observable;
- when experimental observables enter a likelihood, maximum-entropy
  restraint, or population-weight fit, the result is a prior-conditional
  posterior/refined ensemble. The defensible claim is reproduction or
  prediction of observables under the declared prior, not recovery of a
  unique microscopic "true ensemble";
- when an observable is held out from inference, it is an external model
  diagnostic. Monte Carlo error bars quantify estimator uncertainty under
  the target; they do not by themselves certify force-field correctness,
  experimental error, or back-calculator bias.

**Key reframe (already worked out in discussion): we are NOT trying to
estimate a partition function `Z` as an end in itself** (unlike K*, where the
ratio of partition functions IS the deliverable). `Z`'s absolute value has no
independent physical meaning (arbitrary energy zero) and is not
experimentally accessible. **What we actually want to estimate and certify is
`<O>` for each observable `O`** — a normalized ratio, directly comparable to
one real experimental number — or, one level up, a **ratio of partition
functions between two macrostates** (bound/ordered vs. unbound/disordered
ensemble) if/when this feeds a ddG calculation (see §6).

## 3. Historical anchor: Wang & Donald 2006 ([10])

"A data-driven, systematic search algorithm for structure determination of
denatured or disordered proteins" (Comput Syst Bioinformatics Conf 2006).
Formulated disordered-protein structure determination as computing an
*ensemble* directly from sparse restraints (mainly RDC): local backbone
conformations consistent with the data were solved as low-degree polynomial
systems in closed form (same algebraic-geometry machinery as the folded-
protein backbone paper, [15]), combined with systematic search + pruning.
Validated on 2 real proteins (ACBP, eglin C). Never became a widely-used tool
(CSB conference paper, no lineage of follow-on software the way NVR or BWM*/
K* had) — but the *goal* (ensemble from sparse data, no MD, exact/systematic
search rather than stochastic sampling) is exactly right and exactly the gap
that 2025-2026 methods (X-EISD, CALVADOS, bAIes, IDPFold2) still have not
filled with anything that is both MD-free AND carries a rigorous finite-
sample guarantee.

## 4. Current SOTA (2025-2026) and where each falls short of "MD-free +
   certified"

| Method | Candidate generation | Uses per-target expt. data? | MD? | Rigorous finite-sample guarantee? |
|---|---|---|---|---|
| X-EISD | none (reweights an external pool) | yes, 8-9 data types, Bayesian MLE | depends on pool source | no (Bayesian point estimate + uncertainty, not a PAC/(ε,δ) certificate) |
| CALVADOS / AF-CALVADOS | coarse-grained MD (Langevin dynamics) | no (blind, transferable force field) | **yes** — excluded by Bruce's constraint | no |
| bAIes | AF2 prior + MD force field | yes (chem. shift + SAXS reweighting) | yes (MD-based generation step) | no |
| IDPFold / IDPFold2 | end-to-end diffusion/flow-matching on AF3 arch. | no (blind, trained network) | trained on ~MD-derived / AF-DB data | no |
| STARLING (Nature 2026, added after user pushback on 2026-07-03) | VAE + DDPM (diffusion), trained on coarse-grained MD | yes, via Bayesian maximum-entropy (BME) reweighting post-hoc (SAXS/smFRET) | yes (trained on CG-MD) | no (BME is a Bayesian point-estimate/credible-interval refinement, not a PAC/(ε,δ) certificate) |
| ENSEMBLE / Flexible-Meccano+ASTEROIDS / EOM | statistical coil sampling (NOT dynamics) | yes | **no** — MD-free | no (genetic-algorithm / MC selection, heuristic) |
| **Gap** | | | | **nobody found so far gives a certified (ε,δ) bound on ensemble-vs-experiment agreement — but see §8 caveat: web search is not a literature review, verify with Scholar + Bruce before committing** |

The old statistical-coil-sampling family (ENSEMBLE, Flexible-Meccano, EOM) is
already MD-free, which is encouraging, but methodologically pre-dates any
rigorous concentration-inequality treatment — closest in spirit to what we
should upgrade.

## 5. Proposed approach: COHERE-IDP proposal portfolio + certified
   observable estimation

### 5.1 Why this is a natural fit, not a stretch

PACK*'s core trick: a branch-decomposition DP over a pairwise-decomposable,
bounded-branch-width energy graph is *simultaneously* an exact normalizer and
an exact i.i.d. ancestral sampler (Lemma sample-correct in pac_draft). This
fact is old (Koller & Friedman 2009 for PGMs generally; McCaskill 1990 +
Ding & Lawrence 2003 stochastic traceback for RNA) — PACK*'s actual novelty
is using it to do an *exact-reference* free-energy-perturbation-style
correction from a cheap surrogate energy model to an expensive true energy
model, which sidesteps the unknown-normalizer problem of ordinary FEP/BAR/
MBAR. Nothing about that trick is specific to *side-chain* rotamers on a
*fixed backbone* — it applies equally to *backbone dihedral* "rotamers" on a
flexible/disordered chain, provided the interaction graph has bounded
branch-width.

Backbone dihedral rotamer library + short-range (local-in-sequence) pairwise
statistical potential gives bounded branch-width close to "for free": a
1-D chain with only nearby interactions is close to a path graph
(branch-width ~ interaction range), unlike a folded protein's side-chain
packing graph which needs the dual-cutoff sparsification of [jain2017critical]
to get bounded branch-width.

**2026-07-07 update:** that short-range/local proposal should no longer be
treated as the main proposal by itself. A COHERE-IDP toy overlap audit showed
that local/statistical-coil proposals collapse on blocky charge-patterned
chains, exactly the long-range IDP physics this project must handle. The main
route is now a **known-logQ proposal portfolio**: branch-DP when the graph is
low-width, WMB-IS when long-range pairwise terms would blow up DP tables, and
tempered/compactness proposal mixtures chosen by pilot diagnostics. The
essential requirement is not "exact DP everywhere"; it is i.i.d. samples from
a proposal `q(c)` with computable `log q(c)` and adequate overlap with the
target ensemble.

### 5.2 Architecture (revised after toy overlap audit)

1. **Discrete IDR conformational state space.** Discretize backbone φ/ψ (and
   relevant χ) into a rotamer-like library per residue. The first target
   should be short functional IDRs / MoRF-like segments, roughly 30-70
   residues, not proteome-scale long-IDP generation.
2. **Known-logQ proposal portfolio.** Build several MD-free proposals over
   the same discrete space:
   - local/statistical-coil branch-DP proposal when only short-range terms
     are used and branch-width stays small;
   - WMB-IS proposal when long-range charge/contact pairwise terms are needed
     and exact branch-DP tables would explode;
   - independent tempered WMB mixtures, e.g. `β = 0, 0.25, 0.5, 0.75, 1`,
     so sampling does not jump from a local proposal directly to the full
     long-range target;
   - optional compactness/contact surrogate components, included only when
     pilot ESS/tail diagnostics show they help.
3. **Pilot diagnostics / proposal gating.** Do not hard-code all proposal
   tricks together. Use small pilot batches to measure ESS/N, largest-weight
   share, weight-range diagnostics, and certified-interval tightness. Keep
   proposal components that improve overlap; drop components that dilute
   probability mass or loosen weight bounds.
4. **CCD-style local refinement and correction.** For retained proposal
   samples, refine with OSPREY CCD/EPIC (deterministic coordinate descent, not
   dynamics), score with the best available all-atom energy, and optionally
   learn η/residual corrections as in PACK*. The correction machinery is not
   DP-specific; WMB can run on corrected energy matrices too.
5. **Certified observable estimation.** For each experimental observable
   `O(c)`, estimate either `E_q[w(c) O(c)] / E_q[w(c)]` or the explicitly
   scoped expectation under the chosen corrected proposal. The certification
   target must be written honestly as a ratio/self-normalized-IS problem when
   residual true-energy reweighting is used; it is not automatically the
   simpler PACK* single-mean case.

### 5.3 What this actually delivers, and why it's differentiated

Output: for a given IDR/IDP sequence and explicitly declared target, a
**certified `(1-δ)` interval for each target-ensemble expectation
`E_p[O]`**, when the required deterministic bounds are valid. Comparing that
interval with a real measurement is a separate model-checking step that must
also account for experimental uncertainty and back-calculator error. The
certificate controls sampling/estimation error; it does not certify that the
target is the unique physical ensemble. This is the precise differentiated
capability relative to point-estimate reweighting or un-diagnosed black-box
samples.

## 6. Downstream application / bridge back to the main thesis project

SKEMPI2 contains a minority of coupled-folding-binding entries (extreme case:
ACTR/NCBD, mutual synergistic folding — both partners disordered alone,
fold on binding). Standard rigid-backbone ΔΔG pipelines (including our
current multi-backbone K* / PackStar path) cannot represent the entropic
cost of ordering an unbound-state ensemble, because they only ever see the
bound, folded structure. **Action item before committing engineering time to
this whole plan:** screen our existing ~870 SKEMPI hard/failing targets with
a disorder predictor (IUPred2A/3, MobiDB, or cross-reference against DisProt
and/or the new IBPC-Kd atlas of IDR-ordered complexes, 2026) to check whether
failing cases are enriched for disorder / coupled-folding-binding character.
If yes, this becomes a concrete Milestone 4: use the certified `<O>`-matched
ensemble to estimate the free-energy cost of the disorder→order transition
(this is the one place a partition-function *ratio* — `Z_bound` vs.
`Z_unbound-ensemble` — legitimately re-enters, see §2) as an additive
correction term feeding directly into the existing multi-backbone K*
aggregation contract in `restart_plan.md`. This is the highest-value, lowest-
narrative-risk way to justify the whole direction to Bruce: it's not a new
side project, it's a fix for a diagnosed failure mode of the main thesis
project.

## 7. Milestones

- **M0 — Proposal-overlap go/no-go.** Mostly done in `/home/users/lz280/
  COHERE-IDP/` as of 2026-07-07. The toy matrices diagnosed the main risk:
  local-only proposals under-cover long-range charge/contact-driven compact
  states, while tempered WMB mixtures are the best non-cheating direction.
- **M1 — Fast real-observable benchmark, not OSPREY-first.** Use a real IDP/IDR
  sequence and real SAXS/NMR data to establish the scientific loop outside the
  OSPREY bottleneck. First staged case: SASBDB `SASDXC6`, human NHE6
  C-terminal residues G586-A701, real q/I/sigma curve and FASTA. Immediate
  deliverable: backbone-diverse ensemble PDBs -> external SAXS back-calculator
  -> chi2/likelihood -> reweighting/ESS/stability report. **Status
  2026-07-19:** the CA-level ensemble-average loop, matched-QC three-seed
  SASDXC6 analysis, and a preregistered independent SASDNV6/E1A FoXS case are
  complete. Unrestricted MaxEnt failed held-out prediction in both proteins
  and did not pass the advancement rule under any completed
  forward-model/pool diagnostic. The E1A `saxs_md` sensitivity is complete;
  production all-atom/solvent validation and broader biological-case coverage
  remain missing.
- **M2 — IDP WMB adapter.** Reuse WMB as a generic discrete factor-graph
  sampler/certifier, but do not force the IDP state space through OSPREY
  rotamer RCs. Build an adapter for backbone bins, local Ramachandran factors,
  long-range charge/contact factors, exact/evaluable `log q(c)`, and
  state-to-PDB materialization. **Status 2026-07-18:** the generic sparse
  log-factor `WmbModel`, tempered-portfolio overload, discrete backbone-bin
  mapping, deterministic backbone PDB materializer, hard-support pilot job,
  and first SASDXC6 learned steric-surrogate experiment have landed with
  targeted tests. The current target is still a coarse three-state
  Ramachandran base measure with CA hard support; physical target-energy and
  SAXS validation of these WMB samples remain open.
- **M3 — Certification + production observable wrappers.** Add deterministic
  observable bounds and empirical/PAC ratio intervals after the real-observable
  loop is stable. Harden SAXS with FoXS/CRYSOL/Pepsi-SAXS or validated
  `saxs_md` solvent/blank handling; add SPARTA+/SHIFTX2 only after SAXS works.
  The OSPREY bridge remains useful for regression/support tests and possible
  fixed-backbone subproblems, but is no longer the fastest path to the first
  real IDP benchmark.
- **M4 (stretch, conditional on the §6 screen coming back positive)** — plug
  the certified ensemble into the ΔΔG pipeline for disorder-enriched SKEMPI
  hard cases; show improvement over the current rigid-backbone baseline.

## 8. Open questions to resolve with Bruce before committing

- **Exact scope of "no MD."** Does it forbid only *us* running dynamics
  (atomistic or coarse-grained), or also forbid *using* black-box tools
  whose training data included MD (CALVADOS as a comparison baseline is
  presumably fine to cite/benchmark against; using a CALVADOS-family model
  as a component of *our own* pipeline is the ambiguous case; ditto for
  treating IDPFold2 samples as an input prior)? This changes whether M3's
  baseline comparisons are "compare against" only, or whether some hybrid
  use is acceptable.
- **Novelty positioning risk (carried over from PACK* itself, will apply
  here too):** must explicitly cite and distinguish from multi-fidelity
  Monte Carlo / control-variate literature (cheap correlated surrogate +
  few expensive oracle evals + variance reduction) — same caveat flagged
  for the PACK* draft's related-work section, doubly relevant here since
  the "surrogate vs. true energy model" framing is even more exposed to
  this comparison in an IDP context.
- **PAC terminology precedent** — same resolution as PACK*: cite ApproxMC
  ("provides PAC (probably approximately correct) guarantees",
  Chakraborty-Meel-Vardi lineage) and PAC-MDP/PAC-bandit literature as
  precedent for using "PAC" outside strict Valiant-style hypothesis
  generalization; standard PGM-inference papers don't use the word but that
  reflects community habit, not a terminological rule.
- **Naming. Resolved for now:** use **COHERE-IDP**. "PACK*" is already a
  homophone of our own internal nickname "PackStar" for the same underlying
  system, so avoid minting a third confusable "pack*"-family name. The new
  name foregrounds the capability: certified observable-guided heterogeneous
  ensemble refinement for IDPs.

## 9. Key references to build from

- Wang & Donald 2006, "A data-driven, systematic search algorithm for
  structure determination of denatured or disordered proteins" — direct
  historical anchor, [10] in the Donald-lab reading list this plan grew out
  of.
- Koller & Friedman, *Probabilistic Graphical Models* (2009) — exact
  ancestral sampling from a calibrated clique tree, the general fact PACK*
  already cites.
- McCaskill 1990; Ding & Lawrence 2003 — RNA secondary-structure DP +
  stochastic traceback, the closest prior "DP = normalizer + sampler"
  precedent.
- Ermon et al., "Discrete Integration by Hashing and Optimization" (WISH,
  ICML 2013) and Chakraborty-Meel-Vardi ApproxMC lineage — PAC-style
  guarantees for sums/counts via hashing + a cheap oracle; opposite regime
  from PACK*/this plan (their sum is the hard part; ours is the per-
  conformation oracle).
- X-EISD (Bayesian multi-data-type reweighting for IDP ensembles).
- Tesei et al. 2024 (CALVADOS, human IDP proteome ensembles); AF-CALVADOS
  2025 (multi-domain extension).
- bAIes (AlphaFold2-informed Bayesian ensemble reweighting, Nat Comm).
- IDPFold / IDPFold2 (MoE-in-AlphaFold3-diffusion, flow matching, unified
  folded+disordered+multidomain generation).
- STARLING (Novak, Lotthammer, Emenecker & Holehouse, Nature 652:240-250,
  2026; Holehouse lab, WashU) — VAE+DDPM trained on
  coarse-grained MD, seconds-scale IDR ensemble generation from sequence
  alone, refined against SAXS/smFRET via Bayesian maximum-entropy (BME)
  reweighting. Newest, highest-profile competitor found so far; still in the
  "MD-trained + Bayesian point-estimate" bucket, not PAC-certified — but
  found only on a second search pass (2026-07-03), so treat the "gap" claim
  in §4 as provisional until independently verified via Scholar/Bruce, not
  as an established fact.
- IDPEnsembleTools (IDPET, 2026) — standardized ensemble comparison
  (Jensen-Shannon divergence), signals the field still lacks a shared
  ground truth / benchmark convention.
- BioEmu (Microsoft Research AI4Science, Science 2025) — industry precedent
  that "generative ensembles for proteins" is a real, funded, top-venue
  industrial research direction, adjacent but not identical (folded-protein
  dynamics / cryptic pockets, not disorder per se).
- SKEMPI2 (Jankauskaite et al. 2019); IBPC-Kd atlas of IDR-ordered protein
  complexes (2026) — data sources for §6.

## 10. External review + engineering feasibility audit (2026-07-06)

Status: informal review from a planning conversation (literature spot-check,
codebase archaeology, scientific-risk critique), not yet shown to Bruce.
Captured here so it isn't lost.

### 10.1 Literature spot-check: holds up

Every citation checked (STARLING incl. exact *Nature* 652:240-250 (2026)
volume/pages, Wang & Donald 2006 CSB abstract, CALVADOS/AF-CALVADOS, bAIes,
BioEmu, Isomorphic IsoDDE, X-EISD, IDPFold2, IDPEnsembleTools, IBPC-Kd atlas)
is real and accurately characterized (MD-or-not, per-target-data-or-not,
certified-or-not). No fabricated or mischaracterized reference found — a
meaningfully strong signal for the §1/§4 landscape read.

### 10.2 Core scientific risks not resolved by the current draft

1. **Short-range-only surrogate vs. long-range-driven IDP physics.** Phase 1's
   bounded-branch-width trick works *because* the surrogate graph is
   restricted to sequence-local pairwise terms — but real IDP dimensions/
   compaction are frequently long-range-electrostatics-driven (charge
   patterning, Das-Pappu/CALVADOS lineage), and a genuinely flexible backbone
   can't be given a single static long-range interaction graph the way
   dual-cutoff sparsification does for a fixed folded backbone (which pairs
   are spatially close is conformation-dependent). §2's own critique of
   competitor methods — "no amount of reweighting recovers information never
   in the candidate pool" — applies to this surrogate too: if `p_m`
   structurally can't propose long-range-contact conformations, no per-term
   `η` can conjure them. **M3's own validation target, ACTR/NCBD, is likely
   the case most exposed to this gap** (chosen precisely for mutual
   synergistic folding / residual long-range structure in the free state).
   Action: pressure-test before M1's DP engineering. This was started on
   2026-07-07; see §10.6. Result so far: the short-range-only proposal fails
   on blocky charge-patterned chains and should not be the main proposal.
2. **Phase 5's `<O>` estimator semantics are ambiguous, and "simpler than
   PACK*" may not hold.** Two readings of §5.2.4 give different, both-
   important problems: (a) if `<O>` is reweighted back to `E_t` via
   `w=exp(-ξ/RT)` (self-normalized importance sampling / ratio-of-two-
   expectations), plain empirical-Bernstein does not directly apply — needs
   weighted-IS / ratio-estimator concentration machinery (cf.
   Chatterjee-Diaconis on effective-sample-size collapse), which is *harder*
   than PACK*'s single-mean case, not simpler; (b) if `<O>` is just
   `E_{p_η}[O]` with no residual reweighting, the empirical-Bernstein bound
   only certifies **Monte-Carlo/sampling precision**, not **model
   accuracy** — i.e. it does not certify "matches real measurements" as
   worded in §5.3, only "we know our estimate under `p_η` to within ±ε."
   Needs to be nailed down on paper — which reading, and is the resulting
   bound honestly scoped — before M3.
3. **Back-calculator model error is unaccounted for.** SPARTA+/SHIFTX2 are
   known to be less reliable on disordered/coil regions (trained mostly on
   folded-protein data) — X-EISD explicitly models this as a nuisance term;
   the current Phase-5 design doesn't, which understates the uncertainty
   that actually matters for the "match to experiment" claim.
4. **Backbone+side-chain joint discretization is a real added complexity,
   not a parenthetical.** Side chains are not dropped (Phase 1's library
   explicitly includes χ; Phase 2's CCD is all-atom; back-calculators need
   atomic detail) — but combining backbone bins × side-chain rotamers
   multiplies per-node state count (inflates TESS/compute directly,
   independent of branch-width), and χ distributions are backbone-dependent
   (Dunbrack-style), so backbone-bin and χ choice are coupled within a
   residue, not separable — and the existing backbone-dependent rotamer
   statistics (`BackboneDependentRotamerLibrary.java`) were fit on
   folded-protein data; transfer to disordered-range φ/ψ is untested.
5. **The deliverable is a certified scalar, not a structure ensemble** —
   worth stating explicitly so reviewers don't benchmark this against
   STARLING/CALVADOS on "ensemble quality" (a comparison it would likely
   lose, given #1). The actual product is `<O>` (or a ΔΔG-feeding
   partition-function ratio for M4) with a rigorous interval — narrower and
   more specialized than a general-purpose structure generator.

### 10.3 Code-reuse audit (OSPREY3 checkout, 2026-07-06)

More is already built than the plan assumes needs building:

| Component | Status | Location | Size |
|---|---|---|---|
| Branch decomposition + exact tree DP + GPU ancestral sampler | exists, mature | `branchdp/` (19 files) | ~12,100 lines, incl. GPU kernels (`DPGpuFullDP`, `SamplingGpuPhase1`) |
| PackStar partition-function layer | exists | `packstar/` | ~3,600 lines |
| Weighted mini-bucket proposal + WMB-IS state bounder | exists; generic sparse IDP adapter started | `wmb/`, `cohere/IdpBackboneFactorGraph.java`, `spackstar/bound/WmbImportanceSamplingStateBounder.java`, `docs/PackStar_WMB_Importance_Sampling_Plan.md` | `WmbModel.fromLogPotentials(...)` now accepts sparse non-OSPREY unary/pair factors, retains exact WMB proposal `logQ`, and supports tempered scaling. `IdpBackboneFactorGraph` maps discrete φ/ψ states to a deterministic N/CA/C/O PDB. The remaining hard part is a physically justified SASDXC6 factor construction, not sampler plumbing. |
| η/residual correction learning (ridge regression on CCD residuals — literally PACK* Phase 2-3) | exists | `energy/approximation/branch/` | ~1,900 lines |
| Backbone-dependent rotamer library data structure | exists (folded-protein-fit statistics) | `restypes/BackboneDependentRotamerLibrary.java` | — |
| CCD minimization, EPIC | exists, production | `minimization/`, `ematrix/epic/` | — |
| Empirical-Bernstein / bounded-ratio certificate (Phase 5) | partial, fail-closed assumptions implemented | `cohere/RatioConfidenceInterval.java`, `cohere/WmbProposalPortfolio.java`, `cohere/WmbPilotRunner.java` | A bounded self-normalized-ratio interval is available only when deterministic observable bounds and a valid log-weight cap are supplied. This is not yet a useful real-SAXS certificate: physical observable bounds and a defensible target/proposal weight cap are still missing. |
| Back-calculator wrappers (SPARTA+/SHIFTX2, SAXS tools) | partial | `slurm/wrappers/cohere_saxs_md_batch_observable.py`, `slurm/wrappers/cohere_foxs_batch_observable.py`, `slurm/wrappers/cohere_multifoxs_baseline.py` | `saxs_md` and IMP FoXS 2.24.0 now share one q-grid/chi2 convention; FoXS completed both the original 768-frame cross-check and the independent 1,536-frame matched-QC holdout. A QC-filtered MultiFoXS sparse/NNLS baseline is recorded. SPARTA+/SHIFTX2/CRYSOL/Pepsi-SAXS and production all-atom solvent/hydration validation are still missing. |

Git history on `branchdp`/`packstar`/`spackstar` shows only 2 commits — this
was almost certainly imported from lab/collaborator work rather than built in
this checkout, so it predates and is independent of this plan; it is
production code, actively used by the sibling multi-backbone K* project on
grisman A5000s.

### 10.4 Rough code-volume / timeline estimate

This was the original v1 estimate before the 2026-07-07 proposal-overlap
audit. Treat it as an upper-level engineering estimate only; it is **not** a
recommendation to start M1 immediately. M0 proposal-overlap validation now
gates the Java/OSPREY build.

New code is concentrated in the backbone-specific graph/library layer and the
back-calculator/certification layer, not in the DP/WMB engine itself:

- **M1** (backbone interaction-graph builder + discretized backbone/χ
  library seeded from statistical-coil data + smoke test): ~1,500-3,000 new
  lines (Java) + data-prep scripts; **2-4 weeks**.
- **M2** (wire CCD + reuse/extend existing η ridge-regression pipeline):
  ~500-1,500 lines; **1-2 weeks**.
- **M3** (back-calculator wrappers + PAC certificate + benchmark vs X-EISD/
  CALVADOS): ~1,500-2,500 lines + benchmark wall-clock; **3-5 weeks**.
- **M4** (stretch, ΔΔG bridge into existing multi-backbone K* aggregation
  contract): ~500-1,000 lines, gated on §10.5's disorder screen.

Order of magnitude: **~5,000-8,000 new lines total** to reach M3's "minimum
publishable unit," on top of ~17,500 lines of already-existing reusable
DP/PackStar/correction infrastructure — roughly **2-3 calendar months** for
one person, assuming risks #1-2 in §10.2 don't force a redesign. This is a
rough order-of-magnitude estimate, not a committed schedule.

### 10.5 Fast go/no-go status

Two cheap, decoupled checks were recommended before writing any M1 code:

- **Disorder screen** (§6's own action item): run IUPred2A/MobiDB over the
  existing ~870 SKEMPI failing targets. Pure scripting against data already
  on disk (`/usr/xtmp/lz280/packstar_skempi/`). **~1 day.** Decides whether
  M4's thesis-bridge narrative (the strongest pitch to Bruce) has any real
  support.
- **Toy physics-adequacy prototype** (new, motivated by risk #1): outside
  OSPREY entirely, in Python — sample from a short-range-only local
  statistical potential (or even just independent per-residue
  statistical-coil sampling), importance-reweight against a known
  long-range-driven benchmark (Das-Pappu-style charge-patterning Rg series,
  and/or ACTR/NCBD's known transient long-range contacts), and check whether
  reweighting alone can recover the known behavior. **~3-7 days.** This
  directly falsifies/validates the single riskiest design bet (short-range
  surrogate + reweighting ⇒ adequate for real IDPs) before any DP
  infrastructure is built. **Status 2026-07-07: done for the first toy
  model; short-range/local proposal failed on blocky charge-patterned chains.
  Follow-up WMB/mixture matrices also completed; they support tempered
  mixtures at L32/L48 but not a healthy L64/L70 regime (§10.6).**

Current recommendation, updated by §10.24: do **not** make the Java/OSPREY
build the M1 critical path. Start with a real-observable benchmark outside
OSPREY, then bring WMB back through a purpose-built IDP backbone factor-graph
adapter once the SAXS/NMR loop is scientifically grounded.

### 10.6 COHERE-IDP proposal-overlap audit (2026-07-07)

All work in this section is toy-level and outside OSPREY, under
`/home/users/lz280/COHERE-IDP/`. It is designed to test proposal overlap, not
to produce biologically realistic ensembles.

**Project naming / scope.** The idea is now called **COHERE-IDP**. The intended
first biological scope is short functional IDRs / MoRF-like segments, roughly
30-70 residues, not a general long-IDP ensemble generator.

**Continuous bead-chain stress test.** Script: `overlap_probe.py`; Slurm array
`cohere_idp_stress_array.sbatch`, job `12009618`; results aggregated in
`results_slurm_aggregate/`.

- Alternating charge patterns behaved well: at N=100000, median ESS/N stayed
  high even as length increased (`L32≈0.85`, `L48≈0.76`, `L64≈0.68`,
  `L96≈0.54`).
- Blocky charge patterns collapsed badly under local/random-chain proposals:
  at N=100000, median ESS/N was `L32≈2.7e-3`, `L48≈2.0e-4`,
  `L64≈5.1e-5`, `L96≈1.6e-5`; median largest-weight share rose from
  `~0.027` at L32 to `~0.78` at L96, with p90 `~0.98` at L96.
- Interpretation: the original short-range/local proposal + long-range
  reweighting route is a red light for charge-patterned IDRs. Important
  compact/contact-rich states are too rare under the local proposal.

**Discrete pairwise WMB probe.** Script: `wmb_pairwise_probe.py`. This is a
finite-domain toy where WMB can include dense long-range charge pairwise terms
and return sample `logQ`.

- Single WMB helped in some short cases but was not enough for L48+.
  `iBound=1` often underfit the coupling and could be worse than local.
- Tempered WMB mixtures were the strongest signal. The completed aggregate
  matrix in `results_matrix_latest/` shows:
  - L32: best median ESS/N `~0.144`, median max-weight share `~0.031`,
    using 5-temperature WMB mixture (`β=0,0.25,0.5,0.75,1`),
    interaction-aware partition, block-interleave ordering, `iBound=3`.
  - L48: best median ESS/N `~0.029`, median max-weight share `~0.111`,
    same best profile. This clears the provisional minimum line
    (ESS/N ≳ 0.01) but is not yet a comfortable regime.
  - L64: best median ESS/N is `~0.00886`, with median max-weight share
    `~0.257`; this is below the provisional line.
  - L70: best median ESS/N is only `~0.00601`, with median max-weight share
    `~0.278`; increasing chain length did not establish a usable regime.
- Interaction-aware partition and block-interleave ordering can help, but
  neither is monotone. They should be treated as portfolio components selected
  by pilot diagnostics, not guaranteed improvements.
- Compactness surrogate mixtures can help some configurations but can also
  hurt. Example from partial results: at L32 with interaction partition,
  block-interleave ordering, `iBound=3`, adding compactness to a strong
  tempered mixture reduced median ESS/N from `~0.123` to `~0.061`. Do not
  blindly stack all proposal tricks.

**Explicit latent / stratified proposal tests.** Script:
`stratified_latent_probe.py`; Slurm job `12009798`. Smoke testing at L16 gave
`local ESS/N≈0.158`, latent-mixture `≈0.389`, stratified-span `≈0.372`.
That improvement did not survive the completed four-seed L48 aggregate:
local median ESS/N was `0.0119`, the best latent profile was `0.00879`, and
the best explicitly stratified profile was `0.00674`. The stratified sampler
also required about `4.3` base draws per retained sample. No L64/L70 summaries
were produced, so those lengths remain untested rather than pending.

**Current algorithmic conclusion.**

- Red light: local/statistical-coil proposal alone.
- Red/yellow: single WMB alone for L48+.
- Yellow/green: independent tempered WMB mixture, especially 5-temperature
  mixture + interaction-aware partition + block-interleave + `iBound=3`.
- Red/yellow: the tested 15-component compactness latent and four-bin
  pilot-mass stratification did not beat local sampling at L48.
- Practical target: L32 is feasible in this toy; L48 is plausible only with
  the best tempered portfolio; L64/L70 remain below the overlap gate.

### 10.7 Proposal ideas under test / to keep in the matrix

These ideas are not mutually exclusive by default. The correct framing is a
**proposal portfolio with pilot gating**: keep components that improve overlap
and certificate tightness; drop components that dilute probability mass or
make the weight tail worse. Several combinations have already shown negative
interactions, so no component should be enabled unconditionally.

| Idea | What it tests | Compatibility with COHERE-IDP | Current status |
|---|---|---|---|
| **1. Tempered WMB mixture** | Independent mixture of WMB proposals at different long-range-energy strengths, e.g. `β=0,0.25,0.5,0.75,1`; avoids one-shot reweighting from local to full long-range target. | Fully compatible if each component gives exact `log q_beta(c)` and the mixture `log q(c)` is evaluated by log-sum-exp. Avoid MCMC/AIS transitions unless mixing can be certified. | Implemented in `wmb_pairwise_probe.py`; strongest completed toy signal. Best median ESS/N is L32 `~0.144`, L48 `~0.029`, L64 `~0.0089`, and L70 `~0.0060`. |
| **2. Compactness/contact surrogate mixture** | Adds proposal components biased toward compact/contact-rich states, a toy proxy for `Rg`/contact coverage. | Compatible as a known-logQ proposal component when expressed as pairwise/low-dimensional factors. Must be gated because it can over-bias into wrong compact modes. | Implemented as `--compact-bias-strengths`; mixed results. Helps some configurations, hurts others. Example: adding compactness to a strong L32 tempered profile reduced ESS/N from `~0.123` to `~0.061`. |
| **3. Interaction-aware WMB partition** | Keeps strong pairwise interactions together in mini-buckets where possible, instead of only scope-size-greedy partitioning. | Fully compatible; it changes WMB proposal quality/cost, not the estimator identity. | Implemented as `--partition-modes interaction`; helps best L32/L48 profiles but is not monotone. Needs better strength/graph heuristics before production. |
| **4. Better WMB variable ordering** | Tests whether elimination/sample order should reflect charge blocks or long-range couplings. | Fully compatible; standard graphical-model heuristic layer. | Implemented simple `natural` and `block-interleave` orders. Helps in best L32/L48 profile; not uniformly best at L64/L70. Weighted min-fill / charge-aware ordering not yet implemented. |
| **5. True compactness/contact stratification** | Forces equal sampling across compactness/contact bins, then corrects with stratified or multiple-IS weights. This directly attacks missing compact states. | Compatible in principle, but proof is harder because bin probabilities require exact calculation or certified estimates. Toy feasibility can use pilot-estimated bin masses; production certificate needs bin-mass uncertainty accounted for. | `stratified-span` helped the L16 smoke but not four-seed L48: best median ESS/N `0.00674` versus local `0.0119`, at about `4.3` base draws per retained sample. Do not promote this pilot-mass version. |
| **6. Explicit low-rank/global compactness latent proposal** | Introduces a small latent state `z` (expanded/medium/compact) controlling global compactness/contact bias, to cover long-range modes without dense exact DP. | Compatible if `q(c)=Σ_z q(z)q_z(c)` has computable `log q(c)`. This is cleaner for finite-sample certification than rejection stratification. | `latent-mixture` helped the L16 smoke but not four-seed L48: best median ESS/N `0.00879` versus local `0.0119`. A future latent model needs better factors; adding more generic compact components is not supported. |

Important terminology note after the SASDXC6 CA-only work: the
`cohere_sasdxc6_real_benchmark.py` mode named `compact` is not the formal WMB
implementation of ideas 2/5/6.  It is a fast CA-level prototype of the same
scientific hypothesis: compact/contact-biased proposal mass may be needed to
cover SAXS-relevant IDP states.  The 96 x 3 compact-only SASDXC6 result in
§10.37 supports carrying compactness forward into the WMB proposal portfolio,
but also shows the over-compact failure mode that ideas 2/5/6 must control
with mixture weights, stratification, or an explicit latent state rather than
by always sampling the pure compact component.

Current working ranking:

1. Tempered WMB mixture is the main proposal direction.
2. Interaction-aware partition + selected ordering are useful modifiers, not
   guaranteed wins.
3. Compactness surrogate is optional and must be pilot-gated.
4. The tested generic stratification/latent construction is deprioritized;
   the unresolved problem is a physically justified, steric-aware factor
   construction that improves real-target overlap.

### 10.8 First OSPREY implementation landing (2026-07-07)

Decision: start the Java/OSPREY build at the smallest reusable contract, not
with a bespoke back-calculator or a new sampler. The existing high-performance
component to reuse is `WeightedMiniBucket.Proposal`, which already provides
independent samples and exact per-sample `logQ`.

Added package: `src/main/java/edu/duke/cs/osprey/cohere/`.

- `KnownLogQSample<C>`: one proposal sample with `state`, exact `logQ`, and
  unnormalized target score `logTarget`.
- `ObservableBackCalculator<C>`: narrow scalar observable interface. This is
  the intended adapter point for SPARTA+/SHIFTX2/CRYSOL/Debye/etc.; do not
  write a low-performance in-repo replacement unless no suitable component is
  available.
- `ObservableEstimator`: self-normalized IS estimator for `<O>` with ESS/N,
  max-weight share, log-weight range, observed observable range, and an
  explicitly labeled diagnostic bootstrap interval. This intentionally does
  **not** claim a rigorous ratio certificate yet; the future certificate layer
  must add the needed numerator/denominator assumptions.
- `WmbKnownLogQSampler`: bridge from existing
  `WeightedMiniBucket.Proposal.sample()` to `KnownLogQSample<int[]>`, preserving
  exact `logQ` and accepting any target scorer, e.g. `WmbModel::logValue` or a
  future all-atom/CCD residual scorer.

Added tests: `src/test/java/edu/duke/cs/osprey/cohere/TestObservableEstimator.java`.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=1g --max-workers=1 --no-daemon
```

Result: build successful; 5/5 targeted COHERE tests passed. The default test
heap in this checkout is `8g`; in the current environment that caused the test
executor to be killed with exit 137, so COHERE targeted tests should use the
lower heap override unless running on a larger node.

### 10.9 Tempered-mixture and bounded-ratio implementation landing (2026-07-07)

Added the next two reusable Java components under
`src/main/java/edu/duke/cs/osprey/cohere/`.

- `WmbProposalMixture`: exact-logQ mixture over existing
  `WeightedMiniBucket.Proposal` components. This is the production-side form of
  the tempered WMB mixture idea from the toy experiments. It samples by first
  drawing a component, then reuses that component's high-performance WMB
  sampler, and finally evaluates the full mixture proposal probability by
  log-sum-exp across all components.
- `RatioConfidenceInterval`: empirical-Bernstein confidence interval for
  bounded self-normalized importance ratios `E_q[w O] / E_q[w]`. This is only
  exposed when the caller provides deterministic assumptions: observable bounds
  `O in [O_min,O_max]` and a log-weight cap `log w <= logW_max`. Internally it
  rescales weights by `logW_max`, so the certificate works with bounded
  variables and avoids exponent overflow.

Extended `TestObservableEstimator` coverage from 5 to 8 tests:

- Mixture samples preserve exact mixture `logQ`.
- Bounded ratio interval contains a hand-checkable uniform mean.
- Ratio certificate rejects samples that violate the declared weight cap.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 8/8 targeted COHERE tests passed. A prior rerun with
`-DtestMaxHeap=1g` reached the test task but the executor was killed with exit
137 in the current environment; `512m` was stable for this target class.

### 10.10 Proposal portfolio, pilot gate, and Slurm test landing (2026-07-07)

Added the M2 proposal-pipeline skeleton without replacing the existing WMB
sampler.

- `ImportanceWeightDiagnostics`: standalone overlap diagnostics from
  known-logQ samples or raw log weights: ESS/N, max-weight share, and log-weight
  range. This is the cheap pilot signal used to accept/reject proposal
  portfolios before expensive observable back-calculation.
- `WmbProposalPortfolio`: builds exact-logQ mixtures from component specs
  containing weight, i-bound, elimination order, and inverse temperature. The
  tempered path reuses `WmbModel(emat, rcs, assignments, rt / beta)` and
  `WeightedMiniBucket.proposalForModel(..., maxTableCells)`; no custom
  proposal sampler is introduced.
- `WmbPilotRunner`: draws pilot samples from a `WmbProposalMixture`, runs a
  scalar observable through the narrow `ObservableBackCalculator<int[]>`
  interface, and gates the portfolio on ESS/N, max-weight share, and
  log-weight range.
- `slurm/scripts/run_cohere_unit_tests.slurm`: cluster-side targeted unit test
  script using the same `512m`, single-worker Gradle settings that are stable
  in this checkout.

Extended `TestObservableEstimator` coverage from 8 to 16 tests:

- Tempered portfolio samples use exact full-mixture `logQ`.
- Mixture construction rejects same-variable-count proposals with mismatched
  domain sizes.
- Pilot runner accepts a well-overlapped exact WMB portfolio.
- Pilot options reject NaN thresholds.
- Diagnostics report the first failing threshold.
- Portfolio construction rejects malformed elimination orders.
- External-tool observable adapter parses successful process output and rejects
  non-zero tool exits.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 16/16 targeted COHERE tests passed locally.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12010803` ran on `fennario-01` and completed successfully in 15s;
Slurm stdout ended with `BUILD SUCCESSFUL`, stderr was empty.

Back-calculator inventory note: the original environment audit found no
`crysol`, `foxs`, `sparta+`, `shiftx2`, or `pepsi-saxs` executable on PATH,
while the `confdiff` environment contained `saxs_md` and `saxs_rism`. FoXS was
subsequently installed in the isolated `cohere-saxs` environment and integrated
in §10.39; the rest of this paragraph records the original plumbing decision.
`ExternalToolObservableBackCalculator` landed first as generic process
plumbing: it materializes each state into a temporary work directory, executes a
configured command with timeout handling, captures stdout/stderr, and delegates
scalar parsing to a caller-supplied parser. The concrete `saxs_md` q/I wrapper
landed later in §10.22. Actual SPARTA+/SHIFTX2/CRYSOL/FoXS/Pepsi-SAXS
calculation should still be delegated to those high-performance binaries when
they are installed; COHERE should not reimplement those physics internally.

### 10.11 Pilot scan/report landing and stricter proposal compatibility (2026-07-07)

Extended the M2 pilot layer from a single portfolio gate to a reusable
candidate-scan/report path:

- `WmbPilotScan`: runs a list of named WMB proposal-portfolio candidates against
  the same target, observable back-calculator, and pilot options. Each candidate
  owns its component specs and `maxTableCells`, so Slurm arrays can sweep
  i-bounds, orders, inverse temperatures, and mixture weights without parsing
  test logs.
- `WmbPilotReportWriter`: writes pilot scan results as CSV or JSON. The report
  includes pass/fail, failure reason, sample count, observable estimate,
  diagnostic interval metadata, ESS/N, max-weight share, log-weight range,
  observed observable range, and component metadata.
- `WmbProposalMixture`: now rejects proposal mixtures whose WMB models have
  mismatched position or rotamer maps, not just mismatched variable counts or
  domain sizes. This prevents exact-logQ mixture accounting from silently mixing
  semantically different variables.
- `slurm/scripts/run_cohere_unit_tests.slurm`: changed the Gradle command to
  `cleanTest test` so cluster validation actually executes the targeted COHERE
  tests instead of reporting `test UP-TO-DATE` after a local run.

Extended `TestObservableEstimator` coverage from 16 to 19 tests:

- Mixture construction rejects same-domain proposals with mismatched WMB
  position maps.
- Pilot scan runs multiple named candidate portfolios and emits CSV/JSON report
  fields with stable escaping and null JSON interval bounds for absent
  diagnostic intervals.
- Pilot scan rejects duplicate candidate names.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 19/19 targeted COHERE tests passed locally.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12010839` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 19 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.12 Slurm-ready pilot scan job landing (2026-07-07)

Promoted the M2 pilot scan path from library-only code to a runnable job entry:

- `WmbPilotScanJob`: Java `main` entry point driven by a `.properties` file.
  It reads an existing serialized OSPREY `EnergyMatrix`, an optional plain-text
  RC list, optional partial assignments, pilot thresholds, WMB candidate
  portfolios, and CSV/JSON output paths.
- Candidate specs are generated from conservative sweep knobs:
  `candidate.names`, per-candidate/global `iBounds`, `inverseTemperatures`,
  `orders`, `componentWeight`, and `maxTableCells`. The job still builds
  proposals through `WmbProposalPortfolio` and `WeightedMiniBucket`; it does not
  introduce a custom sampler.
- Built-in pilot observables now cover `energy`, `domain_sum`, `domain_value`,
  and `constant` for cheap overlap/debug sweeps.
- External observable mode delegates to the existing
  `ExternalToolObservableBackCalculator`. The job writes
  `domain-values.txt`, `rotamer-values.txt`, and `full-assignments.txt` in each
  work directory and launches a configured command. This is only process
  plumbing; SAXS/NMR/etc. physics should still be delegated to installed
  high-performance back-calculator binaries.
- `slurm/scripts/run_cohere_pilot_scan.slurm`: Slurm wrapper that accepts a
  pilot `.properties` path as an sbatch argument or `COHERE_CONFIG`, builds
  `installDist`, and runs
  `edu.duke.cs.osprey.cohere.WmbPilotScanJob` from the distribution classpath.

Extended `TestObservableEstimator` coverage from 19 to 21 tests:

- Pilot scan job reads relative config paths, serialized `EnergyMatrix`,
  plain-text RCs, candidate sweeps, and writes CSV/JSON reports.
- Pilot scan job delegates observable evaluation to an external command and
  parses the scalar result.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 21/21 targeted COHERE tests passed locally.

```bash
./gradlew installDist --max-workers=1 --no-daemon
java -cp "build/install/osprey/lib/*" \
  edu.duke.cs.osprey.cohere.WmbPilotScanJob
```

Result: `installDist` completed successfully and the distribution classpath
contains `WmbPilotScanJob`; invoking it without a config prints the expected
usage line.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12010863` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 21 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.13 Bounded-ratio certificate/report landing (2026-07-07)

Closed the next M2 gap: pilot scans no longer stop at diagnostic ESS/tail
signals. They can now emit a bounded-ratio confidence certificate when the
observable has deterministic bounds.

- `WmbProposalPortfolio.logWeightUpperBound(...)` derives a conservative global
  cap for the exact-logQ mixture using existing WMB machinery. For target-
  temperature components it reuses
  `WeightedMiniBucket.Proposal.logWeightUpperBound(...)`; for tempered
  components it falls back to target `logZUpper - log q_lower`. The mixture cap
  uses `q_mix >= alpha_i q_i` and takes the tightest component-derived cap.
- `WmbPilotRunner` now evaluates observable values once per sample, then feeds
  the same arrays into both the self-normalized estimate and, when enabled,
  `RatioConfidenceInterval`. This avoids double-calling external
  back-calculators.
- `WmbPilotReportWriter` now writes certificate fields to CSV and JSON:
  estimate, confidence level, lower/upper interval, log-weight cap, denominator
  positivity, cap reason, observable bounds, and the empirical-Bernstein
  numerator/denominator diagnostics.
- `WmbPilotScanJob` supports certificate properties:
  `certificate.enabled`, `certificate.observableLower`,
  `certificate.observableUpper`, `certificate.confidenceLevel`,
  `certificate.logWeightUpper`, `certificate.targetLogZUpper`,
  `certificate.logZIBound`, `certificate.logZOrder`,
  `certificate.logZMaxTableCells`, and `certificate.capMaxAssignments`.
  If no explicit `certificate.logWeightUpper` or
  `certificate.targetLogZUpper` is provided, the job computes target
  `logZUpper` with `WeightedMiniBucket.boundsForModel(...)`.
- `slurm/configs/cohere_pilot_example.properties` documents the Slurm pilot
  config surface, including the external-observable path that delegates SAXS,
  NMR, or other physics to installed high-performance back-calculator binaries.

This still does not implement any SAXS/NMR physics internally. The production
path remains: generate or provide the calculator input from an assignment, call
an existing high-performance tool/wrapper, parse one scalar, and let the WMB
proposal/certificate layer handle the sampling and statistical accounting.

Verification:

```bash
./gradlew compileJava --max-workers=1 --no-daemon
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
./gradlew installDist --max-workers=1 --no-daemon
```

Result: build successful; 22/22 targeted COHERE tests passed locally, and
`installDist` completed successfully.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12010898` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 22 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.14 External observable PDB materialization landing (2026-07-07)

Closed another M3 plumbing gap without implementing SAXS/NMR physics inside
COHERE. External observables can now receive a real OSPREY-generated PDB input
when the pilot config provides a serialized `SimpleConfSpace`.

- `WmbPilotScanJob` now accepts `observable.confSpace` (aliases:
  `observable.confspace`, `confSpace`, `confspace`). The file is read with
  `ObjectIO.read(..., SimpleConfSpace.class)`, so production runs reuse the
  same serialized OSPREY conformational space that produced the energy matrix.
  The job validates that the confspace position count and every assigned/WMB
  rotamer index are compatible before launching the expensive external command.
- For each external-observable sample, the job still writes
  `domain-values.txt`, `rotamer-values.txt`, and `full-assignments.txt`. If a
  confspace is configured, it additionally writes `sample.pdb` by calling
  `SimpleConfSpace.makeMolecule(fullAssignments)` and `PDBIO.writeFile(...)`.
- External commands may now use `{pdbFile}` in addition to the existing
  assignment placeholders. If `{pdbFile}` appears without `observable.confSpace`,
  the job fails early with an explicit configuration error.
- `slurm/configs/cohere_pilot_example.properties` documents the intended
  production path: COHERE materializes the OSPREY sample state, then delegates
  the expensive observable calculation to an installed high-performance binary
  or wrapper such as CRYSOL/FoXS/Pepsi-SAXS or SPARTA+/SHIFTX2.

This is intentionally a bridge, not a new back-calculator. The remaining
production work is to prepare serialized confspace/emat pairs for the IDP
benchmark cases and write thin lab-specific wrappers around the chosen external
tools.

Verification:

```bash
./gradlew compileJava --max-workers=1 --no-daemon
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
./gradlew installDist --max-workers=1 --no-daemon
```

Result: build successful; 23/23 targeted COHERE tests passed locally, and
`installDist` completed successfully.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12010966` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 23 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.15 Batch external observable landing (2026-07-07)

Closed the next M3 throughput gap: external observables no longer have to launch
one process per sample. The single-sample path remains available for simple
wrappers and debugging, but production SAXS/NMR wrappers can now consume an
entire pilot batch in one external call.

- `ObservableBackCalculator` now has a default `values(List<C>)` method. Existing
  lambdas and single-state adapters keep working, while high-throughput
  calculators can override the batch path.
- `ObservableEstimator` and `WmbPilotRunner` now collect sampled states and call
  `calculator.values(...)` once per pilot run. This lets batch-aware calculators
  avoid per-sample process startup and repeated tool initialization.
- Added `ExternalBatchToolObservableBackCalculator`, a generic process adapter
  for external tools that score many states at once. It creates one batch work
  directory, one deterministic `sample-000000` style subdirectory per state,
  runs one external command, parses an ordered vector of scalar values, and
  validates count/finite-value consistency before returning results.
- `WmbPilotScanJob` now supports `observable.batchCommand`. For each batch it
  writes:
  - `manifest.tsv`
  - `sample-dirs.txt`
  - `domain-files.txt`
  - `rotamer-files.txt`
  - `full-assignment-files.txt`
  - `pdb-files.txt` when `observable.confSpace` is configured
- Batch commands may use `{workDir}`, `{manifestFile}`, `{sampleDirsFile}`,
  `{domainFilesFile}`, `{rotamerFilesFile}`, `{fullAssignmentFilesFile}`, and
  `{pdbFilesFile}`. `{pdbFile}` remains single-sample only, and the job now
  rejects it in `observable.batchCommand` with an explicit error.
- `observable.batchResult` accepts `stdout` or `file:<path>` and expects one
  numeric scalar per non-empty line, in sample order. If the number of parsed
  values does not match the sample count, the job fails before reporting an
  estimate.
- `slurm/configs/cohere_pilot_example.properties` now documents the preferred
  production path: a thin batch wrapper around an installed high-performance
  back-calculator, typically driven by `{pdbFilesFile}` or `manifest.tsv`.

This still avoids implementing SAXS/NMR physics in COHERE. The intended
production shape is now: WMB samples assignments, OSPREY materializes assignment
files and optional PDBs, an installed high-performance calculator scores the
batch, and COHERE consumes the returned scalar vector for IS estimates and
certificates.

Verification:

```bash
./gradlew compileJava --max-workers=1 --no-daemon
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 26/26 targeted COHERE tests passed locally. The new
coverage checks the generic batch override path, direct batch external adapter,
and `.properties`-driven `observable.batchCommand` path.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12010986` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 26 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.16 Indexed/table batch-result parsing landing (2026-07-07)

Closed a practical wrapper-integration gap in the batch external-observable
path. Real high-performance back-calculators often emit CSV/TSV tables, headers,
sample ids, or values in tool-defined order rather than a bare ordered vector.
COHERE now handles those formats directly while still delegating the physics to
the external tool.

- `WmbPilotScanJob` batch output parsing now supports:
  - `observable.batchResultDelimiter=auto|whitespace|comma|tab|tsv`
  - `observable.batchResultHeader=true|false`
  - `observable.batchResultValueColumn=<1-based column>`
  - `observable.batchResultIndexColumn=<1-based column>`
  - `observable.batchResultIndexBase=0|1`
- The default remains backward-compatible: `stdout` or `file:<path>` with one
  scalar per non-empty line in sample order.
- If an index column is configured, COHERE validates that each sample index is
  present exactly once, reorders values into sampler order, and fails before
  estimation on duplicates, missing indices, count mismatches, malformed
  columns, or out-of-range indices.
- `slurm/configs/cohere_pilot_example.properties` now documents the table
  parsing options next to the preferred `observable.batchCommand` production
  path.
- `TestObservableEstimator` now includes an end-to-end `.properties` test where
  the same tempered WMB sample batch is scored once in ordered CSV form and once
  in reversed `index,value` CSV form with a header. The two importance-sampling
  estimates must match after parser reordering.

This makes the external-wrapper contract less brittle for tools such as FoXS,
Pepsi-SAXS, CRYSOL, SPARTA+, SHIFTX2, or lab-specific scripts that naturally
return tabular output. COHERE still does not reimplement those calculators; it
only consumes their batch scalar results safely.

Verification:

```bash
./gradlew compileJava --max-workers=1 --no-daemon
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
./gradlew installDist --max-workers=1 --no-daemon
```

Result: build successful; 27/27 targeted COHERE tests passed locally, and
`installDist` completed successfully.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12011033` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 27 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.17 Ratio-certificate numerical-bound tightening (2026-07-07)

Tightened a proof-critical numerical edge in `RatioConfidenceInterval`. The
validator already rejected material violations of the declared deterministic
weight and observable bounds, with a small tolerance for floating-point
roundoff. The interval computation now also clamps in-tolerance roundoff back to
the declared ranges before applying empirical-Bernstein bounds:

- scaled weights use `exp(min(0, logWeight - logWeightUpper))`, so every
  certificate sample stays in `[0, 1]`;
- shifted observable values are clamped to `[0, observableUpper -
  observableLower]`, so the numerator range used by the concentration bound is
  respected;
- clear violations still fail before reporting a certificate.

Added regression coverage with a sample whose log weight and observable are
only `1e-13` above the declared caps. The certificate now accepts that
roundoff-level case and keeps the reported estimate and bounded moments inside
the declared observable/weight ranges.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
./gradlew installDist --max-workers=1 --no-daemon
```

Result: build successful; 28/28 targeted COHERE tests passed locally, and
`installDist` completed successfully.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12011045` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
27s; all 28 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.18 Batch external chunking landing (2026-07-07)

Closed another production-wrapper integration gap: COHERE can now keep the
high-throughput external-observable path while splitting very large pilot
batches across repeated wrapper calls. This is important for tools that are
fast in batch mode but have practical memory, file-count, or wall-time limits
per invocation.

- `ExternalBatchToolObservableBackCalculator` now accepts a `maxBatchSize`
  limit. The default remains one external process for the whole COHERE sample
  batch, preserving the existing high-throughput behavior.
- When a limit is set, the adapter runs consecutive chunks, creates a fresh
  batch work directory for each chunk, validates each returned chunk, and
  concatenates values back into sampler order before importance estimation.
- `WmbPilotScanJob` exposes this through:

```properties
observable.batchSize=256
```

  Aliases `observable.maxBatchSize` and `observable.batchMaxSize` are accepted.
  Indices in `manifest.tsv` and indexed result parsing are chunk-local, which
  matches the external command contract for each invocation.
- `slurm/configs/cohere_pilot_example.properties` now documents the option
  next to the preferred `observable.batchCommand` path.
- Added direct adapter coverage for a 5-sample batch split as `2,2,1`, plus an
  end-to-end `.properties` test proving the `WmbPilotScanJob` config path
  splits an 8-sample external batch as `3,3,2` and preserves the expected
  scalar estimate.

This does not reimplement SAXS/NMR back-calculators. It makes the wrapper layer
more compatible with installed high-performance tools by letting those tools
operate at their preferred batch size.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
./gradlew installDist --max-workers=1 --no-daemon
```

Result: build successful; 30/30 targeted COHERE tests passed locally, and
`installDist` completed successfully.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12011058` ran on `fennario-01`, executed `cleanTest` followed by
the targeted `TestObservableEstimator` suite, and completed successfully in
26s; all 30 COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.19 First real-artifact pilot plumbing (2026-07-07)

Started the first end-to-end pilot that uses materialized OSPREY artifacts
rather than toy matrices:

- Added `CoherePilotArtifactBuilder`, a small CLI that builds a
  `SimpleConfSpace`, computes a real OSPREY `EnergyMatrix` with the existing
  `SimplerEnergyMatrixCalculator`, writes an all-RC `rcs.txt`, and emits
  metadata for a `.properties` pilot scan. This keeps energy and conformation
  handling inside existing OSPREY components.
- Added `slurm/wrappers/cohere_rg_batch_observable.py`, a dependency-free
  batch external observable wrapper. It consumes COHERE-generated sample PDBs,
  computes mass-weighted radius of gyration, compares to a configured
  experimental/reference Rg, and returns indexed CSV (`index,rg,target_rg,chi2`)
  for the existing batch parser.
- Made the OSPREY objects needed by serialized `SimpleConfSpace` artifacts
  serializable: AMBER forcefield parser records and `MutAlignmentCache`.
  `MutAlignmentCache` now treats stored alignments as a transient cache and
  rebuilds it after deserialization.
- Staged a 1I50 chain-L IDR-like pilot config, but marked it blocked. Raw 1I50
  chain L residues were rejected by OSPREY template matching, including
  protonation/template issues around `HIS L53` and full deletion of `L54-L70`
  even with atom-name matching. This needs structure cleaning/protonation before
  it can be a valid IDP/IDR pilot artifact.
- Generated a runnable fallback real-artifact plumbing pilot from the existing
  OSPREY-reduced `examples/2RL0.kstar/2RL0.min.reduce.pdb`, residues A154-A164,
  with six WT-flexible positions:

```text
position.0=A154,28
position.1=A156,5
position.2=A158,9
position.3=A160,2
position.4=A162,19
position.5=A164,9
totalConformations=430920
```

Artifact paths:

```text
slurm/artifacts/cohere_2rl0_rg/2rl0_A153_A164.emat
slurm/artifacts/cohere_2rl0_rg/2rl0_A153_A164.simple.confspace
slurm/artifacts/cohere_2rl0_rg/2rl0_A153_A164.rcs.txt
slurm/artifacts/cohere_2rl0_rg/2rl0_A153_A164.metadata.properties
```

Pilot config:

```text
slurm/configs/cohere_2rl0_rg_pilot.properties
```

Local validation:

```bash
java -cp "build/install/osprey/lib/*" \
  edu.duke.cs.osprey.cohere.WmbPilotScanJob \
  --config /home/users/lz280/IdeaProjects/OSPREY3/slurm/configs/cohere_2rl0_rg_pilot.properties
```

Result: completed successfully; both candidate portfolios passed the current
overlap gates. The local CSV reported ESS fractions of `0.9648` for baseline
and `0.8207` for tempered, with Rg-chi2 estimates around `1.56e-3` and
`1.36e-3`.

Submitted Slurm pilot:

```bash
sbatch slurm/scripts/run_cohere_pilot_scan.slurm \
  /home/users/lz280/IdeaProjects/OSPREY3/slurm/configs/cohere_2rl0_rg_pilot.properties
```

Result: job `12011081` ran on `fennario-01` and completed successfully in
about 18 seconds after the distribution build. Stdout ended with
`COHERE WMB pilot scan finished`, `Candidates: 2`, `Passed: 2`; stderr was
empty. Output paths:

```text
slurm/results/cohere_2rl0_rg_pilot.csv
slurm/results/cohere_2rl0_rg_pilot.json
slurm/logs/cohere_pilot_12011081.out
slurm/logs/cohere_pilot_12011081.err
```

Regression verification after the real-artifact serialization changes:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 30/30 targeted COHERE tests passed locally.

Remaining work before this becomes a true IDP/IDR benchmark, updated after
§10.22-§10.24:

- First prove the real-observable loop outside the OSPREY bottleneck: real
  IDP/IDR sequence + real SAXS/NMR data + backbone-diverse ensemble + external
  back-calculator + chi2/likelihood/reweighting/stability report.
- Keep the current OSPREY/WMB pilot as a certification-engine prototype, not as
  the fastest path to a biologically meaningful IDP ensemble. The current
  OSPREY artifacts mainly explore fixed-backbone side-chain RCs; they do not
  yet represent the backbone conformational diversity that drives IDP SAXS.
- Add a new IDP backbone factor-graph/state-to-PDB adapter before claiming WMB
  as the main IDP engine. WMB is reusable mathematically, but the current
  implementation is wired to OSPREY `SimpleConfSpace`/`EnergyMatrix` objects.
- Replace smoke-mode SAXS handling with production tooling when available:
  FoXS/CRYSOL/Pepsi-SAXS or a validated `saxs_md` solvent/blank setup; add
  SPARTA+/SHIFTX2 or another chemical-shift wrapper only after the SAXS loop is
  stable.
- Add deterministic observable bounds and certificate settings after the real
  observable/reweighting loop is numerically stable.

### 10.20 Artifact-builder regression landing (2026-07-08)

Closed the remaining regression-coverage gap from §10.19 for the real-artifact
plumbing path.

- `CoherePilotArtifactBuilder.run(...)` is now package-visible so tests can
  exercise the builder without going through `main(...)` and risking
  `System.exit(...)` on failure.
- `CoherePilotArtifactBuilder` now accepts user-facing `WT`, `wildtype`, and
  `wild-type` aliases in `--res-types`, normalizing them to OSPREY's internal
  `Strand.WildType` marker. This fixes a real CLI mismatch: the usage text
  advertised `--res-types WT,ALA`, but the implementation previously passed
  literal `WT` into the residue-template lookup.
- Added `artifactBuilderRoundTripsConfSpaceAndPilotWritesSamplePdb()` to
  `TestObservableEstimator`. The test builds a tiny 2RL0-derived pilot artifact
  in a temp directory (`A154-A156`, flexible `A154,A156`), reads the serialized
  `SimpleConfSpace` back with `ObjectIO`, checks the emitted `EnergyMatrix` and
  RC file, then runs `WmbPilotScanJob` with `observable.confSpace` and an
  external `{pdbFile}` command that requires a non-empty PDB containing `ATOM`.
  This directly locks down the `SimpleConfSpace` round-trip and sample-PDB
  materialization path that production SAXS/NMR wrappers will depend on.

Verification:

```bash
./gradlew test --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon
```

Result: build successful; 31/31 targeted COHERE tests passed. The first attempt
found the `WT` alias bug described above; after the normalization fix the new
artifact/PDB regression passed and the full targeted class completed
successfully.

Submitted Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_unit_tests.slurm
```

Result: job `12012247` ran on `fennario-01` and completed successfully in 29s;
all 31 targeted COHERE tests passed, stdout ended with `BUILD SUCCESSFUL`, and
stderr only contained the existing incubator-module warning.

### 10.21 First true 1I50 chain-L artifact and pilot landing (2026-07-08)

Closed the concrete blocker called out in §10.19: a real IDP/IDR-like input now
has OSPREY-compatible cleaned coordinates, a minimized `EnergyMatrix`, and a
passing WMB/Rg pilot.

Root cause of the earlier block:

- Raw `/home/users/lz280/Downloads/1i50.pdb` is heavy-atom-only. OSPREY's
  current template matcher rejected `HIS L53` because unrenamed `HIS` requires
  existing HD/HE hydrogens to infer HID/HIE/HIP. After hand-renaming `HIS L53`
  to `HIE`, all 29 residues in L42-L70 were still deleted because the whole
  segment lacked the hydrogens expected by OSPREY's residue templates.
- The fix was not another template-matching flag. It required structure prep:
  trim to the observed chain-L segment, add missing atoms/hydrogens, and repair
  the N-terminal hydrogen naming mismatch (`H` -> `H1`) produced by PDBFixer.

Added reproducible prep wrapper:

```text
slurm/wrappers/prepare_cohere_pdbfixer_segment.py
```

It uses PDBFixer/OpenMM only for structure preparation, not for dynamics or
sampling. It keeps a specified chain/residue interval, removes heterogens and
waters, clears PDBFixer's missing-residue list so SEQRES gaps are not rebuilt,
adds missing atoms/hydrogens at pH 7, and rewrites the first residue's terminal
`H` atom to OSPREY's expected `H1` name.

Prepared coordinate artifact:

```text
slurm/artifacts/cohere_1i50_chainL_clean/1i50_chainL_42_70_pdbfixer_ph7_osprey.pdb
```

Prep result: 29 residues (`L42-L70`), 501 atoms, 261 hydrogens.

Generated minimized COHERE/OSPREY artifacts:

```text
slurm/artifacts/cohere_1i50_chainL/1i50_chainL.emat
slurm/artifacts/cohere_1i50_chainL/1i50_chainL.simple.confspace
slurm/artifacts/cohere_1i50_chainL/1i50_chainL.rcs.txt
slurm/artifacts/cohere_1i50_chainL/1i50_chainL.assignments.txt
slurm/artifacts/cohere_1i50_chainL/1i50_chainL.metadata.properties
```

Builder command used `--minimize true`, `--template-matching atom-names`, and
six WT-flexible positions:

```text
L42, L48, L53, L58, L64, L70
```

Metadata:

```text
positions=6
totalConformations=3745440
position.0=L42,34
position.1=L48,3
position.2=L53,8
position.3=L58,27
position.4=L64,5
position.5=L70,34
```

Local validation:

```bash
java -cp "build/install/osprey/lib/*" \
  edu.duke.cs.osprey.cohere.WmbPilotScanJob \
  --config slurm/configs/cohere_1i50_chainL_rg_pilot.properties
```

Result: completed successfully; both candidates passed. The local report gave:

```text
baseline: ESS/N=0.999984, max-weight-share=0.01050, log-weight-range=0.0198,
          Rg-chi2 estimate=3.2579e-4
tempered: ESS/N=0.862949, max-weight-share=0.01875, log-weight-range=2.9666,
          Rg-chi2 estimate=3.1359e-4
```

Submitted Slurm pilot:

```bash
sbatch slurm/scripts/run_cohere_pilot_scan.slurm \
  /home/users/lz280/IdeaProjects/OSPREY3/slurm/configs/cohere_1i50_chainL_rg_pilot.properties
```

Result: job `12012276` ran on `fennario-01`, rebuilt `installDist`, completed
successfully in about 19s, and reported `Candidates: 2`, `Passed: 2`; stderr
was empty. Output paths:

```text
slurm/results/cohere_1i50_chainL_rg_pilot.csv
slurm/results/cohere_1i50_chainL_rg_pilot.json
slurm/logs/cohere_pilot_12012276.out
slurm/logs/cohere_pilot_12012276.err
```

Interpretation: the previously missing "real IDP/IDR input -> cleaned PDB ->
minimized OSPREY artifact -> WMB pilot with materialized PDB observables" path
now exists. The current Rg observable is still a plumbing/debug observable, not
a publishable experimental back-calculator benchmark.

### 10.22 First SAXS back-calculator wrapper smoke (2026-07-07)

Moved the 1I50 pilot one layer above the Rg plumbing observable by wiring a
real installed SAXS curve calculator into the external batch-observable path.

Local inventory found:

```text
/home/users/lz280/miniconda3/envs/confdiff/bin/saxs_md
/home/users/lz280/miniconda3/envs/confdiff/bin/saxs_rism
```

`saxs_rism` requires 3D-RISM grid files, so it is not the shortest next step.
`saxs_md` can calculate a SAXS q/I curve from a solute PDB plus a solvent PDB.
For smoke testing only, an empty solvent PDB is accepted by the binary and
produces a deterministic curve; this is not a final solvent-corrected SAXS
model.

Added:

```text
slurm/wrappers/cohere_saxs_md_batch_observable.py
slurm/configs/cohere_1i50_chainL_saxs_md_synthetic_pilot.properties
```

The wrapper:

- consumes COHERE materialized sample PDBs from `{pdbFilesFile}`;
- runs the installed `saxs_md` binary for each sample;
- accepts either `--reference-profile` q/I(/sigma) data or a synthetic
  `--reference-pdb`;
- requires an explicit `--solvent-pdb` or `--empty-solvent`, so smoke-mode
  solvent handling is never implicit;
- fits a scalar intensity scale by weighted least squares unless
  `--no-fit-scale` is set;
- compares predicted curves to the reference profile using q interpolation and
  writes indexed CSV:

```text
index,chi2,scale,n_points,rms_relative,profile
```

Generated the current synthetic reference artifact:

```text
slurm/artifacts/cohere_1i50_chainL_saxs_md/empty-solvent.pdb
slurm/artifacts/cohere_1i50_chainL_saxs_md/1i50_chainL_saxs_md_empty_solvent.dat
```

The profile was generated from:

```text
slurm/artifacts/cohere_1i50_chainL_clean/1i50_chainL_42_70_pdbfixer_ph7_osprey.pdb
```

with `q=0.00..0.50 A^-1` at `dq=0.02 A^-1`. This target is deliberately named
`synthetic`; it is a back-calculator integration smoke, not experimental SAXS.

Verification:

```bash
python3 -m py_compile slurm/wrappers/cohere_saxs_md_batch_observable.py
```

passed. A direct self-reference wrapper call on the cleaned 1I50 PDB returned:

```text
chi2=1.706367099e-08
scale=0.9999992803
n_points=26
rms_relative=6.53139935e-06
```

Local end-to-end pilot:

```bash
java -cp "build/install/osprey/lib/*" \
  edu.duke.cs.osprey.cohere.WmbPilotScanJob \
  --config slurm/configs/cohere_1i50_chainL_saxs_md_synthetic_pilot.properties
```

Result: completed successfully; both candidates passed. The 32-sample report:

```text
baseline: ESS/N=0.999979, max-weight-share=0.03149, log-weight-range=0.0198,
          SAXS-md synthetic chi2 estimate=0.0689331,
          observed chi2 range=[0.0116773, 0.238096]
tempered: ESS/N=0.838003, max-weight-share=0.05324, log-weight-range=2.9041,
          SAXS-md synthetic chi2 estimate=0.0822404,
          observed chi2 range=[0.00558755, 0.291924]
```

Slurm validation:

```bash
sbatch slurm/scripts/run_cohere_pilot_scan.slurm \
  /home/users/lz280/IdeaProjects/OSPREY3/slurm/configs/cohere_1i50_chainL_saxs_md_synthetic_pilot.properties
```

Result: job `12012382` ran on `fennario-01`, rebuilt `installDist`, completed
successfully in about 24s, and reported `Candidates: 2`, `Passed: 2`; stderr
was empty. Output paths:

```text
slurm/results/cohere_1i50_chainL_saxs_md_synthetic_pilot.csv
slurm/results/cohere_1i50_chainL_saxs_md_synthetic_pilot.json
slurm/logs/cohere_pilot_12012382.out
slurm/logs/cohere_pilot_12012382.err
```

Interpretation: the previous claim "we have no back-calculator" is now too
broad. We still do not have a publishable experimental SAXS/NMR benchmark, but
we do have a working external SAXS curve-calculator path:

```text
COHERE sample PDB -> saxs_md q/I curve -> scaled chi2 -> WMB observable report
```

Historical next-step note, superseded by §10.23-§10.24: the synthetic
reference has now been replaced by a staged real SASBDB profile (`SASDXC6`),
and the plan no longer treats the OSPREY fixed-backbone side-chain pilot as the
shortest path to a scientifically meaningful IDP SAXS benchmark.

### 10.23 First real experimental SAXS artifact: SASBDB SASDXC6 (2026-07-07)

Moved beyond synthetic SAXS by staging a real public SASBDB curve and matching
sequence.

Source:

```text
SASDXC6 - Sodium/hydrogen exchanger 6 protein (HsNHE6) C-terminal end (G586-A701)
Sample: human NHE6 C-terminal construct, monomer, expected MW 13 kDa
SAXS: EMBL P12/PETRA III, buffer-subtracted, Rg_Guinier 3.1 nm
UniProt: Q92581 residues 586-701
```

Downloaded/staged artifacts:

```text
slurm/artifacts/cohere_sasdxc6_real_saxs/SASDXC6.dat
slurm/artifacts/cohere_sasdxc6_real_saxs/SASDXC6.zip
slurm/artifacts/cohere_sasdxc6_real_saxs/SASDXC6_6012.fasta
```

Observed file facts:

```text
SASDXC6.dat: 2651 lines; 2596 numeric q/I/sigma points
q range: 0.02253841..7.27662 nm^-1 = 0.002253841..0.727662 A^-1
positive intensities: 2580; nonpositive intensities: 16
FASTA length: 116 residues, sequence corresponds to G586-A701
SASDXC6.zip: contains only the experimental data file, no atomic model
```

The SASBDB file uses `s` in nm^-1 and columns `s, I(s), sigma`, while
`saxs_md` uses A^-1. The SAXS wrapper therefore needs explicit unit conversion
before real-data scoring:

```text
--reference-q-scale 0.1
```

The first few low-q experimental points include negative intensities after
background subtraction. Real-data scoring should either fit only a physically
safe q window or use a documented filter such as `--positive-reference-only`;
negative points should not silently dominate chi2.

Structure-prep artifacts for a deterministic starting model:

```text
slurm/wrappers/fasta_to_extended_backbone_pdb.py
slurm/artifacts/cohere_sasdxc6_real_saxs/sasdxc6_A586_701_extended_backbone.pdb
slurm/artifacts/cohere_sasdxc6_real_saxs/sasdxc6_A586_701_pdbfixer_ph7_osprey.pdb
```

The FASTA-to-PDB script generated an extended-chain backbone with 116 residues
and 464 backbone atoms. PDBFixer then completed side-chain atoms and hydrogens:

```text
residues=116 first=586:GLY last=701:ALA
atoms=1710 hydrogens=835
```

Important caveat: this is a deterministic extended-chain starting structure,
not a real IDP ensemble generator. Its value is making the real SAXS
back-calculator/reweighting loop testable. A publishable benchmark still needs
backbone-diverse conformers from a defensible IDP conformer generator, ensemble
library, or explicit backbone-state model.

Attempted OSPREY artifact build:

```text
java -cp "build/install/osprey/lib/*" \
  edu.duke.cs.osprey.cohere.CoherePilotArtifactBuilder \
  --pdb slurm/artifacts/cohere_sasdxc6_real_saxs/sasdxc6_A586_701_pdbfixer_ph7_osprey.pdb \
  --out slurm/artifacts/cohere_sasdxc6_real_saxs \
  --first A586 --last A701 \
  --flex A602,A632,A650,A667,A678,A692 \
  --prefix sasdxc6_A586_701 \
  --res-types WT --add-wild-type false --minimize true --cpus 4 \
  --template-matching atom-names
```

The build completed but emitted:

```text
WARNING: 2 Strand residue(s) could not be matched to templates and were automatically deleted:
DXH A 670
DXH A 698
```

This is a concrete example of why OSPREY should not stay on the critical path
for the fastest real IDP/SAXS benchmark. The issue is not central science; it is
template/protonation/interface friction. The prep script was updated with a
`--his-name` option so PDBFixer `HIS` residues can be normalized to OSPREY
`HID/HIE/HIP` when the OSPREY path is needed, but this should be treated as
supporting infrastructure, not the main benchmark route.

### 10.24 Route decision: fast real benchmark first, OSPREY/WMB second (2026-07-07)

Decision after the SASDXC6 work and user pushback:

**If the immediate goal is the fastest credible real IDP/IDR result, do not put
OSPREY on the critical path.** Use OSPREY/WMB as the second-stage
certification/differentiation engine.

Reasoning:

- OSPREY's current reusable path is strongest for fixed-backbone
  side-chain/rotamer conformation spaces. IDP SAXS is dominated by backbone
  conformational diversity and global dimensions.
- The existing WMB pilot implementation is wired to OSPREY `SimpleConfSpace`,
  `EnergyMatrix`, RC assignments, and OSPREY PDB materialization. This is useful
  for validating orchestration, but not yet the right state space for IDPs.
- Back-calculation physics is already delegated to external tools. The fastest
  real benchmark loop is simply:

```text
real sequence/profile -> backbone-diverse ensemble PDBs -> SAXS/NMR
back-calculator -> chi2/likelihood -> reweighting/ESS/stability report
```

- WMB itself is not lost. It is a generic discrete factor-graph method. To use
  it correctly for the IDP route, build a new adapter around backbone bins /
  local Ramachandran factors / charge-contact factors / state-to-PDB
  materialization. That is a real engineering task, but smaller and cleaner than
  forcing IDP backbone diversity through OSPREY rotamer RCs.

New priority order:

1. **Fast non-OSPREY real SAXS benchmark.** Use SASDXC6 first: real FASTA, real
   SASBDB q/I/sigma, backbone-diverse conformers, `saxs_md` or a production
   SAXS tool, chi2/reweighting/ESS/multi-seed stability.
2. **Production back-calculator hardening.** Add q-unit handling, q-window
   filtering, positive-intensity handling, and solvent/blank strategy; swap to
   FoXS/CRYSOL/Pepsi-SAXS if available.
3. **IDP WMB adapter.** Build a discrete backbone factor-graph interface that
   can provide `logQ(c)`, sample states, and materialize PDBs without depending
   on OSPREY `SimpleConfSpace`.
4. **Certificate.** Add deterministic observable bounds and empirical/PAC
   intervals after the real-observable loop is stable.
5. **OSPREY bridge retained.** Keep the 1I50 and SASDXC6 OSPREY artifacts as
   regression/support tests for the external-observable and report plumbing,
   not as the main claim for IDP ensemble quality.

### 10.25 SASDXC6 non-OSPREY benchmark scaffold (2026-07-08)

Implemented the first executable scaffold for M1's fast real-observable route:

```text
slurm/wrappers/cohere_sasdxc6_real_benchmark.py
slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

The driver has two deliberately separate modes:

- **Backbone smoke ensemble:** generate lightweight backbone-only PDBs from the
  SASDXC6 FASTA using simple `coil/expanded/compact` phi/psi proposal modes,
  then report Rg, end-to-end distance, Rg chi2 against the SASBDB Guinier Rg
  target (`31 A`), weighted mean Rg, and ESS. This is a proposal/materialization
  smoke, not a final IDP ensemble generator.
- **External all-atom PDB scoring:** consume `--score-pdb` or
  `--score-pdb-list`, run the real SASBDB q/I/sigma profile through
  `cohere_saxs_md_batch_observable.py`, and report SAXS chi2, scale, n_points,
  likelihood-reweighting ESS, and weighted mean Rg. This is the intended route
  once a defensible all-atom conformer ensemble exists.

Hardened the SAXS wrapper for real SASBDB curves:

```text
slurm/wrappers/cohere_saxs_md_batch_observable.py
```

Additions:

- `--reference-stride` and `--max-reference-points` to thin dense experimental
  curves after q-unit conversion / q-window filtering / positive-intensity
  filtering.
- sample `qcut = max(reference_q) + dq` so `saxs_md`'s rounded output grid
  covers the last scoring point.

Important boundary found during smoke testing:

```text
Generated backbone-only PDB -> saxs_md
```

failed with:

```text
Unable to recognize atom
```

The same wrapper path works on the PDBFixer all-atom SASDXC6 artifact. The
benchmark driver now refuses `--run-saxs` on generated backbone-only PDBs unless
`--allow-backbone-only-saxs` is explicitly set. For real SAXS scoring, pass
all-atom conformers through `--score-pdb-list` or `--score-pdb`.

Local backbone-only Rg smoke:

```bash
python3 slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  --samples 3 \
  --seeds 20260707 \
  --out-dir slurm/work/cohere_sasdxc6_real_saxs_benchmark_smoke_no_saxs2
```

Result:

```text
rg_min=14.72530563
rg_mean=35.51387973
rg_max=62.89103287
target_rg=31
rg_best_chi2=0.1721750877
rg_ess_fraction=0.3369698807
rg_weighted_mean=28.84826029
```

Local all-atom real SAXS smoke against the staged PDBFixer structure:

```bash
python3 slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  --score-pdb slurm/artifacts/cohere_sasdxc6_real_saxs/sasdxc6_A586_701_pdbfixer_ph7_osprey.pdb \
  --seeds 20260707 \
  --reference-q-max 0.05 \
  --reference-stride 20 \
  --max-reference-points 12 \
  --timeout-seconds 60 \
  --out-dir slurm/work/cohere_sasdxc6_real_saxs_benchmark_smoke_all_atom2 \
  --run-saxs
```

Result:

```text
SAXS status=ok
Rg=114.4630278 A
SAXS chi2=75.29430832
scale=2.150819093e-09
n_points=10
rms_relative=1.361941762
```

Interpretation: the current all-atom structure is still the deterministic
extended-chain PDBFixer artifact, so the poor Rg/SAXS fit is expected. The
useful result is that the real profile, q-unit conversion, q-windowing,
positive-intensity filtering, q-grid coverage, `saxs_md`, CSV merge, and
reweighting summary now run outside OSPREY.

Slurm entry point:

```bash
sbatch slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

By default it scores the single staged all-atom PDB as a plumbing smoke. For a
real ensemble, provide:

```bash
COHERE_SASDXC6_PDB_LIST=/path/to/all_atom_conformers.list \
sbatch slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

Submitted the default Slurm smoke:

```bash
sbatch slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

Result: job `12012680` ran on `fennario-01`, completed in about 31s, and stderr
was empty. Output paths:

```text
slurm/logs/cohere_sasdxc6_saxs_12012680.out
slurm/logs/cohere_sasdxc6_saxs_12012680.err
slurm/work/cohere_sasdxc6_real_saxs_benchmark/summary.tsv
slurm/work/cohere_sasdxc6_real_saxs_benchmark/seed_20260707/saxs_values.csv
```

The default q-window/stride smoke (`q_max=0.5 A^-1`, `reference_stride=10`,
`max_reference_points=200`) scored 179 real SAXS points:

```text
SAXS status=ok
Rg=114.4630278 A
SAXS chi2=32.1610611
scale=1.840107609e-09
n_points=179
rms_relative=0.6634608611
```

Verification:

```bash
python3 -m py_compile \
  slurm/wrappers/cohere_saxs_md_batch_observable.py \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  slurm/wrappers/cohere_rg_batch_observable.py \
  slurm/wrappers/fasta_to_extended_backbone_pdb.py

bash -n slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

Both passed locally.

### 10.26 SASDXC6 PDBFixer all-atom ensemble benchmark (2026-07-08)

The first real SAXS benchmark has now moved beyond a single deterministic
extended-chain artifact.  Added a batch PDBFixer completion helper and an
end-to-end Slurm pipeline:

```text
slurm/wrappers/cohere_complete_pdbfixer_ensemble.py
slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm
```

Pipeline:

```text
SASDXC6 FASTA
  -> backbone-diverse coil/expanded/compact proposal PDBs
  -> PDBFixer all-atom completion for chain A586-A701
  -> completed PDB list
  -> real SASDXC6 saxs_md scoring
  -> chi2 / Rg / ESS summary
```

`cohere_complete_pdbfixer_ensemble.py` keeps the requested chain and residue
range, removes heterogens/waters, clears PDBFixer missing-residue rebuilds,
adds missing heavy atoms and hydrogens at the requested pH, writes
`completed-pdb-list.txt`, and records per-sample prep metadata.  This preserves
the current project's accepted structure-prep route while avoiding OSPREY for
the fast real-observable benchmark.

Local two-sample smoke:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  --samples 2 \
  --seeds 20260707 \
  --modes coil,expanded \
  --out-dir slurm/work/cohere_sasdxc6_pdbfixer_ensemble_smoke/backbone

/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_complete_pdbfixer_ensemble.py \
  --pdb-list slurm/work/cohere_sasdxc6_pdbfixer_ensemble_smoke/backbone/seed_20260707/pdb-list.txt \
  --out-dir slurm/work/cohere_sasdxc6_pdbfixer_ensemble_smoke/all_atom/seed_20260707 \
  --chain A \
  --first-res 586 \
  --last-res 701 \
  --ph 7.0

/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  --score-pdb-list slurm/work/cohere_sasdxc6_pdbfixer_ensemble_smoke/all_atom/seed_20260707/completed-pdb-list.txt \
  --seeds 20260707 \
  --reference-q-max 0.5 \
  --reference-stride 10 \
  --max-reference-points 200 \
  --timeout-seconds 120 \
  --out-dir slurm/work/cohere_sasdxc6_pdbfixer_ensemble_smoke/score_q05 \
  --run-saxs
```

Result:

```text
completed PDBs: 2/2
residues per PDB: 116
atoms per PDB: 1710
hydrogens per PDB: 835
q_max=0.05 smoke best chi2: 6.190933129 over 10 points
q_max=0.5 smoke best chi2: 14.9311131 over 179 points
```

Submitted the default six-sample Slurm pipeline:

```bash
sbatch slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm
```

Result: job `12012733` ran on `fennario-01`, completed in about 4m43s, and
stderr was empty. Output paths:

```text
slurm/logs/cohere_sasdxc6_pdbfixer_12012733.out
slurm/logs/cohere_sasdxc6_pdbfixer_12012733.err
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/backbone/summary.tsv
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/all_atom/seed_20260707/pdbfixer_metadata.tsv
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/all_atom/seed_20260707/completed-pdb-list.txt
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/score/summary.tsv
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/score/seed_20260707/saxs_values.csv
```

Backbone proposal summary:

```text
samples=6
modes=coil,expanded,compact
Rg min/mean/max=14.72530563 / 38.7206787 / 83.53637664 A
Rg best chi2=0.1721750877
Rg ESS fraction=0.2904754367
Rg weighted mean=27.09390557 A
```

All six generated PDBs completed with PDBFixer:

```text
status=ok for 6/6
residues=116
atoms=1710
hydrogens=835
first=586:GLY
last=701:ALA
```

Real SAXS scoring summary (`q_max=0.5 A^-1`, `reference_stride=10`,
`max_reference_points=200`, 179 positive-intensity points):

```text
samples=6
Rg min/mean/max=15.07878921 / 38.40611774 / 82.78859432 A
SAXS status=ok
SAXS best chi2=14.96086822
SAXS best index=0
SAXS ESS fraction=0.166666669
SAXS weighted mean Rg=28.33891391 A
```

Per-conformer SAXS chi2:

```text
index  chi2         rms_relative
0      14.96086822  0.4726334854
1      22.18285653  0.3723869409
2      34.70826407  0.5740875341
3      18.96154718  0.5373913087
4      17.06497408  0.4500921389
5      15.17049     0.4559434924
```

Interpretation:

- This is the first reproducible SASDXC6 path with real experimental q/I/sigma,
  generated backbone-diverse conformers, all-atom PDBFixer completion, `saxs_md`
  scoring, and reweighting summaries running outside OSPREY.
- The best all-atom ensemble member improves over the deterministic extended
  artifact default benchmark (`chi2=14.96` vs `32.16` over the same 179-point
  q-window/stride setup).
- The SAXS ESS fraction is effectively one conformer out of six.  That is a
  useful discrimination signal, but not yet a healthy ensemble posterior.
- The conformers are still scaffold proposals, not a defensible final IDP
  generator: the first geometry QC pass shows compact/coil proposals can carry
  severe nonlocal clashes, PDBFixer is still atom completion rather than
  side-chain packing, OpenMM relaxation is too slow for the default fast path,
  and there is still no MD/coil-library generator or solvent-model cross-check
  against FoXS/CRYSOL/Pepsi-SAXS.

Immediate next technical steps, after the first 48-conformer raw run:

1. Finish the in-flight QC-passed SAXS scoring for the same 48-conformer job and
   compare it against the raw all-conformer fit.
2. Keep reporting QC pass fraction by proposal mode; do not accept lower raw
   SAXS chi2 when the winning conformer fails steric QC.
3. Validate the new compact backbone self-avoidance retry with full PDBFixer
   heavy-atom QC.  The CA-level filter improves compact proposals but does not
   prove side-chain physicality.
4. Keep OpenMM minimization as a separate experiment until the full 116-residue
   runtime is acceptable.
5. Use q-window sensitivity (`q_max=0.2`, `0.3`, `0.5`) as a required report
   field because the raw winner can change with q-window.

Verification:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python -m py_compile \
  slurm/wrappers/cohere_complete_pdbfixer_ensemble.py \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  slurm/wrappers/cohere_saxs_md_batch_observable.py

bash -n \
  slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm \
  slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

Both passed locally.

### 10.27 SASDXC6 geometry QC and filtered SAXS scoring (2026-07-08)

The PDBFixer ensemble helper now records geometry QC in the same
`pdbfixer_metadata.tsv` file and writes a second list:

```text
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/all_atom/seed_20260707/qc-passed-pdb-list.txt
```

Default QC-passed criteria:

```text
status in {ok, skipped}
residue_count_delta = 0
residue_number_gap_count = 0
ca_adjacent_outlier_count = 0
nonlocal_heavy_clash_count_lt_1p5_a = 0
```

The metadata now includes atom/residue counts, expected residue delta, CA
count, residue-number gaps, heavy-atom and CA Rg, CA end-to-end distance,
adjacent CA distance min/mean/max, nonlocal heavy-atom minimum distance, and
counts of nonlocal heavy-atom pairs below 1.5 A and 2.0 A.

Six-sample QC result:

```text
all samples: residues=116, residue_count_delta=0, ca_count=116
all samples: residue_number_gap_count=0, ca_adjacent_outlier_count=0

index  mode      rg_ca       min_nonlocal_heavy  clashes_lt_1p5  close_lt_2p0
0      coil      29.036199   0.476530            14              41
1      expanded  62.967581   1.519248            0               1
2      compact   14.829266   0.294635            97              238
3      coil      19.100770   0.450063            29              72
4      expanded  83.553473   1.976362            0               1
5      compact   23.154142   0.450455            41              121
```

Only 2/6 structures pass the strict no-1.5-A-clash filter, and both are
expanded-mode conformers:

```text
000001_sample_000001_expanded_pdbfixer.pdb
000004_sample_000004_expanded_pdbfixer.pdb
```

Filtered SAXS scoring command:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  --score-pdb-list slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/all_atom/seed_20260707/qc-passed-pdb-list.txt \
  --seeds 20260707 \
  --reference-q-max 0.5 \
  --reference-stride 10 \
  --max-reference-points 200 \
  --timeout-seconds 120 \
  --out-dir slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/score_qc_passed \
  --run-saxs
```

Filtered result:

```text
samples=2
Rg min/mean/max=61.82638294 / 72.30748863 / 82.78859432 A
SAXS status=ok
SAXS best chi2=17.06497408
SAXS best index=1
SAXS ESS fraction=0.5
SAXS weighted mean Rg=82.78859432 A
SAXS n_points=179
```

Per-filtered-conformer SAXS chi2:

```text
filtered index  original index  chi2         rms_relative
0               1               22.18285653  0.3723869409
1               4               17.06497408  0.4500921389
```

Interpretation:

- The previous all-conformer best (`chi2=14.96086822`) came from original
  index 0, which has 14 nonlocal heavy-atom pairs below 1.5 A. It is useful as
  a SAXS sensitivity signal, but not as a physically acceptable benchmark hit.
- The strict QC-passed subset has worse best SAXS fit (`17.06497408`) and much
  larger Rg values. This means the current scaffold generator is not producing
  enough compact, clash-free conformers near the SASDXC6 SAXS optimum.
- Scaling to 48-200 conformers should now report both raw best chi2 and
  QC-passed best chi2, plus QC pass rate by proposal mode. Otherwise the
  benchmark can reward geometrically invalid compact structures.

Optional OpenMM minimization is also wired in, but remains off by default:

```text
cohere_complete_pdbfixer_ensemble.py:
  --minimize
  --minimize-max-iterations
  --minimize-tolerance
  --minimize-backbone-restraint-kcal-mol-a2
  --minimize-platform

run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm:
  COHERE_SASDXC6_MINIMIZE=1
  COHERE_SASDXC6_MINIMIZE_MAX_ITERATIONS=...
  COHERE_SASDXC6_MINIMIZE_BACKBONE_RESTRAINT=...
  COHERE_SASDXC6_SCORE_QC_PASSED=1
```

Tiny minimization smoke, using only residues A586-A590, passed:

```text
residues=5
atoms=71
min_nonlocal_heavy_distance_a=2.914398
nonlocal_heavy_clash_count_lt_1p5_a=0
initial_potential=937.168173 kJ/mol
final_potential=-496.140761 kJ/mol
minimize_seconds=0.302308
```

Full 116-residue single-structure minimization was interrupted after exceeding
90 seconds even with a 10-iteration budget, so it should not be part of the
default fast path yet. Treat minimization as a separate Slurm experiment after
the generator/QC pass-rate problem is understood.

Updated route:

1. Generate more backbone proposals, but score and report the QC-passed subset
   separately from all generated structures.
2. Improve the proposal generator or add a fast side-chain/clash repair stage
   until compact conformers can pass geometry QC.
3. Only then run larger real SAXS ensembles and q-window sensitivity
   (`q_max=0.2`, `0.3`, `0.5`).

Verification:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python -m py_compile \
  slurm/wrappers/cohere_complete_pdbfixer_ensemble.py \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  slurm/wrappers/cohere_saxs_md_batch_observable.py

bash -n \
  slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm \
  slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm
```

Both passed locally.

### 10.28 SASDXC6 QC summary and dual SAXS scoring switch (2026-07-08)

The PDBFixer ensemble helper now writes a machine-readable QC aggregate:

```text
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark/all_atom/seed_20260707/qc_summary.tsv
```

The metadata rows also include a proposal `mode` inferred from the input PDB
name (`coil`, `expanded`, `compact`, or `unknown`).  The summary reports total
samples, completed/failed counts, QC-passed counts, QC pass fractions, CA Rg
range, nonlocal heavy-atom minimum-distance statistics, total/mean nonlocal
1.5 A clashes, total/mean nonlocal 2.0 A close contacts, and mean adjacent-CA
outlier count for both all samples and each proposal mode.

Six-sample QC summary:

```text
group          samples  completed  failed  qc_passed  qc_pass_fraction
all            6        6          0       2          0.333333
mode:coil      2        2          0       0          0.000000
mode:expanded  2        2          0       2          1.000000
mode:compact   2        2          0       0          0.000000
```

The Slurm benchmark driver now supports two additional controls:

```text
COHERE_SASDXC6_SCORE_BOTH=1
COHERE_SASDXC6_SKIP_EXISTING=1
```

`SCORE_BOTH=1` scores `completed-pdb-list.txt` into `score_all/` and, when the
QC-passed list is nonempty, scores `qc-passed-pdb-list.txt` into
`score_qc_passed/` in the same job.  `SKIP_EXISTING=1` reuses existing PDBFixer
outputs and recomputes metadata/QC summaries, which makes reruns cheap when the
wrapper logic changes.

Local six-sample smoke:

```bash
COHERE_SASDXC6_SCORE_BOTH=1 \
COHERE_SASDXC6_SKIP_EXISTING=1 \
bash slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm
```

Outputs:

```text
score_all/summary.tsv:
  samples=6
  saxs_best_chi2=14.96086822
  saxs_best_index=0
  saxs_ess_fraction=0.166666669
  saxs_weighted_mean_rg=28.33891391

score_qc_passed/summary.tsv:
  samples=2
  saxs_best_chi2=17.06497408
  saxs_best_index=1
  saxs_ess_fraction=0.5
  saxs_weighted_mean_rg=82.78859432
```

This made the first 48-conformer benchmark straightforward: run with
`COHERE_SASDXC6_SCORE_BOTH=1`, report raw and QC-passed SAXS fits side by side,
and use `qc_summary.tsv` to decide whether the generator or the repair stage
needs to be fixed first.

### 10.29 SASDXC6 48-conformer raw SAXS and generator cleanup (2026-07-08)

The first 48-conformer SASDXC6 benchmark was submitted as Slurm job `12012798`:

```bash
sbatch --export=ALL,COHERE_SASDXC6_SAMPLES=48,COHERE_SASDXC6_SCORE_BOTH=1,COHERE_SASDXC6_SKIP_EXISTING=1,COHERE_SASDXC6_OUT_ROOT=/home/users/lz280/IdeaProjects/OSPREY3/slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark_s48 \
  slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm
```

PDBFixer completed all 48 structures.  Geometry QC remains the limiting issue:

```text
group          samples  completed  qc_passed  qc_pass_fraction
all            48       48         21         0.437500
mode:coil      16       16         5          0.312500
mode:expanded  16       16         16         1.000000
mode:compact   16       16         0          0.000000
```

Raw all-conformer SAXS scoring finished:

```text
samples=48
saxs_best_chi2=3.951338931
saxs_best_index=15
saxs_ess_fraction=0.02083333333
saxs_weighted_mean_rg=35.4594177 A
saxs_n_points=179
```

The raw best conformer is original index 15, mode `coil`, with
`rg_ca_a=35.454446 A`, but it fails strict geometry QC:

```text
nonlocal_heavy_clash_count_lt_1p5_a=37
nonlocal_heavy_close_count_lt_2p0_a=89
```

This is the strongest evidence so far that raw SAXS chi2 cannot be treated as a
valid benchmark success criterion without a geometry filter.

The q-window rescore now runs from existing `saxs_md` profiles, without
re-running `saxs_md`:

```text
slurm/wrappers/cohere_saxs_profile_rescore.py
slurm/wrappers/cohere_sasdxc6_q_sweep_report.py
```

Raw 48-conformer q-window result:

```text
qmax  n_points  best_chi2    original_index  mode  qc_passed
0.2   71        3.523611516  18              coil  1
0.3   107       3.328667787  15              coil  0
0.5   179       3.951338931  15              coil  0
```

Interpretation:

- The low-q window (`qmax=0.2 A^-1`) chooses a QC-passed coil conformer.
- Adding higher-q data switches the raw winner to a QC-failed coil conformer.
- q-window sensitivity must therefore be part of the benchmark report, not an
  optional post-hoc check.

The report helper now maps local indices from filtered subsets back to original
conformer indices and modes:

```text
slurm/wrappers/cohere_sasdxc6_benchmark_report.py
```

Job `12012798` finished at `Wed Jul 8 00:55:27 EDT 2026`, and the final
QC-passed SAXS summary is now available:

```text
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark_s48/score_qc_passed/summary.tsv
samples=21
saxs_best_chi2=5.23218989
saxs_best_local_index=16
saxs_best_original_index=36
saxs_best_mode=coil
saxs_best_qc_passed=1
saxs_best_rg_ca_a=37.662971 A
saxs_ess_fraction=0.0476560963
saxs_weighted_mean_rg=37.41417031 A
saxs_n_points=179
```

The combined machine-readable report is:

```text
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark_s48/benchmark_report.tsv
```

QC-passed q-window result:

```text
qmax  n_points  best_chi2    local_index  original_index  mode  qc_passed
0.2   71        3.523611516  8            18              coil  1
0.3   107       4.315527051  8            18              coil  1
0.5   179       5.23218989   16           36              coil  1
```

The combined q-window report is:

```text
slurm/work/cohere_sasdxc6_pdbfixer_ensemble_benchmark_s48/q_sweep/q_sweep_report.tsv
```

The final interpretation is unchanged but sharper: the raw all-conformer score
is numerically better (`3.951338931`) only because it can select a QC-failed
coil conformer.  The best valid QC-passed full-window hit is original index 36
with chi2 `5.23218989`; under lower-q windows, original index 18 is the best
QC-passed conformer.

Generator cleanup:

`cohere_sasdxc6_real_benchmark.py` now records backbone proposal QC in
`metadata.csv` and retries compact proposals with a CA-level self-avoidance
filter.  New fields:

```text
proposal_attempts
ca_min_nonlocal_distance
ca_nonlocal_clash_count
ca_self_avoidance_passed
```

Default generation now applies self-avoidance to `compact` mode with:

```text
self_avoidance_attempts=256
self_avoidance_min_ca_distance=3.0 A
self_avoidance_max_ca_clashes=0
self_avoidance_local_separation=2
```

Backbone-only 48-sample smoke with the new defaults:

```text
mode      count  CA-pass  CA-pass-fraction  mean_CA_clashes  mean_attempts
coil      16     8        0.500000          2.750000         1.000000
expanded  16     16       1.000000          0.000000         1.000000
compact   16     9        0.562500          0.500000         143.375000
```

Older compact behavior, with no self-avoidance, had 0/16 CA-level passes and a
mean of 50.0625 nonlocal CA clashes.  The retry filter is therefore useful, but
it is still only a backbone-level prefilter.  The next confirmation must be a
fresh PDBFixer/heavy-atom QC run using the new generator defaults, not another
raw SAXS-only run.

### 10.30 Fresh self-avoidance-default PDBFixer confirmation (2026-07-08)

The follow-up 48-conformer PDBFixer/heavy-atom QC confirmation was submitted as
Slurm job `12017072`:

```bash
sbatch --export=ALL,COHERE_SASDXC6_SAMPLES=48,COHERE_SASDXC6_SCORE_BOTH=1,COHERE_SASDXC6_SKIP_EXISTING=0,COHERE_SASDXC6_OUT_ROOT=/home/users/lz280/IdeaProjects/OSPREY3/slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48 \
  slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm
```

This run uses the new generator defaults from §10.29 and writes to a fresh
output root:

```text
slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48
```

Backbone generation completed immediately with the expected self-avoidance
default statistics:

```text
mode      count  CA-pass  CA-pass-fraction  mean_CA_clashes  mean_attempts
coil      16     8        0.500000          2.750000         1.000000
expanded  16     16       1.000000          0.000000         1.000000
compact   16     9        0.562500          0.500000         143.375000
```

Job `12017072` completed successfully on `fennario-01` at
2026-07-08 13:58:39 EDT.  The final output files are:

```text
slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48/all_atom/seed_20260707/qc_summary.tsv
slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48/score_all/summary.tsv
slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48/score_qc_passed/summary.tsv
slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48/benchmark_report.tsv
slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48/q_sweep/q_sweep_report.tsv
```

Heavy-atom QC:

```text
group          samples  completed  failed  qc_passed  qc_pass_fraction
all            48       48         0       24         0.500000
mode:coil      16       16         0       6          0.375000
mode:expanded  16       16         0       16         1.000000
mode:compact   16       16         0       2          0.125000
```

SAXS scoring:

```text
score_set   samples  best_chi2    local_index  original_index  mode  QC-passed  heavy_clashes_lt_1p5_A
all         48       3.839304434  15           15              coil  no         36
qc_passed   24       5.231221717  17           36              coil  yes        0
```

The all-conformer best score is still a QC-failed coil conformation, so it
should not be treated as a physical hit.  The best full-window QC-passed hit is
again original index `36`, with `chi2=5.231221717`; this is essentially tied
with the previous QC-passed full-window result (`5.23218989`) while improving
the heavy-atom QC pass count from `21/48` to `24/48`.

The q-window sweep for this run is:

```text
score_set   qmax  best_chi2    original_index  mode     QC-passed
all         0.2   2.990986884  23              compact  no
all         0.3   2.716290981  23              compact  no
all         0.5   3.839304434  15              coil     no
qc_passed   0.2   3.523611516  18              coil     yes
qc_passed   0.3   4.315527051  18              coil     yes
qc_passed   0.5   5.231221717  36              coil     yes
```

Verification:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python -m py_compile \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  slurm/wrappers/cohere_saxs_profile_rescore.py \
  slurm/wrappers/cohere_sasdxc6_q_sweep_report.py \
  slurm/wrappers/cohere_sasdxc6_benchmark_report.py
```

The wrapper compile checks passed locally.

### 10.31 PDBFixer benchmark wrapper automation update (2026-07-08)

The Slurm benchmark wrapper now removes two manual steps from future SASDXC6
PDBFixer runs:

- `run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm` exposes the generator
  self-avoidance controls through environment variables:
  `COHERE_SASDXC6_SELF_AVOIDANCE_MODES`,
  `COHERE_SASDXC6_SELF_AVOIDANCE_ATTEMPTS`,
  `COHERE_SASDXC6_SELF_AVOIDANCE_MIN_CA_DISTANCE`,
  `COHERE_SASDXC6_SELF_AVOIDANCE_MAX_CA_CLASHES`, and
  `COHERE_SASDXC6_SELF_AVOIDANCE_LOCAL_SEPARATION`.
- When `COHERE_SASDXC6_SCORE_BOTH=1`, the wrapper now automatically writes:
  `benchmark_report.tsv` and `q_sweep/q_sweep_report.tsv`.
- Report generation can be disabled with `COHERE_SASDXC6_GENERATE_REPORTS=0`.

A read-only status helper was added:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_run_status.py \
  --out-root slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48
```

It reports backbone completion, PDBFixer completed/failed/QC-passed counts,
`score_all` and `score_qc_passed` completion, report existence, and current
best raw/QC-passed SAXS `chi2` values when summaries exist.  On the in-flight
`selfavoidance_default_s48` root, the helper currently sees 48 backbone PDBs,
33 CA self-avoidance-passed proposals, and no completed PDBFixer/SAXS summary
yet.

Verification:

```bash
bash -n slurm/scripts/run_cohere_sasdxc6_pdbfixer_ensemble_benchmark.slurm

/home/users/lz280/miniconda3/envs/confdiff/bin/python -m py_compile \
  slurm/wrappers/cohere_sasdxc6_run_status.py \
  slurm/wrappers/cohere_sasdxc6_benchmark_report.py \
  slurm/wrappers/cohere_sasdxc6_q_sweep_report.py
```

Both checks passed locally.  The status helper was also run against the finished
`cohere_sasdxc6_pdbfixer_ensemble_benchmark_s48` root and reported existing
`benchmark_report.tsv` and `q_sweep/q_sweep_report.tsv` artifacts.

### 10.32 Self-avoidance sweep and stricter report fields (2026-07-08)

Two additional code cleanups were completed while job `12017072` was scoring:

- `cohere_sasdxc6_benchmark_report.py` and
  `cohere_sasdxc6_q_sweep_report.py` now include best-conformer QC detail
  fields in their TSV/JSON outputs: minimum nonlocal heavy distance, nonlocal
  heavy clash/close-contact counts, CA self-avoidance minimum distance, CA
  clash count, CA pass flag, and proposal attempts.
- `cohere_sasdxc6_selfavoidance_sweep.py` was added as a conservative sweep
  driver for attempts/min-CA-distance/max-CA-clashes/local-separation grids.
  By default it writes a manifest and summary only; it submits jobs only when
  `--submit` is passed explicitly.

Example dry-run:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_selfavoidance_sweep.py \
  --work-root /tmp/cohere_sasdxc6_selfavoidance_sweep_test \
  --label test \
  --samples 4 \
  --attempts 64,128 \
  --min-ca-distance 3.0 \
  --max-ca-clashes 0 \
  --local-separation 2
```

The dry-run wrote:

```text
/tmp/cohere_sasdxc6_selfavoidance_sweep_test/selfavoidance_sweep_manifest.tsv
/tmp/cohere_sasdxc6_selfavoidance_sweep_test/selfavoidance_sweep_summary.tsv
```

Verification:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python -m py_compile \
  slurm/wrappers/cohere_sasdxc6_selfavoidance_sweep.py \
  slurm/wrappers/cohere_sasdxc6_benchmark_report.py \
  slurm/wrappers/cohere_sasdxc6_q_sweep_report.py \
  slurm/wrappers/cohere_sasdxc6_run_status.py
```

The compile and dry-run checks passed.  The report-field changes were also
tested on the finished `cohere_sasdxc6_pdbfixer_ensemble_benchmark_s48` root
using `/tmp` output paths.

### 10.33 STARLING SASDXC6 baseline smoke (2026-07-08)

The STARLING repository is public:

```text
https://github.com/idptools/starling/
```

STARLING was installed into an isolated environment under `/usr/xtmp`:

```text
/usr/xtmp/lz280/conda_envs/starling
```

The installed CLI reports STARLING `2.0.2`, using the v2.0.0 VAE/DDPM model
weights from the GitHub release.  Direct login-node runs were killed after DDIM
sampling, so a Slurm wrapper was added:

```text
slurm/scripts/run_starling_sasdxc6.slurm
```

Slurm job `12017223` completed successfully on `fennario-01`:

```text
/usr/xtmp/lz280/starling_runs/sasdxc6_s48_slurm/sasdxc6_starling_s48.starling
/usr/xtmp/lz280/starling_runs/sasdxc6_s48_slurm/sasdxc6_starling_s48_STARLING.pdb
/usr/xtmp/lz280/starling_runs/sasdxc6_s48_slurm/sasdxc6_starling_s48_STARLING.xtc
```

STARLING generated a 48-conformer CA-only coarse-grained ensemble for the same
SASDXC6 sequence.  `starling2info` reports:

```text
Number of conformations: 48
Average radius of gyration: 32.04404190319428 A
Average end-to-end distance: 78.01904296875 A
Structures: yes
```

The XTC trajectory was split into one PDB per frame with:

```text
slurm/wrappers/starling_xtc_to_pdb_list.py
```

CA-level diagnostics on the 48 STARLING frames:

```text
Rg min/mean/max: 19.93103183 / 32.04144579 / 53.19649451 A
CA self-avoidance passed: 36/48
mean CA nonlocal clash count: 0.5208333333
max CA nonlocal clash count: 3
```

CA-only SAXS scoring with the existing SASDXC6 `saxs_md` wrapper completed:

```text
score_set     samples  best_chi2    best_index  best_Rg_A    CA-passed
all           48       2.292908052  27          39.86256685  yes
CA-passed     36       2.292908052  21          39.86256685  yes
```

The q-window sweep on the CA-only STARLING profiles was:

```text
score_set  qmax  best_chi2    local_index
all        0.2   1.736743651  42
all        0.3   2.374461125  42
all        0.5   2.292908052  27
CA-passed  0.2   1.736743651  31
CA-passed  0.3   2.374461125  31
CA-passed  0.5   2.292908052  21
```

Important caveat: this is not directly comparable to the all-atom PDBFixer
benchmarks above.  STARLING's structure output is CA-only.  A two-frame
PDBFixer smoke test showed that PDBFixer can complete STARLING CA-only frames,
but strict all-atom QC failed for both tested frames because of adjacent CA
distance outliers and heavy-atom clashes:

```text
PDBFixer smoke: 2/2 completed, 0/2 QC-passed
```

Therefore STARLING is clearly stronger than the toy generator on coarse-grained
Rg/SAXS diagnostics, but it has not yet been shown to beat the current pipeline
under the strict all-atom heavy-QC metric.  A full 48-frame STARLING-to-PDBFixer
Slurm run would be the next direct heavy-QC comparison, but the smoke result
suggests strict peptide-geometry QC will likely reject many STARLING frames
unless a CA-to-all-atom reconstruction/relaxation step is added.

### 10.34 CG/CA-first baseline comparison pivot (2026-07-08)

After the STARLING baseline completed, the immediate project direction was
shifted from all-atom-first generation to CA/coarse-grained ensemble quality
first, with all-atom reconstruction retained as a downstream validation layer.

Two reporting changes were added:

```text
slurm/wrappers/cohere_sasdxc6_q_sweep_report.py
slurm/wrappers/cohere_sasdxc6_compare_baselines.py
```

`cohere_sasdxc6_q_sweep_report.py` now falls back to score-level
`metadata.csv` when no PDBFixer all-atom metadata exists.  This fixes the
STARLING/external CA-only case: q-window reports now include best-frame Rg,
CA nonlocal distance, CA clash count, CA self-avoidance pass status, and the
original frame index inferred from the PDB path.  Strict heavy-QC fields are
left blank for CA-only inputs instead of being treated as failed all-atom QC.

The new comparison wrapper reads existing benchmark artifacts and writes:

```text
slurm/work/cohere_sasdxc6_baseline_comparison.tsv
slurm/work/cohere_sasdxc6_baseline_comparison.json
```

Default compared cases:

```text
ours_selfavoidance = slurm/work/cohere_sasdxc6_pdbfixer_selfavoidance_default_s48
starling_ca       = slurm/work/starling_sasdxc6_ca_s48_benchmark
```

Key comparison rows:

```text
case                 score_set   samples  best_chi2    best_original_index  CA_passed/total  q-window best original indices
ours_selfavoidance   all         48       3.839304434  15                   33/48            23,23,15
ours_selfavoidance   qc_passed   24       5.231221717  36                   22/24            18,18,36
starling_ca          all         48       2.292908052  27                   36/48            42,42,27
starling_ca          qc_passed   36       2.292908052  27                   36/36            42,42,27
```

The comparison confirms the current state:

- STARLING is substantially stronger on CA-only SAXS and Rg diagnostics.
- The current toy generator improves from raw to self-avoidance-heavy QC only
  modestly; strict QC-passed best SAXS remains around chi2 5.23.
- The raw all-atom best in the toy generator remains physically invalid under
  heavy-clash QC.
- STARLING's heavy-QC status is intentionally blank in the comparison table
  because STARLING has only been scored as a CA-only ensemble so far.

Verification:

```bash
/home/users/lz280/miniconda3/envs/confdiff/bin/python -m py_compile \
  slurm/wrappers/cohere_sasdxc6_q_sweep_report.py \
  slurm/wrappers/cohere_sasdxc6_compare_baselines.py

/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_q_sweep_report.py \
  --out-root slurm/work/starling_sasdxc6_ca_s48_benchmark \
  --allow-missing-score-sets

/home/users/lz280/miniconda3/envs/confdiff/bin/python \
  slurm/wrappers/cohere_sasdxc6_compare_baselines.py
```

All checks completed successfully.

### 10.35 CA-only scoring fix and compact sweep result (2026-07-08)

The first CA-first generator sweeps initially failed in `saxs_md` before any
useful SAXS comparison could be made.  The failure was:

```text
saxs_md did not write a nonempty profile
Unable to recognize atom
```

The root cause was the generated PDB atom-name fixed-column formatting.  The
wrapper wrote CA atoms as `  CA`; STARLING/standard PDB CA-only files use
` CA `.  `saxs_md` is sensitive to this field.  The fix was added in:

```text
slurm/wrappers/cohere_sasdxc6_real_benchmark.py
```

Changes:

- corrected PDB atom-name alignment for one-letter elements;
- added `--saxs-ca-only-input`;
- when enabled, the wrapper writes `saxs_ca_pdbs/` and
  `saxs-ca-pdb-list.txt`, and scores that temporary CA-only list while leaving
  the original backbone PDB list and metadata unchanged.
- extended `slurm/wrappers/cohere_sasdxc6_run_status.py` to recognize these
  root-level CA sweep outputs in addition to the PDBFixer benchmark layout.

Verification:

```bash
python -m py_compile slurm/wrappers/cohere_sasdxc6_real_benchmark.py
python -m py_compile slurm/wrappers/cohere_sasdxc6_run_status.py

python slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  --out-dir /tmp/cohere_sasdxc6_ca_only_smoke \
  --samples 2 \
  --modes coil,compact \
  --self-avoidance-modes coil,compact \
  --self-avoidance-attempts 8 \
  --run-saxs \
  --saxs-ca-only-input \
  --timeout-seconds 60

python slurm/wrappers/cohere_sasdxc6_run_status.py \
  --out-root slurm/work/cohere_sasdxc6_ca_sweep_compactonly_s48
```

The smoke run completed and wrote `/tmp/cohere_sasdxc6_ca_only_smoke/summary.tsv`.

After the fix, three 48-conformer CA-only sweeps were rerun:

```text
case             modes                 attempts  best_chi2   best_index  best_mode  best_Rg_A    best_CA_pass
noexpanded       coil,compact          256       2.2420141   6           coil       29.69522212  yes
compactonly      compact               512       2.020890157 18          compact    29.71557931  yes
compact2_coil1   compact,compact,coil  256       2.541070172 26          coil       39.12370395  yes
STARLING_CA      external              n/a       2.292908052 27          external   39.86256685  yes
```

The strongest result from this small run is the pure compact proposal:

```text
out_root: slurm/work/cohere_sasdxc6_ca_sweep_compactonly_s48
samples: 48
Rg min/mean/max: 17.84543862 / 27.96327724 / 48.9189812 A
SAXS best chi2: 2.020890157
SAXS best index: 18
best proposal attempts: 382
best CA min nonlocal distance: 3.362397918 A
best CA nonlocal clash count: 0
```

This is the first local generator result that beats the 48-frame STARLING CA
baseline on the current `saxs_md`/SASDXC6 chi2 metric.  It should not yet be
treated as a final scientific result because it is one seed and only 48 samples,
but it changes the engineering direction: the immediate next target should be
multi-seed CA/CG benchmarking of compact-biased proposals, not more PDBFixer
side-chain rescue.

CA shape diagnostics were also generated:

```text
slurm/work/cohere_sasdxc6_ca_sweep_shape_diagnostics/ca_shape_summary.tsv
slurm/work/cohere_sasdxc6_ca_sweep_shape_diagnostics/ca_shape_by_seqsep.tsv
slurm/work/cohere_sasdxc6_ca_sweep_shape_diagnostics/ca_shape_delta_vs_reference.tsv
```

Key shape comparison against STARLING CA:

```text
case           CA-adjacent outliers  nonlocal pair mean A  Rg mean A    CA pass fraction
STARLING_CA    73.875/frame          41.46464618           32.04144579  0.75
noexpanded     0/frame               40.3816968            31.34196669  0.8541666667
compactonly    0/frame               36.02061341           27.96327724  0.875
compact2_coil1 0/frame               38.66429288           29.9953268   0.75
```

Interpretation:

- STARLING remains the stronger established prior-work baseline, but its local
  CA geometry is not peptide-like under the adjacent-distance diagnostic.
- The local generator now has valid CA geometry and can match or beat STARLING
  on this seed's SAXS chi2 when compact-biased.
- Pure compact is probably over-compact relative to STARLING in global pair
  distances, so the next sweep should vary compact strength and seed count
  rather than lock in this exact setting.
- All-atom should remain a downstream validation step for selected CA frames.

### 10.36 Compact-only multi-seed CA sweep launch (2026-07-08)

The immediate follow-up to the single-seed compact result is now launched:
verify that the `compactonly` chi2 2.020890157 result is not a lucky
single-seed/48-conformer hit.

Code changes:

- added `slurm/wrappers/cohere_sasdxc6_ca_sweep_grid.py`
  - plans/runs CA-only SAXS sweep cases;
  - writes `ca_sweep_grid_manifest.tsv`;
  - writes per-seed `ca_sweep_grid_summary.tsv`;
  - writes per-case `ca_sweep_grid_case_summary.tsv`;
  - reports best chi2, best index, best mode/Rg, proposal attempts, CA
    nonlocal clash stats, and profile counts.
- extended `slurm/scripts/run_cohere_sasdxc6_real_saxs_benchmark.slurm`
  to run generated-backbone CA-only SAXS jobs via environment variables:
  - `COHERE_SASDXC6_GENERATE_BACKBONE`
  - `COHERE_SASDXC6_RUN_SAXS`
  - `COHERE_SASDXC6_SAXS_CA_ONLY_INPUT`
  - `COHERE_SASDXC6_MODES`
  - `COHERE_SASDXC6_SELF_AVOIDANCE_MODES`
  - `COHERE_SASDXC6_SELF_AVOIDANCE_ATTEMPTS`
  - `COHERE_SASDXC6_SELF_AVOIDANCE_MIN_CA_DISTANCE`
  - `COHERE_SASDXC6_SELF_AVOIDANCE_MAX_CA_CLASHES`
  - `COHERE_SASDXC6_SELF_AVOIDANCE_LOCAL_SEPARATION`
  - `COHERE_SASDXC6_TIMEOUT_SECONDS`
- updated `cohere_sasdxc6_real_benchmark.py` seed parsing to accept either
  comma-separated or colon-separated seeds.  Colon-separated seeds are needed
  for direct Slurm `--export` usage because commas delimit exported variables.

Verification:

```bash
python -m py_compile \
  slurm/wrappers/cohere_sasdxc6_real_benchmark.py \
  slurm/wrappers/cohere_sasdxc6_ca_sweep_grid.py

python slurm/wrappers/cohere_sasdxc6_ca_sweep_grid.py \
  --work-root slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2 \
  --samples 96 \
  --seeds 20260707:20260708:20260709 \
  --case compactonly512:compact:compact:512:3.0:0:2 \
  --dry-run
```

Historical launch record (superseded by the final three-seed result in
§10.37):

```text
job_id: 12017457
case: compactonly512
out_root: slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/compactonly512
samples: 96 per seed
seeds: 20260707, 20260708, 20260709
modes: compact
self-avoidance: compact, attempts=512, min CA distance=3.0 A, max clashes=0
SAXS: saxs_md via --saxs-ca-only-input, qmax=0.5, stride=10
final status: completed; see §10.37
```

Historical in-flight checkpoint (kept only for provenance; no longer the
current status):

```text
job_id: 12017457
status at that historical checkpoint: still running on grisman/fennario-01
root summary.tsv: not written yet because the job was launched before per-seed incremental summary writes were added
completed seed: 20260707
profiles for seed 20260707: 96/96
seed 20260707 best chi2: 1.810623035
seed 20260707 best index: 84
seed 20260707 best Rg: 31.93665646 A
seed 20260707 best end-to-end: 88.10926126 A
seed 20260707 best proposal attempts: 108
seed 20260707 best CA min nonlocal distance: 3.116043812 A
seed 20260707 best CA nonlocal clash count: 0
seed 20260707 CA self-avoidance pass: yes
seed 20260707 CA-pass count: 85/96
```

This is encouraging because the first larger-N seed beats the previous
single-seed/48 compact result (`2.020890157`) and the STARLING 48-frame CA
baseline (`2.292908052`) on the current `saxs_md` chi2 metric.  It is still not
the final multi-seed answer: the acceptance criterion for this check is the
three-seed case summary after seeds `20260708` and `20260709` complete.

There is also a superseded job `12017453` from an initial direct `sbatch`
attempt where the comma-separated seed list was split by Slurm export parsing.
That job should be treated only as an incidental 1-seed/96 run, not as the
multi-seed validation result.

### 10.37 CA sweep grid tooling and compact-only multi-seed result (2026-07-08)

The compact-only 96 x 3 validation run finished.  The single-seed/48 result was
not a lucky one-off on the current `saxs_md`/SASDXC6 CA-only metric: all three
96-conformer compact seeds beat the STARLING 48-frame CA baseline chi2
`2.292908052`.

Final grid reports:

```text
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_manifest.tsv
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_summary.tsv
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_case_summary.tsv
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_status.tsv
```

Final per-seed SAXS summary:

```text
seed      profiles  best_chi2   best_index  CA_passed
20260707  96/96     1.810623035 84          85/96
20260708  96/96     1.694473488 17          87/96
20260709  96/96     2.036518114 24          80/96
```

Final case summary:

```text
case: compactonly512
seeds done/ok/running/failed/missing: 3/3/0/0/0
best chi2 min/mean/median/max: 1.694473488 / 1.8472 / 1.810623035 / 2.036518114
best seed/index: 20260708 / 17
mean Rg mean: 27.4329 A
mean CA-pass fraction: 0.875
mean proposal attempts: 186.542
```

New code work completed:

- `cohere_sasdxc6_ca_sweep_grid.py`
  - supports `--parallel-seeds`, where each case x seed is an independent
    Slurm/local run rooted under `CASE/seed_runs/seed_SEED`;
  - manifest rows now record real seed/job granularity, command, q settings,
    timeout, and Slurm job id after submit;
  - supports `--rerun-missing` and `--rerun-failed` for seed-level recovery;
  - writes `ca_sweep_grid_status.tsv` with profiles done/expected, seed status,
    best chi2/index, summary existence, and stale-summary flag;
  - recovers summaries directly from `metadata.csv` + `saxs_values.csv` when
    root `summary.tsv` is absent or incomplete;
  - optionally runs CA shape diagnostics with `--shape-diagnostics`.
- `cohere_sasdxc6_recover_summary.py`
  - standalone recovery tool for rebuilding `summary.tsv/json` from existing
    per-seed `metadata.csv` and `saxs_values.csv`.
- `cohere_sasdxc6_export_best_frames.py`
  - exports best/top-k PDBs per case/seed;
  - writes metadata sidecars and optional `pdbfixer-input-list.txt`;
  - supports `--passed-ca-only` for all-atom validation candidate lists.

Generated shape diagnostics:

```text
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_shape_diagnostics/ca_shape_summary.tsv
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_shape_diagnostics/ca_shape_by_seqsep.tsv
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/ca_sweep_grid_shape_diagnostics/ca_shape_delta_vs_reference.tsv
```

Key shape result versus STARLING CA-passed reference:

```text
case                          nonlocal_pair_mean_A  lt8_contact_fraction  lt12_contact_fraction  Rg_mean_A  CA_pass_fraction
starling_ca_passed            43.51761445           0.006611292242        0.0478790388           33.69784479 1
compactonly512_seed_20260707  36.26907894           0.02718424934         0.0903748124           28.21519408 0.8854166667
compactonly512_seed_20260708  34.96523962           0.02823707499         0.09384703462          27.13814713 0.90625
compactonly512_seed_20260709  34.78109547           0.02775837085         0.09268908813          26.94550465 0.8333333333
```

Interpretation:

- The compact proposal is robust on this metric across three seeds.
- It is also systematically more compact than STARLING by CA pair-distance
  diagnostics and has more short nonlocal contacts.  The current best chi2
  should therefore be treated as a SAXS-fit lead, not as a final ensemble prior.
- The next scientific sweep should vary compactness/contact strength and keep
  reporting both chi2 and CA shape diagnostics.

Top-frame exports:

```text
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/best_frame_export_top3/
slurm/work/cohere_sasdxc6_ca_compactonly_ms96_v2/best_frame_export_top3_ca_passed/
```

The unrestricted top-3 export includes the raw best frame for seed `20260709`
(`index 24`, chi2 `2.036518114`), but that frame has one CA nonlocal clash.  The
`best_frame_export_top3_ca_passed` directory is the safer input list for
downstream all-atom validation.

### 10.38 Mixed compactness CA sweep result (2026-07-08)

The mixed compactness sweep finished all 15 case x seed jobs.  Its purpose was
to test whether the pure-compact SAXS win from §10.37 can be retained while
reducing the over-compact CA shape signature.  This is also the bridge from the
quick CA `compact` mode back to WMB ideas 2/5/6: the CA modes are only
prototype proposal components, but the comparison tells which
compactness/contact mixture should be promoted into the WMB portfolio.

Final grid:

```text
work_root: slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1
samples: 256 per seed
seeds: 20260707, 20260708, 20260709
layout: one Slurm job per case x seed
SAXS: saxs_md CA-only input, qmax=0.5, stride=10, max points=200
```

Cases:

```text
case                 modes                 self_avoidance_modes  attempts  min_CA_A  max_CA_clashes
compact512           compact               compact               512       3.0       0
compact1024          compact               compact               1024      3.0       0
coil_compact256      coil,compact          coil,compact          256       3.0       0
compact2_coil1_512   compact,compact,coil  compact,coil          512       3.0       0
expanded_compact256  expanded,compact      compact               256       3.0       0
```

Submitted Slurm jobs:

```text
compact512:          12017472, 12017473, 12017474
compact1024:         12017475, 12017476, 12017477
coil_compact256:     12017478, 12017479, 12017480
compact2_coil1_512:  12017481, 12017482, 12017483
expanded_compact256: 12017484, 12017485, 12017486
```

Final reports:

```text
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/ca_sweep_grid_manifest.tsv
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/ca_sweep_grid_summary.tsv
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/ca_sweep_grid_status.tsv
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/ca_sweep_grid_case_summary.tsv
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/ca_sweep_grid_shape_diagnostics/ca_shape_summary.tsv
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/ca_sweep_grid_shape_diagnostics/ca_shape_delta_vs_reference.tsv
```

Final case-level SAXS summary:

```text
case                 seeds ok/running/failed/missing  best chi2 min/mean/median/max                best seed/index  mean Rg A  mean CA-pass fraction  mean proposal attempts
compact512           3/0/0/0                           1.673747235 / 1.7185 / 1.694473488 / 1.787283568  20260707 / 98   27.3997    0.865886               186.576
compact1024          3/0/0/0                           1.673747235 / 1.7185 / 1.694473488 / 1.787283568  20260707 / 98   27.3785    0.988281               244.781
coil_compact256      3/0/0/0                           1.694473488 / 1.76527 / 1.787283568 / 1.814059846 20260708 / 17   31.2384    0.825521               62.2031
compact2_coil1_512   3/0/0/0                           1.787283568 / 1.81965 / 1.810623035 / 1.861028604 20260709 / 97   30.0857    0.915365               127.04
expanded_compact256  3/0/0/0                           1.694473488 / 1.8349 / 1.787283568 / 2.022945752  20260708 / 17   50.2404    0.822917               60.3646
```

Shape diagnostics, averaged by case and compared to the STARLING CA-passed
reference:

```text
case                 n  nonlocal_pair_mean_A  lt8_contact_fraction  lt12_contact_fraction  metadata_Rg_mean_A  CA_pass_fraction
starling_ca_passed   1  43.5176145            0.006611292           0.047879039            33.6978448          1.000000
compact512           3  35.3030053            0.027606754           0.092340976            27.3997209          0.865885
compact1024          3  35.2858780            0.027310597           0.091739160            27.3784773          0.988281
coil_compact256      3  40.2336139            0.019057398           0.072571992            31.2383804          0.825521
compact2_coil1_512   3  38.7425667            0.021634880           0.078574189            30.0857135          0.915365
expanded_compact256  3  63.7993635            0.014074872           0.055789247            50.2403633          0.822917
```

Top CA-passed frames were exported for downstream all-atom validation:

```text
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/best_frame_export_top3_ca_passed/
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/best_frame_export_top3_ca_passed/best_frame_export.tsv
slurm/work/cohere_sasdxc6_ca_mixed_compact_ms256_v1/best_frame_export_top3_ca_passed/pdbfixer-input-list.txt
```

Interpretation:

- Pure compact remains the numerical best on the current `saxs_md` CA-only
  metric.  Increasing compact self-avoidance attempts from 512 to 1024 improves
  CA-pass fraction from `0.865886` to `0.988281`, but it does not change the
  best chi2 distribution or solve the over-compact shape signature.
- `coil_compact256` is the best scientific compromise in this sweep: it keeps
  STARLING-beating chi2 (`min=1.694473488`, `mean=1.76527`) while moving the
  nonlocal pair mean from about `35.3 A` toward STARLING's `43.5 A`, and lowering
  short-contact inflation relative to pure compact.
- `compact2_coil1_512` also improves shape relative to pure compact, but its
  chi2 is worse than `coil_compact256`.
- `expanded_compact256` is too expanded by Rg/pair-distance diagnostics.  It can
  still hit good chi2 for one seed because that seed selects a compact frame,
  but its ensemble-level shape is not a good next prior.

The next CA-level sweep should refine `coil_compact256` rather than pure
compact: vary the coil/compact mixture ratio, compact self-avoidance strength,
and optional Rg/contact-window constraints, while retaining the same
case-level chi2 + CA shape diagnostics.  Pure compact should remain a control
and an upper bound for this SAXS metric, not the primary ensemble prior.

### 10.39 Ensemble-average correction, FoXS cross-check, and clean holdout (2026-07-18)

The interpretation in §10.38 used the best individual conformer chi2 as its
main SAXS ranking. That is useful for coverage screening, but it is not the
quantity measured by SAXS. The experimental curve is a population-average
intensity, so the primary generator score must first average calculated
intensities and then fit/score that one ensemble curve:

```text
wrong primary ranking:
min_c chi2(I_c(q), I_exp(q))

ensemble ranking:
chi2(sum_c w_c I_c(q), I_exp(q))
```

The new backend-independent report is:

```text
slurm/wrappers/cohere_saxs_ensemble_report.py
slurm/wrappers/tests/test_cohere_saxs_ensemble_report.py
```

It:

- aligns each calculated profile to the experimental q grid;
- averages raw intensities with explicit normalized frame weights;
- fits one intensity scale after averaging;
- reports ensemble chi2, RMS relative error, ESS, Rg, CA pair distances,
  contact fractions, and CA-pass fraction;
- supports a named profile subdirectory, so the same metadata and conformers
  can be compared under `saxs_md`, FoXS, or a future backend;
- distinguishes proposal-mixture ESS from a PAC certificate or fitted BME
  posterior.

Re-analysis of the existing 3 x 256 `coil_compact256` library changed the
scientific picture. Uniform 1:1 ensemble scores were:

```text
backend: saxs_md CA-only, qmax=0.5 A^-1

seed       frames  uniform ensemble chi2
20260707   256     2.462481331
20260708   256     2.452188119
20260709   256     2.355892962
pooled     768     2.420421304
```

Scanning only the total coil/compact weight, while keeping frames uniform
within each mode, gave:

```text
seed       best compact weight  ensemble chi2  ESS fraction
20260707   0.100                1.773400802    0.609756
20260708   0.125                1.763038317    0.640000
20260709   0.125                1.771643519    0.640000
pooled     0.100                1.763416765    0.609756
```

The pooled 10% compact ensemble has:

```text
weighted Rg:                 34.5340 A
weighted nonlocal pair mean: 44.4442 A
weighted contact <8 A:       0.0118272
weighted contact <12 A:      0.0560604
weighted CA-pass fraction:   0.965104
```

This is substantially closer to the STARLING CA-passed shape scale
(`Rg=33.6978 A`, nonlocal pair mean `43.5176 A`) than the uniform 1:1 library.
The optimum is not invariant to q-window: pooled `qmax=0.2` and `0.3 A^-1`
selected 15% compact, while `qmax=0.5 A^-1` selected 10%. Therefore 10-15%
is a screening range, not yet a fixed physical population estimate.

#### Independent FoXS backend

IMP FoXS 2.24.0 was installed in:

```text
/home/users/lz280/miniconda3/envs/cohere-saxs
```

The new wrapper:

```text
slurm/wrappers/cohere_foxs_batch_observable.py
slurm/wrappers/tests/test_cohere_foxs_batch_observable.py
```

isolates FoXS's input-adjacent output in temporary directories, normalizes
profiles to `sample-NNNNNN.dat`, records backend/version/q/scoring metadata,
and supports both CA/residue (`-r`) and all-atom modes. It intentionally
calculates fixed per-frame profiles and fits only the common intensity scale;
it does not optimize hydration parameters independently for every frame.

All 768 existing frames were recalculated with FoXS CA mode. The cross-check
selected the same coil-rich regime:

```text
seed       best compact weight  FoXS ensemble chi2  ESS fraction
20260707   0.050                1.501381378         0.552486
20260708   0.075                1.475285013         0.580552
20260709   0.100                1.521444886         0.609756
pooled     0.075                1.493162325         0.580552
```

Absolute chi2 values should not be compared across calculators because their
forward models differ. The useful cross-check is that both independently
reject 1:1/pure-compact weighting and select a small compact component:
approximately 7.5% with FoXS and 10% with `saxs_md`.

Reproducible outputs:

```text
slurm/work/cohere_sasdxc6_ensemble_mode_fraction_v1/
slurm/work/cohere_sasdxc6_ensemble_qrobust_v1/
slurm/work/cohere_sasdxc6_foxs_ca_ensemble_v1/
```

#### First attempted holdout and detected confound

An independent 6-ratio x 3-seed scan completed all 18 jobs and 4,320 SAXS
profiles:

```text
slurm/work/cohere_sasdxc6_ca_coilrich_ms240_v1/
jobs: 12139949-12139966
```

Its raw fixed-mixture result favored pure coil. That result must not be used as
a mixture validation, because the scan accidentally changed two variables:

```text
old coil_compact256:
  self-avoidance modes = coil,compact
  coil CA-pass fraction = 1.0

new ratio scan:
  self-avoidance modes = compact only
  coil CA-pass fraction ~= 0.239
```

The new scan's compact frames had approximately 94-99% CA pass while its coil
frames had only approximately 23-24% pass. Compact fraction was therefore
confounded with geometry QC and random-stream changes. Keeping this failed
validation in the record is important: its completion does not make it valid
evidence.

The grid summarizer now also writes a true pooled ensemble score across all
completed seeds, rather than only min/mean/median/max of per-seed scores:

```text
ensemble_pooled_frames
ensemble_pooled_chi2
ensemble_pooled_scale
ensemble_pooled_rms_relative
```

#### Corrected high-QC sample holdout

A cleaner validation completed as one balanced proposal library per seed, with
both modes under identical self-avoidance settings. Mixture fractions were
evaluated on the same frames, removing case-to-case Monte Carlo and QC
confounding:

```text
work root: slurm/work/cohere_sasdxc6_ca_balanced_highqc_ms512_v1
samples: 512 per seed (coil,compact alternating)
seeds: 20260713, 20260714, 20260715
self-avoidance modes: coil,compact
attempts: 1024
min nonlocal CA distance: 3.0 A
max CA clashes: 0
jobs: 12141318, 12141319, 12141320
final status: all COMPLETED, exit code 0
profiles: 512/512 per seed, 1,536/1,536 total
```

The retry policy produced a nearly matched valid library:

```text
mode      CA passes / frames  pass fraction  mean attempts/frame
coil      768 / 768           1.000000       4.21
compact   761 / 768           0.990885       251.80
all       1529 / 1536         0.995443       -
```

The high compact retry cost is itself a warning: this is a quality-controlled
holdout library, not a tractable exact-logQ sampler. The rejection normalizer
is unknown.

The predeclared 0/5/7.5/10/12.5/15/20% compact scan gave:

```text
backend: saxs_md CA-only, qmax=0.5 A^-1

case            best compact weight  ensemble chi2  ESS fraction
seed_20260713   0.100                1.714932152    0.609756
seed_20260714   0.150                1.668687076    0.671141
seed_20260715   0.150                1.739700040    0.671141
pooled          0.125                1.701353805    0.640000
```

The pooled uniform 1:1 library scored `2.302395717`, so the predeclared
coil-rich reweighting improves the independent pooled ensemble substantially.
The pooled 12.5% compact ensemble has `Rg=35.0788 A`, nonlocal CA pair mean
`45.0806 A`, contact `<8 A=0.0118886`, and weighted CA-pass fraction
`0.998861`.

FoXS 2.24.0 independently recalculated all 1,536 frames. A later file-format
audit found that the original high-QC run retained seven-column FoXS
`-p` partial profiles and the generic two-column reader incorrectly treated
their first partial term as the final intensity. Those original FoXS chi2
values are superseded. The conformers were recalculated without `-p`, and the
ensemble report now rejects partial profiles rather than silently reading
them as total curves. MultiFoXS itself consumes the partial terms correctly,
so this correction does not invalidate the separate MultiFoXS baseline.

```text
backend: FoXS CA/residue mode, qmax=0.5 A^-1

case            best compact weight  ensemble chi2  ESS fraction
seed_20260713   0.075                1.450767115    0.580552
seed_20260714   0.125                1.411667843    0.640000
seed_20260715   0.125                1.485332762    0.640000
pooled          0.100                1.442460776    0.609756
```

The corrected FoXS pooled uniform 1:1 library scored `2.170575066`. Absolute chi2 values
are not compared across calculators; the cross-backend evidence is the stable
regime: pooled 10% compact under FoXS and 12.5% under `saxs_md`, with every
seed choosing 7.5-15%.

The low-q sensitivity remains real:

```text
backend   qmax A^-1  pooled best compact weight  pooled chi2
saxs_md   0.2        0.20 (grid endpoint)        1.216009109
saxs_md   0.3        0.20 (grid endpoint)        1.327575708
FoXS      0.2        0.15                        1.152737706
FoXS      0.3        0.15                        1.196343223
```

For `saxs_md`, 20% is a lower-bound statement within the registered grid, not
an interior optimum. The defensible conclusion is therefore not a physical
population estimate. It is: a small but nonzero compact proposal component is
reproducibly required, and the preferred weight depends on forward model and
q-window.

Reproducible holdout reports:

```text
slurm/work/cohere_sasdxc6_ca_balanced_highqc_ms512_v1/
  ensemble_saxs_md_qmax0p2_v1/
  ensemble_saxs_md_qmax0p3_v1/
  ensemble_saxs_md_qmax0p5_v1/
  ensemble_foxs_fixed_v2_qmax0p2/
  ensemble_foxs_fixed_v2_qmax0p3/
  ensemble_foxs_fixed_v2_qmax0p5/
```

The older `ensemble_foxs_qmax*_v1` directories remain only as a provenance
record of the detected partial-profile misuse and must not be cited.

Current focused Python regression suite:

```text
python3 -m unittest discover -s slurm/wrappers/tests -p 'test_cohere_*.py' -v
Ran 19 tests
OK
```

### 10.40 MultiFoXS baseline, honest known-logQ CA proposal, and generic M2 adapter (2026-07-18)

#### QC-filtered MultiFoXS comparison

`cohere_multifoxs_baseline.py` now provides two comparison modes over FoXS
partial profiles:

- sparse MultiFoXS ensembles with 1-4 states;
- the bundled NNLS fit, with parsed nonzero weights, ESS, and maximum weight.

On the older 3 x 256 `coil_compact256` pool, requiring the metadata CA-pass
flag retained 634/768 input profiles. MultiFoXS clustering exposed 404
representatives to NNLS. Results were:

```text
sparse states     1     2     3     4
best chi          1.37  1.12  1.11  1.10

NNLS chi:         1.12
positive weights: 8 / 404 clustered profiles
ESS:              3.8506
max weight:       0.455
weight sum:       0.998 (printed-weight precision)
```

All eight positive NNLS frames pass CA QC. This is a strong pure-fitting
baseline, but its effective ensemble is very sparse and selected on the same
experimental curve. It is neither an out-of-sample population estimate nor a
thermodynamic/PAC ensemble.

Reproducible output:

```text
slurm/work/cohere_sasdxc6_multifoxs_qc_nnls_768_v1/
```

#### Exact/evaluable proposal accounting without rejection

The original CA generator used deterministic mode cycling and, for selected
modes, finite rejection/retry followed by a "best failed" fallback. That
process does not have the iid mixture density assumed by
`KnownLogQSample`, and the rejection acceptance normalizer is unknown.

`cohere_sasdxc6_real_benchmark.py` now has an opt-in `--known-logq` mode that:

- draws one global `coil/expanded/compact` mode iid from declared
  `--mode-weights`;
- draws residue-local Ramachandran regions conditionally and independently;
- evaluates the full mixture density
  `q(c)=sum_m alpha_m product_i q_m(phi_i,psi_i | aa_i)` by log-sum-exp;
- uses wrapped-normal densities on `(-180,180]`, with the numerical image
  cutoff recorded in `proposal.json`;
- writes full φ/ψ arrays and full-precision proposal values to
  `torsions.jsonl`, plus audit columns in `metadata.csv`;
- refuses self-avoidance rejection (`--self-avoidance-modes none` is
  required), while still reporting CA geometry as post-hoc QC.

The legacy default RNG stream was checked byte-for-byte against the former
torsion sampler for 300 fixed-seed draws, so old SAXS artifacts remain
reproducible.

A 3 x 512, 1:1 coil/compact SASDXC6 pilot quantified the cost of honest
unconditioned sampling:

```text
mode      frames  CA-pass fraction  mean Rg A   mean proposal log q
coil      752     0.2660            31.591      -1204.503
compact   784     0.0051            20.856      -1156.537
all       1536    0.1328            26.111      -1180.021
```

This route is now an auditable known-logQ local/latent baseline, not the final
proposal. If CA validity is treated as target hard support, approximately 87%
of draws have zero target weight; compact draws are especially poor. This
empirically confirms that steric/contact information must enter the WMB
proposal rather than be repaired by undocumented rejection.

The Java importance-sampling path now represents that hard support directly:
`KnownLogQSample`, the observable estimator, diagnostics, and the ratio
certificate accept `logTarget=-Infinity` / `logWeight=-Infinity` as an exact
zero weight. They reject NaN, positive infinity, and an all-zero-weight sample
set. Thus QC-failed conformers can be retained in the proposal audit without
silently dropping proposal failures or assigning them artificial finite
target mass.

The diagnostic pairs bootstrap now redraws a replicate if it contains only
zero-weight samples, because that replicate has an undefined
self-normalized-ratio denominator. This conditioning is explicitly
diagnostic; it is not promoted to a finite-sample certificate.

Reproducible output:

```text
slurm/work/cohere_sasdxc6_known_logq_ms512_v1/
```

#### First non-OSPREY M2 adapter

The previous M2 text said a purpose-built IDP factor-graph adapter was needed;
that statement is now partially resolved:

```text
src/main/java/edu/duke/cs/osprey/wmb/WmbModel.java
src/main/java/edu/duke/cs/osprey/cohere/WmbProposalPortfolio.java
src/main/java/edu/duke/cs/osprey/cohere/WmbPilotScan.java
src/main/java/edu/duke/cs/osprey/cohere/IdpBackboneFactorGraph.java
src/main/java/edu/duke/cs/osprey/cohere/IdpBackboneProposalJob.java
src/test/java/edu/duke/cs/osprey/cohere/TestIdpBackboneFactorGraph.java
```

The new path:

- constructs a sparse generic WMB graph directly from unary and pairwise log
  potentials, with identity state maps and true graph non-edges;
- scales any such model into exact-logQ tempered portfolio components without
  an `EnergyMatrix` or `RCs`;
- represents per-residue discrete named `(phi,psi)` states;
- deterministically materializes a sampled assignment as an N/CA/C/O PDB for
  external observable wrappers.

This was initially sampler/materializer plumbing rather than a physical
SASDXC6 model. Section 10.41 resolves the first empirical proposal gate with
an explicit coarse target and a held-out steric surrogate, but not the
physical-target gate. A pairwise surrogate must not be described as an exact
all-atom torsion energy: global Cartesian contacts are not generally pairwise
functions of individual torsion bins.

Targeted verification after the adapter and hard-support changes:

```text
./gradlew test \
  --tests edu.duke.cs.osprey.cohere.TestIdpBackboneFactorGraph \
  --tests edu.duke.cs.osprey.cohere.TestObservableEstimator \
  --tests edu.duke.cs.osprey.wmb.TestWeightedMiniBucket \
  -DtestMaxHeap=512m --max-workers=1 --no-daemon

TestIdpBackboneFactorGraph: 3 passed
TestObservableEstimator:    34 passed
TestWeightedMiniBucket:     18 passed
total:                      55 passed, 0 failed
BUILD SUCCESSFUL
```

### 10.41 First real-sequence steric-aware WMB proposal pilot (2026-07-18)

The first non-OSPREY adapter is now exercised on the full 116-residue
SASDXC6 sequence rather than only on a tiny factor graph:

```text
src/main/java/edu/duke/cs/osprey/cohere/IdpBackboneProposalJob.java
slurm/wrappers/cohere_sasdxc6_discrete_steric_factors.py
slurm/wrappers/tests/test_cohere_sasdxc6_discrete_steric_factors.py
```

The experiment deliberately separates target and proposal:

- the target is an iid three-state `(PPII, beta, alpha)` Ramachandran base
  measure with probabilities `(0.34, 0.285, 0.375)`, multiplied by exact zero
  when deterministic CA QC finds any nonlocal distance below `3 A`;
- those probabilities are a coarse coil/compact-average stress-test prior,
  not a claimed physical SASDXC6 population;
- the proposal adds shared state-count and lag-1-through-8 state-pair
  logistic factors trained to predict CA validity;
- the learned factors affect proposal sampling only. Every sample is scored
  against the separate base target and real deterministic geometry, with
  exact mixture `logQ`;
- the WMB portfolio contains `iBound={2,3}`,
  `beta={0,0.5,1}`, and natural/reverse orders. The newly supported `beta=0`
  component is the exact uniform infinite-temperature endpoint.

The steric surrogate used 12,000 training conformers and an independent
4,000-conformer validation set. CA-pass fractions were `0.0539` and `0.0605`;
classifier AUC was `0.749` in training and `0.727` in validation. This is
moderate held-out signal, not a solved steric model.

All proposal strengths below used the same independent 1,024-sample pilot
seed:

```text
proposal scale  CA-pass fraction  ESS/N       max weight share  weighted CA Rg A
baseline        0.05469           0.05469     0.01786           32.328
0.25            0.08789           0.06264     0.02817           32.517
0.50            0.08789           0.05486     0.05078           31.964
1.00            0.10059           0.04859     0.06579           32.622
2.00            0.12109           0.03205     0.08029           32.797
```

The `0.25` surrogate is the first real-sequence WMB result that improves both
coverage and total-sample overlap: relative to the exact local baseline, CA
pass rises by about `61%` and ESS/N by about `14.5%`. Stronger factors keep
raising the pass fraction but lower ESS and worsen tail concentration. This
is direct evidence for pilot gating rather than maximizing acceptance.

This is still a method-development result, not a paper claim:

- scale selection examined the same pilot seed, so a predeclared `0.25`
  multi-seed replication is required;
- the target is a coarse discrete prior plus CA hard support, not an
  all-atom/thermodynamic ensemble;
- these WMB samples have not yet been scored as ensemble-average SAXS curves;
- exact zero weights make the observed log-weight range infinite, and a
  useful deterministic weight cap/certificate remains missing.

Reproducible training, configs, sample tables, and summaries:

```text
slurm/work/cohere_sasdxc6_discrete_steric_factor_v1/
```

### 10.42 Locked eight-seed WMB replication result (2026-07-18)

The `0.25` selection result above is now separated from its evaluation before
any replication output is generated. The machine-readable plan is:

```text
slurm/wrappers/cohere_discrete_wmb_replication.py
slurm/wrappers/tests/test_cohere_discrete_wmb_replication.py
slurm/work/cohere_sasdxc6_discrete_steric_factor_v1/
  replication_scale_0p25_v1/replication_plan.json
```

The locked protocol excludes selection seed `20260721` and evaluates
baseline versus the already selected `scale=0.25` proposal at seeds
`20260722` through `20260729`, with `1,024` samples per arm and seed. Config
hashes, target/QC equality, sample counts, and output paths are recorded in
the plan. Arms are matched by replicate seed identifier; this is not a claim
of samplewise common random numbers.

The sole primary endpoint is within-seed `ESS/N` difference
(`scale=0.25 - baseline`). The engineering replication gate was fixed as:

```text
median paired ESS/N difference > 0
and at least 7 of 8 paired differences > 0
```

Exact directional and two-sided sign-test p-values will be reported, but this
gate is an engineering reproducibility rule rather than a formal statistical
certificate. QC pass fraction, maximum weight share, weighted CA Rg, pooled
importance diagnostics, and cross-seed dispersion are secondary. No seed or
scale will be changed after results are seen.

QC-pass PDBs are written during the locked runs so their exact importance
weights can feed a subsequent ensemble-average SAXS/FoXS analysis. That SAXS
analysis is explicitly exploratory and is not part of this proposal gate.
The immutable plan records `predeclared_not_run` because it is the pre-result
snapshot; results are kept separately in `replication_result.json`.

All 16 planned runs subsequently completed without config-hash or sample-count
drift:

```text
seed      baseline ESS/N  scale=.25 ESS/N  paired delta
20260722  0.05078         0.05439           +0.00361
20260723  0.05469         0.06225           +0.00756
20260724  0.04883         0.06049           +0.01166
20260725  0.04785         0.05569           +0.00784
20260726  0.05957         0.05104           -0.00853
20260727  0.07129         0.05511           -0.01617
20260728  0.06543         0.06225           -0.00318
20260729  0.05664         0.07046           +0.01382
```

The median paired difference was `+0.00559`, but only `5/8` pairs were
positive. The exact directional sign-test value was `p=0.3633` (two-sided
`p=0.7266`). Therefore the locked `7/8` engineering gate **failed**. The
single-seed claim that `scale=0.25` reproducibly improves total-sample overlap
is rejected and must not be carried forward.

The secondary diagnostics explain the failure rather than hiding it:

```text
pooled over 8,192 samples/arm  baseline       scale=.25
CA-pass fraction               0.05688        0.08752
ESS/N                          0.05688        0.05842
maximum weight share           0.00215        0.00800
log target-normalizer estimate -2.86673       -2.82902
relative SE of normalizer      0.04499        0.04436
weighted CA Rg (A)             31.918         32.073
```

The learned proposal raises feasible-state coverage in every seed and by
about `53.9%` in the pooled sample, but converts most of that gain into
unequal weights. Its pooled ESS/N gain is only about `2.7%`, while its largest
pooled weight share is about `3.7x` baseline. Weighted Rg and target
normalizer estimates remain mutually compatible at this Monte Carlo
resolution. The supported conclusion is therefore narrower: the sparse
steric factors contain reusable feasibility signal, but the fixed tempered
mixture does not yet turn that signal into a robust efficiency gain.

Reproducible paired rows and machine-readable diagnostics:

```text
slurm/work/cohere_sasdxc6_discrete_steric_factor_v1/
  replication_scale_0p25_v1/paired_results.tsv
  replication_scale_0p25_v1/replication_result.json
```

### 10.43 Exact-weight ensemble-average FoXS follow-up (2026-07-18)

The locked runs wrote `1,183` finite-weight/QC-pass PDBs: `466` from the
baseline arm and `717` from `scale=0.25`. A new auditable adapter maps each
full-precision `logTarget-logQ` value to exactly one PDB and FoXS profile:

```text
slurm/wrappers/cohere_wmb_saxs_report.py
slurm/wrappers/tests/test_cohere_wmb_saxs_report.py
slurm/work/cohere_sasdxc6_discrete_steric_factor_v1/
  replication_scale_0p25_v1/saxs_foxs_ca_v1/
```

IMP FoXS 2.24.0 ran in residue/CA mode once per PDB through a shared
`q<=0.5 A^-1` profile protocol. For each arm and seed, intensities were first
averaged with the already fixed self-normalized importance weights. Only
after that ensemble curve was fixed was one global intensity scale fitted
against SASDXC6. No per-conformer chi2 entered a weight, and SAXS did not
change the target, proposal, or failed replication decision.

Results across three q-range diagnostics:

```text
q max  pooled chi2          cross-seed normalized-curve RMS relative SD
       baseline  scale=.25  baseline  scale=.25
0.2    1.6886    1.6462     0.02786   0.03015
0.3    1.6692    1.6110     0.02531   0.02627
0.5    1.8942    1.8460     0.02165   0.02180
```

The pooled baseline and WMB curve shapes differ by only `0.48%`, `0.59%`,
and `0.51%` RMS at q maxima `0.2`, `0.3`, and `0.5 A^-1`, respectively. This
is useful end-to-end evidence that the two exact-logQ routes estimate the
same declared target observable after reweighting.

It is not evidence that WMB improved the likelihood or recovered a physical
ensemble. The WMB curve has lower chi2 in only `5/8` seeds at every reported
q range. Its across-seed chi2 standard deviations are larger:

```text
q max  baseline chi2 SD  scale=.25 chi2 SD
0.2    0.2728            0.3501
0.3    0.2250            0.3321
0.5    0.1770            0.2706
```

The small pooled chi2 reduction is below the observed replicate variability
and was not an optimized endpoint. Together with §10.42, the supported
conclusion is that the surrogate adds feasibility coverage but does not
deliver a reproducible ESS or observable-stability gain. This rules out the
current fixed `scale=0.25` tempered mixture as the paper proposal. Further M2
work must change the proposal construction and pass a new independent gate;
rerunning or relabeling this result is not an option.

Machine-readable weighted curves, paired seed rows, and semantics:

```text
slurm/work/cohere_sasdxc6_discrete_steric_factor_v1/
  replication_scale_0p25_v1/saxs_foxs_ca_v1/per_seed_saxs.tsv
  replication_scale_0p25_v1/saxs_foxs_ca_v1/pooled_saxs.tsv
  replication_scale_0p25_v1/saxs_foxs_ca_v1/saxs_report.json
```

### 10.44 Session closeout and MaxEnt/MAP handoff (2026-07-18)

**Historical handoff:** this sequence was executed on 2026-07-19 and is
superseded by the result in §10.45. It is retained to show that the objective
and convergence corrections were specified before the blocked-CV output was
generated.

The project is not paper-ready and the active research goal remains open.
The present evidence supports a clean separation between two questions:

1. **Sampling efficiency.** The exact-logQ WMB machinery works, but the
   specific lag-1-through-8 steric surrogate, `scale=0.25`, and fixed
   12-component tempered mixture failed its locked replication gate. It is
   no longer the paper proposal.
2. **Observable-conditioned inference.** The next defensible question is
   whether KL-regularized MaxEnt/MAP refinement predicts held-out q regions
   better than matched fixed-library baselines. This estimates a
   prior-conditional ensemble; it does not identify a unique physical
   ensemble.

A retrospective four-block nested-CV implementation has been staged but
**has not been run or validated**:

```text
slurm/wrappers/cohere_saxs_maxent_cv.py
slurm/wrappers/tests/test_cohere_saxs_maxent_cv.py
slurm/configs/cohere_sasdxc6_maxent_blocked_cv_v1.json
slurm/scripts/run_cohere_sasdxc6_maxent_cv.slurm
```

It uses the corrected fixed-total FoXS profiles from §10.39, filters to the
1,529 CA-QC-passing candidates, assigns equal prior mass to coil and compact
modes, freezes the training-fitted intensity scale on each held-out block,
and compares:

- the balanced prior with scale fitting only;
- a two-mode coil/compact maximum-likelihood fit;
- unrestricted NNLS as a high-variance finite-pool baseline;
- KL-regularized MaxEnt/MAP with tau selected only inside each outer
  training split.

This remains method development, not an independent biological test: the
candidate generator and mixture family were historically inspected against
the same SASDXC6 curve, and neighboring SAXS q points are correlated.

Static review identified two required corrections before execution:

- define the regularized objective as
  `0.5 * mean_q(standardized residual^2) + tau * KL(w || w0)`, so the meaning
  of `tau` does not change with the number of inner/outer training points;
  in the current sum-form dual this requires `tau_sum = n_train * tau`;
- describe and diagnose the scale/weight loop as block-coordinate
  optimization, relax its overly strict scale convergence threshold, and
  avoid any claim of a global joint optimum.

Do not interpret or publish output from the staged version before those
corrections and a fresh compute-node validation. Slurm job `12149186`
completed on `fennario-01` with 26 existing Python tests and all 55 targeted
Java tests passing; it began before the three new MaxEnt tests existed, so
those tests are still unvalidated.

The exact next-session sequence is:

1. apply the two MaxEnt corrections above and statically inspect the diff;
2. submit `slurm/scripts/run_cohere_unit_tests.slurm` through `sbatch`;
3. only after the new MaxEnt tests pass, submit
   `slurm/scripts/run_cohere_sasdxc6_maxent_cv.slurm`;
4. report outer-fold predictive chi2, selected tau, ESS/max weight, mode
   mass, and limitations, then update this section and remove superseded
   claims.

**Compute policy:** the login node is for lightweight inspection and editing
only. Python tests, Gradle, FoXS, optimization, and all experiments must run
through Slurm on a compute node.

### 10.45 Corrected MaxEnt/MAP blocked-CV result and environment migration (2026-07-19)

The §10.44 handoff was executed without running computation on the login
node. The MaxEnt/MAP implementation now optimizes the declared objective

```text
0.5 * mean_q(standardized residual^2) + tau * KL(w || w0)
```

by passing `tau_sum = n_train * tau` to the existing sum-form dual. Thus a
fixed `tau` has the same meaning in inner, outer, and full-data fits even when
the number of fitted q points changes. A regression test duplicates every q
point three times and requires the fitted scale and weights to remain
unchanged.

The scale/weight loop is now explicitly a block-coordinate optimization, not
a claimed global joint optimizer. Each iteration records the normalized
primal objective before the update, after the convex fixed-scale weight
block, and after the closed-form scale block. It fails if either block
materially increases the objective and reports scale relative change, weight
L1 change, objective relative change, and the dual gradient. Convergence
requires:

```text
scale relative change     <= 1e-6
weight L1 change          <= 1e-5
objective relative change <= 1e-8
```

The fixed-scale convex dual retains L-BFGS-B as its first solver. A small-tau
unit case exposed an L-BFGS line-search termination, so a trust-region Newton
fallback with an analytic Hessian-vector product was added. A failed solver
status is accepted only when the final maximum absolute dual gradient is at
most `1e-6`; primal monotonicity is still checked independently.

Compute-node validation history:

```text
job 12168571: 30/31 Python tests passed; the small-tau line-search case failed
job 12168602: 5/5 targeted MaxEnt tests passed after the trust-region fallback
job 12168639: 31/31 focused Python tests and all 55 targeted Java tests passed
               BUILD SUCCESSFUL
```

The corrected retrospective blocked nested-CV then ran as Slurm job
`12168648` on `fennario-01`. It used the pre-staged protocol without changing
the candidate library or q blocks:

```text
candidates:       1,529 CA-QC-passing frames
                  768 coil, 761 compact
prior:            equal total coil/compact mass, uniform within each mode
q points:         179 fixed-total FoXS points, q=0.002817..0.499929 A^-1
outer split:      four contiguous blocks (45,45,45,44 points)
tau grid:         0.1, 1, 10, 100, 1000
tau selection:    inner blocked CV excluding the outer block
held-out policy:  fit scale on training q points and freeze it on the test block
```

The config SHA-256 recorded in the result is:

```text
b6742e67d197c2249364ae4eaf2b95ca0f5e4c1196fb456a2cc456ac47e76a06
```

Primary outer-held-out result:

```text
method          pooled chi2  fold-mean chi2  fold SD  mean ESS/N  mean compact mass
balanced prior  2.997430     2.992239         2.153729 0.999979    0.500000
two-mode MLE    2.652953     2.647462         2.118919 0.599293    0.086864
NNLS            2.611670     2.605936         2.460899 0.004083    0.311622
KL-MaxEnt/MAP   4.384891     4.371389         5.398443 0.787662    0.349921
```

The MaxEnt fold details explain the pooled failure:

```text
outer fold  held-out q A^-1    tau  train chi2  test chi2  ESS/N     max weight  compact mass
0           0.0028..0.1267     0.1  1.2291      12.4626    0.46535   0.007206    0.20503
1           0.1296..0.2533     1    1.4016       1.4286    0.84146   0.001406    0.37292
2           0.2561..0.3795     1    1.3567       1.6400    0.85884   0.001444    0.36025
3           0.3823..0.4999    10    1.8194       1.9545    0.98499   0.000823    0.46148
```

All four reported MaxEnt fits converged monotonically in 14-50 block
iterations. Final scale changes were below `9.6e-7`, weight L1 changes below
`5.6e-6`, objective relative changes below `4e-10`, and final maximum dual
gradients were approximately `7e-9` to `1.1e-6`. The negative held-out result
therefore cannot be dismissed as an unconverged optimization.

The predeclared pooled endpoint does **not** support unrestricted
per-conformer KL-MaxEnt/MAP as the paper method. It is worse than the
balanced prior because the low-q extrapolation fold is catastrophic, even
though it improves over that prior on the other three folds. A post-result
edge/interior decomposition is diagnostic only:

```text
method          interior folds 1-2  edge folds 0,3
balanced prior  1.85446             4.15324
two-mode MLE    1.55351             3.76475
NNLS            1.27703             3.96131
KL-MaxEnt/MAP   1.53427             7.26754
```

The two-mode MLE is the current low-complexity baseline lead: it beats the
balanced prior in all four folds, retains mean `ESS/N=0.599`, and is only
slightly worse in pooled chi2 than NNLS, whose mean `ESS/N=0.00408`
corresponds to an effective ensemble of only about six frames. This does not
promote the two-mode fit to a biological or paper claim. The generator and
mixture family were historically inspected on SASDXC6, only four correlated
q blocks exist, and the low-q fold selects essentially pure coil while the
other folds select roughly 10-12.5% compact mass.

The full-data MaxEnt fit selected `tau=1`, training chi2 `1.40697`,
`ESS/N=0.85370`, maximum weight `0.001335`, and compact mass `0.36914`.
Those are descriptive same-data values and must not be substituted for the
outer-held-out result or interpreted as a physical population estimate.

Reproducible outputs:

```text
slurm/work/cohere_sasdxc6_ca_balanced_highqc_ms512_v1/
  maxent_blocked_cv_v1/
    result.json
    outer_fold_results.tsv
    inner_tau_selection.tsv
    summary.tsv
    final_maxent_weights.tsv
```

The result JSON, outer-fold TSV, and summary TSV have SHA-256 values
`34b286c6...c00655`, `92659e39...39382`, and
`abfb5fd8...2a54b`, respectively. Preserve this v1 directory as the
retrospective result; do not overwrite it while developing a successor.

The isolated FoXS environment was also moved off the home filesystem through
Slurm job `12168676`. Its physical prefix is now:

```text
/usr/xtmp/lz280/conda_envs/cohere-saxs
```

The historical prefix
`/home/users/lz280/miniconda3/envs/cohere-saxs` is a compatibility symlink to
that directory. A byte-level rsync checksum dry-run, conda package-list
comparison, Python 3.14.6 launch, and FoXS 2.24.0 version check all passed on
the compute node before the 2.5 GB home copy was deleted.

Evidence-constrained next route:

1. freeze the negative unrestricted-MaxEnt result and retain two-mode MLE and
   sparse NNLS as baselines, not discoveries;
2. quantify candidate-pool, forward-model, and fold-edge sensitivity without
   using those diagnostics to relabel the v1 endpoint;
3. choose and preregister a second biological SAXS case before claiming any
   observable-conditioned method generalizes;
4. if a hierarchical/mode-aware refinement is developed from this failure,
   label it exploratory on SASDXC6 and evaluate it only on frozen or new
   data;
5. continue to report sampling efficiency, observable-conditioned
   prediction, and ensemble identifiability as separate axes. The failed WMB
   efficiency gate in §10.42 and the failed unrestricted-MaxEnt predictive
   gate here cannot be combined into a positive method claim.

### 10.46 `saxs_md` forward-model sensitivity and second-case preregistration (2026-07-19)

The same-library `saxs_md` sensitivity was registered only after the primary
FoXS result in §10.45. It is therefore diagnostic and cannot replace or rescue
that endpoint. The first Slurm attempt, job `12168809`, stopped before any
optimization because the partial-profile gate mistakenly counted the words
in the un-commented `saxs_md` program banner as profile columns. The shared
`first_numeric_column_count()` helper now skips nonnumeric text and counts the
first actual q/I row. A regression test accepts a banner followed by a
two-column curve while the existing seven-column FoXS partial-profile test
still requires rejection.

Slurm job `12168924` validated that correction:

```text
32/32 focused Python tests passed
55/55 targeted Java tests passed
BUILD SUCCESSFUL
```

The original sensitivity config has SHA-256
`caa4d5a143ace3c9d606a994b84b3270811d37c9fa715f0110707c9509e1589b`.
After the format fix, job `12168954` reached an inner MaxEnt fit but stopped
fail-closed at its 100-iteration ceiling. Its final weight and objective
changes met their thresholds; the scale relative change was
`1.14786e-6`, narrowly above the unchanged `1e-6` threshold. It wrote no
fold result. Rather than relax convergence, a new auditable v2 config changed
only the maximum iterations from 100 to 250 and used a new output directory:

```text
slurm/configs/cohere_sasdxc6_maxent_blocked_cv_saxs_md_sensitivity_v2.json
SHA-256: 1a95baa0b2776eff37b5b912fc4b02017c6f2ee77dd874a0d62f6e30f22d0f70
```

Slurm job `12169004` completed v2 with the same 1,529 candidates, reference
curve, q points, four outer blocks, tau grid, baselines, and objective as the
FoXS analysis. Only the pre-existing calculated conformer profiles changed
from fixed-total FoXS to `saxs_md`.

```text
method          pooled chi2  fold-mean chi2  fold SD   mean ESS/N  mean compact mass
balanced prior  3.717399     3.710696         2.996750  0.999979    0.500000
two-mode MLE    3.919969     3.909706         4.019578  0.628960    0.110753
NNLS            6.726127     6.698605        10.009713  0.003753    0.469925
KL-MaxEnt/MAP   7.735309     7.705612        11.262124  0.790207    0.369604
```

MaxEnt again improved over the balanced prior on the three non-low-q folds
but failed catastrophically when extrapolating into the lowest-q block:

```text
outer fold  held-out q A^-1    tau  prior chi2  MaxEnt chi2  ESS/N    max weight  compact mass
0           0.0028..0.1267     0.1  8.1524      24.5849       0.34015  0.012251    0.15298
1           0.1296..0.2533     1    1.5809       1.4339       0.84492  0.001528    0.39279
2           0.2561..0.3795    10    2.5986       2.4137       0.98875  0.000804    0.46579
3           0.3823..0.4999    10    2.5108       2.3899       0.98700  0.000804    0.46685
```

All four outer fits converged monotonically in 14-52 block iterations under
the unchanged thresholds. The post-result interior/edge decomposition is:

```text
method          interior folds 1-2  edge folds 0,3
balanced prior  2.08977             5.36331
two-mode MLE    1.82250             6.04100
NNLS            1.65506            11.85417
KL-MaxEnt/MAP   1.92379            13.61213
```

Thus two conclusions are now supported. First, unrestricted per-conformer
MaxEnt is negative under both forward models; its low-q extrapolation failure
is not attributable only to FoXS. Second, the baseline ranking is itself
forward-model-sensitive: two-mode MLE led the FoXS comparison but the
balanced prior leads under `saxs_md`. Neither result identifies a physical
coil/compact population. The full-data `saxs_md` MaxEnt fit selected
`tau=10`, training chi2 `2.04812`, `ESS/N=0.98985`, maximum weight
`0.000784`, and compact mass `0.47021`; these remain descriptive same-data
values.

Frozen sensitivity outputs:

```text
slurm/work/cohere_sasdxc6_ca_balanced_highqc_ms512_v1/
  maxent_blocked_cv_saxs_md_sensitivity_v2/
    result.json
    outer_fold_results.tsv
    inner_tau_selection.tsv
    summary.tsv
    final_maxent_weights.tsv
```

The result, outer-fold, and summary files have SHA-256 values
`6c137262...579d5`, `b0277cd1...6305`, and `e3e20c52...0537`,
respectively. Preserve both failed v1 job logs and the completed v2 output;
the iteration-only amendment is part of the audit trail.

#### Independent biological case locked before data access

The second case is now preregistered as
[SASDNV6 Early E1A SEC-SAXS](https://www.sasbdb.org/data/SASDNV6/), linked
by SASBDB to the 2022 study
[Conformational buffering underlies functional selection in intrinsically disordered protein regions](https://www.sasbdb.org/project/1609/).
The official entry describes a monomeric approximately 13 kDa Early E1A
protein, inline SEC-SAXS collection at SWING/SOLEIL, a downloadable curve and
FASTA, and no deposited structural model. This choice was recorded after
viewing only catalog metadata and before downloading numeric curve rows or
the FASTA, generating E1A candidates, or calculating any profile.

The machine-readable registration is:

```text
slurm/configs/cohere_sasdnv6_second_case_preregistration_v1.json
SHA-256: 91375a4613cfae3257187319f48aac6f49db0ec27ba12923dd462c473643ec0e
registered: 2026-07-19T11:41:13-04:00
```

It freezes the following primary protocol:

```text
seeds:                 20260720, 20260721, 20260722
candidates per seed:   512, alternating 256 coil / 256 compact
CA QC:                 both modes; 1024 attempts; min nonlocal CA 3.0 A;
                       zero clashes; local separation <=2 excluded
library adequacy:      >=90% CA pass in every seed x mode cell, before SAXS
prior:                 equal total coil/compact mass, uniform within mode
primary forward model: fixed-total IMP FoXS 2.24.0 residue/CA profiles
reference:             deposited sigma where valid; deterministic <=200-point
                       q<=0.5 A^-1 grid; at least 40 retained points
CV/objective:          the corrected four-block nested-CV and mean-q MaxEnt
tau grid:              0.1, 1, 10, 100, 1000
optimizer:             unchanged tolerances, maximum 250 block iterations
```

The generator and QC may not use the catalog Rg/Dmax or curve intensity. The
related non-SEC E1A concentration series (`SASDNN6`, `SASDNP6`, `SASDNQ6`)
is excluded from selection and tuning and may only be examined after the
primary `SASDNV6` result is frozen.

Unrestricted MaxEnt advances only if its pooled held-out chi2 is lower than
both the balanced prior and two-mode MLE, it beats the prior in at least
three of four folds, no fold exceeds twice the matched prior chi2, mean
outer-fit `ESS/N>=0.10`, and the ranking survives preregistered
candidate-pool diagnostics. Passing would justify another frozen test, not a
unique-ensemble claim. If two-mode MLE alone wins, a later mode-aware or
hierarchical method remains exploratory and must be tested on a new case; it
cannot be developed on E1A and then called an E1A validation.

### 10.47 Candidate-pool sensitivity and frozen SASDNV6 build (2026-07-19)

#### SASDXC6 candidate-pool and curve-level identifiability

Before running any pool perturbation, the following registry locked one full
replay, three leave-one-seed-library-out variants, and two disjoint balanced
half-pools:

```text
slurm/configs/cohere_sasdxc6_candidate_pool_sensitivity_v1.json
SHA-256: a400933d5269fe5122a8b4df7f9ed03308b8d170754305b66c01794d58590487
```

The half-pools use candidate indices modulo four: remainders `0,1` versus
`2,3`. Because the original library alternates coil/compact candidates, each
half retains both modes within every seed. The updated runner also writes the
scaled full-data MaxEnt prediction at every q point. This allows comparison
of observable curves between disjoint pools without pretending that their
microscopic weights correspond.

Slurm job `12169350` completed all six fits but encountered a jq reserved-word
error while constructing the aggregate table. No optimization was rerun.
Report-only jobs `12169454` and `12169471` reused the completed result
directories, fixed the variable name, and added fold diagnostics. The full
replay reproduced all four primary pooled chi2 values with exactly zero
difference at the preregistered `1e-8` tolerance.

```text
variant                  N     prior     two-mode  NNLS      MaxEnt    tau  ESS/N    compact
full replay              1529  2.99743   2.65295   2.61167   4.38489   1    0.85370  0.36914
leave out seed 20260713  1018  2.85160   2.61255   2.79802   4.27445   1    0.86181  0.37275
leave out seed 20260714  1020  3.20852   2.91152   4.85617   5.05864  10    0.98798  0.46631
leave out seed 20260715  1020  2.94861   2.46493   2.55499   3.79064   1    0.85144  0.36483
balanced half A           765  2.77873   2.70169   2.73568   4.48988  10    0.99040  0.46972
balanced half B           764  3.23720   2.60951   4.99621   4.19089   1    0.83936  0.35887
```

Unrestricted MaxEnt remains worse than both the prior and two-mode MLE in
every perturbation. The same fold pattern also persists without exception:
MaxEnt improves over the prior on folds 1-3, while its lowest-q fold is
`1.70-2.36x` worse than the matched prior. This makes the failure more
specific than a single unlucky candidate pool: fitting higher-q blocks does
not reliably extrapolate global scale/shape into the low-q block.

Two-mode MLE has lower pooled chi2 than the prior in all six FoXS libraries,
so its FoXS result is candidate-pool-stable. It is still not a method claim:
§10.46 showed that `saxs_md` reverses that ordering. NNLS is unstable under
both seed removal and half-pooling.

The full-data MaxEnt descriptive quantities are also pool-sensitive:

```text
selected tau:        1 or 10
compact mass range:  0.35887..0.46972
ESS/N range:         0.83936..0.99040
curve RMS relative
  to full replay:    0.00119..0.04393
curve max relative
  difference:        0.00319..0.06164
```

Removing seed `20260714` changes the full-data scaled curve by `4.39%` RMS
and `6.16%` maximum; balanced half A changes it by `3.90%` RMS. Thus even the
observable curve is only moderately pool-stable, while mode mass and
microscopic weights are plainly not identified.

Reproducible reports:

```text
slurm/work/cohere_sasdxc6_ca_balanced_highqc_ms512_v1/
  candidate_pool_sensitivity_v1/
    variant_summary.tsv
    curve_stability.tsv
    fold_diagnostics.tsv
    full_replay_check.tsv
    materialized_configs/
    results/
```

Their four aggregate SHA-256 values are `a3930c43...d35f0`,
`7397fe27...adf86`, `d0c0acd9...b64f1`, and `e30e84ca...9ea3`.

#### SASDNV6 artifact audit and q-unit amendment

The preregistered SASDNV6 artifacts were downloaded only through Slurm job
`12169173` and atomically published under:

```text
slurm/artifacts/cohere_sasdnv6_sec_saxs_v1/
```

The downloaded curve exactly matches the copy inside the full-entry ZIP. Raw
artifact SHA-256 values are:

```text
SASDNV6.dat          b6cbec8f...41279
SASDNV6.zip          2820b751...de58c
SASDNV6_3629.fasta   6dc9ad96...72471
```

The first ingestion validator reported zero positive sigma rows even though
the file visibly had a third numeric column. This was an audit-script defect:
the SASBDB curve uses CRLF, and the awk numeric check had not removed the
carriage return attached to the third field. It did not affect the original
bytes, hashes, q ordering, archive comparison, or the Python profile reader,
which strips lines before parsing. CRLF-aware Slurm audit job `12169292`
confirmed:

```text
numeric curve rows:          1,234
positive deposited sigmas:  1,234
strictly increasing q:      yes
FASTA length:               114 aa
archive curve match:        exact
```

`ingest_metadata_v2.json` has SHA-256 `c052e4ef...f187e`; the v1 metadata is
retained with an explicit supersession record. The analysis will use every
positive deposited sigma, not the 5% fallback.

Artifact inspection also falsified the provisional q-unit assumption in the
base preregistration. The deposited curve is already in `A^-1`:

```text
deposited q_min = 0.01778791 A^-1
official Dmax   = 175 A
q_min*Dmax/pi   = 0.990862
official SASBDB validation value = 1.0
```

Multiplying q by `0.1` would instead give `0.0991` and contradict the
official validation metric. Before any E1A profile calculation or fit, this
was corrected in a frozen amendment:

```text
slurm/configs/cohere_sasdnv6_second_case_preregistration_amendment_v2.json
SHA-256: 2aa4a355a44c0ec03b0ead5588557732f4f61d8d7ea04ddee97e7e598fcb8a65
registered: 2026-07-19T12:05:01-04:00
```

The only protocol change is `q_scale=1.0`. The mechanical reference rule
retains 1,058 positive rows through `q<=0.5 A^-1`, selects `stride=6`, and is
expected to leave 178 points.

The complete primary E1A FoXS/CV config was then frozen, still before any
profile calculation or inference:

```text
slurm/configs/cohere_sasdnv6_maxent_blocked_cv_foxs_v1.json
SHA-256: 8e12ba9140778e3dc169afb52f7ca8641d62fbe10998a38b94a9bc349c68a238
registered: 2026-07-19T12:07:57-04:00
```

Generator support for `--skip-target-rg-summary` prevents the curve-blind
library build from even computing a legacy target-Rg reweighting diagnostic.
Slurm job `12169427` passed 35 focused Python tests, including this policy
and the pool/curve additions, plus all 55 targeted Java tests. Candidate
array job `12169481` completed the three registered 512-frame libraries at
seeds `20260720-20260722`. Large outputs physically reside under
`/usr/xtmp/lz280/cohere-idp/`; the project work path is a symlink. The
registered CA-QC gates and the subsequent frozen primary analysis are
reported in §10.48.

### 10.48 Independent E1A primary result and preregistered decision (2026-07-19)

#### Curve-blind candidate-library gate

Candidate array job `12169481` generated exactly 512 structures per seed,
alternating 256 coil and 256 compact attempts. The registered CA-QC gate was
evaluated before any FoXS inference:

```text
seed       coil passed / 256   compact passed / 256   gate
20260720   256                 254                    pass
20260721   256                 253                    pass
20260722   256                 250                    pass
```

All six seed-mode cells exceed the frozen `231/256` minimum. The inference
library therefore contains 1,525 candidates: 768 coil and 757 compact.
Concurrent array tasks briefly exposed a symlink-creation race in the
project work path. The accidental self-referential link was removed before
scoring, the creation command was hardened with `ln -sT`, and inspection
confirmed that no candidate bytes or gate counts were affected.

#### FoXS execution-only amendment

The first profile array, Slurm job `12170956`, stopped after only 16-51
profiles per seed. FoXS itself had returned profiles, but Python
`TemporaryDirectory` cleanup on shared `/usr/xtmp` encountered delayed
directory metadata and raised `Directory not empty`. No baseline, MaxEnt, or
other inference ran, and no partial profile values were used to alter the
scientific protocol. The incomplete `foxs_fixed_profiles_v1` directories and
logs remain failed-execution provenance.

The FoXS wrapper was changed only to create isolated work directories on
node-local scratch and to atomically publish completed profiles. Slurm job
`12171034` then passed 35 focused Python tests and all 55 targeted Java tests.
The validated wrapper SHA-256 is:

```text
c28e4c521ac4aabf398f852c2d255d78b3e0c6160d00d83aa9ebf44c8328a3cf
```

Before inference, the execution-only amendment was frozen as:

```text
slurm/configs/cohere_sasdnv6_maxent_blocked_cv_foxs_v2.json
SHA-256: a485742f7a0b902ccbee5710bd95939fc06e20705716186f5a56a634b43f39f6
registered: 2026-07-19T13:20:49-04:00
```

Protein, candidates, QC filter, prior, FoXS version and options, q grid,
uncertainties, folds, tau grid, objective, optimizer tolerances, baselines,
decision rule, and claim limits are unchanged from v1. Clean profile array
job `12171068` completed all 512 profiles for every seed. The per-seed
value-file SHA-256 values are `637e0d...80c64`, `3c66df...f4c2e`, and
`4692df...9b36`; the corresponding sorted profile-manifest SHA-256 values
are `8c7a77...3712`, `3bac6f...24e6`, and `da1c23...e5b2`.

#### Frozen FoXS blocked nested-CV result

Slurm job `12171084` completed the preregistered primary comparison on 178 q
points in four contiguous blocks of sizes 45, 45, 44, and 44. Pooled
held-out standardized residual chi2 is:

```text
method             pooled held-out chi2   mean outer ESS/N   mean compact mass
balanced prior     7.2053152615            0.999948           0.500000
two-mode MLE       5.2539959060            0.575943           0.060226
sparse NNLS        6.3337792607            0.004246           0.081712
KL-MaxEnt/MAP      7.2870985147            0.812550           0.376675
```

The mode-mass and ESS columns are descriptive diagnostics, not physical
population estimates. Fold-level scores and the nested-CV-selected MaxEnt
state are:

```text
fold   q range (A^-1)       prior       two-mode    NNLS        MaxEnt      tau    ESS/N
0      0.01779..0.13817     24.32413    16.96465    18.76595    24.32298    1000   0.999948
1      0.14091..0.26117      2.53694     1.92073     3.93691     2.54053    1000   0.999967
2      0.26390..0.38119      0.55833     0.61559     0.78010     0.55854    1000   0.999969
3      0.38392..0.49998      1.11890     1.32463     1.62409     1.44704       0.1  0.250317
```

All four MaxEnt outer fits satisfied the frozen block-coordinate convergence
diagnostics, in 2, 3, 3, and 67 iterations. Folds 0-2 select `tau=1000` and
remain essentially at the prior. Training on the three lower-q blocks instead
selects `tau=0.1`, strongly favors the coil proposal, and predicts the held-out
highest-q block worse than the prior. Thus the negative endpoint is not an
optimizer non-convergence result.

Unrestricted MaxEnt fails the preregistered predictive gate:

1. its pooled score is worse than both the prior and two-mode MLE;
2. it is lower than the prior in only one fold, and that fold-0 difference is
   negligible, rather than in at least three folds;
3. although the ESS, convergence, and no-`>2x`-prior safeguards pass, those
   safeguards cannot compensate for failure of the predictive criteria.

The independent E1A result therefore confirms the no-go for unrestricted
per-conformer MaxEnt under the current candidate generators and SAXS
forward-model protocol. It does not show that observable conditioning is
useless, and it does not test sampling efficiency.

Two-mode MLE is the best pooled FoXS predictor on E1A, independently repeating
its FoXS lead on SASDXC6. However, it loses to the prior in folds 2 and 3, and
§10.46 already showed that the SASDXC6 two-mode ordering reverses under
`saxs_md`. It is consequently a low-complexity mode-level hypothesis, not yet
a backend-stable method claim. Any hierarchical successor designed after
seeing this E1A result requires a third frozen biological case; E1A cannot be
relabelled as its validation.

The full-data MaxEnt fit is retained only as a descriptive observable
projection:

```text
selected tau:       0.1
training chi2:      1.2929094126
ESS/N:              0.3176799363
maximum weight:     0.0063110923
compact mass:       0.0167202051
iterations:         59
```

These quantities are not held-out evidence and do not identify a physical
compact population. Reproducible output:

```text
slurm/work/cohere_sasdnv6_ca_balanced_highqc_ms512_v1/
  maxent_blocked_cv_foxs_v2/
    result.json
    outer_fold_results.tsv
    summary.tsv
    final_maxent_curve.tsv
```

Their SHA-256 values are, respectively,
`580576556cfdf60ce8cc33e0cdb88db8295f48ff3ab497d5d4238dcbf7b5b5d1`,
`7fa15a98d1a9c7f32e485d2fd99c95281e324d83776f5a0cd22e001711248880`,
`6844478865fb641afe3c9cc3d34433f4ce41377dfada745a362309c1987ce62c`,
and `9fa744badd04363ad6d9bd614c5de0261d9630110c1336f6a948d825991aca90`.

Before any E1A pool perturbation, the mandatory full replay, three
leave-one-seed-library-out variants, and two disjoint balanced half-pools were
materialized in the following registry:

```text
slurm/configs/cohere_sasdnv6_candidate_pool_sensitivity_foxs_v1.json
SHA-256: 281153d4b0bf1212426374f12c321e8971b7a139f27daf85deab320e8dd9b775
registered: 2026-07-19T13:28:13-04:00
```

The half-pools use the same within-seed alternating-mode-preserving
`index mod 4` rule as §10.47. These candidate-pool perturbations and the
secondary `saxs_md` forward-model sensitivity remain mandatory diagnostics;
they cannot reverse or rescue the frozen primary endpoint.

After the FoXS result was frozen, but before calculating any E1A `saxs_md`
profile or fit, the preregistered secondary forward-model sensitivity was
materialized as:

```text
slurm/configs/cohere_sasdnv6_maxent_blocked_cv_saxs_md_sensitivity_v1.json
SHA-256: b35eaae96d0f8d1d6e295b2bd0cada34d2f167e6d3f7697869a235987c8d85d0
registered: 2026-07-19T13:32:00-04:00
```

It reuses the identical 1,525-candidate CA-QC filter, balanced prior,
178-point reference grid, deposited uncertainties, blocked nested-CV folds,
tau grid, objective, optimizer, and baselines. Only the forward profiles
change to the existing CA-only `saxs_md`/empty-solvent pipeline. This is a
backend-sensitivity diagnostic, not a second primary endpoint and not an
opportunity to retune the E1A analysis.

### 10.49 E1A FoXS candidate-pool stability (2026-07-19)

Slurm job `12171199` completed all six variants in the pre-run registry from
§10.48. The full replay reproduced every pooled primary score with exactly
zero difference at the registered `1e-8` tolerance:

```text
variant                  N     prior     two-mode   NNLS       MaxEnt     tau   ESS/N    compact
full replay              1525  7.20532   5.25400    6.33378    7.28710    0.1   0.31768  0.01672
leave out seed 20260720  1015  7.32528   5.07750    4.88287    7.40750    0.1   0.31483  0.01569
leave out seed 20260721  1016  7.29391   5.14883   12.64337    7.37607    0.1   0.32866  0.01662
leave out seed 20260722  1019  7.00958   5.53200    6.89281    7.09007    0.1   0.30789  0.01798
balanced half A           763  7.54335   5.45420    9.69227    7.61980    0.1   0.32104  0.01505
balanced half B           762  6.90191   5.06808    4.76674    6.98823    0.1   0.31438  0.01856
```

Unrestricted MaxEnt is worse than both the balanced prior and two-mode MLE
in all six pools. Its fold pattern is also invariant: the nominal fold-0
improvement over the prior is only `0.0042-0.0052%`, it loses folds 1-3, and
its fold-3 score is `1.27-1.31x` the matched prior. Candidate-pool
perturbation therefore neither explains nor rescues the negative independent
endpoint.

Two-mode MLE has the lowest pooled score among prior, two-mode, and MaxEnt in
all six variants. Its q-block heterogeneity is equally stable: it beats the
prior in folds 0 and 1 but loses in folds 2 and 3 every time. This supports a
pool-stable FoXS mode-level signal, not a uniformly extrapolating or
backend-stable method. Sparse NNLS remains a finite-library fitting bound:
its pooled score ranges from `4.76674` to `12.64337`.

Unlike the earlier SASDXC6 pool perturbations, the E1A full-data MaxEnt
observable projection is numerically stable:

```text
selected tau:                  0.1 in all variants
ESS/N range:                   0.30789..0.32866
descriptive compact mass:      0.01505..0.01856
curve RMS vs full replay:      0.00135..0.00173
curve maximum vs full replay:  0.00284..0.00468
```

This is stability of a training-fit observable curve that failed held-out
prediction, not validation of its microscopic weights or compact population.
Reproducible reports:

```text
slurm/work/cohere_sasdnv6_ca_balanced_highqc_ms512_v1/
  candidate_pool_sensitivity_foxs_v1/
    variant_summary.tsv
    curve_stability.tsv
    fold_diagnostics.tsv
    full_replay_check.tsv
    materialized_configs/
    results/
```

The four aggregate SHA-256 values are
`abf013568fde3521f38f7100f1cc9237e5943ddfbc08c40cb20ae0b5f8647a7b`,
`5331f41c3f6009868bf6462a594225572103f8e0980118184465d0223850f01d`,
`e8bf02261bc3ec4c3f0ec31c9e07ec155404b625c5e3ab360bf41b47189be5f7`,
and `6cf539ff4be5213684eb8e36fc69c5af71015dba9768f5152671cbc5623328d1`.

### 10.50 Third-case data embargo before successor-method development (2026-07-19)

Metadata-only screening selected
[SASDU37, the ANAC013 residues 161-274 IDR](https://www.sasbdb.org/data/SASDU37/),
as the next untouched biological case. The official entry reports a
114-residue monomer, 13 kDa expected versus 14 kDa experimental molecular
weight, inline SEC-SAXS collection, downloadable curve/archive/FASTA, and no
model related to the curve. This gives a matched-length case from a different
organism and protein family without importing a deposited EOM ensemble into
the prior.

The main alternatives inspected at catalog level were the TRPV4 IDR
`SASDQM8` and FATZ-1 `SASDJJ6`, both of which have multiple deposited EOM
models, plus two urea-condition ruler constructs with larger
experimental-versus-expected molecular-weight discrepancies. No numeric
SASDU37 curve row, FASTA residue, archive member, conformer, profile, or fit
was inspected before registration.

```text
slurm/configs/cohere_sasdu37_third_case_preregistration_v1.json
SHA-256: 7aae7af87b91a61175e9ab6495375a5d22c2661a76ce57997e43570ef86501aa
registered: 2026-07-19T13:47:31-04:00
```

This is deliberately stronger than the E1A timing rule. The FASTA alone may
be downloaded through Slurm for curve-blind generation of three new
`3 x 512` libraries at seeds `20260723-20260725`, using the unchanged
per-seed/mode QC gate. The numeric curve and full archive remain embargoed
until all of the following are frozen and hash-registered:

1. the exact successor objective and population interpretation;
2. hyperparameter grids, tie breaks, optimizer and convergence diagnostics;
3. blocked nested-CV endpoint, baselines and predictive/stability gates;
4. matched-model and off-library synthetic calibration;
5. exact inference code and primary config.

If no successor method survives development calibration, SASDU37 remains a
third preregistered no-go/benchmark case. Either way, it is one manually
screened external case, not a representative cohort.

FASTA-only ingestion then ran through Slurm job `12171453`. The compute job
resolved the FASTA link from a retained official entry-page snapshot, verified
the expected 114-residue alphabet and length, and downloaded neither the
numeric curve nor the full archive:

```text
artifact root: slurm/artifacts/cohere_sasdu37_fasta_v1/
FASTA URL:     https://www.sasbdb.org/molecule/SASDU37/5176.fasta
FASTA SHA-256: c31b098c935106b6fb8f7baae203a3681e55822db71364cc9d25afe0c0e81682
metadata hash: 5ec84561ec0427d2e28178a66c0a1a113d2dfb77809b2b6154ea0aaf59592a9f
curve fetched: no
archive fetched: no
```

The artifact bytes reside under `/usr/xtmp/lz280/cohere-idp/artifacts/`; the
project artifact path is a compatibility symlink.

### 10.51 E1A secondary forward-model result and method decision (2026-07-19)

Slurm array `12171364` generated exactly 512 CA-only `saxs_md` profiles for
each frozen E1A seed. Each task checked the secondary config, frozen FoXS
primary result, SASDNV6 curve, `saxs_md` executable, profile wrapper, and
CA-input helper hashes before execution. The three profile-manifest SHA-256
values are `b923f8...2d50`, `4346e5...2428`, and `188309...0964`; their
completion-record SHA-256 values are `d3c73f...5b73`,
`a861ce...5b00`, and `bacb10...180`.

Dependent Slurm job `12171368` then completed the frozen secondary blocked
nested-CV comparison on the identical 1,525 candidates and 178 reference
points:

```text
method             pooled held-out chi2   mean outer ESS/N   mean compact mass
balanced prior     5.1372741708            0.999948           0.500000
two-mode MLE       2.4596654578            0.533751           0.027108
sparse NNLS       14.2524777062            0.004470           0.029241
KL-MaxEnt/MAP      3.9733549352            0.600661           0.237178
```

Fold-level results are:

```text
fold   q range (A^-1)       prior      two-mode   NNLS       MaxEnt     tau   ESS/N
0      0.01779..0.13817     16.74284    6.49225   50.62411   11.92207    0.1  0.810972
1      0.14091..0.26117      1.96654    1.57455    3.70269    1.99350  100    0.999181
2      0.26390..0.38119      0.59050    0.55565    0.61582    0.57648    0.1  0.305679
3      0.38392..0.49998      1.05752    1.14467    1.48043    1.26571    0.1  0.286813
```

All MaxEnt outer fits converged under the frozen diagnostics in 46, 6, 84,
and 84 iterations. MaxEnt improves the pooled score relative to the prior in
this secondary backend, but it still fails the preregistered advancement
logic: it is worse than two-mode MLE and beats the prior in only two of four
folds. Because this is a secondary sensitivity, it could not have rescued
the failed FoXS primary endpoint in any event.

Two-mode MLE is the best pooled predictor and beats the prior in three of four
`saxs_md` folds. Thus the E1A mode-level pooled lead survives both FoXS and
`saxs_md`, although its q-local pattern is not identical: it wins only folds
0-1 under FoXS and folds 0-2 under `saxs_md`. SASDXC6 still reverses the
two-mode/prior ordering under `saxs_md`, so a universal backend-stability
claim remains unsupported. The correct conclusion is a promising
low-complexity successor hypothesis, not a validated physical compact
population.

Sparse NNLS again behaves as an unstable finite-pool bound: its lowest-q
held-out score is `50.62`, producing the worst pooled result. The full-data
MaxEnt fit remains descriptive only:

```text
selected tau:       0.1
training chi2:      1.2090482534
ESS/N:              0.3469975486
maximum weight:     0.0071073360
compact mass:       0.0306088326
iterations:         73
```

Reproducible outputs:

```text
slurm/work/cohere_sasdnv6_ca_balanced_highqc_ms512_v1/
  maxent_blocked_cv_saxs_md_sensitivity_v1/
    result.json
    outer_fold_results.tsv
    summary.tsv
    final_maxent_curve.tsv
```

Their SHA-256 values are
`43445fafe94dba016b58c8f8ee165b547ab8f57fe7d4fef753da36119b3abf78`,
`4f7b2532c5bf37712178eab10f9e083b373e92ae24c7779a3eaabaf65baa0cf4`,
`4bdfe5e937a12f767a2e872190c858b0d60fbbec910f25c560c9dda3be35856f`,
and `8da3767b9e433daf954298bcaf567549cce423feaa2850201a2040f5c2c25159`.

#### Decision after the complete second-case evidence

Unrestricted per-conformer MaxEnt is a no-go as the paper method under the
current protocol. It fails the SASDXC6 FoXS primary, SASDXC6 `saxs_md`
sensitivity, E1A FoXS primary, all ten non-replay pool perturbations across
the two proteins, and the E1A advancement rule despite improving its
secondary pooled prior score.

The next development route is therefore complexity-controlled and explicitly
mode-aware, with three safeguards:

1. demonstrate matched-model positive controls and off-library failure
   behavior before another experimental endpoint;
2. develop and choose all successor details using only SASDXC6 and SASDNV6;
3. keep the registered SASDU37 numeric curve embargoed until the method,
   hyperparameters, code hashes and decision gate are frozen.

This route tests whether the repeatable low-dimensional signal can be retained
without assigning unrestricted microscopic weights. It does not presume that
the current coil/compact proposal labels are uniquely physical.

### 10.52 Session closeout and successor-method handoff (2026-07-19)

> Same-day continuation note: this stopping point was superseded by §10.53.
> In particular, the curve-blind SASDU37 library was subsequently submitted
> and completed, and the synthetic-calibration registry plus an unvalidated
> runner draft were added. Preserve this section as the historical
> pre-submission record, but use §10.53 for the current handoff.

#### Scientific stopping point

The current evidence closes two routes without closing the COHERE-IDP
research program:

1. the locked steric-surrogate WMB proposal does not have a stable sampling
   efficiency advantage;
2. unrestricted per-conformer MaxEnt is not the paper method.

The second conclusion now rests on two biological cases, two forward-model
implementations where available, blocked nested-CV, ten non-replay
candidate-pool perturbations, explicit scale/weight convergence diagnostics,
and full-data-versus-held-out separation. It is stronger than a single
negative run, but it is not a claim that all maximum-entropy ideas or all
observable conditioning fail.

The repeatable positive lead is narrower: two-mode MLE is the best pooled
FoXS prior/regularized comparison for both proteins and is also best under
E1A `saxs_md`. SASDXC6 `saxs_md` reverses the two-mode/prior ordering, and
q-block wins are heterogeneous. Therefore the evidence supports developing a
complexity-controlled mode-aware successor; it does not yet support calling
two-mode MLE a backend-stable paper method or interpreting its fitted mass as
a physical population.

#### Final completed Slurm actions

```text
12171199   E1A FoXS candidate-pool sensitivity       completed
12171364   E1A CA-only saxs_md profiles, 3 x 512     completed
12171368   E1A saxs_md blocked nested-CV             completed
12171453   SASDU37 FASTA-only ingestion              completed
12171465   focused Python + targeted Java validation completed
```

The final queue inspection found no active or pending COHERE-IDP job. Other
userspace jobs in the queue are unrelated and were left untouched.

Job `12171465` passed 36 focused Python tests and all 55 targeted Java tests.
It validates the policy that curve-blind candidate generation with
`--run-saxs` disabled no longer requires an experimental reference-profile
file. Current relevant code hashes are:

```text
candidate generator:
  a0b54ed78c286a545dfafb38b8ea613ef67c938b65ec567b2919cab8d3a1ea80
candidate-generator test:
  d978fe16bd82ddbfd57f5e816756e14aa74475acea797b0fb261bdbec1f93ae6
MaxEnt/CV wrapper:
  1364b69fe1146edd0fe5ca5f1af92530c17a611dfdb76511371cc2cd0c3fd236
```

#### SASDU37 embargo and prepared next action

The third case remains protected by the registry
`cohere_sasdu37_third_case_preregistration_v1.json` with SHA-256
`7aae7af8...01aa`. Only the official 114-residue FASTA and catalog-page
snapshot exist locally. FASTA SHA-256 is `c31b098c...1682`; the numeric curve
and full archive have not been downloaded.

The following curve-blind candidate-library script is prepared and passed
static `bash -n`, but was deliberately **not submitted** at session close:

```text
slurm/scripts/run_cohere_sasdu37_candidate_library.slurm
SHA-256: 76197a783cd98c92b958ccb29bf3a59b8a76d5caf57d6217580110a9672aab35
```

It locks seeds `20260723-20260725`, 512 candidates per seed, alternating
coil/compact modes, the unchanged `>=231/256` seed-mode QC gate, the
Slurm-validated generator hash, and an explicit check that the embargoed
curve path does not exist. Submitting this job is safe before method
development because it reads only the FASTA and computes no target-Rg or SAXS
score.

#### Next-session order

1. Re-read §§10.45-10.52 and verify that the SASDU37 numeric curve/archive
   remain absent.
2. If continuing candidate preparation, submit
   `run_cohere_sasdu37_candidate_library.slurm` only through `sbatch`; report
   the six seed-mode QC cells without accessing SAXS data.
3. Before implementing a successor, freeze a machine-readable synthetic
   calibration registry covering matched-library positive controls,
   off-library targets, correlated-noise seeds, mode-only signals, and
   within-mode signals. Run all calibration and tests only through Slurm.
4. Develop and compare the complexity-controlled mode-aware successor using
   only SASDXC6 and SASDNV6. Keep prior, two-mode MLE, NNLS and unrestricted
   MaxEnt as fixed baselines, and retain blocked nested-CV and convergence
   diagnostics.
5. Freeze the exact successor objective, hyperparameter grid, code/config
   hashes and primary decision gate. Only then may a Slurm job download and
   audit the SASDU37 numeric curve/archive.
6. Use the untouched SASDU37 endpoint to decide the paper claim. If the
   successor fails, retain the result in the transparent no-go/benchmark
   framing rather than tuning against it.

Hard constraint for every continuation: the login node is for lightweight
inspection and editing only. Python tests, Gradle, candidate generation,
FoXS, `saxs_md`, optimization, calibration and experiments must run through
Slurm.

Copy-paste continuation prompt:

```text
请恢复并继续 COHERE-IDP 项目。

工作目录：
/home/users/lz280/IdeaProjects/OSPREY3

首先完整读取：
docs/PACKStar_IDP_Ensemble_Extension_Plan.md
尤其是 §10.45-§10.52，并以 §10.52 的 session handoff 为接力点。

当前科学判定：
unrestricted per-conformer MaxEnt 已被两蛋白、blocked nested-CV、
两种前向模型和候选池扰动判为 no-go；two-mode 只支持开发
complexity-controlled mode-aware successor，尚不是论文 claim。

下一步：
1. 严格保持 SASDU37 数值曲线和完整 archive 封存；
2. 可先仅通过 sbatch 提交已准备的
   slurm/scripts/run_cohere_sasdu37_candidate_library.slurm，
   只检查曲线盲候选库的六个 seed-mode QC gate；
3. 在实现 successor 前，先冻结 matched-model/off-library/correlated-
   noise synthetic calibration registry；
4. 所有测试和 calibration 仅用 Slurm；
5. 只用 SASDXC6 与 SASDNV6 开发并锁定 successor；
6. 方法、超参数、代码哈希和决策门槛完全冻结后，才允许通过 Slurm
   下载 SASDU37 数值曲线并作 untouched external validation。

Goal：
持续推进 COHERE-IDP，直到方法、基线、消融、稳定性、可复现产物和
科学 claim 达到可以开始撰写 research paper 的成熟度。严格区分
observable-conditioned inference、sampling efficiency 和 ensemble
identifiability；只作证据支持的 claim。

硬性约束：
绝不在 login 节点运行 Python 测试、Gradle、候选生成、FoXS、
saxs_md、优化、calibration 或实验；所有计算任务必须通过 sbatch
提交到计算节点。
```

### 10.53 Curve-blind third-case library and current session handoff (2026-07-19)

> Same-day continuation note: the runner/testing status and next-session
> order in this section were superseded by §10.54. The scientific decision
> and SASDU37 embargo remain unchanged.

This section supersedes §10.52 as the exact stopping point.

#### Completed curve-blind SASDU37 candidate preparation

After the §10.52 closeout was written, the already registered curve-blind
candidate-library script was submitted as Slurm array `12171467`. All three
tasks completed with exit code `0:0` on a compute node. The official FASTA was
the only SASDU37 scientific input. Each provenance record independently says:

```text
curve_blind:                    true
numeric_curve_available_to_job: false
target_rg_summary_skipped:      true
saxs_scoring_run:               false
```

All six preregistered seed-mode CA-QC cells exceeded the unchanged
`>=231/256` gate:

```text
seed       coil passed   compact passed   gate
20260723   256/256       250/256          pass
20260724   256/256       254/256          pass
20260725   256/256       253/256          pass
```

The output root is:

```text
slurm/work/cohere_sasdu37_ca_balanced_highqc_ms512_v1/
```

Per-seed reproducibility hashes are:

```text
seed       metadata.csv SHA-256
20260723   9e7b8ba967c0f5f38ddde9593e0f234289c289a5c0c7e4a66e1847f9ba62b9db
20260724   bb5df627e7f71fab99835fc38aac98f38b38c0e9d33ee8e011798b970ebbdc7c
20260725   dc08b4f6edb0f16b9e6acd51a0405723f4a2cf06769c8cb7bf22f37d1d5e0484

seed       candidate_gate.tsv SHA-256
20260723   a564cdfff02eb85e0b35a1c3304fc4d3bf4bddd4bb63886e1c8297baf79ba81d
20260724   733ee8e15fbb293b617282dde4b6a9287efc07a4ed6987d07f3654558070a556
20260725   580e5de2622c34af53dc1876c2250dcf6eeebac0b700a4647eb9f883834d1d4c

seed       generation_provenance.json SHA-256
20260723   8d52625f4b89e43d69c56d2c05dbd3007508d9483a9dd416930e8699bab706a7
20260724   12a1f877be2212c312c53692c44c6b21a51a8db275279f1976b7123dff72eddc
20260725   335f4ebab6bc09b6d12212e09dd65a629936b86b54b6fde74e6f55166d91a5a9
```

One array task lost the harmless race to create the shared logical-root
symlink and emitted `File exists` on stderr. The script then verified that
the existing symlink resolved to the registered physical root, continued,
and completed successfully; all three outputs and provenance records are
present.

The embargoed path
`slurm/artifacts/cohere_sasdu37_curve_embargoed/SASDU37.dat` remains absent.
No SASDU37 curve, archive, SAXS profile, fit, or target-Rg statistic has been
opened or generated. Candidate preparation therefore does not spend the
untouched external endpoint.

#### Synthetic calibration freeze and unfinished implementation

Before implementing a successor method or inspecting the third-case curve, a
machine-readable calibration design was registered:

```text
slurm/configs/cohere_saxs_synthetic_calibration_v1.json
SHA-256: 0d647a1593fbdfbeeff3974bc6a29b46ac111f1af7955f075f69458be72da0a9
registered: 2026-07-19T14:08:00-04:00
```

It freezes 280 scenarios over the two already opened development cases,
FoXS/`saxs_md`, matched full-pool, leave-one-seed off-library and
cross-backend settings, five known truth families, and noise-free, iid and
strongly q-correlated noise. The frozen comparators are balanced prior,
two-mode MLE, sparse NNLS and unrestricted KL-MaxEnt/MAP with the corrected
mean-per-q objective and blocked nested-CV. Synthetic recovery is explicitly
an implementation/inferential-capacity calibration, not biological
validation.

A first implementation draft exists at:

```text
slurm/wrappers/cohere_saxs_synthetic_calibration.py
current draft SHA-256:
ccbf80233168da319d71f1954cb301c5af52a5d6b178a05353eddd87e2d2aac1
```

This runner is **not validated or frozen code**. It has received only
login-safe text inspection; no Python compile, import, unit test, scenario
expansion, optimizer run, or calibration has been executed. No dedicated
test file, prepare/array/aggregate Slurm scripts, scenario manifest, or
calibration result exists yet. The registry itself passes `jq` syntax and
structural checks, and the candidate-library script passes `bash -n`; these
checks do not validate the Python runner.

At final queue inspection there was no active or pending COHERE-IDP job.
Other user jobs were unrelated and were left untouched.

#### Why the next step is calibration rather than more target-driven sampling

The project is still directed toward a paper, but the current bottleneck is
not raw candidate count. Both opened cases already have multi-seed,
mode-balanced libraries, pool perturbations and two forward-model analyses.
The locked WMB proposal did not improve sampling efficiency stably, while
unrestricted per-conformer MaxEnt did not improve held-out prediction
reliably. Repeating either proposal without a new hypothesis would add
correlated candidates without resolving support mismatch, model complexity,
or ensemble non-identifiability.

Sampling has therefore not been abandoned. The untouched SASDU37 candidate
library was just generated without its curve, and the synthetic registry
contains explicit off-library tests that ask when a fixed candidate pool is
insufficient. A new round of sampling becomes justified if those tests
identify a concrete missing-support failure or if a new proposal has a
predeclared efficiency gate. It should not be triggered merely because an
inference method lost.

#### Repository state at handoff

The OSPREY3 worktree is intentionally dirty and contains pre-existing and
current COHERE-IDP work. In particular, the main plan is modified and the
SASDU37 preregistration, SASDU37 candidate script, synthetic registry and
synthetic runner currently appear as untracked files. They are research
artifacts, not disposable scratch. Preserve unrelated user changes and do
not use `git clean`, `git reset --hard` or bulk checkout to make the tree
look clean. The legacy `/home/users/lz280/COHERE-IDP/README.md` was updated
but that standalone directory is not a Git worktree.

#### Exact next-session order

1. Read §§10.51-10.53 and verify again that the SASDU37 curve/archive remain
   absent. Do not resubmit `12171467` or overwrite its candidate library.
2. Review the full synthetic runner against the registered design. Add
   focused tests for exactly 280 unique scenarios, truth-weight construction,
   deterministic noise, candidate identity/q-grid rules and tiny-library
   behavior.
3. Run those Python tests only through
   `slurm/scripts/run_cohere_unit_tests.slurm`. Fix failures and repeat only
   through Slurm.
4. After tests pass, freeze the runner hash in a separate execution record
   and add hash-checking prepare, bounded-concurrency array and aggregate
   Slurm scripts. Run a compute-node smoke check before the 280-scenario
   calibration.
5. Use the calibration to define falsifiable gates for a
   complexity-controlled, mode-aware successor. Develop it only with
   SASDXC6/SASDNV6 and retain all four frozen baselines, blocked nested-CV,
   forward-model sensitivity and candidate-pool diagnostics.
6. Freeze the successor objective, interpretation, hyperparameters, code and
   config hashes, convergence criteria and primary decision gate. Only then
   may a Slurm job download and audit the SASDU37 numeric curve/archive.
7. Treat the untouched SASDU37 result as external validation. A failure is a
   reportable no-go/benchmark outcome, not permission to tune against it.

Hard constraint: the login node is for lightweight inspection and editing
only. Python/Java tests, candidate generation, FoXS, `saxs_md`, optimization,
synthetic calibration and all experiments must use `sbatch`.

Copy-paste continuation prompt:

```text
请恢复并继续 COHERE-IDP 项目。

工作目录：
/home/users/lz280/IdeaProjects/OSPREY3

首先完整读取：
docs/PACKStar_IDP_Ensemble_Extension_Plan.md
尤其是 §10.51-§10.53，并以 §10.53 为唯一当前 handoff；§10.52 已被
同日后续工作 supersede。

论文方向与当前科学判定：
项目仍以达到可撰写 research paper 的证据成熟度为目标。锁定 WMB
proposal 不具备稳定 sampling-efficiency 优势；unrestricted
per-conformer KL-MaxEnt/MAP 在两蛋白、blocked nested-CV、候选池
扰动和两种前向模型证据下是当前 protocol 的 no-go。two-mode 的
重复 pooled lead 只支持开发 complexity-controlled mode-aware
successor，尚不支持 backend-stable 或 physical-population claim。

已完成：
- SASDU37 曲线盲候选库 Slurm array 12171467 已完成，六个
  seed-mode CA-QC gate 全部通过；
- SASDU37 数值曲线与完整 archive 仍未下载，embargo 必须保持；
- synthetic calibration registry 已冻结为
  slurm/configs/cohere_saxs_synthetic_calibration_v1.json，
  SHA-256 为
  0d647a1593fbdfbeeff3974bc6a29b46ac111f1af7955f075f69458be72da0a9；
- slurm/wrappers/cohere_saxs_synthetic_calibration.py 只是未验证草稿，
  没有 Python/Slurm 测试结果，也没有 calibration 输出。

下一步：
1. 保持 SASDU37 curve/archive 封存，不要重复提交或覆盖候选库；
2. 审查 synthetic runner，并补 exactly-280-scenarios、truth weights、
   deterministic noise、identity/q-grid 和 tiny-library tests；
3. 仅通过 sbatch 跑测试，绝不在 login 节点运行 Python；
4. 测试通过后冻结 runner/execution hashes，准备 hash-checking
   prepare/array/aggregate Slurm scripts，先 smoke 再跑 280 scenarios；
5. 只用 SASDXC6/SASDNV6 calibration 与真实 endpoint 开发并冻结
   mode-aware successor；
6. 方法、超参数、代码/config hashes、收敛和 decision gate 全部冻结
   后，才允许通过 Slurm 下载 SASDU37 数值曲线做 untouched external
   validation。

Goal：
持续推进 COHERE-IDP，直到方法、基线、消融、稳定性、可复现产物和
科学 claim 达到可以开始撰写 research paper 的成熟度。严格区分
observable-conditioned inference、sampling efficiency 和 ensemble
identifiability；只作证据支持的 claim。

环境：
cohere-saxs 的物理环境位于
/usr/xtmp/lz280/conda_envs/cohere-saxs；
/home/users/lz280/miniconda3/envs/cohere-saxs 仅为兼容 symlink。

硬性约束：
绝不在 login 节点运行 Python 测试、Gradle、候选生成、FoXS、
saxs_md、优化、calibration 或实验；所有计算任务必须通过 sbatch
提交到计算节点。
```

### 10.54 Synthetic-calibration core validation and session handoff (2026-07-19)

This section supersedes §10.53 as the exact stopping point. The session
advanced only the frozen synthetic-calibration implementation and its tests;
it did not start a calibration scenario or access the SASDU37 endpoint.

#### Integrity fixes made before execution

Review of the first runner draft found several implementation-integrity gaps
that did not alter the registered 280-scenario scientific design:

1. scenario expansion previously used every globally registered noise model
   rather than the noise-model list declared by each scenario submatrix;
2. only the total scenario count was checked, not the frozen
   `140/84/56` family counts;
3. AR(1) generation used a default `rho` instead of reading the registered
   value;
4. truth scale, finite/positive curve values, mode labels, compact mass, Rg
   length/finite values, attainable per-mode ESS and candidate-name
   uniqueness needed explicit validation;
5. sparse-eight selection now implements candidates nearest the four
   registered within-mode Rg quantiles with deterministic name tie-breaking;
6. aggregate manifest labels now use the same stable comma-separated
   serialization as the scenario manifest.

The frozen scientific registry was not changed:

```text
slurm/configs/cohere_saxs_synthetic_calibration_v1.json
SHA-256: 0d647a1593fbdfbeeff3974bc6a29b46ac111f1af7955f075f69458be72da0a9
```

Eight focused tests were added in
`slurm/wrappers/tests/test_cohere_saxs_synthetic_calibration.py`. They lock:

- exactly 280 unique deterministic scenario IDs and the expected
  `140/84/56` family, truth-family and noise counts;
- submatrix-specific noise expansion and family-count drift rejection;
- balanced/mode-shift masses and per-mode Rg-tilt `ESS/N=0.30`;
- deterministic sparse-eight support and weights;
- invalid truth-weight input rejection;
- deterministic iid/AR(1) noise and strong lag-one correlation;
- matched, off-library and cross-backend identity/q-grid rules;
- stable manifest label rendering.

#### Compute-node validation

All validation ran through Slurm job `12171668`; nothing was executed on the
login node. The job completed on `fennario-01` with exit code `0:0` in
34 seconds:

```text
focused Python tests: 44 passed (8 new synthetic-calibration tests)
targeted Java tests:  55 passed
Gradle result:        BUILD SUCCESSFUL
```

Hashes of the exact validated inputs are:

```text
MaxEnt/CV wrapper:
1364b69fe1146edd0fe5ca5f1af92530c17a611dfdb76511371cc2cd0c3fd236

synthetic-calibration runner:
700e45dd0621859cdd251359ebb297afc757790a8405625fe910883a38eb6f45

synthetic-calibration tests:
50c7d0c6709a51e1963519d1472d758f06a3ac88c4302a71c6e911e59427bb61

Slurm validation script:
fe6c2060cc809ec03e74938e585529fc439d0c7abd24f6aa6fcacc8a10d4ff8d
```

Validation log SHA-256 values are:

```text
slurm/logs/cohere_validate_12171668.out
d471d8140026aa4548841784a86ce27ec2ae3f47b6a05da6452441b74356981e

slurm/logs/cohere_validate_12171668.err
c06687b3f26f28d5d7b0619d9a2e400bc9eeace96ea6035a43f1699e29023433
```

This is core unit-validation evidence, not end-to-end calibration evidence.
The tests import the runner and exercise scenario construction, synthetic
truth/noise logic and integrity rules on tiny libraries. They do not load all
four real profile libraries, run blocked nested-CV, exercise prepare/atomic
result/aggregate paths, or demonstrate optimizer convergence on any of the
280 scenarios.

#### Protected state at closeout

The SASDU37 numeric curve/archive remain absent. The curve-blind candidate
library from `12171467` is unchanged. No
`slurm/work/cohere_saxs_synthetic_calibration_v1` output root, scenario
manifest, optimizer result or aggregate exists. No active or pending
COHERE-IDP Slurm job remained at queue inspection; other user jobs were left
untouched.

The runner and its test are still untracked research artifacts in the dirty
OSPREY3 worktree. Do not clean or reset them. Hash `700e45...6f45` identifies
the tested implementation candidate; it is not yet the final execution hash,
because an end-to-end smoke may expose a required correction. Any correction
must produce a new hash and be revalidated through Slurm rather than silently
reusing this one.

#### Exact next-session order

1. Read §§10.51-10.54, treat §10.54 as the only current handoff, and verify
   again that the SASDU37 curve/archive and calibration output root are
   absent.
2. Create a separate machine-readable execution/smoke registration that
   binds the scientific registry hash, runner hash, MaxEnt wrapper hash,
   test hash, real-library config hashes, output isolation policy and
   expected smoke assertions.
3. Add hash-checking Slurm scripts for an isolated prepare/integration smoke.
   The smoke must load all four registered libraries and verify q grids,
   candidate identities, Rg metadata, 280-row manifest serialization and one
   small end-to-end scenario without using SASDU37.
4. If the smoke changes code, version the runner hash and rerun the 44
   Python plus 55 Java validation through Slurm. Do not overwrite or
   retrospectively edit the registered scientific design.
5. Only after a clean smoke, freeze bounded-concurrency scenario-array and
   aggregate scripts and launch the 280-scenario calibration. Interpret its
   gates before implementing the mode-aware successor.
6. Keep SASDU37 sealed until calibration, successor objective,
   hyperparameters, code/config hashes, convergence rules and primary
   decision gate are frozen.

Hard constraint: the login node is for lightweight inspection and editing
only. Python/Java tests, profile/library loading, candidate generation, FoXS,
`saxs_md`, optimization, smoke runs, calibration and experiments must use
`sbatch`.

Copy-paste continuation prompt:

```text
请恢复并继续 COHERE-IDP 项目。

工作目录：
/home/users/lz280/IdeaProjects/OSPREY3

首先完整读取：
docs/PACKStar_IDP_Ensemble_Extension_Plan.md
尤其是 §10.51-§10.54，并以 §10.54 为唯一当前 handoff；§10.52 和
§10.53 的旧 prompt 均已被同日后续工作 supersede。

当前科学判定：
项目继续朝 research paper 的证据成熟度推进。锁定 WMB proposal
不具备稳定 sampling-efficiency 优势；unrestricted per-conformer
KL-MaxEnt/MAP 是当前 protocol 的 no-go；two-mode 只支持开发
complexity-controlled mode-aware successor，尚不是 backend-stable
或 physical-population claim。

已完成：
- SASDU37 曲线盲候选库 array 12171467 已完成且六个 QC gate 通过；
- SASDU37 数值曲线/archive 仍未下载；
- 280-scenario scientific registry 保持冻结，SHA-256 为
  0d647a1593fbdfbeeff3974bc6a29b46ac111f1af7955f075f69458be72da0a9；
- synthetic runner 的 scenario/noise/truth/identity 完整性问题已修正；
- Slurm job 12171668 已通过 44 个 Python 测试和 55 个 Java 测试；
- 当前经测试的 runner SHA-256 为
  700e45dd0621859cdd251359ebb297afc757790a8405625fe910883a38eb6f45；
- 尚无 prepare manifest、end-to-end smoke、optimizer/calibration
  result 或 aggregate。

下一步：
1. 保持 SASDU37 curve/archive 封存，不要覆盖候选库；
2. 先冻结独立 execution/smoke registration，绑定 registry、runner、
   MaxEnt wrapper、test 和四个 real-library config hashes；
3. 添加 hash-checking Slurm prepare/integration-smoke scripts，在隔离
   output root 上验证四个 library、q/candidate identity、Rg metadata、
   280-row manifest 和一个小型 end-to-end scenario；
4. 若 smoke 要求改代码，生成新 hash 并仅通过 sbatch 重跑 44+55
   tests；
5. smoke 通过后才冻结并启动 bounded-concurrency 280-scenario
   calibration，随后再开发 mode-aware successor；
6. successor 的方法、超参数、代码/config hashes、收敛和 decision
   gate 全冻结后，才允许通过 Slurm 下载 SASDU37 曲线作 untouched
   external validation。

Goal：
持续推进 COHERE-IDP，直到方法、基线、消融、稳定性、可复现产物和
科学 claim 达到可以开始撰写 research paper 的成熟度。严格区分
observable-conditioned inference、sampling efficiency 和 ensemble
identifiability；只作证据支持的 claim。

环境：
cohere-saxs 物理环境位于
/usr/xtmp/lz280/conda_envs/cohere-saxs；
/home/users/lz280/miniconda3/envs/cohere-saxs 仅为兼容 symlink。

硬性约束：
绝不在 login 节点运行 Python/Java 测试、library/profile loading、
候选生成、FoXS、saxs_md、优化、smoke、calibration 或实验；所有
计算任务必须通过 sbatch 提交到计算节点。
```

### 10.55 Frozen-calibration v1 failure and v2 prelaunch registration (2026-07-19)

This section supersedes §10.54 as the current handoff. Work remained directed
toward paper-level evidence and did not access SASDU37 numeric data. It
completed the integration smoke, launched the frozen v1 calibration,
classified its failure, and registered a strictly failure-preserving v2
execution before any v2 output existed.

#### Clean integration smoke and formal v1 launch

The isolated integration smoke was frozen by:

```text
slurm/configs/cohere_saxs_synthetic_calibration_smoke_execution_v1.json
SHA-256:
8501ea9d67e49782747aaa1f2958757c3f0de6ac44eeabe55bdc9c8c05fcc9b2

slurm/wrappers/cohere_saxs_synthetic_smoke.py
SHA-256:
955c14eed86d1615849baa3666b5ee74dbd1722100f7ee58d82627505c65cf5b

slurm/scripts/run_cohere_saxs_synthetic_calibration_smoke.slurm
SHA-256:
409e4a238389779675900f7187f4000cb2ac715a4db60776f8b2de2d212d855b
```

Slurm job `12171735` completed with exit code `0:0` in 21 seconds. It
loaded all four development libraries without SASDU37:

- SASDXC6: 1,529 candidates and 179 retained q points for both FoXS and
  `saxs_md`;
- SASDNV6: 1,525 candidates and 178 retained q points for both backends;
- candidate order and q grids matched across backends within each case;
- the frozen registry expanded to 280 scenarios and reproduced manifest TSV
  SHA-256
  `c917f8dbed1cf82d77614389e032d5f348da16646ca3b378f7f25f77d6d51fdb`;
- scenario zero, a noiseless matched balanced-prior truth, gave pooled
  prior held-out chi2 `1.8340835502335667e-29`, and its outer and full
  MaxEnt fits converged.

The smoke report is:

```text
slurm/work/cohere_saxs_synthetic_calibration_smoke_v1/smoke_report.json
SHA-256:
0d99492bb1121506024d1708cc7ea3d6207bbf71d4d80ee037e5586351a7b7a5
```

The formal v1 execution record was then frozen as:

```text
slurm/configs/cohere_saxs_synthetic_calibration_execution_v1.json
SHA-256:
98606e664e14012f17681d7194a4e41c79da8aa9c669d1142ceb529cb204aade
```

Prepare job `12171740` completed and published the exact 280-row manifest.
Scenario array `12171741` then produced a scientifically important execution
failure:

```text
scenario tasks COMPLETED: 79
scenario tasks FAILED:    201
published result.json:     79

197 failures: MaxEnt block-coordinate optimization did not converge
  4 failures: MaxEnt dual optimization failed
  0 failures: other causes
```

The partial 79-result-set composite SHA-256 is
`7df379d2ddfd1e13ba3af1f2b1e811502932a4d5c7369d169ec58befed550c2f`.
Aggregate job `12171762` became `DependencyNeverSatisfied` and was cancelled
rather than left pending. The entire v1 output root and logs must be retained
unchanged as the failed first formal execution.

This was not a scheduler, memory, walltime, file-integrity or SASDU37-embargo
failure. It demonstrated that the current unrestricted per-conformer MaxEnt
implementation does not satisfy its frozen convergence gate robustly at the
registered 250-iteration budget. It does **not** prove that every possible
MaxEnt formulation is invalid. Combined with the completed real-data
blocked-CV no-go, however, it is strong evidence against using this
high-dimensional formulation as the paper's primary method.

#### Why v1 could not simply be aggregated

The v1 runner allowed a `RuntimeError` from any nested or final MaxEnt fit to
abort the entire scenario. Consequently the 201 failed tasks also discarded
already well-defined balanced-prior, two-mode and NNLS rows. Increasing the
iteration budget, relaxing the gradient threshold or weakening the hard gate
after seeing these failures would be post hoc tuning and is prohibited.
Repeating v1 would also erase the distinction between the registered failed
execution and a repaired execution layer.

The permitted correction is therefore mechanical and versioned:

1. retain the exact scientific registry, tau grid, q folds, objective,
   tolerances, 250-iteration limit, scenarios, baselines and gates;
2. attempt each outer MaxEnt pipeline and the full-data pipeline under those
   unchanged settings;
3. record a MaxEnt `RuntimeError` by stage, outer fold, selected tau when
   available, exception type and exact message;
4. publish all three baseline CV summaries and full fits regardless of that
   MaxEnt outcome;
5. publish MaxEnt summaries only when all required fits for that summary are
   complete;
6. let captured MaxEnt nonconvergence fail the original hard convergence gate
   at aggregate time rather than fail the scenario task;
7. continue to abort on integrity, input and non-MaxEnt programming errors.

#### Failure-preserving v2 implementation and compute-node validation

The new execution layer and targeted tests are:

```text
slurm/wrappers/cohere_saxs_synthetic_calibration_v2.py
SHA-256:
2162ee255742e81f08174d95d1106039df0a02b3ff3da8411861a60016bf1838

slurm/wrappers/tests/test_cohere_saxs_synthetic_calibration_v2.py
SHA-256:
d1a6ade8ec28b62be74da3e76185319c97742010938b8a7020ae6919d96d626d
```

The three new tests require that:

- all four outer folds and the full fit for each baseline survive captured
  outer/full MaxEnt nonconvergence;
- MaxEnt summary/full rows are published when and only when their fits are
  complete;
- non-`RuntimeError` input/integrity failures are not hidden.

Only Slurm job `12172042` executed them. It completed on `fennario-01` with
exit code `0:0`:

```text
focused Python tests: 50 passed (3 new v2 tests)
targeted Java tests:  55 passed
Gradle result:        BUILD SUCCESSFUL

stdout SHA-256:
93c48206c6a72001522ee153d40ca5bca08c4b8d8868cadb79acea7b79d4dd70
stderr SHA-256:
5710f5c4634739a9a6f449498172a37572c27fd5e639c241517f4f7fdfed562f
```

The unchanged code/design anchors remain:

```text
scientific registry:
0d647a1593fbdfbeeff3974bc6a29b46ac111f1af7955f075f69458be72da0a9

v1 calibration runner:
700e45dd0621859cdd251359ebb297afc757790a8405625fe910883a38eb6f45

MaxEnt/CV wrapper:
1364b69fe1146edd0fe5ca5f1af92530c17a611dfdb76511371cc2cd0c3fd236

unit-validation Slurm script:
fe6c2060cc809ec03e74938e585529fc439d0c7abd24f6aa6fcacc8a10d4ff8d
```

#### v2 execution freeze before launch

The v2 execution record was written after validation and before either its
logical or physical output root existed:

```text
slurm/configs/cohere_saxs_synthetic_calibration_execution_v2.json
SHA-256:
5fcf1190558bc287595ff6afe385f93cf01f10a891fa92fce5d1c1850da579b6
registered: 2026-07-19T18:20:44-04:00
```

It explicitly states `scientific_change: none`, binds all source and smoke
hashes, records the `79/201` v1 incident, and prohibits threshold or iteration
changes. Its new output roots are:

```text
logical:
/home/users/lz280/IdeaProjects/OSPREY3/slurm/work/cohere_saxs_synthetic_calibration_v2

physical:
/usr/xtmp/lz280/cohere-idp/cohere_saxs_synthetic_calibration_v2
```

Both roots and both SASDU37 embargo paths were absent at registration. The
three launch scripts passed login-safe `bash -n` inspection and are frozen
before submission:

```text
slurm/scripts/run_cohere_saxs_synthetic_calibration_v2_prepare.slurm
SHA-256:
91091f45378e7fe80edd6ba276eab87e05e301173929af5093bc33aa120e2554

slurm/scripts/run_cohere_saxs_synthetic_calibration_v2_array.slurm
SHA-256:
90e1e97a80410c348e61cbead45ab90e68a24c77d12d115e2df4485960ecdee8

slurm/scripts/run_cohere_saxs_synthetic_calibration_v2_aggregate.slurm
SHA-256:
d5034c3d11b318b85908b1ac6476c33bba7f0480d2c443d88759dbd9e68ac892
```

Every script rechecks the execution, scientific-registry, v1/v2 runner,
MaxEnt wrapper, v2-test and smoke-report hashes and aborts if a SASDU37
embargo path exists. Every scenario result must contain three complete
baseline methods and either a complete or explicitly failed MaxEnt status.
The aggregate requires exactly 280 manifest-matched result files, writes
baseline and failure tables before gate evaluation, and is expected to exit
nonzero if the unchanged `maxent_convergence` hard gate fails. Such a nonzero
aggregate is a scientific gate result, not permission to rerun or tune.

#### Exact continuation order

1. Treat this §10.55 registration and all listed hashes as immutable. Verify
   the v2 roots and SASDU37 embargo paths are still absent.
2. Submit the v2 prepare script. Require a 280-row v2 manifest with the same
   frozen TSV hash and a runtime registry whose scientific change is `none`.
3. Only after prepare passes, submit the `0-279%20` v2 array; then submit the
   v2 aggregate with `afterok` dependency on the array.
4. Require every scenario task to exit successfully and publish its baseline
   results. Do not interpret MaxEnt success/failure counts until all 280
   results are present.
5. Retain aggregate output even if the aggregate job exits nonzero. Interpret
   every original hard/diagnostic gate, stratify MaxEnt failures by scenario
   family, truth family, noise and backend, and compare complete baseline
   held-out performance.
6. Use this evidence to decide whether to freeze a low-complexity two-mode or
   other mode-aware successor, or to stop inference-method development and
   frame the paper as a calibrated no-go/benchmark result. Do not access
   SASDU37 numeric data yet.

Hard constraint: the login node is for lightweight inspection, hashing and
editing only. Python/Java tests, profile loading, optimization, calibration,
FoXS, `saxs_md` and all experiments must use `sbatch`.

### 10.56 Complete frozen calibration and session handoff (2026-07-19)

This section supersedes §10.55 as the only current handoff. The
failure-preserving v2 execution is complete. It retained the frozen v1
scientific design, produced all 280 scenario results and closed the current
unrestricted per-conformer MaxEnt/MAP route as a candidate paper method. It
did not access SASDU37 numeric data or start successor development.

#### Slurm execution outcome

Prepare job `12172079` completed on `fennario-01` with exit code `0:0`. It
published:

```text
runtime registry SHA-256:
c617df3c0273b823e8448e0e38b19441d94da76354eadbfd01585d4d80c9a51a

scenario manifest JSON SHA-256:
37218bcc090e1462f6085d44077e82422992deaafeea26de87f6706e131302c3

scenario manifest TSV SHA-256:
c917f8dbed1cf82d77614389e032d5f348da16646ca3b378f7f25f77d6d51fdb
```

The manifest TSV is byte-identical to the v1 manifest, and the runtime
registry states `scientific_change: none`. Prepare stdout SHA-256 is
`166690d827174092227e94504c0ce23d917be4c6488a197b96b05ec896eb026b`;
stderr is empty with SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

All 280 tasks in scenario array `12177766` completed with exit code `0:0`.
Every task atomically published its three complete baseline methods and a
complete or explicitly failed MaxEnt status. All 280 array stderr logs are
empty. The deterministic composites are:

```text
280 result.json files:
6834e306261aab8175333346a63103cb43c8409f8d8ca1608d7d313b4af82a0c

280 scenario stdout logs:
39d600b6a5a4df7e05b6b3304ab932348e380394415cb3ce8897017f20961011

280 empty scenario stderr logs:
cfcfb4a15707acb471bc4374bae0a9d04850c7812ecf0cc682a0d98645d0a4ce
```

Dependent aggregate job `12177811` ran only after all array tasks passed. It
wrote the complete aggregate and then exited `1:0`, exactly as registered,
because one preregistered hard scientific gate failed. This is not an
execution or serialization failure. Aggregate stdout and stderr SHA-256
values are:

```text
stdout:
3d4f42cef54548b8da0385ac3850962aa742c0e38a94eb0da710cd48e561c829

stderr:
d000318b057dcc7995c2c53f4d7c3501bb5e7c9024dfddcfa58be3fbe4052b10
```

#### Frozen gates

The aggregate contains 280 scenarios, 921 complete CV method-summary rows,
1,118 complete full-fit rows and 422 structured MaxEnt failure records.
Every gate except MaxEnt convergence passed:

| gate | severity | passed/evaluated | result |
|---|---:|---:|---:|
| manifest count | hard | 1/1 | pass |
| exact result set | hard | 280/280 | pass |
| three-baseline completeness | hard | 280/280 | pass |
| matched balanced-prior oracle | hard | 4/4 | pass |
| matched two-mode oracle | hard | 4/4 | pass |
| matched sparse NNLS extrapolation | diagnostic | 4/4 | pass |
| MaxEnt convergence | hard | 79/280 | **fail** |
| identity and q-grid integrity | hard | 280/280 | pass |

The aggregate anchor is:

```text
slurm/work/cohere_saxs_synthetic_calibration_v2/aggregate.json
SHA-256:
241c517e66fb94c585ca58a86bcec78f6eda5613d304b8308b8dfa6cc8ccaa91
```

Its referenced artifact hashes are:

```text
aggregate_method_summary.tsv
63ed474efea6bafabd1320627056e08d4213a0305981194111525e3c739d19aa

aggregate_full_fit_metrics.tsv
f3f7f48fad7ead3380bd46641e28429e59364b63b4907f6ef0da4bbcc74c1b63

scenario_status.tsv
db84df27739f01e87e847df2eba7bab1921bd49bb2e55d350aa5881ad6b9b65d

maxent_failures.tsv
52330715ad9153d4832ba6d3790a5ba235d01f2ace5697f11ed4da7c5b9b5e2f

integrity_gates.tsv
63f4a52b3e15141434ee9f1cb79ec7b29353b700c13f10862cf21e22ad461393
```

These v2 outputs and all v1 outputs must now be treated as immutable
evidence. Do not rerun, complete, delete or overwrite either root.

#### Where MaxEnt failed

The 79 v2 scenarios marked complete are exactly the same 79 scenario IDs
that published v1 results. This is strong execution evidence that v2
preserved the old convergence behavior rather than improving it; it only
prevented the other 201 scenarios from losing their baseline results.

Failure is not confined to correlated noise, candidate-support mismatch or
cross-backend mismatch:

```text
scenario family                     complete   failed
matched full pool                         50       90
off-library leave-one-seed                15       69
cross-backend full pool                   14       42

noise model                          complete   failed
none                                       8       32
iid Gaussian                              18       62
AR(1), rho=0.95                           53      107
```

Within the cleanest noiseless matched condition itself, only `8/20`
scenarios completed and `12/20` failed. Both noiseless off-library
(`0/12`) and noiseless cross-backend (`0/8`) sets had no complete MaxEnt
scenario. Failure by truth family was:

```text
truth family                         complete   failed
balanced prior                            17       11
mode shift 0.10                           14       84
mode shift plus within-mode Rg tilt       20       78
sparse eight                              10       18
within-mode Rg tilt                       18       10
```

Of 422 failure records, 420 occurred during nested outer tau selection and
two during the final full fit; 415 were frozen block-coordinate
nonconvergence and seven were frozen dual-optimization failures. No other
error type was captured.

#### Scientific decision at this stopping point

This calibration does **not** establish that the general maximum-entropy
principle is invalid. It establishes the narrower but actionable result that
the current unrestricted per-conformer KL-MaxEnt/MAP formulation and frozen
block-coordinate implementation are not robust enough for the paper's
primary method:

- the real SASDXC6 and SASDNV6 blocked-CV endpoints already lacked a stable
  held-out advantage;
- candidate-pool and forward-model sensitivities did not rescue that real
  endpoint;
- now the synthetic calibration fails its hard convergence gate in
  `201/280` scenarios, including `12/20` noiseless matched cases.

Therefore no more threshold tuning, iteration extension or unrestricted
MaxEnt rerun belongs to this frozen route. A future low-dimensional,
mode-aware method may still use entropy or a prior internally, but it must be
registered as a new method rather than described as a repair of the failed
baseline.

The four prior oracle and four two-mode oracle checks, plus all four sparse
NNLS diagnostic checks, show that the calibration harness and the three
baseline model families can recover their own noiseless matched targets.
They do **not** yet show that two-mode is generally best, biologically
correct, backend-stable or microscopically identifiable. The 280-scenario
comparative held-out and full-fit tables have not yet received a frozen
descriptive analysis. No method-ranking claim should be extracted ad hoc
from individual rows.

This result also does not justify blind resampling. MaxEnt fails even when
the truth curve is generated from the identical candidate pool, so missing
candidate support cannot explain the main convergence failure. Sampling
should resume only if the later baseline analysis isolates an off-library
support deficit or a newly registered proposal has a falsifiable efficiency
gate.

#### Protected state and next-session order

The SASDU37 numeric curve and archive remain absent after every job-level
embargo check. Its curve-blind `3 x 512` candidate library is unchanged. At
final queue inspection no COHERE-IDP job from this execution remained active
or pending; unrelated user jobs were left untouched.

Next session:

1. Read §§10.55-10.56 and treat §10.56 as the only current handoff.
2. Register, before comparative interpretation, a deterministic descriptive
   analysis plan for the frozen aggregate. It should cover all 280 paired
   prior/two-mode/NNLS held-out results, full-fit curve and compact-mass
   errors, failure strata, and the selection-biased 79-scenario complete
   MaxEnt subset labeled diagnostic only.
3. Implement and validate that report only through Slurm. Do not modify or
   rerun either calibration root.
4. Use the paired/stratified results to decide between:
   - a newly frozen low-complexity mode-aware successor with explicit null
     and matched-signal calibration gates; or
   - a calibrated benchmark/no-go paper framing if no robust successor
     signal exists.
5. Do not tune unrestricted per-conformer MaxEnt further. Do not start new
   sampling without support-deficit evidence or a registered efficiency
   hypothesis.
6. Keep SASDU37 sealed until the successor objective, hyperparameters, code
   hashes, convergence rule and decision gate are all frozen. Only then may
   its numeric curve be downloaded through Slurm for one untouched external
   validation.

Copy-paste continuation prompt:

```text
请恢复并继续 COHERE-IDP 项目。

工作目录：
/home/users/lz280/IdeaProjects/OSPREY3

首先完整读取：
docs/PACKStar_IDP_Ensemble_Extension_Plan.md
尤其是 §10.55-§10.56，并以 §10.56 为唯一当前 handoff。

当前科学判定：
项目继续朝 research paper 的证据成熟度推进。锁定 WMB proposal
不具备稳定 sampling-efficiency 优势。当前 unrestricted
per-conformer KL-MaxEnt/MAP 路线已经同时在真实 blocked-CV 和冻结
synthetic calibration 下失败，不再作为论文主方法；这不是对所有
MaxEnt 思想的否定。two-mode 仍只是待检验的 low-complexity
successor 假设，不是 biological population 或 backend-stable claim。

刚完成：
- v2 prepare 12172079 COMPLETED 0:0；
- v2 scenario array 12177766 的 280/280 tasks 全部 COMPLETED 0:0；
- aggregate 12177811 写完全部产物后按设计因 hard gate FAILED 1:0；
- manifest/result/baseline/prior-oracle/two-mode-oracle/NNLS/identity gates
  全通过；
- MaxEnt convergence 只有 79/280，通过 hard gate 失败；201/280
  场景失败，且 noiseless matched 也有 12/20 失败；
- aggregate.json SHA-256 为
  241c517e66fb94c585ca58a86bcec78f6eda5613d304b8308b8dfa6cc8ccaa91；
- SASDU37 数值 curve/archive 仍未下载。

下一步：
1. 不修改、不重跑 v1/v2 calibration roots；
2. 在阅读 method-ranking 结果前，先冻结 aggregate descriptive
   analysis 计划；
3. 仅通过 Slurm 生成 paired/stratified prior、two-mode、NNLS 报告；
   79 个完整 MaxEnt 场景只能作为 selection-biased diagnostic；
4. 根据结果决定是否冻结新的 low-complexity mode-aware successor，
   或将论文收束为 calibrated benchmark/no-go；
5. 不再给 unrestricted MaxEnt 放宽阈值或增加迭代，也不在没有
   support-deficit/efficiency 假设时盲目重采样；
6. successor 方法、超参数、code/config hashes、convergence 和
   decision gate 全冻结后，才允许通过 Slurm 下载 SASDU37 curve
   做一次 untouched external validation。

Goal：
持续推进 COHERE-IDP，直到方法、基线、消融、稳定性、可复现产物和
科学 claim 达到可以开始撰写 research paper 的成熟度。严格区分
observable-conditioned inference、sampling efficiency 和 ensemble
identifiability；只作证据支持的 claim。

环境：
cohere-saxs 物理环境位于
/usr/xtmp/lz280/conda_envs/cohere-saxs；
/home/users/lz280/miniconda3/envs/cohere-saxs 仅为兼容 symlink。

硬性约束：
绝不在 login 节点运行 Python/Java 测试、library/profile loading、
FoXS、saxs_md、优化、calibration、report computation 或实验；所有
计算任务必须通过 sbatch。
```

### 10.57 Frozen baseline-analysis plan and restored certificate branch (2026-07-21)

Before parsing or interpreting any numeric method value from the frozen v2
aggregate tables, the complete baseline-analysis protocol was registered at:

```text
slurm/configs/cohere_saxs_synthetic_baseline_analysis_plan_v1.json
SHA-256:
1cbded543008427ff540a0dc17925694f1e9621b4e32d55a06b5466a0039f175
```

A source-schema audit then identified two intentional structural blanks before
any aggregate metric was parsed: baseline `selected_tau` is not applicable,
and weight total variation exists only for identical candidate identities.
The original plan remains byte-frozen; the schema-only amendment is:

```text
slurm/configs/cohere_saxs_synthetic_baseline_analysis_amendment_v1.json
SHA-256:
1458d5ed509b34cee04ee1cb0fa7f9db53f8b220a438c7d7863b09d97064a132
```

The amendment changes no endpoint, comparison, stratum, threshold or decision
rule. It requires explicit structural-NA counts and rejects every other
missing or nonfinite baseline value.

Fixture-only Slurm validation job `12189609` passed all `7/7` focused tests on
`fennario-01`; it parsed no frozen aggregate metric. The resulting one-shot
execution lock is:

```text
slurm/configs/cohere_saxs_synthetic_baseline_analysis_execution_v1.json
SHA-256:
895ef8928d2327e28ad23facd59822c3112e4f1902bf8bb70f382ea0d15314d7
```

It pins the plan, amendment, runner, tests, Slurm test script and both test-log
hashes before the first comparative report run.

The first locked report attempt, Slurm job `12189618`, failed closed during
input-schema validation and published zero report files. The v1 amendment had
incorrectly assumed weight total variation was undefined for every
cross-backend scenario. In fact, cross-backend full-pool scenarios retain the
same candidate identities, so microscopic weight total variation remains
defined even though the forward-model profiles differ. The failed execution,
stdout and stderr are preserved and pinned.

No endpoint, comparison, observed value or decision threshold was changed.
The post-validation-failure schema correction is explicit:

```text
slurm/configs/cohere_saxs_synthetic_baseline_analysis_amendment_v2.json
SHA-256:
b628a007ac2ecb864030198b54fb6a12e1ccd3bb3ccd94d4fafb14c54abb10bb
```

The frozen v1 runner was not edited. A small v2 shim changes only the
candidate-identity structural-NA rule. Slurm job `12189668` passed all `11/11`
v1+v2 fixture tests without recomputing the aggregate or overwriting the failed
run. The schema-corrected execution lock is:

```text
slurm/configs/cohere_saxs_synthetic_baseline_analysis_execution_v2.json
SHA-256:
2f741861f127f146bc8ff3a39176480cb3dda9f783f9826fd12edadf26ce87ca
```

Only the already reported aggregate counts, schemas and frozen file hashes
were inspected before registration. No baseline ranking, paired difference,
full-fit metric or MaxEnt complete-case metric was read. The plan pins all
eight calibration input hashes and requires hash verification before any TSV
metric is parsed.

The primary analysis unit is one `scenario_id`; q folds are not treated as
independent replicates. The primary endpoint is frozen as pooled outer
held-out chi2. All 280 complete balanced-prior, two-mode and NNLS scenarios
are compared by deterministic paired summaries and registered strata. Full-fit
curve error, compact-mass error, ESS, maximum weight and total variation are
secondary diagnostics. MaxEnt is excluded from every ranking and decision;
its available rows are reported only as a selection-biased complete-case
diagnostic.

The mode-aware development decision is also frozen before result access. It
uses two matched difference-in-differences contrasts: mode shift versus the
balanced null, and mode shift plus within-mode tilt versus tilt alone. Both
must pass the registered direction, win-fraction, Holm-adjusted exact sign-test
and cross-scenario-family replication requirements. A pass authorizes a newly
named, complexity-controlled successor, not a two-mode biological claim. A
failure routes the mode-aware branch to calibrated benchmark/no-go.

The original certificate contribution is now restored as an independent
required branch. After the frozen baseline report, a separate plan must
preregister a raw CA-Rg observable, deterministic chain-geometry bounds,
known-logQ target/proposal semantics, automatic conservative weight cap,
sample counts, confidence level and a quantitative non-vacuity gate. That
experiment proceeds regardless of the mode-aware decision. The fixed
`scale=0.25` WMB efficiency failure and unrestricted-MaxEnt failure do not
constitute a certificate result.

Generated analysis outputs will live physically under
`/usr/xtmp/lz280/cohere-idp/`; the repository retains only source, configs,
documentation and stable links. SASDU37 numeric artifacts remain absent.

#### Frozen baseline-analysis result

The schema-corrected formal report job `12189685` generated and atomically
published all 11 registered files. Its launch wrapper then exited nonzero only
because `readlink -f` canonicalized `/usr/xtmp` as `/usr/project/xtmp` in a
post-publication literal-path check. It did not alter or invalidate any output.
Independent read-only verifier job `12189688` subsequently completed on
`fennario-01` and checked every output hash, the exact 11-file set, total byte
count `2,802,128`, decision, embargo declaration and physical/project-link
mapping. The frozen report is:

```text
slurm/work/cohere_saxs_synthetic_baseline_analysis_v1/
provenance.json SHA-256:
8ce13c8a07eff50c57c068876fe6f2ba0b2e79487781c8c32324b75a8092f937
successor_decision.json SHA-256:
be9843e9b70c5927638b6d29ae6dadf287944cfc50165b28362e7a5fb1092f72
```

The registered low-complexity mode-aware gate is **no-go**. Two-mode strongly
improves held-out prediction for pure mode shift: the matched
mode-shift-versus-balanced contrast has all `28/28` differences negative,
median `-1.9403`, and Holm-adjusted one-sided sign-test
`p=7.45e-9`. The second required contrast reverses direction under combined
mode shift plus within-mode Rg tilt: only `2/28` differences are negative,
median `+0.8898`, and adjusted `p=0.9999999`. This failure replicates across
matched, off-library and cross-backend strata. Overall, two-mode beats the
prior in only `108/280` scenarios (four ties) with median held-out delta
`+0.1468`; NNLS is also worse than the prior overall despite favorable
full-fit diagnostics in some strata. Full-fit reconstruction therefore cannot
rescue the preregistered held-out gate.

This is a no-go for developing the proposed two-mode successor from this
calibration, not a universal no-go for every low-dimensional or entropy-based
method. MaxEnt complete cases remain selection-biased diagnostics and did not
enter any rank or decision. The correct next step is the independently
required certificate information-content experiment.

### 10.58 Preregistered certified known-logQ Rg information test (2026-07-21)

Before implementing the generic-IDP certificate bridge or reading any
certificate interval, the full protocol was frozen at:

```text
slurm/configs/cohere_known_logq_rg_certificate_plan_v1.json
SHA-256:
6b7f62c5717a911ae5b36d05cc021856b697e236a318c7f34e49ca258a545e27
```

The experiment uses the already declared full 116-residue SASDXC6 coarse
three-state base target multiplied by deterministic CA hard support. It uses
the exact local baseline proposal, not the learned steric proposal that failed
§10.42. One fixed seed generates `65,536` samples; literal prefixes of
`1,024`, `4,096`, and `16,384` are diagnostic only, while the final row is the
sole primary 95% certificate.

Raw CA Rg has lower bound zero. The deterministic materializer uses a maximum
adjacent-CA path length of
`1.525 + 1.329 + 1.458 = 4.312 A`. Together with
`Rg^2=N^-2 sum_{i<j} d(i,j)^2` and
`d(i,j)<=4.312(j-i)`, this gives the registered upper bound
`4.312 sqrt((116^2-1)/12)`, approximately `144.38763 A` (implementation
rounding may only increase it).

The log-weight cap is automatic. A target-temperature component may reuse its
local WMB cap only after exact factor equality establishes that the proposal
base model pointwise dominates the hard-support target. If target and proposal
factors differ, every component must instead use the universally safe
`targetLogZUpper - log q_lower` bound before applying
`q_mix >= alpha_i q_i`; no target-temperature shortcut is allowed.

The primary result is `informative_go` only if the certified denominator is
positive and interval width is at most 25% of the deterministic physical Rg
range. A positive-denominator interval strictly narrower than the physical
range but wider than 25% is `finite_but_wide`: the route remains alive but
needs a separately registered efficiency/MoRF-scale follow-up. A nonpositive
certified denominator or full-range interval is `vacuous_no_go` for this
full-chain baseline. Thresholds, sample budget and confidence level cannot be
revised after output access.

The estimated retained output is at most 8 files and 50 MiB, with no PDBs.
Physical output is registered under `/usr/xtmp/lz280/cohere-idp/`; only a
stable project link and small source/config/audit files belong under the home
workspace. Focused implementation-validation job `12189785` completed `0:0`
on `fennario-01` and passed all 59 targeted Java tests. The tests include a
constructed counterexample proving that a beta=1 local cap can be invalid when
proposal and target factor models differ, coverage by the new generic fallback,
the deterministic Rg bound, prefix certificates, and fresh atomic publication.

The one-shot execution lock pins all relevant source/test/script and validation
log hashes:

```text
slurm/configs/cohere_known_logq_rg_certificate_execution_v1.json
SHA-256:
db5b75229c02d5e3438a37f7c7ea154ad7e94f6fcc63c5b3b94adbfef6bdc7b7
```

Formal certificate job `12189812` has been submitted. No certificate result
had been read when this execution lock and job were created. SASDU37 numeric
data remain sealed.

#### Preserved pre-sampling failure and mechanical correction

Job `12189812` failed `1:0` on `fennario-01` after 17 seconds, before any
sample or certificate row was generated. `Files.createDirectories` had been
called on `slurm/work`, which is itself an existing directory symlink, and
raised `FileAlreadyExistsException`. The final output, staging output and
project link were all absent. This was a publication-path implementation
failure, not a certificate result. The failed execution and logs remain
unchanged; their stdout and stderr SHA-256 values are respectively
`fe2f3adaf2cd5a863b851dfb873d828234462f3ed1ba34dacc978082baf375f4`
and
`60e2a01e384832a046a7382325b93e04f5d0b8a205978d95b9b3033a40cacf3b`.

The versioned amendment was frozen as:

```text
slurm/configs/cohere_known_logq_rg_certificate_execution_amendment_v1.json
SHA-256:
6b64f7ad527d48ecacf6275b463ac175316c152bf09893f9722bf394fecc1c09
```

It permits only treating an existing `Files.isDirectory` path, including a
followed directory symlink, as an already satisfied parent. No scientific
field changed. Regression job `12189813` then completed `0:0` on
`fennario-01` in 49 seconds and again passed all 59 targeted Java tests,
including atomic publication through a symlink parent. Its stdout SHA-256 is
`0f5f289691c6a53336289558664e025dd6342f23189d41a94f34b9c0544db539`;
stderr retained the expected Gradle warning log with SHA-256
`d5af0caca94a4cd5ef66ff9e472cfa31f2c6bba65e0360ac8e60ded9c086f894`.

The corrected one-shot lock is:

```text
slurm/configs/cohere_known_logq_rg_certificate_execution_v2.json
SHA-256:
58d383de06c378895514d1625695c177a9a43a63acca0f26b8c5b2f704484cf2
```

It retains the original plan, target, proposal, support, seed, bounds, cap
rules, sample count, confidence level and decision threshold byte-for-byte.

#### Formal certificate result

Formal job `12189826` completed `0:0` on `fennario-01` in 25 seconds.
Independent read-only verifier job `12189827` completed `0:0` on the same
node in four seconds. The verifier checked the exact eight-file set, every
registered source/config hash, all `65,536` sample rows, proposal/target
factor equality, support semantics, state counts, log target, exact logQ,
log weights, cap compliance, interval arithmetic, provenance and the absence
of SASDU37 numeric access. The frozen project link is:

```text
slurm/work/cohere_known_logq_rg_certificate_v1
-> /usr/xtmp/lz280/cohere-idp/cohere_known_logq_rg_certificate_v1

files: 8
bytes: 22,296,152

certificate.json SHA-256:
8b16615059a37c4a64eca53f765abddfda6714152c48a2afae2e17fca5ddc960

provenance.json SHA-256:
6ef081e0fc294941f6324a0c71e5466a8283845d7516f539c4e632d7e0539984
```

Of `65,536` proposals, `3,788` passed deterministic hard support. Exact
target/proposal factor equality justified the
`pointwise_dominated_same_factor_model` cap with
`targetLogZUpper=0` and `logWeightUpper=0`. The maximum observed log weight
was `2.13e-13`, floating-point roundoff within the frozen `1e-12` audit
tolerance, and the cap was conservatively respected.

The registered checkpoints were:

| samples | role | estimate (A) | 95% certified interval (A) | relative physical width | decision |
|---:|---|---:|---:|---:|---|
| 1,024 | diagnostic | 32.5740 | `[0, 144.3876]` | 1.0000 | diagnostic only |
| 4,096 | diagnostic | 32.5394 | `[15.5194, 60.4291]` | 0.3110 | diagnostic only |
| 16,384 | diagnostic | 32.0732 | `[24.8203, 41.0248]` | 0.1122 | diagnostic only |
| 65,536 | primary | 31.9669 | `[28.6376, 35.6422]` | 0.04851 | **informative_go** |

Thus the omitted information-content experiment has a decisive positive
answer at its fixed full-chain budget: the certificate is not merely finite
but approximately 4.85% of the deterministic physical Rg range, far below
the preregistered 25% gate. The small prefixes also show why the result is an
information-efficiency statement rather than a formal-existence claim: the
1,024-sample certificate is vacuous even though the exact-logQ machinery is
valid.

#### Preregistered independent replication

The development seed was excluded before independent validation. Four new
seeds, with no replacement or pooling, were frozen at:

```text
slurm/configs/cohere_known_logq_rg_certificate_validation_plan_v1.json
SHA-256:
873f9cb5ddf86f3fdce12bc584db7cfa6c64f29f64c9a05f3979d841e97b0145

slurm/configs/cohere_known_logq_rg_certificate_validation_execution_v1.json
SHA-256:
05170e0213c2aee074ffb3c13050027bac0be76c45d9d4419a9d8f91dc4747ff
```

Every seed retained the same target, proposal, hard support, automatic cap,
`65,536`-sample budget, checkpoints, 95% confidence level and 25% threshold.
The replication gate required all four runs to pass independently, all caps
and denominators to certify, a nonempty four-interval common intersection,
and overlap of every validation interval with the development interval.

Replication job `12189854` completed `0:0` on `fennario-01` in 87 seconds.
The independent verifier `12189856` completed `0:0` on `fennario-01` in 15
seconds after scanning all `262,144` validation sample rows and verifying all
32 files (`89,176,525` bytes):

| seed | support-pass samples | estimate (A) | 95% certified interval (A) | relative physical width | decision |
|---:|---:|---:|---:|---:|---|
| 20260801 | 3,709 | 32.1025 | `[28.7209, 35.8399]` | 0.04930 | informative_go |
| 20260802 | 3,674 | 32.1416 | `[28.7373, 35.9059]` | 0.04965 | informative_go |
| 20260803 | 3,723 | 32.2266 | `[28.8406, 35.9681]` | 0.04936 | informative_go |
| 20260804 | 3,656 | 32.3500 | `[28.9167, 36.1475]` | 0.05008 | informative_go |

The common validation-interval intersection is
`[28.9167175584, 35.8399244185] A`; the estimate range is
`[32.1025356306, 32.3499742795] A`. Every validation interval overlaps the
development interval. The frozen gate therefore returns **pass**. Run stdout
and verifier stdout SHA-256 values are respectively
`cf1418fed7e4e5c5df970726af0de25d92c7a3a6c8950a94f2badeea5c627535`
and
`ecb042e5f97320c4cd22954f61bfa58a02e70a28f745d2d8247409492c7dfb17`;
both stderr logs are empty.

#### Final main-method decision and paper boundary

The post-result decision record is frozen at:

```text
slurm/configs/cohere_idp_main_method_decision_v1.json
SHA-256:
b9e21af95543c2a1c3e7e84e24d7bdaf0e5befabd55345d6b4bbc652495657d0
```

The current method-development branch is complete and routes to
**calibrated benchmark/no-go with a replicated certificate method**:

- unrestricted per-conformer KL-MaxEnt/MAP is not the paper primary method;
- the proposed two-mode successor is a no-go under its frozen matched-signal
  gate, because its pure mode-shift success reverses under mode-plus-tilt;
- NNLS remains an oracle/integrity diagnostic, not a validated successor;
- certified known-logQ raw CA-Rg estimation is a positive, independently
  replicated methods contribution for this declared 116-residue coarse
  hard-support target at `N=65,536`;
- no observable-conditioned SAXS successor passed, so SASDU37 numeric curve
  and archive remain sealed.

The supported positive claim is intentionally narrow: under the declared
target/proposal pair, the implemented ratio certificate reproducibly returns
an informative finite-sample 95% interval for raw CA-Rg. It does not establish
a physical or unique microscopic ensemble, a successful SAXS-conditioned
refinement method, learned-proposal efficiency, two biological populations,
or informativeness for arbitrary observables, targets, chain lengths or
budgets. Conversely, the benchmark is a no-go for the tested unrestricted
MaxEnt and proposed two-mode routes, not for all entropy-based or
low-dimensional methods.

Paper writing and reproducibility packaging can now begin from these frozen
claims. Any later successor is a separately preregistered research branch;
it is not a post hoc repair or unfinished requirement of this decision.

### 10.59 New active branch: generic certificate-based observable inference (2026-07-21)

The user explicitly opened a new algorithm-development branch after the
§10.58 benchmark/no-go decision. This does not modify or reinterpret that
frozen result. The new objective is stronger: build a genuinely
observable-conditioned, certificate-based inference algorithm that supports
bounded observables in general and beats balanced prior, two-mode and NNLS on
the identical held-out task. A narrow Rg-information result alone is no longer
sufficient.

Before COPPER source implementation or any COPPER result, the complete v1
protocol was frozen at:

```text
slurm/configs/cohere_copper_observable_inference_plan_v1.json
SHA-256:
ee1558cfe11b32ac5c49d7f6a68812847c225baf7659896451900ee4e990ecb3
```

COPPER expands to **Certified Observable-Partition Prior Ensemble
Reweighting**. Its core accepts state IDs, a positive known proposal mass,
scalar or vector per-state observable values, deterministic coordinate bounds
and an observation-model plugin. SAXS is the first full adapter because its
280 scenarios and three complete baselines are already frozen; Rg,
end-to-end distance and compact occupancy form the preregistered non-SAXS
genericity suite.

For a fixed training set, COPPER recursively partitions whitened observable
space by deterministic proposal-weighted principal-component median splits.
Within each leaf it preserves the proposal conditional distribution. It fits
only the leaf masses, with nested-CV selection over
`K in {1,2,4,8,16}` and shrinkage
`gamma in {0,0.25,0.5,0.75,1}`. For leaf prior mass `Q_k` and fitted mass
`m_k`, the target is
`p_i=m_k q_i/Q_k`, so `log(p_i/q_i)=log(m_k)-log(Q_k)` and the exact global
cap is the maximum of those leaf ratios. This avoids an unproved cap and is
not unrestricted per-conformer MaxEnt.

The primary prediction uses `65,536` proposal samples per outer fold and the
same four contiguous held-out blocks. Simultaneous coordinate certificates
use deterministic support bounds and propagate to a held-out chi2 interval.
Success requires both classes of gate:

- all 280 scenarios and 1,120 folds complete;
- COPPER has negative median paired held-out-chi2 delta, at least 60% non-tie
  wins and Holm-adjusted one-sided `p<=0.05` against **each** of prior,
  two-mode and NNLS;
- direction replicates in all three scenario families, mode-plus-tilt beats
  both prior and two-mode, and the balanced-prior null is protected;
- every cap and denominator certifies, at least 95% of scenarios have complete
  simultaneous exact-target coverage, and the registered interval-width
  gates pass;
- the non-SAXS scalar/vector adapters pass through the identical generic API.

A failed v1 candidate is preserved and followed by a separately registered
mathematical v2; its gates may not be weakened after results. The active goal
does not permit an early return to benchmark/no-go. All generated bytes remain
under `/usr/xtmp/lz280`, all computation remains Slurm-only, unrestricted
MaxEnt remains closed, and SASDU37 numeric artifacts remain sealed until both
development and independent validation pass.
