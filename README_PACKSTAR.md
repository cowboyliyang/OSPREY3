# PACK*

Probably Approximately Correct partition-function estimation for provable
ensemble-based protein design, implemented in OSPREY 3.

PACK* fits a tractable proposal over conformational assignments, normalizes
and samples it with branch-decomposition dynamic programming, and corrects
to the fixed continuously minimized target through importance sampling
(classical free-energy perturbation). State partition-function intervals
propagate to K* binding-score intervals.

The original [OSPREY README](README.md) remains unchanged. This fork retains
OSPREY's [GPLv2 license](LICENSE.txt) and [citation information](CITING_OSPREY.txt).

## Method and assumptions

Proposal fitting changes the proposal energy E_eta, not the target CCD energy
E_t. With xi = E_t - E_eta and w = exp(-xi/RT), the partition-function identity
is q_t = q_eta * E_p_eta[w]. Training, proposal refits, and sample-size planning
precede the independent final batch.

Using frozen scales mu and kappa, PACK* separates the weight into the clipped
component y = min(w/(mu*kappa),1) and excess u = max(w/(mu*kappa)-1,0).
It combines an empirical-Bernstein interval for E[y] with a Clopper-Pearson
upper bound for the exceedance probability pi = P(u>0).

The upper endpoint assumes E[u | u>0] <= m_u. Its excess contribution is
m_u * piUpper. This is a conditional mean-excess bound, not a maximum-weight
bound. The final-batch diagnostic can reject the assumption; failing to reject
does not establish it. Coverage requires the stated assumption, supported
proposal, fixed pre-final choices, and conditionally independent final draws.
It concerns the specified finite conformational space and minimization model,
not force-field accuracy or experimentally measured binding affinity.

## Build and test

Use JDK 17 and the included Gradle wrapper. CPU execution is supported; GPU DP
and sampling additionally require a compatible CUDA runtime/device. See the
original OSPREY documentation linked from [README.md](README.md) for its native
dependencies.

From the repository root, define BUILD_DIR as an absolute directory for build
outputs, caches, and temporary files. Choose resource allocations and Java heap
sizes appropriate to your environment.

```bash
#!/bin/bash
set -euo pipefail
: "${BUILD_DIR:?Set BUILD_DIR to an absolute build directory}"
mkdir -p "$BUILD_DIR/tmp"
export TMPDIR="$BUILD_DIR/tmp"
export GRADLE_USER_HOME="$BUILD_DIR/gradle-cache"
export OSPREY_GRADLE_XTMP_BUILD_ROOT="$BUILD_DIR/build"
export PACKSTAR_CLASSPATH_FILE="$BUILD_DIR/test_classpath.txt"
./gradlew test writeExperimentClasspath --no-daemon --max-workers=1 \
  --project-cache-dir "$BUILD_DIR/project-cache" \
  -I slurm/scripts/gradle_xtmp.init.gradle \
  -I slurm/h200/classpath.init.gradle \
  -Dorg.gradle.jvmargs="-Xmx8g -Djava.io.tmpdir=$TMPDIR" -DtestMaxHeap=8g \
  --tests edu.duke.cs.osprey.kstar.TestKStarScoreLog10 \
  --tests edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC \
  --tests edu.duke.cs.osprey.packstar.TestPackStarAutoPreflight \
  --tests edu.duke.cs.osprey.packstar.TestPackStarAdmissionDecision \
  --tests edu.duke.cs.osprey.packstar.TestPackStarTripleDecompositionCosts \
  --tests edu.duke.cs.osprey.branchdp.TestBranchDpAdmission \
  --tests edu.duke.cs.osprey.branchdp.TestBranchDecompositionStrategies
```

The generated test classpath is written to BUILD_DIR/test_classpath.txt. The
bundled Gradle init scripts redirect build outputs and export that classpath.

## Run a design

The configurable design driver is
`edu.duke.cs.osprey.markstar.bench.GenericPDBBench`; it lives in the test source
set, so use the exported **test** classpath. Supply a prepared PDB with residue
IDs matching your mutable/flexible lists. Explicit chain selections avoid
ambiguity in multichain systems. Mutable specifications use entries such as
`A42=ALA,SER`, separated by semicolons; flexible lists contain residue IDs.

Define BUILD_DIR, PDB_PATH, and OUTPUT_DIR as absolute paths. BUILD_DIR is the
directory used in the build step. Set NUM_CPUS to the number of available CPU
workers. The residue and chain selections below are placeholders for your own
system.

```bash
#!/bin/bash
set -euo pipefail
export TMPDIR="$OUTPUT_DIR/tmp"
mkdir -p "$OUTPUT_DIR" "$TMPDIR"
java --add-opens java.base/java.util=ALL-UNNAMED \
  --add-opens java.base/java.lang=ALL-UNNAMED \
  --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
  -Xmx16g -Djava.io.tmpdir="$TMPDIR" \
  -Dosprey.bench.method=packstar \
  -Dosprey.bench.designId=my_design \
  -Dosprey.bench.pdbPath="$PDB_PATH" \
  -Dosprey.bench.proteinChains=A -Dosprey.bench.ligandChains=B \
  '-Dosprey.bench.mutable=A42=ALA,SER' \
  '-Dosprey.bench.flexible=A43;B15' \
  -Dosprey.bench.outputDir="$OUTPUT_DIR" \
  -Dosprey.bench.numCPUs="${NUM_CPUS:-4}" \
  -Dpackstar.cutoff.strategy=RESIDUAL_BUDGET \
  -Dpackstar.cutoff.residualBudget=1 \
  -Dpackstar.dp.gpu=false -Dpackstar.pac.sampling.gpu=false \
  -Dpackstar.dp.mmap.dir="$OUTPUT_DIR/dp_mmap" \
  -Dpackstar.pac.frequencySeverity.tripleEta=false \
  -Dpackstar.pac.frequencySeverity.outputDir="$OUTPUT_DIR/diagnostics" \
  -cp "$(<"$BUILD_DIR/test_classpath.txt")" \
  edu.duke.cs.osprey.markstar.bench.GenericPDBBench
```

The residual budget controls interaction sparsification, not statistical
confidence. Sparse-target and full-interaction runtimes/scores should be
interpreted with that model difference in mind. Memory and runtime depend on
rotamer counts and decomposition table sizes, not only on position count.

Important settings:

| JVM property | Meaning |
| --- | --- |
| `packstar.pac.confidence` | Per-state failure-probability budget delta, not coverage probability. |
| `packstar.pac.targetEpsilon` | Requested state interval width, measured as 1 - lower/upper. |
| `packstar.pac.randomSeed` | Base sampling seed. |
| `packstar.pac.maxEstSamples` | Maximum final-batch size. |
| `packstar.pac.frequencySeverity.severityCap` | Conditional mean-excess bound m_u; default 20 is an assumption. |
| `packstar.pac.frequencySeverity.tripleEta` | Optional higher-order proposal corrections; false selects the pair-only route. |
| `packstar.dp.gpu` / `packstar.pac.sampling.gpu` | GPU DP / proposal sampling switches. |

Historical property names, diagnostic column names, and the statistical helper
class name are retained for compatibility. In the code and manuscript, their
statistical meanings are conditional mean excess, clipped component, and
exceedance probability. For simultaneous sequence/state guarantees, allocate
failure probability across every interval being used; a per-state delta is
not automatically a design-wide guarantee.

## Results and source map

The design driver writes `my_design_packstar.csv`: 29 columns covering sequence
rank, log10 score/interval endpoints, seven fields for each of Protein, Ligand,
and Complex, and elapsed time. Each state's fields are partition-function
lower/upper bounds, status, epsilon, conformation count, scoring count, and
partial-minimization count. Historical 41-column files require their original
reader or explicit conversion; do not mix schemas in a shard merge.

`result.numConfs` and `getNumConfsEvaluated()` report actual full-conformation
CCD calls, including work completed before a certificate abort. Result stats
also expose `fullCCD`, `ccdSampleRecords` (draws including duplicates and cache
hits), and `sampleCallbacks` (successfully delivered listener calls).
Registering `setSampleListener` delivers one trace per draw, in draw order on
the compute thread, after each CCD batch and before certificate checks. Stages
are `TRAIN`, `DISCOVERY`, `VALIDATION`, and `ESTIMATION`; the sample index is
zero-based across the run. Listener exceptions propagate to the caller.
Paths that skip sampling emit no traces. The small-system coverage tool requires
trace counts to match the recorded draws; zero draws are reported as
`energy_check=NOT_SAMPLED` with `max_energy_error=NaN`, never as a passed energy
comparison. Its `covered` column only reports interval containment; interpret
it together with `status`, especially for aborted runs with unbounded intervals.

A normally exiting workload can contain `INCOMPLETE_ESTIMATES`: all configured
rows can be present even when some state estimates do not meet the requested
precision. Report workload completion and estimated-sequence counts separately.

| Source | Role |
| --- | --- |
| `src/main/java/edu/duke/cs/osprey/packstar/` | Public entry points, adaptation, final estimation, and observables. |
| `PackStarEstimator.java` | Proposal fitting/refits, fixed final batch, and target correction. |
| `PackStarFrequencySeverityPAC.java` | Clipped-component interval, exceedance bound, and rejection diagnostic. |
| `src/main/java/edu/duke/cs/osprey/branchdp/` | Shared graph decomposition, normalization, and sampling machinery. |
| `slurm/wrappers/` | Sequence-shard validation, result merging, and comparison. |

For an existing Java OSPREY conformational space, the primary partition-function
entry point is `PackStarPartitionFunction`, which implements
`PartitionFunction.WithConfDB`. MARK* remains available as a baseline through
the same design driver with `-Dosprey.bench.method=markstar`.
