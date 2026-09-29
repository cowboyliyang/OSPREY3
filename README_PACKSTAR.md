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

Run compilation and tests through Slurm. Save the following as a build job;
submit it from the repository root with your site's partition/account settings.
At Duke, use `--account=grisman`. Place Slurm stdout/stderr under
`/usr/xtmp/$USER`, for example with
`--output=/usr/xtmp/$USER/packstar-build-%j.out` and a corresponding `--error`.
Allow approximately 2 GiB and 20,000 files for an initial source/dependency
cache; build outputs can require additional space.

```bash
#!/bin/bash
set -euo pipefail
RUN_DIR="/usr/xtmp/$USER/packstar-build/$SLURM_JOB_ID"
mkdir -p "$RUN_DIR/source" "$RUN_DIR/tmp"
cd "$SLURM_SUBMIT_DIR"
git archive HEAD | tar -x -C "$RUN_DIR/source"
export TMPDIR="$RUN_DIR/tmp"
export GRADLE_USER_HOME="/usr/xtmp/$USER/packstar-gradle-cache"
export PACKSTAR_CLASSPATH_FILE="$RUN_DIR/test_classpath.txt"
cd "$RUN_DIR/source"
./gradlew test writeExperimentClasspath --no-daemon --max-workers=1 \
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

The source snapshot and generated classpath remain under RUN_DIR. Allocate
enough memory for the Java heaps and compilation overhead; the repository's
local validation job uses 4 CPUs and 24 GiB. Its Slurm scripts contain Duke
paths and allocation settings and should be adapted before use elsewhere.

## Run a design

The configurable design driver is
`edu.duke.cs.osprey.markstar.bench.GenericPDBBench`; it lives in the test source
set, so use the exported **test** classpath. Supply a prepared PDB with residue
IDs matching your mutable/flexible lists. Explicit chain selections avoid
ambiguity in multichain systems. Mutable specifications use entries such as
`A42=ALA,SER`, separated by semicolons; flexible lists contain residue IDs.

Run the following inside a Slurm job after defining BUILD_DIR, PDB_PATH, and
OUTPUT_DIR as absolute paths. BUILD_DIR is the preceding RUN_DIR. The residue
and chain selections below are placeholders for your own system. Keep outputs
and temporary files under `/usr/xtmp/$USER`.

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
  -Dosprey.bench.numCPUs="${SLURM_CPUS_PER_TASK:-4}" \
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
