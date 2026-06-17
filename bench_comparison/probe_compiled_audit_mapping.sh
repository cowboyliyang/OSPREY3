#!/bin/bash
# Probe whether S11 AuditLeafCCD CSV assignments can be mapped to a compiled ConfSpace.
set -euo pipefail
cd /home/users/lz280/IdeaProjects/OSPREY3

OUTDIR=/usr/xtmp/lz280/bench_comparison
PDBDIR=/usr/xtmp/lz280/dance_bench/pdbs_prepped
JAVA=${JAVA:-/home/users/lz280/java/jdk-17.0.2+8/bin/java}
MAIN=edu.duke.cs.osprey.markstar.bench.CompiledAuditMappingProbe

DID=${DID:-3u7y}
STATE=${STATE:-Complex}
PDB=${PDB:-$PDBDIR/$DID/${DID}.min.reduce.renum.pdb}
INPUT=${INPUT:-$OUTDIR/audit_leaves/gnn_s11_leafonly/3u7y/Complex/seq_00034.csv}
MUT=${MUT:-"G384;G382"}
FLEX=${FLEX:-"G385;L587;L649"}
COMPILE=${COMPILE:-true}
PREVIEW_CONFS=${PREVIEW_CONFS:-4}
COMPILE_THREADS=${COMPILE_THREADS:-1}
JAVA_XMX=${JAVA_XMX:-32g}
RUN_BUILD=${RUN_BUILD:-true}
CUDA_SWEEP=${CUDA_SWEEP:-false}
MAX_CONFS=${MAX_CONFS:-512}
WARMUP_CONFS=${WARMUP_CONFS:-64}
REPEATS=${REPEATS:-2}
STREAMS_LIST=${STREAMS_LIST:-"8,16,32,64,128"}
BATCH_SIZES=${BATCH_SIZES:-"256,512,1024,2048,4096"}
PRECISION=${PRECISION:-Float32}
COMPARE=${COMPARE:-false}
COMPARE_CUDA=${COMPARE_CUDA:-false}
CCD_RESULTS=${CCD_RESULTS:-}
COMPARE_OUT=${COMPARE_OUT:-}
COMPARE_PRECISION=${COMPARE_PRECISION:-Float64}
TOLERANCE=${TOLERANCE:-0.1}
INCLUDE_HYDROXYLS=${INCLUDE_HYDROXYLS:-true}
INCLUDE_NONHYDROXYL_HGROUPS=${INCLUDE_NONHYDROXYL_HGROUPS:-true}
if [ -z "${SERVICE_PORT+x}" ]; then
  if [ -n "${SLURM_JOB_ID:-}" ]; then
    SERVICE_PORT=$((44342 + (SLURM_JOB_ID % 10000)))
  else
    SERVICE_PORT=44342
  fi
fi

if [ "$RUN_BUILD" = "true" ]; then
  ./gradlew testClasses
fi

CP="build/classes/java/main:build/classes/java/test:build/classes/kotlin/main:build/classes/kotlin/test:build/resources/main:build/resources/test"
for jar in lib/*.jar; do CP="$CP:$jar"; done
for jar in $(find /home/users/lz280/.gradle/caches/modules-2 -name "*.jar" -path "*/files-2.1/*" 2>/dev/null | sort -u); do CP="$CP:$jar"; done

"$JAVA" \
  --add-modules jdk.incubator.foreign \
  --add-opens java.base/java.util=ALL-UNNAMED \
  --add-opens java.base/java.lang=ALL-UNNAMED \
  --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
  -Xmx"$JAVA_XMX" \
  -Dosprey.bench.pdbPath="$PDB" \
  -Dosprey.bench.mutable="$MUT" \
  -Dosprey.bench.flexible="$FLEX" \
  -Dosprey.bench.designId="$DID" \
  -Dosprey.audit.state="$STATE" \
  -Dosprey.audit.input="$INPUT" \
	  -Dosprey.service.port="$SERVICE_PORT" \
	  -Dosprey.compiledAudit.previewConfs="$PREVIEW_CONFS" \
	  -Dosprey.compiledAudit.compile="$COMPILE" \
	  -Dosprey.compiledAudit.compileThreads="$COMPILE_THREADS" \
	  -Dosprey.compiledAudit.cudaSweep="$CUDA_SWEEP" \
	  -Dosprey.compiledAudit.maxConfs="$MAX_CONFS" \
	  -Dosprey.compiledAudit.warmupConfs="$WARMUP_CONFS" \
	  -Dosprey.compiledAudit.repeats="$REPEATS" \
	  -Dosprey.compiledAudit.streamsList="$STREAMS_LIST" \
	  -Dosprey.compiledAudit.batchSizes="$BATCH_SIZES" \
	  -Dosprey.compiledAudit.precision="$PRECISION" \
	  -Dosprey.compiledAudit.compare="$COMPARE" \
	  -Dosprey.compiledAudit.compareCuda="$COMPARE_CUDA" \
	  -Dosprey.compiledAudit.ccdResults="$CCD_RESULTS" \
	  -Dosprey.compiledAudit.compareOut="$COMPARE_OUT" \
	  -Dosprey.compiledAudit.comparePrecision="$COMPARE_PRECISION" \
	  -Dosprey.compiledAudit.tolerance="$TOLERANCE" \
	  -Dosprey.compiledAudit.includeHydroxyls="$INCLUDE_HYDROXYLS" \
	  -Dosprey.compiledAudit.includeNonHydroxylHGroups="$INCLUDE_NONHYDROXYL_HGROUPS" \
	  -cp "$CP" \
	  "$MAIN"
