# Benchmark Guide: MARK* and BranchMARK* with GNN Integration

## Overview

The benchmark framework compares partition function computation strategies
for protein design using MARK* and BranchMARK* algorithms, with optional
GNN (Graph Neural Network) energy surrogates for accelerated convergence.

All benchmarks are run through a single entry point:
`edu.duke.cs.osprey.markstar.RunBenchmark` (for MARK* strategies) or
`edu.duke.cs.osprey.markstar.RunBranchBenchmark` (for BranchMARK* strategies).

## Available Strategies

| Strategy | Description | GPU Required |
|----------|-------------|:---:|
| `ccd` | MARK* with CCD minimization (baseline) | No |
| `single` | MARK* with single-conf GNN inference | Yes |
| `strategy4` | A* scan + batch GNN + selective CCD | Yes |
| `strategy5` | Hybrid subtree GNN expansion | Yes |
| `strategy6` | GNN + Conformal Prediction pool | Yes |
| `strategy7` | Decoupled GNN pool (best throughput) | Yes |
| `branch_ccd` | BranchMARK* with CCD only (no GNN) | No |
| `branch_gnn` | BranchMARK* with decoupled GNN pool | Yes |

## System Properties

### Required

| Property | Description |
|----------|-------------|
| `osprey.gnn.benchStrategy` | Strategy to run (see table above) |

### Optional

| Property | Default | Description |
|----------|---------|-------------|
| `osprey.gnn.confSpace` | `medium` | Conformation space: `medium` (108 sequences) or `highrot` (625 sequences, large rotamer count) |
| `osprey.gnn.numSeqs` | `20` | Number of sequences to compute K* scores for |
| `osprey.gnn.numCPUs` | `4` | Number of CPU threads for minimization |
| `osprey.gnn.gpuBatchSize` | `1000` | Accumulate this many conformations before firing one ONNX GPU batch |
| `osprey.gnn.ematCache` | auto | Energy matrix cache file pattern (e.g., `emat_cache/highrot_4pos.*.dat`) |
| `osprey.gnn.eval.proteinModelPath` | `gnn_data/.../protein/model/gnn_model.onnx` | Path to protein GNN ONNX model |
| `osprey.gnn.eval.complexModelPath` | `gnn_data/.../complex/model/gnn_model.onnx` | Path to complex GNN ONNX model |
| `osprey.gnn.trainingConfSpace` | (empty) | Training confspace for RC mapping when inference confspace differs: `all20`, `highrot`, or `medium` |

## SLURM Submission

### Prerequisites

1. Build the project:
   ```bash
   cd /path/to/OSPREY3
   ./gradlew testClasses
   ```

2. Ensure GNN models are trained and available under `gnn_data/`.

3. ONNX Runtime GPU library is at `lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib`.

### Using bench_s7.sh

The script `bench_s7.sh` submits SLURM jobs for `ccd` baseline and `strategy7`
with varying GPU batch sizes on the grisman partition (A5000 GPUs):

```bash
# Default batch sizes (500, 1000, 2000):
bash bench_s7.sh

# Custom batch sizes:
bash bench_s7.sh 256 512 1024
```

### Manual Submission

For any strategy, submit directly with `sbatch`:

```bash
JAVA=/path/to/java
ORT_LIB=lib/ort-cuda12/onnxruntime-linux-x64-gpu-1.20.0/lib
CUDA_LIB=/usr/local/cuda-12.8/lib64
JARGS="--add-opens java.base/java.util=ALL-UNNAMED \
       --add-opens java.base/java.lang=ALL-UNNAMED \
       --add-opens java.base/java.lang.invoke=ALL-UNNAMED \
       -Xmx32g -Xms2g \
       -Djava.library.path=$ORT_LIB:$CUDA_LIB"

# Example: BranchMARK* + GNN on grisman A5000
sbatch --partition=grisman --account=grisman \
       --gres=gpu:a5000:1 --cpus-per-task=8 --mem=32G --time=02:00:00 \
       --job-name=bench_branch_gnn \
       --output=bench_logs/branch_gnn_%j.out \
       --error=bench_logs/branch_gnn_%j.err \
       --wrap "cd /path/to/OSPREY3 \
         && export LD_LIBRARY_PATH=$ORT_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH \
         && $JAVA $JARGS \
            -Dosprey.gnn.benchStrategy=branch_gnn \
            -Dosprey.gnn.confSpace=highrot \
            -Dosprey.gnn.numSeqs=5 \
            -Dosprey.gnn.gpuBatchSize=1000 \
            -Dosprey.gnn.eval.proteinModelPath=gnn_data/2RL0_all20_4pos_merged/protein/model/gnn_model.onnx \
            -Dosprey.gnn.eval.complexModelPath=gnn_data/2RL0_all20_4pos_merged/complex/model/gnn_model.onnx \
            -Dosprey.gnn.trainingConfSpace=all20 \
            -cp \"\$(cat bench_logs/.classpath.txt)\" \
            edu.duke.cs.osprey.markstar.RunBranchBenchmark"

# Example: BranchMARK* CCD baseline (no GPU needed, but same resources for fair comparison)
sbatch --partition=grisman --account=grisman \
       --gres=gpu:a5000:1 --cpus-per-task=8 --mem=32G --time=02:00:00 \
       --job-name=bench_branch_ccd \
       --output=bench_logs/branch_ccd_%j.out \
       --error=bench_logs/branch_ccd_%j.err \
       --wrap "cd /path/to/OSPREY3 \
         && export LD_LIBRARY_PATH=$ORT_LIB:$CUDA_LIB:\$LD_LIBRARY_PATH \
         && $JAVA $JARGS \
            -Dosprey.gnn.benchStrategy=branch_ccd \
            -Dosprey.gnn.confSpace=highrot \
            -Dosprey.gnn.numSeqs=5 \
            -cp \"\$(cat bench_logs/.classpath.txt)\" \
            edu.duke.cs.osprey.markstar.RunBranchBenchmark"
```

## Output

Each run prints:
1. Energy matrix computation time
2. Per-sequence MARK*/BranchMARK* progress (epsilon convergence, leaf/internal round stats)
3. GNN pool statistics (if GNN enabled): batch count, bounded nodes, budget usage
4. Summary table: method, wall time, GNN confs evaluated, CCD confs minimized, number of sequences
5. Top K* scores ranked by lower bound

## GNN Training Pipeline

1. **Export training data** (CPU-only, no GPU):
   ```bash
   ./gradlew test --tests "...TestBranchMARKStar.exportGNNDataAllMutable" \
       -Dosprey.gnn.numSamples=200000 \
       -Dosprey.gnn.outputDir=gnn_data/2RL0_all20_4pos
   ```

2. **Train GNN model** (GPU recommended):
   ```bash
   python gnn/train.py --data-dir gnn_data/2RL0_all20_4pos_merged/protein --epochs 200
   python gnn/train.py --data-dir gnn_data/2RL0_all20_4pos_merged/complex --epochs 200
   ```

3. **Calibrate Conformal Prediction** bounds:
   ```bash
   python gnn/calibrate_cp.py --data-dir gnn_data/2RL0_all20_4pos_merged/protein
   ```
   This outputs the `cpQ` quantile value used by the GNN pool strategies.
