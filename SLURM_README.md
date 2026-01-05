# Running Tests on Duke CS Cluster with SLURM

This directory contains SLURM batch scripts for running OSPREY tests on the Duke CS compute cluster.

## Available Scripts

### 1. **submit_comets_performance_tests.sh** - All Performance Comparison Tests
Runs all COMETS vs COMETS-Z performance comparison tests (small + full + bounded memory variants).
- **Time**: 12 hours
- **Memory**: 16GB
- **CPUs**: 8
- **Tests**: All tests in `TestCometsVsCometsZPerformance`

### 2. **submit_comets_small_test.sh** - Small Tests Only (RECOMMENDED TO START)
Runs only the small test cases (4 positions, 10 sequences).
- **Time**: 2 hours
- **Memory**: 8GB
- **CPUs**: 4
- **Tests**: `compareSmall2RL0` and `compareSmall2RL0BoundedMemory`

### 3. **submit_comets_full_test.sh** - Full 2RL0 Test (LONG RUNNING)
Runs the full 2RL0 test (8 positions, 25 sequences).
- **Time**: 24 hours
- **Memory**: 32GB
- **CPUs**: 8
- **Tests**: `compare2RL0Full`
- **Warning**: This will take many hours!

### 4. **submit_mskstar_tests.sh** - Original MSKStar Tests
Runs all original MSKStar tests (kept for backward compatibility).
- **Time**: 12 hours
- **Memory**: 16GB
- **CPUs**: 8

### 5. **submit_cometsz_tests.sh** - COMETS-Z Tests
Runs all COMETS-Z tests (the renamed version of MSKStar).
- **Time**: 12 hours
- **Memory**: 16GB
- **CPUs**: 8

## How to Use

### Make scripts executable (first time only):
```bash
chmod +x submit_*.sh
```

### Submit a job:
```bash
# Start with small tests to verify everything works
sbatch submit_comets_small_test.sh

# Or run all performance tests
sbatch submit_comets_performance_tests.sh

# Or run a specific script
sbatch submit_comets_full_test.sh
```

### Check job status:
```bash
# View your jobs
squeue -u lz280

# View all jobs in the queue
squeue

# View detailed info about a specific job
scontrol show job <job_id>
```

### Cancel a job:
```bash
scancel <job_id>

# Cancel all your jobs
scancel -u lz280
```

### View output:
```bash
# Output files are created in the current directory
# Format: <job_name>_<job_id>.out and <job_name>_<job_id>.err

# View live output (updates as job runs)
tail -f comets_perf_12345.out

# View errors
tail -f comets_perf_12345.err

# Search for specific test results
grep "Speedup:" comets_perf_12345.out
grep "PASSED\|FAILED" comets_perf_12345.out
```

## Important Notes

### Adjust Email Notifications
Edit the script and change:
```bash
#SBATCH --mail-user=lz280@duke.edu
```
to your email address, or remove these lines if you don't want emails.

### Adjust Java Version
If the cluster has a different Java module, update:
```bash
module load Java/17
```

Check available modules:
```bash
module avail Java
```

### Adjust Partition
The scripts use `--partition=compsci`. Check available partitions:
```bash
sinfo
```

Common partitions on Duke CS cluster:
- `compsci` - General compute
- `compsci-gpu` - GPU nodes
- `scavenger` - Lower priority, preemptible

### Resource Requirements

If you get out-of-memory errors, increase memory:
```bash
#SBATCH --mem=32G  # or 64G
```

If tests are slow, increase CPUs:
```bash
#SBATCH --cpus-per-task=16
```

### Test Results Location

Test results are saved to:
```
build/test-results/test/TEST-*.xml
build/reports/tests/test/index.html
```

To copy results after job completes:
```bash
# Find the job output file
ls -lt comets_perf_*.out | head -1

# Copy test reports to a safe location
cp -r build/reports/tests/test ~/test_results_$(date +%Y%m%d_%H%M%S)
```

## Troubleshooting

### Job doesn't start
```bash
# Check job status
squeue -u lz280

# Common reasons:
# - Resources not available (wait in queue)
# - Invalid partition
# - Requesting too much memory/CPUs
```

### Job fails immediately
```bash
# Check error file
cat comets_perf_<job_id>.err

# Common issues:
# - Java module not found
# - JAVA_HOME not set correctly
# - Gradle daemon issues
```

### Out of memory
Increase `--mem` in the script:
```bash
#SBATCH --mem=32G  # or higher
```

### Gradle daemon issues
Add to script before running gradle:
```bash
# Stop any existing daemons
./gradlew --stop

# Run with daemon disabled
./gradlew --no-daemon test --tests "..."
```

## Example Workflow

### 1. Test locally first (quick sanity check):
```bash
# On login node (don't run long jobs here!)
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsVsCometsZPerformance.compareSmall2RL0" --dry-run
```

### 2. Submit small test to cluster:
```bash
sbatch submit_comets_small_test.sh
```

### 3. Monitor progress:
```bash
# Watch job queue
watch -n 10 squeue -u lz280

# Watch output
tail -f comets_small_<job_id>.out
```

### 4. Once small test passes, run full suite:
```bash
sbatch submit_comets_performance_tests.sh
```

### 5. Collect results:
```bash
# After job completes
grep "Speedup:" comets_perf_<job_id>.out
grep "PASSED" comets_perf_<job_id>.out

# Save results
cp comets_perf_<job_id>.out ~/results/
```

## Advanced: Running Multiple Tests in Parallel

You can submit multiple jobs to run different tests simultaneously:

```bash
# Submit different test suites
sbatch submit_comets_small_test.sh
sbatch submit_mskstar_tests.sh
sbatch submit_cometsz_tests.sh

# Check all your jobs
squeue -u lz280
```

## Questions?

Check SLURM documentation:
```bash
man sbatch
man squeue
man scancel
```

Or Duke CS cluster documentation (if available).
