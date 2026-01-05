# Phase 2 Subtree Caching Test Summary

## Test Configuration Changes

### Previous Configuration (Why Cache Failed)
- **Scale**: 5 flexible residues only
- **Problem**: Too small for meaningful subtree caching
  - Only 4 possible subtrees after branch decomposition
  - High RCs variation → low cache hit rate
  - Fragment filtering (`conf.size() < 3`) skipped most calls

### New Configuration
- **Scale**: 1 mutable + 7 flexible residues
- **Rationale**:
  - **1 Mutable position (G648)**: Creates sequence space with 3 mutations (ALA, VAL, LEU)
  - **7 Flexible positions**: Provides larger conformations for better subtree reuse
  - **Larger conf space**: More opportunities for subtree sharing across conformations

## Branch Decomposition Strategy

### Current Algorithm (Simple Binary Split)
```
Algorithm: Greedy Balanced Partitioning
1. Sort positions
2. Split at mid = positions.size() / 2
3. Find separator (positions with cross-partition edges)
4. Recursively build subtrees
```

### Example for 8 positions (1 mutable + 7 flexible):
```
                    Root [0-7]
                   sep={3,4}
                  /          \
           [0,1,2,3]        [4,5,6,7]
           sep={1,2}        sep={5,6}
           /      \          /       \
        [0,1]    [2,3]    [4,5]   [6,7]
        sep={1}  sep={3}  sep={5}  sep={7}
        /  \     /  \     /  \     /  \
       [0] [1]  [2] [3]  [4] [5]  [6] [7]
```

**Valid Subtrees (internal nodes with size > 1)**:
- Root: [0-7] (8 positions) - entire conformation
- [0-3], [4-7] (4 positions each)
- [0,1], [2,3], [4,5], [6,7] (2 positions each)

## Why Larger Scale Should Work Better

### 1. More Subtrees Available
- **5 residues**: Only 4 subtrees
- **8 residues (1+7)**: 7 subtrees (excluding root)

### 2. Better Reuse Opportunities
When minimizing different sequences:
- **Sequence 1**: G648=ALA, others wild-type
- **Sequence 2**: G648=VAL, others wild-type
- **Shared subtrees**: Positions [1-7] have same RCs → Cache hit!

### 3. Fragment Filtering Fixed
The code skips cache for `conf.size() < 3`, but:
- Full conformations now have 8 positions → Always use cache
- Energy matrix fragments (singles/pairs) correctly skip cache

## Expected Results

### Cache Statistics to Watch
```
=== Subtree DOF Cache Statistics (TRUE SUBTREE CACHING) ===
Total conformation queries:   [should be > 0]
Total subtree queries:        [should be > 0]
Full cache hits:              [goal: > 0]
Partial cache hits:           [goal: > 0]
Subtree hit rate:             [goal: > 20%]
```

### Performance Comparison
Testing on scales: 7 and 9 flexible residues

| Scale | Original | Phase 1 (DP) | Phase 1+2 | P2 vs P1 Speedup |
|-------|----------|--------------|-----------|------------------|
| 7 res | ?        | ?            | ?         | Target: 1.1-1.3x |
| 9 res | ?        | ?            | ?         | Target: 1.2-1.5x |

## Key Improvements Made

### 1. Fixed Configuration
- ✅ Changed from 5 flexible to 1 mutable + 7 flexible
- ✅ Commented out small-scale (5 res) test
- ✅ Added mutations to create sequence space

### 2. Confirmed CCD Usage
- ✅ SimpleCCDMinimizer is always used (no threshold check)
- ✅ CCD called for all conformations with DOFs > 0
- ✅ No special logic for 3+ residues

### 3. Test Configuration
```java
// G648: 1 mutable position (ALA, VAL, LEU)
protein.flexibility.get("G648").setLibraryRotamers("ALA", "VAL", "LEU").setContinuous();

// G649-G652, A172, A156, A192: 7 flexible positions (wild-type only)
// This creates conformations where many subtrees are shared across sequences
```

## Files Modified

1. **TestDPvsOriginal.java**
   - Modified `buildConfSpace()` to add 1 mutable + flexible
   - Changed test scales from {5,7,9} to {7,9}
   - Commented out `testPhase2SubtreeCaching()` (5 res test)

2. **run_phase2_1mut7flex_test.sh**
   - New SLURM submission script
   - Tests comprehensive comparison with new config

## Running the Test

```bash
cd /home/users/lz280/IdeaProjects/OSPREY3
sbatch run_phase2_1mut7flex_test.sh
```

Monitor with:
```bash
squeue -u lz280
tail -f phase2_*.out
```

## Debug Information

The code includes debug output:
```java
[CachedMinimizer DEBUG #N] USING CACHE for conf.size()=X
[Phase 2] TRUE Subtree DOF Cache initialized
[Phase 2] Branch decomposition: BranchDecomposition[nodes=X, leaves=Y, internal=Z, branchWidth=W]
```

Look for these messages to confirm cache is being used!

## Next Steps

1. **Run the test** and check cache statistics
2. **If still no cache hits**, investigate:
   - Check if `subtrees.isEmpty()` in `getSubtrees()`
   - Add debug logging in `SubtreeDOFCache.minimizeWithCache()`
   - Print actual subtree configurations and RCs
3. **If cache works**, analyze speedup and iterate on larger scales

## Expected Bottlenecks

Even with cache working, Phase 2 might be slower than Phase 1 due to:
- **Cache lookup overhead**: HashMap lookups for each subtree
- **Boundary refinement**: Extra minimization at subtree boundaries
- **ConstrainedMinimizer overhead**: Setting up constrained optimization

Target is **modest speedup (10-50%)**, not dramatic improvement.
