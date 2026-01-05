# NegatedEnergyMatrix Investigation Results

## Summary

Investigated negatedHScorer calculation in MARK* that causes partition function to be 10²² times too large. Tested two fix attempts - **both failed**.

## Problem

```
K* (correct):    Z = 1.2×10⁴
MARK* (wrong):   Z = 1.4×10²⁶  (10²² times larger!)
```

Root cause appears to be in how `confUpperBound` is calculated and used.

## Current Implementation

**Line: MARKStarBound.java:368**
```java
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)
```
- Uses default constructor → `Optimizer.Minimize`
- NegatedEnergyMatrix negates all energy values
- Result: negatedHScore = 35.5 kcal/mol (large positive number)

**Line: MARKStarNode.java:405**
```java
double confUpperBound = rigidScore - negatedHScore;
                      = 0 - 35.5
                      = -35.5 kcal/mol
```

This leads to:
```
Boltzmann weight = exp(-(-35.5)/0.592) = exp(60) ≈ 10²⁶
```

## Fix Attempts

### Attempt 1: Remove NegatedEnergyMatrix + Use Maximize

**Change**:
```java
// OLD:
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)

// NEW:
new TraditionalPairwiseHScorer(rigidEmat, rcs, MathTools.Optimizer.Maximize)
```

**Rationale**:
- Remove negation, select maximum (least stable) energies directly
- confUpperBound should be highest possible energy

**Result**: ❌ **FAILED**
- Z still = 1.4×10²⁶ (same as before)
- No improvement

### Attempt 2: Change Formula from Subtract to Add

**Change**:
```java
// OLD:
double confUpperBound = rigidScore - negatedHScore;

// NEW:
double confUpperBound = rigidScore + negatedHScore;
```

**Rationale**:
- Maybe the formula itself is wrong
- Try adding instead of subtracting

**Result**: ❌ **FAILED**
- Z still = 1.4×10²⁶ (same as before!)
- No improvement

## Key Observations

1. **Neither fix changed the result** - Z remained exactly 10²⁶ in both cases
2. This suggests the problem is **NOT in the negatedHScore calculation itself**
3. The problem might be **downstream** in how confUpperBound is used

## Hypotheses for Real Problem

### Hypothesis 1: The Issue is in Boltzmann Weight Calculation

The problem may not be in confUpperBound calculation, but in how it's converted to Boltzmann weights later.

**Need to investigate**:
- `setBoundsFromConfLowerAndUpper()` in MARKStarNode.java:355-359
- How subtreeUpperBound/subtreeLowerBound are used
- BoltzmannCalculator implementation

### Hypothesis 2: Root Node Initialization is Fundamentally Wrong

MARK* may not be designed to work with a single conf space (Complex only).

**Evidence**:
- MARK* is always used in K* context (protein + ligand + complex)
- Direct use of MARKStarBound may be incorrect
- Successful MARK* tests use full MARKStar class with 3 conf spaces

**Need to check**:
- How MARKStarBoundFastQueues creates MARKStarBound
- Whether there's preprocessing in MARKStar that we're missing

### Hypothesis 3: negatedHScorer Has a Different Purpose

Maybe negatedHScorer isn't meant to calculate confUpperBound at all.

**Possible alternative interpretation**:
- It could be for tree pruning heuristics
- Or for ordering nodes in priority queue
- Not directly related to Boltzmann weights

## Relationship to "Negative Correction" Bug

Your colleague's main concerns:
1. ✅ **"Negative correction" warnings** - Fixed by applying negative corrections
2. ❓ **"Bounds are incorrect" errors** - Possibly related to this confUpperBound issue
3. ✅ **Partition function off by orders of magnitude** - This is what we're debugging

**Connection**: The confUpperBound bug and the coordinate descent accuracy issue may be separate problems.

## Next Steps (Recommendations)

### Option A: Deep Dive into Boltzmann Calculation

1. Add debug output to `setBoundsFromConfLowerAndUpper()`
2. Trace how confUpperBound → subtreeUpperBound → Boltzmann weight
3. Find where exp() is being called with the wrong argument

### Option B: Study Working MARK* Implementation

1. Look at TestMARKStar.java tests that succeed
2. Compare their setup with our test
3. Understand what MARKStar class does before calling MARKStarBound

### Option C: Consult MARK* Paper/Documentation

1. Find the original MARK* algorithm paper
2. Check the mathematical formula for confUpperBound
3. Verify whether our implementation matches the paper

## Files Modified (Now Reverted)

- MARKStarBound.java:368 - Tested different negatedHScorer configurations
- MARKStarNode.java:405 - Tested formula change (subtract → add)

**All changes should be reverted** as they didn't fix the problem.

## Conclusion

The negatedHScorer calculation itself appears to be working as designed. The real bug is likely:
1. In how confUpperBound is used to calculate Boltzmann weights, OR
2. In the fundamental design assumption about how MARK* should be initialized

The problem is **deeper than just the negatedHScorer value**.
