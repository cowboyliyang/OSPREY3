# MARK* Bug修复最终报告

## 执行摘要

通过分析团队之前的调查（来自Wiki文档）和我的深入测试，成功识别并修复了MARK*算法中的**两个关键bug**：

1. **Negative Correction Bug**: 当triple minimization找到比pairwise能量矩阵更低的能量时，代码只打印警告但没有实际修正
2. **qprime赋值错误**: qprime应该是upper-lower的差值，而不是upper bound本身

这两个bug都已修复，**所有测试现在都通过了**。

---

## Bug #1: Negative Correction处理错误

### 问题描述

**文件**: [MARKStarBound.java:1076](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java#L1076)

**根本原因**:
当triple minimization找到的能量**低于**pairwise能量矩阵的lower bound时，原代码只打印警告，但**没有实际修正**这个错误的bound：

```java
// 原始代码（Line 1071-1077）
double lowerbound = minimizingEmat.getInternalEnergy(tuple);
if (tripleEnergy - lowerbound > 0) {
    double correction = tripleEnergy - lowerbound;
    correctionMatrix.setHigherOrder(tuple, correction);
}
else
    System.err.println("Negative correction for "+tuple.stringListing());
    // ⚠️ 问题：只打印警告，但没有修正！
```

### 为什么会发生？

根据团队7/22/25的分析：

> "In theory, the minimized triplet energy should always be >= than the triplet energy calculated from the pairwise minimized energy matrix."
>
> **但实际上**: Pairwise minimization可能没有找到真正的global minimum，导致：
> - Pairwise能量矩阵的lower bound **过高**（不够optimistic）
> - Triple minimization找到了**更低**的能量（更接近真实的global minimum）
> - 结果：`tripleEnergy < lowerbound` → negative correction

### 影响

这个bug会导致：
1. **Lower bounds过高** → Partition function的bounds不准确
2. **Partition function被高估** → 可能导致算法过早认为epsilon满足条件
3. **K* estimates不准确** → 虽然complex和ligand的错误可能"抵消"，导致K*结果看起来"碰巧正确"

在我的测试中，这导致：
- MARK*的partition function比K*大**10²²倍**
- Free energy差异达到**30 kcal/mol**
- Algorithm不评估任何conformation就错误地认为收敛

### 修复方案

**实施团队提出的"Handler for root cause"策略**：

```java
// 修复后的代码（Line 1071-1084）
double lowerbound = minimizingEmat.getInternalEnergy(tuple);
double correction = tripleEnergy - lowerbound;
if (correction > 0) {
    correctionMatrix.setHigherOrder(tuple, correction);
}
else {
    // ✅ 修复：当triple energy更低时，使用更准确的triple能量
    System.err.println("Negative correction for "+tuple.stringListing() +
                     " (correction=" + correction + "). Using triple energy as correction.");
    // 存储负correction来降低bound
    correctionMatrix.setHigherOrder(tuple, correction);
}
```

**关键改进**:
- 现在**总是**调用`correctionMatrix.setHigherOrder()`，无论correction是正还是负
- 负correction会**降低**能量矩阵的lower bound，使其更接近实际的global minimum
- 这使得bounds更准确，partition function计算正确

---

## Bug #2: qprime赋值错误

### 问题描述

**文件**: [MARKStarBound.java:248](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java#L248)

**根本原因**:
`qprime`应该表示upper和lower bounds之间的**差距(gap)**，但原代码直接赋值为upper bound：

```java
// 原始代码（Line 246-248）
values.pstar = rootNode.getUpperBound();
values.qstar = rootNode.getLowerBound();
values.qprime= rootNode.getUpperBound();  // ⚠️ 错误！
```

### 正确的实现

查看K*的实现 [GradientDescentPfunc.java:486-489](src/main/java/edu/duke/cs/osprey/kstar/pfunc/GradientDescentPfunc.java#L486-L489)：

```java
// K*的正确实现
values.qstar = state.getLowerBound();
values.qprime = bigMath()
    .set(state.getUpperBound())
    .sub(state.getLowerBound())  // ✅ qprime = upper - lower
    .get();
```

### 影响

`qprime`在测试代码中被用于计算`pstar`：

```java
// TestMARKStarVsKStarPartitionFunction.java:260
result.pstar = pfuncResult.values.qprime.add(pfuncResult.values.qstar);
// 即: pstar = qprime + qstar = (upper - lower) + lower = upper
```

如果qprime错误地等于upper：
- `pstar = upper + lower` → **错误的上界**
- 导致partition function bounds计算错误

### 修复方案

```java
// 修复后的代码（Line 246-250）
values.pstar = rootNode.getUpperBound();
values.qstar = rootNode.getLowerBound();
// ✅ 修复：qprime是gap，不是upper bound本身
// 与GradientDescentPfunc.java:486-489保持一致
values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
```

---

## 测试结果

### 测试套件
- **文件**: [TestMARKStarVsKStarPartitionFunction.java](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java)
- **目的**: 对比MARK*和K* (GradientDescentPfunc)的partition function计算结果

### 测试系统
1. **3-flex system**: 3个flexible residues
2. **4-flex system**: 4个flexible residues
3. **5-flex system**: 5个flexible residues
4. **All systems**: 综合测试所有系统

### 修复前 vs 修复后

#### 修复前（有bug）
```
MARK* partition function: 2.6×10²⁶
K* partition function:    1.2×10⁴
差异:                     10²² 倍！

MARK* free energy:        -36.1 kcal/mol
K* free energy:           -5.6 kcal/mol
差异:                     -30.5 kcal/mol

Conformations evaluated:  0 (MARK*算法错误地认为已收敛)
```

#### 修复后（测试通过）
```bash
$ ./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMARKStarVsKStarPartitionFunction"

> Task :test
TestMARKStarVsKStarPartitionFunction > testSmallSystem3Flex() PASSED ✓
TestMARKStarVsKStarPartitionFunction > testSmallSystem4Flex() PASSED ✓
TestMARKStarVsKStarPartitionFunction > testSmallSystem5Flex() PASSED ✓
TestMARKStarVsKStarPartitionFunction > testAllSmallSystems()  PASSED ✓

BUILD SUCCESSFUL
```

所有测试现在都**通过**，说明：
- MARK*和K*的partition function计算结果一致（在epsilon容差内）
- Free energy bounds合理
- Algorithm正确地评估conformations并收敛

---

## 代码修改总结

### 文件修改
- **[MARKStarBound.java](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java)**
  - Line 248-250: 修复qprime赋值
  - Line 1071-1084: 修复negative correction处理

### Git Diff
```diff
@@ -245,7 +245,9 @@ public class MARKStarBound implements PartitionFunction.WithConfDB {
         debugPrint(String.format("Average Z reduction per minimization: %12.6e",averageReduction));
         values.pstar = rootNode.getUpperBound();
         values.qstar = rootNode.getLowerBound();
-        values.qprime= rootNode.getUpperBound();
+        // qprime should be the gap (upper - lower), not the upper bound itself
+        // This matches the implementation in GradientDescentPfunc.java:486-489
+        values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());

@@ -1068,12 +1069,19 @@ public class MARKStarBound implements PartitionFunction.WithConfDB {
             double tripleEnergy = minimizedTuple.energy;

             double lowerbound = minimizingEmat.getInternalEnergy(tuple);
-            if (tripleEnergy - lowerbound > 0) {
-                double correction = tripleEnergy - lowerbound;
+            double correction = tripleEnergy - lowerbound;
+            if (correction > 0) {
+                correctionMatrix.setHigherOrder(tuple, correction);
+            }
+            else {
+                // When triple energy is lower than pairwise bound, the matrix bound is too high.
+                // This can happen when pairwise minimization doesn't find the global minimum
+                // but triple minimization does. We should use the lower (more accurate) triple energy.
+                System.err.println("Negative correction for "+tuple.stringListing() +
+                                 " (correction=" + correction + "). Using triple energy as correction.");
+                // Store a negative correction to lower the bound
                 correctionMatrix.setHigherOrder(tuple, correction);
             }
-            else
-                System.err.println("Negative correction for "+tuple.stringListing());
```

---

## 与团队调查的关联

### 团队Wiki文档发现
从您提供的Wiki文档，团队已经识别了关键问题：

1. **"Negative correction"警告频繁出现** (7/22/25, 7/24-25/25)
2. **Root cause hypothesis**: Pairwise minimization没有找到global minimum
3. **Proposed fixes**:
   - Fix #1: 增加coordinate descent的步数
   - Fix #2: **使用triple能量更新lower bound** ← 我们实施的方案

### 我的贡献
1. **确认了团队的hypothesis**: Triple minimization确实找到了比pairwise更低的能量
2. **实施了团队提出的Fix #2**: 使用negative correction来修正过高的bounds
3. **发现并修复了额外的qprime bug**: 团队文档没有提到，但也是关键问题
4. **创建了自动化测试**: `TestMARKStarVsKStarPartitionFunction.java`验证修复效果

---

## 理论背景

### 为什么Pairwise Minimization不够精确？

**Coordinate Descent的局限性**:
1. 只在每个维度（dihedral angle）上单独优化
2. 可能陷入local minimum
3. 不同residues之间的coupling可能导致pairwise minimization错过更好的conformation

**Triple (Multi-body) Minimization的优势**:
1. 同时优化多个residues
2. 能够探索更大的构象空间
3. 更可能找到global minimum

**结果**: Triple minimization经常找到比pairwise matrix预测的更低的能量。

### MARK* Algorithm的设计意图

MARK*使用A*搜索 + bounds propagation：
1. **Energy Matrix** (pairwise): 提供快速的bounds估计
2. **Triple Minimization**: 在关键位置精炼bounds
3. **Correction Matrix**: 存储额外的higher-order corrections

**关键**: Correction matrix必须能够处理**负correction**，因为triple minimization可能发现matrix bounds过高。

---

## 未来工作建议

### 1. 增强Coordinate Descent精度（团队建议的Fix #1）

虽然Fix #2（使用negative correction）解决了问题，但也可以考虑：
- 增加minimization的步数（SimpleCCDMinimizer）
- 调整line search的tolerance
- 可能减少"Negative correction"事件的频率

**但注意**: 团队7/24/25的测试显示，简单地增加步数没有效果，因为算法通常在improvement变负时终止（不同的stopping condition）。

### 2. 测试其他MARK*变体

需要验证这些修复在其他MARK*实现中是否也需要：
- [MARKStarBoundAsync.java](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBoundAsync.java)
- [MARKStarBoundRigid.java](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBoundRigid.java)
- [GradientDescentMARKStarPfunc.java](src/main/java/edu/duke/cs/osprey/kstar/pfunc/GradientDescentMARKStarPfunc.java)

**注**: 这些文件可能有类似的negative correction处理逻辑。

### 3. 监控Negative Correction频率

添加统计来跟踪：
- Negative correction发生的频率
- Correction的平均大小
- 对最终partition function的影响

这可以帮助评估coordinate descent的质量。

### 4. 完整的系统测试

在实际的K*计算中测试（protein + ligand + complex）：
- 验证修复在完整的K*工作流中有效
- 确保没有引入新的问题
- 对比修复前后的K*结果

---

## 结论

通过分析团队之前的深入调查和我的测试，成功修复了MARK*算法中的两个关键bug：

1. **Negative Correction Bug** (MARKStarBound.java:1076)
   - **问题**: 忽略了triple minimization找到的更低能量
   - **影响**: Partition function被高估10²²倍
   - **修复**: 总是应用correction，包括负correction

2. **qprime Bug** (MARKStarBound.java:248)
   - **问题**: qprime应该是gap (upper-lower)，不是upper本身
   - **影响**: Partition function bounds计算错误
   - **修复**: 使用`upper.subtract(lower)`

**测试结果**: ✅ 所有测试通过

这些修复解决了MARK*算法的核心正确性问题，使其partition function计算与K*一致。

---

## 参考资料

### 内部文档
- [团队Wiki: MARKStar bugs](https://wiki.duke.edu/display/donaldlab/MARKStar+bugs)
- [TEST_RESULTS_SUMMARY.md](TEST_RESULTS_SUMMARY.md)
- [MARKSTAR_PROBLEM_SUMMARY.md](MARKSTAR_PROBLEM_SUMMARY.md)

### 代码文件
- [MARKStarBound.java](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java)
- [GradientDescentPfunc.java](src/main/java/edu/duke/cs/osprey/kstar/pfunc/GradientDescentPfunc.java)
- [TestMARKStarVsKStarPartitionFunction.java](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java)

### 测试命令
```bash
# 运行所有partition function对比测试
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMARKStarVsKStarPartitionFunction"

# 运行单个测试
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMARKStarVsKStarPartitionFunction.testSmallSystem3Flex"
```

---

**报告日期**: 2025-01-10
**修复实施**: Liyang Zhang & Claude
**测试验证**: 所有测试通过 ✅
