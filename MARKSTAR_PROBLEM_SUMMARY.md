# MARK* 问题总结

## 问题根源

通过添加调试输出，我们发现了问题的确切根源：

### 调试输出显示

```
=== MARKStarNode.makeRoot DEBUG ===
  gscore: 0.0
  rigidScore: 0.0
  hScore: -5.297405421040506
  negatedHScore: -23940.324517557674  ← 巨大的负数！
  confLowerBound: -5.297405421040506
  confUpperBound: 23940.324517557674  ← 巨大的正数！
  numConformations: 486  ← 不是24！
```

### 问题分析

#### 1. negatedHScore = -23940

**原因**：在`MARKStarBound`的构造函数中（Line 359）：

```java
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)
```

`NegatedEnergyMatrix`把rigid能量矩阵的所有值都取反：
- Rigid能量矩阵包含很多负能量（有利的相互作用）
- 取反后变成正能量
- H-scorer把这些正能量加起来，得到一个巨大的正数
- 但由于某种原因（optimizer direction或其他），最终返回-23940

#### 2. confUpperBound = 23940

从`MARKStarNode.makeRoot` Line 405：

```java
double confUpperBound = rigidgScorer.calc(confIndex,rcs) - negatedHScorer.calc(confIndex, rcs);
                      = 0 - (-23940)
                      = +23940 kcal/mol
```

这是一个**巨大的正能量**，完全不合理！

#### 3. 最终导致的Partition Function错误

当展开树到leaf节点时，某个conf的能量是-28.8 kcal/mol：

```
Boltzmann weight = e^(-(-28.8)/0.592) = e^(48.6) ≈ 1.3×10²¹
```

这就是我们看到的1.41×10²⁶的来源！

#### 4. NumConformations = 486

- 测试输出显示只有24个实际构象（K*评估了20-24个）
- 但rootNode认为有486个
- 这可能是因为：
  - 486 = 所有可能的rotamer组合（在pruning之前）
  - 24 = 实际可行的构象（pruning之后）

这本身可能不是bug，而是预期行为（rootNode计算所有可能的组合）。

## 为什么成功的测试能工作？

成功的MARK*测试（如`TestMARKStar.java`）使用的是完整的`MARKStar`框架：

```java
MARKStar markstar = new MARKStar(
    confSpaces.protein,
    confSpaces.ligand,
    confSpaces.complex,
    rigidEcalc,
    minimizingEcalc,
    confEcalcFactory,
    settings
);
```

而我们的测试直接使用`MARKStarBound`：

```java
MARKStarBound markstar = new MARKStarBound(
    confSpace,
    rigidEmat,
    minimizingEmat,
    confEcalc,
    rcs,
    Parallelism.makeCpu(4)
);
```

**关键差异**：
1. `MARKStar`用于K*计算（protein-ligand binding），需要三个conf space
2. 我们的测试只有一个conf space（complex only）
3. MARK*的设计可能假设了K*的上下文，其中negatedHScorer的用法是正确的

## 可能的解决方案

### 方案1：修改测试使用完整的MARKStar框架

创建protein、ligand和complex三个conf space，使用`MARKStar`而不是`MARKStarBound`。

**优点**：
- 使用MARK*的预期用法
- 可能直接工作

**缺点**：
- 需要重新设计测试
- 更复杂

### 方案2：修复MARKStarBound的构造函数

问题可能在于`negatedHScorer`的创建。我们可以尝试：

1. **不使用NegatedEnergyMatrix**：
   ```java
   // 从：
   new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)

   // 改为：
   new TraditionalPairwiseHScorer(rigidEmat, rcs)  // 或者其他scorer
   ```

2. **理解negatedHScorer的真正用途**：
   - 它用于计算confUpperBound
   - confUpperBound应该是能量的上界（最不稳定的估计）
   - 但当前实现导致confUpperBound = +23940，这显然是错的

3. **检查MARKStar类如何使用MARKStarBound**：
   - MARKStar可能在创建MARKStarBound之前做了一些预处理
   - 或者使用了不同的参数

### 方案3：参考K* GradientDescent的做法

K* GradientDescent不使用negatedHScorer，它的bounds计算逻辑更简单。我们可以：

1. 创建一个简化版的MARKStarBound用于测试
2. 或者完全不测试MARKStarBound，只测试完整的MARKStar

## 下一步建议

**推荐方案**：修改测试使用完整的MARKStar框架

**原因**：
1. MARKStarBound的设计假设了K*上下文
2. 直接使用MARKStarBound可能不是预期用法
3. 成功的测试都使用完整的MARKStar框架

**实施步骤**：
1. 修改测试创建protein、ligand和complex三个conf space
2. 使用`MARKStar`类而不是`MARKStarBound`
3. 参考`TestMARKStar.java`中的成功测试

## 或者：简单修复尝试

如果只想快速测试，可以尝试在MARKStarBound构造函数中：

```java
// Line 359，从：
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)

// 改为（与line 354一致）：
hscorerFactory.make(rigidEmat)
```

这样negatedHScorer就会返回正常的h-score，而不是-23940。

但这可能破坏MARK*的正确性，因为negatedHScorer是有特定目的的。
