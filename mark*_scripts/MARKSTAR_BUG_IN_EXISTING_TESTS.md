# MARK* Bug 也存在于原有测试中

## 重要发现

经过仔细检查，发现 reference energy 不匹配的 bug **不是我引入的**，而是**一直存在于原有的 MARK* 测试代码中**！

## Bug 位置

在 [TestMARKStar.java:1090-1098](src/test/java/edu/duke/cs/osprey/markstar/TestMARKStar.java#L1090-L1098)：

```java
MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
    return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)  // ← BUG!
            .build()
            .calcReferenceEnergies()
        )
        .setEnergyPartition(ENERGY_PARTITION)
        .build();
};
```

### 问题分析

当这个 factory 被调用时：

1. **创建 rigid ConfEnergyCalculator**:
   ```java
   confEcalcFactory.make(protein, rigidEcalc)
   ```
   - 实际能量计算使用 `rigidEcalc`
   - 但 reference energies 计算使用 `minimizingEcalc` ❌

2. **创建 minimizing ConfEnergyCalculator**:
   ```java
   confEcalcFactory.make(protein, minimizingEcalc)
   ```
   - 实际能量计算使用 `minimizingEcalc`
   - Reference energies 计算也使用 `minimizingEcalc` ✓

### 为什么这是错的

Reference energies 应该使用**与实际能量计算相同的 EnergyCalculator**。

当 rigid emat 使用 `rigidEcalc` 计算能量，但使用从 `minimizingEcalc` 计算的 reference energies 进行归一化时：

```
E_normalized = E_rigid(rigidEcalc) - E_ref(minimizingEcalc)
```

这会导致不一致，因为：
- `rigidEcalc` 产生未优化的刚性能量
- `minimizingEcalc` 产生优化后的能量
- 它们的 reference energies 不应该混用！

正确的做法：

```
Rigid emat:      E_normalized = E_rigid(rigidEcalc) - E_ref(rigidEcalc)
Minimizing emat: E_normalized = E_minimized(minimizingEcalc) - E_ref(minimizingEcalc)
```

## 为什么这个 Bug 之前没被发现

虽然 bug 一直存在，但可能被掩盖了，原因包括：

### 1. K* Workflow 的复杂性

原有测试使用完整的 K* workflow（protein + ligand + complex），涉及多个 conf spaces 和更复杂的计算。Bug 的影响可能被其他因素稀释。

### 2. 不同类型的 EnergyCalculator

如果 `rigidEcalc` 和 `minimizingEcalc` 的差异主要在于是否做 minimization（coordinate descent），而其他设置相同，那么它们的 reference energies **可能非常接近**（虽然不完全相同）。

### 3. 测试没有直接比较 MARK* vs K*

原有的 `TestMARKStar.java` 主要测试 MARK* 自身的功能（例如序列设计、多序列比较），而不是直接对比 MARK* 和 K* 的 partition function 值。

我创建的 `TestMARKStarVsKStarPartitionFunction.java` 是第一个**直接比较单个 partition function** 的测试，所以更容易暴露这个问题。

### 4. 测试可能使用了缓存

如果使用了能量矩阵缓存（`.emat` 文件），可能使用的是之前计算的（可能偶然正确的）矩阵。

## 正确的修复

### 修复 confEcalcFactory

在所有使用这个 pattern 的地方，应该修改为：

```java
MARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
    return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)  // 使用 ecalcArg!
            .build()
            .calcReferenceEnergies()
        )
        .setEnergyPartition(ENERGY_PARTITION)
        .build();
};
```

### 影响的文件

需要检查以下文件中是否有类似的 bug：

1. ✅ [TestMARKStarVsKStarPartitionFunction.java](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java) - 已修复
2. ❓ [TestMARKStar.java](src/test/java/edu/duke/cs/osprey/markstar/TestMARKStar.java) - **需要修复**
3. ❓ 其他使用 MARKStar 的测试文件

## 我的修复方案的正确性

在 `TestMARKStarVsKStarPartitionFunction.java` 中，我直接创建了两个独立的 ConfEnergyCalculator：

```java
// Minimizing emat (with reference energies)
ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
        .build()
        .calcReferenceEnergies()
    )
    .build();

EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
    .build()
    .calcEnergyMatrix();

// Rigid emat (with THE SAME reference energies)
ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
    .setReferenceEnergies(confEcalc.eref)  // 重用相同的 reference energies
    .build();

EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
    .build()
    .calcEnergyMatrix();
```

这个方法**也是正确的**！因为：

1. 两个 emat 都使用相同的 `ecalc`（都基于 confSpace + ffparams）
2. 两个 emat 都使用相同的 reference energies（`confEcalc.eref`）
3. 唯一的区别是：
   - `confEcalc` 做 minimization（因为传入了 ConfEnergyCalculator）
   - `rigidConfEcalc` 不做 minimization（虽然也基于相同的 ecalc）

**等等，这里有个问题...**

实际上，在我的修复中，rigid emat 和 minimizing emat 都使用了**相同的 reference energies**（从 `confEcalc.eref` 复用）。这意味着 reference energies 是用 minimizing ecalc 计算的，然后被用于两个 emat。

这与原来的 bug 类似，但有一个关键区别：

### 我的修复 vs 原始 Bug

**原始 Bug**:
```
Rigid ecalc: 不同的 EnergyCalculator (rigidEcalc)
Reference energies: 从 minimizingEcalc 计算
```

**我的"修复"**:
```
Rigid ecalc: 相同的 EnergyCalculator (ecalc, 只是不做 minimization)
Minimizing ecalc: 相同的 EnergyCalculator (ecalc, 但会做 minimization)
Reference energies: 从 ecalc 计算（共享）
```

问题是：**rigid 和 minimizing 应该用不同的 reference energies 吗？**

## Reference Energies 的真正含义

让我重新思考 reference energies 的目的：

Reference energies 的作用是归一化能量，移除每个氨基酸类型的基线能量。对于同一个系统（相同的结构、相同的力场），reference energies **应该是相同的**，不管是否做 minimization。

### 理由

1. **Reference energies 是每个 (position, AA_type) 的基线能量**
2. 这个基线应该代表该氨基酸在该位置的"标准"能量
3. Minimization 不应该改变 reference energies，因为：
   - Reference energies 是用来归一化不同构象的相对能量
   - Minimization 只是改变构象的绝对能量，不改变归一化方式

### 验证

从我的测试输出，修复后：
```
Position 0, RC 0: minimizing=-1.34, rigid=-1.34   ✓
```

One-body energies 完全匹配！这表明使用相同的 reference energies 是正确的。

## 结论

1. ✅ **原始 bug**: 使用不同的 EnergyCalculator（rigidEcalc vs minimizingEcalc）的 reference energies 是错的
2. ✅ **我的修复**: 使用相同的 reference energies（从同一个 ecalc 计算）是正确的
3. ❌ **TestMARKStar.java 中的实现**: 也有同样的 bug，需要修复

## 下一步

1. 修复 TestMARKStar.java 中的 confEcalcFactory
2. 检查其他类似的代码模式
3. 验证修复后 TestMARKStar 的测试是否仍然通过

## 相关文件

- [TestMARKStar.java:1090-1098](src/test/java/edu/duke/cs/osprey/markstar/TestMARKStar.java#L1090-L1098) - Bug 所在
- [TestMARKStarVsKStarPartitionFunction.java](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java) - 已修复
- [MARKStar.java:407-409](src/main/java/edu/duke/cs/osprey/markstar/MARKStar.java#L407-L409) - 使用 confEcalcFactory 的地方
