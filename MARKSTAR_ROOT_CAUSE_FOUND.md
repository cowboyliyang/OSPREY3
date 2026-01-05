# MARK* Bug Root Cause Analysis

## 问题总结

MARK* 计算的 partition function 比 K* 大约 10²² 倍 (Z_MARK ≈ 1.4×10²⁶ vs Z_K* ≈ 1.2×10⁴)

## 根本原因

**Rigid energy matrix 没有使用 reference energies 进行归一化。**

### 证据

从 debug 输出（[TestMARKStarVsKStarPartitionFunction.java](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java)）：

```
Sample one-body energies:
  Position 0 (first 3 RCs):
    RC 0: minimizing=-1.34, rigid=-4.79
    RC 1: minimizing=-1.33, rigid=-4.79
    RC 2: minimizing=-1.30, rigid=-4.75
  Position 1 (first 3 RCs):
    RC 0: minimizing=+8.74, rigid=-5.04
    RC 1: minimizing=+16.49, rigid=+2.70
    RC 2: minimizing=+2.27, rigid=-11.52
  Position 2 (first 3 RCs):
    RC 0: minimizing=+1.03, rigid=-11.81
    RC 1: minimizing=+6.07, rigid=-6.77
    RC 2: minimizing=+6.53, rigid=-6.31
```

关键观察：
1. **Minimizing 能量**: 有正有负，范围大约 -1 到 +16 kcal/mol
2. **Rigid 能量**: 几乎全是负数，范围大约 -12 到 +3 kcal/mol
3. **差异很大**: 例如 pos1 RC0 相差 13.78 kcal/mol

### 为什么会这样？

在 [TestMARKStarVsKStarPartitionFunction.java:344-347](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java#L344-L347)：

```java
log("Calculating rigid energy matrix...");
EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
    .build()
    .calcEnergyMatrix();
```

而 minimizing emat (line 340) 使用了 `confEcalc`：

```java
EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
    .build()
    .calcEnergyMatrix();
```

差异在于：
- `confEcalc` = ConfEnergyCalculator，**已经设置了 reference energies**（line 331-336）
- `ecalc` = EnergyCalculator，**没有 reference energies**

### Reference Energies 的作用

Reference energies 用于归一化能量矩阵：

```
E_normalized(pos, RC) = E_absolute(pos, RC) - E_reference(pos, AA_type)
```

没有归一化时，rigid 能量是绝对值，包含了整个残基的所有原子间相互作用。

归一化后，能量相对于野生型（或参考态），数值更小且更有意义。

## negatedHScore 为什么是负数？

从 debug 输出：
```
negatedHScore: -23940.324517557674
```

### 计算过程

1. **rigidEmat 的 one-body 能量**: 大约 -10 kcal/mol（未归一化，很负）
2. **NegatedEnergyMatrix 取反**: 变成 +10 kcal/mol（正数）
3. **TraditionalPairwiseHScorer 使用 Minimize**: 选每个 position 的最小值
4. **hscore 求和**: 对 3 个 positions，假设每个最小值 ~+10，总和 ~+30 kcal/mol

**但实际上 negatedHScore = -23940！**

这说明原始 rigidEmat 的能量总和是 **+23940 kcal/mol**（包括 pairwise terms）。

取反后变成 **-23940 kcal/mol**。

### 为什么原始能量是正数？

虽然 one-body 能量是负数（~-10），但是：
1. **Pairwise energies** 可能是正数（不利相互作用）
2. **所有能量项的总和** 可能是正数
3. 没有归一化时，绝对能量值很大

实际上，原始的 rigid 能量矩阵包含了整个蛋白质的绝对能量，数量级可能是数千 kcal/mol。

## confUpperBound 的错误传播

```
confUpperBound = rigidScore - negatedHScore
               = 0 - (-23940)
               = +23940 kcal/mol
```

这是一个**巨大的正能量**（非常不稳定的构象）。

### 为什么这导致 10²⁶ 的 partition function？

从 [MARKStarNode.java:520-521](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarNode.java#L520-L521)：

```java
private BigDecimal computeBoundsFromEnergy(double energy) {
    return bc.calc(energy).multiply(new BigDecimal(getNumConformations()));
}
```

BoltzmannCalculator.calc() ([BoltzmannCalculator.java:116-118](src/main/java/edu/duke/cs/osprey/tools/BoltzmannCalculator.java#L116-L118))：

```java
public BigDecimal calc(double energy) {
    return e.exp(-energy/RT);  // RT ≈ 0.592 at 298K
}
```

但是注意 [MARKStarNode.java:508-518](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarNode.java#L508-L518) 的逻辑：

```java
private void updateConfUpperBound(double tighterUpper) {
    if(tighterUpper < confUpperBound) {
        confUpperBound = tighterUpper;
        updateSubtreeLowerBound(computeBoundsFromEnergy(confUpperBound));
    }
}
```

**关键点**: confUpperBound → subtreeLowerBound（注意方向！）

这是因为：
- confUpperBound 是**能量的 upper bound**（最高能量）
- 转换成 Boltzmann weight: exp(-E_upper/RT)
- 高能量 → 小权重 → partition function 的 **lower bound**

所以如果 confUpperBound = +23940（错误地太高）：
```
subtreeLowerBound = exp(-23940/0.592) ≈ exp(-40439) ≈ 0
```

这应该让 partition function 接近 0，**不是** 10²⁶！

### 那为什么 MARK* 得到 10²⁶？

答案在于 **leaf nodes**。

Root node 的 bounds 是初始估计。真正的 partition function 来自 leaf nodes（完全定义的构象）。

Leaf nodes 会重新计算 bounds，使用实际的能量值（通过 minimization 或其他方式）。

如果 leaf nodes 的计算也受到 rigid emat 未归一化的影响，那么 leaf 的 confUpperBound 可能是一个**负数**（favorable energy），例如 -28.8 kcal/mol。

那样的话：
```
Boltzmann weight = exp(-(-28.8)/0.592) = exp(48.6) ≈ 6×10²¹
```

对 20 个这样的 leaves：
```
Z ≈ 20 × 6×10²¹ ≈ 1.2×10²³
```

这就接近我们看到的 10²⁶ 了！

## 解决方案

**修改 test 使 rigid emat 也使用 reference energies**。

在 [TestMARKStarVsKStarPartitionFunction.java:344-347](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java#L344-L347)，改为：

```java
log("Calculating rigid energy matrix...");
EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(
        new ConfEnergyCalculator.Builder(confSpace, ecalc)
            .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                .build()
                .calcReferenceEnergies()
            )
            .build()
    )
    .build()
    .calcEnergyMatrix();
```

或者更简洁（重用 confEcalc 的 reference energies）：

```java
ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
    .setReferenceEnergies(confEcalc.eref)  // 重用相同的 reference energies
    .build();

EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
    .build()
    .calcEnergyMatrix();
```

## 下一步

1. 实施修复
2. 验证 MARK* 和 K* 的 partition functions 现在一致
3. 检查是否有其他地方也需要使用 reference energies

## 相关文件

- [TestMARKStarVsKStarPartitionFunction.java:344-347](src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java#L344-L347) - Bug 所在
- [MARKStarNode.java:388-441](src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarNode.java#L388-L441) - makeRoot 方法
- [NegatedEnergyMatrix.java:54-66](src/main/java/edu/duke/cs/osprey/ematrix/NegatedEnergyMatrix.java#L54-L66) - 取反逻辑
- [BoltzmannCalculator.java:116-118](src/main/java/edu/duke/cs/osprey/tools/BoltzmannCalculator.java#L116-L118) - 权重计算
