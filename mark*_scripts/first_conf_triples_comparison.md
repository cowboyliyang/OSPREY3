# 第一个构象最小化后Triples对比分析

## 构象信息
- **构象**: `conf=17 13 12 8 27 7`
- **Pairwise Score**: -47.954216
- **Minimized Energy**: -42.961689
- **Gap**: 4.992527

## 版本对比

### 原版（Original）
- **Z Comparison**: `6.447882e+33` vs `1.118824e+29`
  - internalZ (6.447882e+33) > leafZ (1.118824e+29)
  - **决策**: 处理内部节点（因为internalZ更大）
- **处理**: 处理了 **1358个内部节点** 后收敛
- **Triples生成**: 42个
- **Correction应用**: **已应用**
  - Corrected Score: **-44.293191**
  - Correction值: -47.954216 - (-44.293191) = **-3.661025**

### 改进版（Improved）
- **Z Comparison**: `1.104669e+34` vs `1.258240e+35`
  - internalZ (1.104669e+34) < leafZ (1.258240e+35)
  - **决策**: 处理叶子节点（因为leafZ更大）
- **处理**: 处理了 **2个叶子节点**（继续处理更多构象）
- **Triples生成**: 42个（与原版相同）
- **Correction应用**: **未应用**
  - Corrected Score: **-47.954216**（等于pairwise score）
  - 没有应用triple correction

## Triples详细对比

两个版本生成的triples**完全相同**，都是42个，包括以下唯一的triple tuples：

1. `Res 0 RC 17 Res 2 RC 12 Res 5 RC 7` - correction: 0.038580
2. `Res 1 RC 13 Res 2 RC 12 Res 5 RC 7` - correction: 3.654442
3. `Res 1 RC 13 Res 2 RC 12 Res 4 RC 27` - correction: 3.617463
4. `Res 0 RC 17 Res 2 RC 12 Res 3 RC 8` - correction: 0.116694
5. `Res 1 RC 13 Res 2 RC 12 Res 3 RC 8` - correction: 3.661025
6. `Res 0 RC 17 Res 1 RC 13 Res 4 RC 27` - correction: 0.296918
7. `Res 0 RC 17 Res 1 RC 13 Res 3 RC 8` - correction: 0.477396
8. `Res 0 RC 17 Res 1 RC 13 Res 2 RC 12` - correction: 3.426351
9. `Res 2 RC 12 Res 3 RC 8 Res 5 RC 7` - correction: 0.118361
10. `Res 2 RC 12 Res 4 RC 27 Res 5 RC 7` - correction: 0.084430
11. `Res 1 RC 13 Res 3 RC 8 Res 4 RC 27` - correction: 1.135649
12. `Res 0 RC 17 Res 2 RC 12 Res 4 RC 27` - correction: 0.052475
13. `Res 2 RC 12 Res 3 RC 8 Res 4 RC 27` - correction: 0.647610
14. `Res 0 RC 17 Res 3 RC 8 Res 4 RC 27` - correction: 0.596120

**注意**: 某些triples在日志中出现了多次（重复计算），但唯一的triple tuples数量为14个。

## 关键差异分析

### 1. Z Comparison值的差异

**原版**:
- internalZ = 6.447882e+33
- leafZ = 1.118824e+29
- 比例: internalZ/leafZ ≈ 5.76e+4（内部节点Z值远大于叶子节点）

**改进版**:
- internalZ = 1.104669e+34
- leafZ = 1.258240e+35
- 比例: internalZ/leafZ ≈ 0.088（叶子节点Z值大于内部节点）

**差异原因**: 
- 改进版的leafZ值从1.118824e+29增加到1.258240e+35（增加了约6个数量级）
- 改进版的internalZ值从6.447882e+33增加到1.104669e+34（增加了约1.7倍）
- 这导致决策从"处理内部节点"变为"处理叶子节点"

### 2. Correction应用差异

**原版**: 应用了triple correction，corrected score = -44.293191
**改进版**: 未应用triple correction，corrected score = -47.954216

#### Correction应用的条件

根据源代码分析，correction在两个地方可能被应用：

1. **在`populateQueues`中（处理内部节点时）**:
   ```java
   if(!node.isMinimized() && node.getConfLowerBound() < confCorrection
       && node.getConfLowerBound() - confCorrection > 1e-5) {
       // 应用correction
   }
   ```
   条件：`pairwiseLowerBound < confCorrection` 且差值 > 1e-5

2. **在`processFullConfNode`中（处理叶子节点时）**:
   ```java
   if(node.getConfLowerBound() < confCorrection || node.gscore < confCorrection) {
       // 应用correction
   }
   ```
   条件：`pairwiseLowerBound < confCorrection` **或** `gscore < confCorrection`

#### 为什么原版应用了correction而改进版没有？

**关键发现**:
- Pairwise score = -47.954216
- Corrected score = -44.293191（原版）
- Correction值 = -44.293191 - (-47.954216) = **+3.661025**（正值，表示correction使能量变差）

**问题**: 如果correction值（-44.293191）比pairwise值（-47.954216）更差（更高），那么条件`pairwiseLowerBound < confCorrection`应该是**false**（因为-47.954216 < -44.293191是false）。

**原版为什么应用了correction？**

可能的原因：
1. **在`populateQueues`中的特殊逻辑**：原版在处理内部节点时，可能使用了不同的correction计算方式（第616行：`confCorrection = Math.min(correctgscore, node.rigidScore) + hscore`），其中`hscore = node.getConfLowerBound() - node.gscore`。这个计算可能产生了不同的confCorrection值。

2. **Correction值的动态更新**：在处理1358个内部节点的过程中，correctionMatrix可能被更新，导致后续节点的correction值发生变化。

3. **Bound更新后的影响**：在`updateBound()`调用后（第502行），节点的bound可能被更新，使得correction应用条件在后续循环中满足。

**改进版为什么没有应用correction？**

可能的原因：
1. **处理路径不同**：改进版直接处理叶子节点，跳过了`populateQueues`中的correction应用逻辑（第617-633行）。

2. **条件不满足**：在`processFullConfNode`中，条件`pairwiseLowerBound < confCorrection`不满足（因为-47.954216 < -44.293191是false），且`gscore < confCorrection`也不满足。

3. **Correction值未更新**：由于没有处理内部节点，correctionMatrix可能没有机会被更新，导致correction值保持为初始值（可能等于或接近pairwise值）。

#### 根本原因

**核心差异在于处理路径和correction计算方式**：

1. **原版（处理内部节点）**：
   - 在`populateQueues`中（第617-633行），correction的计算方式：
     ```java
     double correctgscore = correctionMatrix.confE(node.assignments);
     double hscore = node.getConfLowerBound() - node.gscore;
     double confCorrection = Math.min(correctgscore, node.rigidScore) + hscore;
     ```
   - 应用条件：`node.getConfLowerBound() < confCorrection && node.getConfLowerBound() - confCorrection > 1e-5`
   - **关键**：`confCorrection`不是简单的`correctgscore`，而是经过调整的值（`Math.min(correctgscore, node.rigidScore) + hscore`）
   - 这个调整可能使得条件更容易满足

2. **改进版（处理叶子节点）**：
   - 在`processFullConfNode`中（第956行），correction的计算方式：
     ```java
     double confCorrection = correctionMatrix.confE(node.assignments);
     ```
   - 应用条件：`node.getConfLowerBound() < confCorrection || node.gscore < confCorrection`
   - **关键**：直接使用`correctionMatrix.confE()`的值，没有调整
   - 由于`confCorrection = -44.293191`比`pairwiseLowerBound = -47.954216`更高（更差），条件`-47.954216 < -44.293191`为**false**，所以correction未应用

**为什么原版应用了correction？**

可能的原因：
1. **`hscore`的调整**：在`populateQueues`中，`confCorrection = Math.min(correctgscore, node.rigidScore) + hscore`，其中`hscore = node.getConfLowerBound() - node.gscore`。这个调整可能使得`confCorrection`值发生变化，满足应用条件。

2. **节点状态不同**：在处理内部节点时，节点的`gscore`、`rigidScore`等状态可能与处理叶子节点时不同，导致条件判断结果不同。

3. **Bound更新时机**：在处理1358个内部节点的过程中，`updateBound()`被多次调用，节点的bound可能被更新，使得后续节点的correction应用条件更容易满足。

**关键点**：correction的应用不仅取决于correction值本身，还取决于：
1. 处理路径（内部节点 vs 叶子节点）
2. Correction的计算方式（是否经过`hscore`调整）
3. 节点的当前状态（bound、gscore、rigidScore等）
4. CorrectionMatrix的更新时机

### 3. 处理路径差异

**原版**: 
- 处理1358个内部节点后收敛
- 不再继续最小化

**改进版**:
- 处理2个叶子节点
- 继续探索更多构象

## Z Comparison的含义

根据源代码（MARKStarBound.java:506），Z Comparison比较的是：
- **internalZ**: 内部节点的分区函数（partition function）值
- **leafZ**: 叶子节点的分区函数值

**决策逻辑**:
```java
if(MathTools.isLessThan(internalZ, leafZ)) {
    // 处理叶子节点
} else {
    // 处理内部节点
}
```

**物理意义**:
- Z值代表该类型节点的"贡献度"或"重要性"
- 较大的Z值意味着该类型节点对整体bound的贡献更大
- 算法优先处理Z值较大的节点类型，以更快地收紧bound

## 总结

1. **Triples完全相同**: 两个版本生成的42个triples完全一致，说明triple生成逻辑没有变化。

2. **主要差异在Z值计算**: 改进版的leafZ值显著增加，导致决策从"处理内部节点"变为"处理叶子节点"。

3. **Correction应用不同**: 原版应用了correction，改进版未应用，这可能与处理路径的改变有关。

4. **处理策略不同**: 原版收敛于内部节点处理，改进版继续探索叶子节点，可能导致不同的搜索路径和最终结果。
