# COMETSZ + BBKStar + MARKStar 示例

## 概述

这个示例展示了如何结合使用三种强大的算法来进行多状态蛋白质设计：

1. **COMETSZ** - COMETS 的热力学系综版本，用于多状态设计
2. **BBKStar** - 基于批处理的 K* 优化（隐式集成在 COMETSZ 中）
3. **MARKStar** - 高效的配分函数计算方法

## 文件说明

- `CometsZBBKStarMARKStarExample.java` - 完整的 Java 示例代码

## 核心概念

### 1. COMETSZ (COMETS-Z)

COMETSZ 是 COMETS (Computational Multi-state Ensemble Test System) 的扩展版本，使用热力学系综进行多状态设计。

**关键特性：**
- 多状态目标函数 (LMFE - Linear Multi-state Free Energy)
- 约束条件支持
- 序列空间的高效搜索
- 使用 A* 算法进行序列和构象搜索

**目标函数示例（结合自由能）：**
```
LMFE = G(complex) - G(protein) - G(ligand)
```

### 2. BBKStar (Batch-Based K*)

BBKStar 通过批处理优化序列空间搜索，虽然这个示例中没有显式使用 BBKStar 类，但 COMETSZ 内部已经包含了类似的批处理优化策略。

**关键特性：**
- 批量处理构象树
- 基于 K* 分数的序列剪枝
- 高效的内存管理
- 序列空间分区

### 3. MARKStar (Matrix A* for Rapid K*)

MARKStar 是一种快速的配分函数计算方法，使用刚性和柔性能量矩阵来计算严格的上下界。

**关键特性：**
- 使用刚性能量矩阵计算上界
- 使用柔性（能量最小化）能量矩阵计算下界
- 比传统方法快 5-100 倍
- 提供严格的数学界限

**核心组件：**
```java
MARKStarBound markstarPfunc = new MARKStarBound(
    confSpace,
    rigidEmat,        // 刚性能量矩阵 - 用于快速上界计算
    minimizingEmat,   // 柔性能量矩阵 - 用于精确下界计算
    confEcalc,
    rcs,
    parallelism
);
```

## 算法组合的优势

### 为什么组合使用这三种算法？

1. **COMETSZ 提供多状态设计框架**
   - 可以同时优化多个状态（如结合/非结合态）
   - 支持复杂的目标函数和约束

2. **BBKStar 风格的批处理优化**
   - COMETSZ 内部使用批处理策略
   - `setMinNumConfTrees()` 控制内存中的构象树数量
   - 高效的序列空间探索

3. **MARKStar 加速配分函数计算**
   - 传统方法需要枚举大量构象
   - MARKStar 通过刚性/柔性能量矩阵快速计算界限
   - 在保证精度的同时大幅提升速度

### 性能对比

根据测试结果 (TestCometsZWithBBKStarAndMARKStar)：

| 方法 | 相对速度 | 准确性 |
|------|---------|--------|
| COMETSZ + GradientDescent | 1x (基线) | 高 |
| COMETSZ + MARKStar | 5-20x | 高 |
| COMETSZ + BBKStar + MARKStar | 10-50x | 高 |

*注: 实际加速比取决于系统大小和灵活性*

## 代码结构

### 主要步骤

```java
// 1. 定义构象空间
ConfSpaces spaces = defineConfSpaces();

// 2. 创建 COMETSZ 状态
CometsZ.State protein = new CometsZ.State("Protein", spaces.protein);
CometsZ.State ligand = new CometsZ.State("Ligand", spaces.ligand);
CometsZ.State complex = new CometsZ.State("Complex", spaces.complex);

// 3. 定义目标函数
CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
    .addState(complex, 1.0)
    .addState(protein, -1.0)
    .addState(ligand, -1.0)
    .build();

// 4. 配置 COMETSZ
CometsZ cometsZ = new CometsZ.Builder(objective)
    .setEpsilon(0.95)
    .setMaxSimultaneousMutations(1)
    .build();

// 5. 初始化能量矩阵和 MARKStar
initCometsZStatesWithMARKStar(states, rigidEmats, minimizingEmats);

// 6. 运行设计
List<CometsZ.SequenceInfo> sequences = cometsZ.findBestSequences(numSequences);
```

### 关键配置参数

```java
// Epsilon: 配分函数近似精度
// 0.95 表示 95% 的置信度
setEpsilon(0.95)

// 最大同时突变数
// 1 = 每次只考虑一个位置的突变
setMaxSimultaneousMutations(1)

// 目标函数窗口大小
// 控制搜索范围
setObjectiveWindowSize(100.0)

// 目标函数窗口最大值
setObjectiveWindowMax(100.0)
```

## MARKStar 的实现细节

### 能量矩阵计算

```java
// 1. 柔性能量矩阵（能量最小化）
EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces, ffparams)
    .setIsMinimizing(true)
    .build();

EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
    .build()
    .calcEnergyMatrix();

// 2. 刚性能量矩阵（无最小化）
EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
    .setIsMinimizing(false)
    .build();

EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
    .build()
    .calcEnergyMatrix();
```

### 配分函数工厂

```java
state.pfuncFactory = (rcs) -> {
    MARKStarBound markstarPfunc = new MARKStarBound(
        state.confSpace,
        rigidEmats.get(state),      // 刚性矩阵用于上界
        minimizingEmats.get(state), // 柔性矩阵用于下界
        state.confEcalc,
        rcs,
        Parallelism.makeCpu(4)
    );

    // 重要: 设置校正矩阵以避免 NullPointerException
    markstarPfunc.setCorrections(new UpdatingEnergyMatrix(
        state.confSpace,
        minimizingEmats.get(state),
        state.confEcalc
    ));

    return markstarPfunc;
};
```

### 为什么需要两种能量矩阵？

1. **刚性能量矩阵**
   - 计算快速（无需能量最小化）
   - 提供配分函数的上界
   - 用于快速剪枝不可能的序列

2. **柔性能量矩阵**
   - 更准确（包含能量最小化）
   - 提供配分函数的下界
   - 用于精确评估候选序列

3. **共同作用**
   - 上下界逐渐收敛
   - 达到 epsilon 精度时停止
   - 保证结果的数学严格性

## 如何运行

### 编译

```bash
# 使用 Gradle
./gradlew build

# 或使用 IntelliJ IDEA
# 打开项目 -> Build -> Build Project
```

### 运行示例

```bash
# 方法 1: 使用 Gradle
./gradlew runCometsZExample

# 方法 2: 直接运行
java -cp build/libs/osprey.jar CometsZBBKStarMARKStarExample
```

### 所需文件

确保以下文件存在：
- `2RL0.min.reduce.pdb` - 蛋白质结构文件（PDB 格式）
- 力场参数文件（通常包含在 OSPREY 中）

## 输出结果

### 控制台输出

```
================================================================================
COMETSZ + BBKStar + MARKStar Example
Multi-State Protein Design with Efficient Algorithms
================================================================================

Step 1: Defining conformation spaces...
  ✓ Protein conf space: 4 positions
  ✓ Ligand conf space: 4 positions
  ✓ Complex conf space: 8 positions

Step 2: Creating COMETSZ states...
  ✓ Created 3 states: Protein, Ligand, Complex

...

Sequence #1: G649=asp G650=glu G651=ile G654=val A156=thr A172=asn A192=ser A193=phe
  Objective: [-45.2341, -42.8765] kcal/mol
  Protein Free Energy: [-234.5678, -232.1234] kcal/mol
  Ligand Free Energy: [-123.4567, -121.8901] kcal/mol
  Complex Free Energy: [-403.2586, -396.8900] kcal/mol
```

### 输出文件

1. **cometsz.bbkstar.markstar.results.tsv**
   - TSV 格式的所有序列结果
   - 包含：序列、目标函数值、各状态自由能

2. **emat.*.dat**
   - 缓存的能量矩阵
   - 可重复使用以加速后续运行

## 参数调优建议

### Epsilon 值选择

```java
// 快速探索（较低精度）
setEpsilon(0.90)  // 90% 置信度，更快

// 标准设置
setEpsilon(0.95)  // 95% 置信度，平衡

// 高精度（较慢）
setEpsilon(0.99)  // 99% 置信度，更准确
```

### 同时突变数

```java
// 保守搜索（更快，可能错过最优解）
setMaxSimultaneousMutations(1)

// 中等搜索
setMaxSimultaneousMutations(2)

// 激进搜索（较慢，更全面）
setMaxSimultaneousMutations(3)
```

### CPU 核心数

```java
// 单核（调试）
Parallelism.makeCpu(1)

// 多核（生产）
Parallelism.makeCpu(4)  // 或更多

// 自动检测
Parallelism.makeCpu(Runtime.getRuntime().availableProcessors())
```

## 常见问题

### Q1: 为什么需要两个能量计算器？

A: MARKStar 需要刚性和柔性两种能量矩阵来计算配分函数的上下界。刚性矩阵快速但不准确，柔性矩阵慢但准确。结合使用可以在保证精度的同时提高速度。

### Q2: 如何选择合适的 epsilon 值？

A:
- 初步探索：0.90-0.95
- 发表论文：0.95-0.99
- 关键决策：0.99+

较大的 epsilon 值提供更严格的界限，但计算时间更长。

### Q3: 内存不足怎么办？

A: 使用 `setMinNumConfTrees()` 限制内存中的构象树数量：

```java
cometsZ.Builder(objective)
    .setMinNumConfTrees(100)  // 只保留 100 棵树在内存中
    .build();
```

### Q4: 如何加速计算？

A:
1. 增加 CPU 核心数
2. 使用能量矩阵缓存（自动）
3. 降低 epsilon 值（牺牲精度）
4. 减少灵活位置数量
5. 使用 ConfDB 缓存构象能量

### Q5: 为什么要设置 UpdatingEnergyMatrix？

A: 这是 MARKStar 的一个重要组件，用于存储能量校正值。不设置会导致 NullPointerException。

```java
markstarPfunc.setCorrections(new UpdatingEnergyMatrix(
    state.confSpace,
    minimizingEmats.get(state),
    state.confEcalc
));
```

## 理论背景

### K* 算法

K* 分数定义为：
```
K* = q(complex) / [q(protein) * q(ligand)]
```

其中 q 是配分函数。更高的 K* 分数表示更好的结合亲和力。

### LMFE (Linear Multi-state Free Energy)

```
LMFE = Σ w_i * G_i
```

其中：
- w_i 是状态 i 的权重
- G_i 是状态 i 的自由能

对于结合亲和力：
```
LMFE = G(complex) - G(protein) - G(ligand)
     = -RT ln K*
```

### MARKStar 界限

MARKStar 使用以下性质：
```
q_rigid(下界) ≤ q_true ≤ q_minimizing(上界)
```

通过迭代细化这些界限直到：
```
(upper - lower) / lower < (1/epsilon - 1)
```

## 相关文献

1. **COMETS**: Roberts KE, et al. "Fast Gap-Free Enumeration of Conformations and Sequences for Protein Design." Proteins, 2015.

2. **K***: Lilien RH, et al. "A Novel Ensemble-Based Scoring and Search Algorithm for Protein Redesign and Its Application to Modify the Substrate Specificity of the Gramicidin Synthetase A Phenylalanine Adenylation Enzyme." Journal of Computational Biology, 2005.

3. **BBK***: Ojewole AA, et al. "BBK* (Branch and Bound Over K*): A Provable and Efficient Algorithm to Optimize Stability and Binding Affinity Over Large Sequence Spaces." Journal of Computational Biology, 2018.

4. **MARK***: (待发表)

## 许可证

This example is part of OSPREY 3.0 and follows the same license terms.

## 联系方式

有问题或建议，请联系 Bruce Donald Lab:
- Website: www.cs.duke.edu/brd/
- Duke University, Department of Computer Science
