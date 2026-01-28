# MARK* vs K* Partition Function 测试结果详细分析

## 测试概述

测试文件：`src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java`

测试系统：3个flexible residues (G649, G650, G651)，epsilon=0.10

运行时间：2024年执行

---

## 1. K* (GradientDescent) 执行结果

### 1.1 基本统计
- **Q* (lower bound)**: 1.22e+04
- **P* (upper bound)**: 1.22e+04
- **Gap (P* - Q*)**: 1.51e+02
- **Free Energy**: [-5.587, -5.580]
- **构象评估**: 24个
- **最小化构象**: 20个
- **执行时间**: 194 ms

### 1.2 执行过程分析

从log输出可以看到K*的执行轨迹：

```
Conf#  Assignments   Score       Energy      Lower Bound  Upper Bound  Delta
1      [8,3,8]      -5.085246   -4.772786   3.495281     Infinity     1.000
2      [8,5,8]      -5.261103   -4.891207   3.841756     Infinity     1.000
...
20     [0,2,8]      -2.647334   -2.162909   4.083185     4.097184     0.032
```

**关键观察**：

1. **能量关系正确**：
   - Score ≤ Energy (所有20个构象)
   - Score是A*给出的下界（最小化前）
   - Energy是实际最小化后的能量
   - 平均差异：~0.4 kcal/mol

2. **Bounds收敛**：
   - Delta从1.0逐渐降到0.032
   - Upper bound从Infinity收敛到4.097
   - Lower bound稳步增长到4.083
   - 最终gap: log10(upper) - log10(lower) ≈ 0.014

3. **分数分布**：
   - 最低能量: -4.891207 (conf#2: [8,5,8])
   - 最高能量: -1.590522 (conf#19: [2,4,8])
   - 范围: ~3.3 kcal/mol

4. **算法行为**：
   - 前7个构象：仅能量计算（scores=0）
   - 第8个构象后：开始打分（scores增长到280）
   - 这是gradient descent的典型行为：先minimize，后score

---

## 2. MARK* 执行结果

### 2.1 基本统计
- **Q* (lower bound)**: 1.406e+26
- **P* (upper bound)**: 1.406e+26
- **Gap (P* - Q*)**: 1.97e+16
- **Free Energy**: -35.706
- **构象评估**: 0个
- **节点展开**: 18个
- **执行时间**: 9 ms

### 2.2 异常现象

```
Running until leaf is found...
MARK* expanded: 18, queued: 0, scored/sec: 27011
Found a leaf!
```

**问题识别**：

1. **过早终止**：
   - MARK*只运行了`runUntilNonZero()`阶段
   - 找到一个leaf后立即停止
   - 没有进入`tightenBoundInPhases()`主循环

2. **Partition function异常**：
   - 值比K*大10^22倍
   - Q* = P* (gap为0意味着没有uncertainty)
   - 但实际gap=1.97e+16显示这不正确

3. **没有最小化**：
   - 0个构象被评估
   - 只做了树的初始展开（18个节点）
   - 没有实际的能量计算

---

## 3. 问题根源分析

### 3.1 MARK*算法流程

从代码分析，MARK*的`compute()`方法流程：

```java
public void compute(int maxNumConfs) {
    int previousConfCount = workDone();

    if(!nonZeroLower) {
        runUntilNonZero();      // 阶段1：找到一个leaf
        updateBound();
    }

    while (epsilonBound > targetEpsilon &&     // 阶段2：收紧bounds
           workDone()-previousConfCount < maxNumConfs
           && isStable(stabilityThreshold)) {
        tightenBoundInPhases();
        ...
    }

    values.pstar = rootNode.getUpperBound();
    values.qstar = rootNode.getLowerBound();
    values.qprime = rootNode.getUpperBound();  // ⚠️ qprime = pstar
}
```

### 3.2 关键问题

**问题1：While循环没有执行**

循环条件：
```java
epsilonBound > targetEpsilon  // 需要：epsilon > 0.10
```

可能原因：
- `epsilonBound`在`updateBound()`后已经 ≤ 0.10
- 或者初始值就是一个很小的数

**问题2：Partition Function Values定义不同**

K*:
```java
values.qstar = state.getLowerBound();              // 真正的下界
values.qprime = state.getUpperBound() - state.getLowerBound();  // gap
// pstar = qstar + qprime
```

MARK*:
```java
values.qstar = rootNode.getLowerBound();
values.pstar = rootNode.getUpperBound();
values.qprime = rootNode.getUpperBound();  // ⚠️ 这里赋值不对
```

MARK*的`qprime`应该是gap，但代码中赋值为upper bound！

**问题3：MARK*使用不同的bounds语义**

MARK*的文档注释：
```java
// Overwrite the computeUpperBound and computeLowerBound methods
public static class Values extends PartitionFunction.Values {
    @Override
    public BigDecimal calcUpperBound() {
        return pstar;  // MARK*中pstar是upper bound
    }

    @Override
    public BigDecimal calcLowerBound() {
        return qstar;  // qstar是lower bound
    }
}
```

这与K*的语义不同！

### 3.3 能量计算的差异

**K*的能量计算**：
- 使用Boltzmann weights: e^(-E/RT)
- Lower bound: Σ e^(-E_minimized/RT)
- Upper bound: Σ e^(-E_score/RT) + unscored * min_weight

**MARK*的能量计算**：
- 使用树结构的bounds
- 节点的upper/lower bounds传播到root
- 没有明确的Boltzmann weight计算（或者在节点级别）

---

## 4. 能量边界关系验证

### 4.1 预期的边界关系

对于单个构象：
```
RigidScore ≤ MinimizingScore ≤ MinimizedEnergy ≤ RigidEnergy
    ↓            ↓                  ↓                 ↓
 最紧下界    A*下界            实际能量          最松上界
```

对于partition function：
```
Q* (lower) ≤ True Z ≤ P* (upper)
```

### 4.2 K*的验证

从20个构象的数据：
```
Score        Energy       Relationship
-5.261103   -4.891207    ✅ Score < Energy
-5.085246   -4.772786    ✅ Score < Energy
-3.597335   -3.131962    ✅ Score < Energy
...
-2.647334   -2.162909    ✅ Score < Energy
```

**所有构象满足 Score ≤ Energy ✅**

### 4.3 MARK*的验证

无法验证，因为：
- 没有构象被实际评估
- 没有score/energy配对数据
- 只有树节点的bounds

---

## 5. 测试系统特征

### 5.1 系统规模

- **Flexible residues**: 3个
- **总构象数**: 通过RCs计算
- **位置和残基类型**:
  - G649: PHE (wild-type) + ALA, VAL (mutations)
  - G650: ASP (wild-type)
  - G651: GLU (wild-type)

### 5.2 能量矩阵

- **Minimizing energy matrix**: 213个entries (计算时间: 318 ms)
- **Rigid energy matrix**: 213个entries (计算时间: 76 ms)
- **Reference energies**: 24个residue confs (计算时间: 76 ms)

### 5.3 构象空间复杂度

根据K*评估了24个构象：
- 低能量构象集中在: [8,*,8], [0,*,8] 区域
- [8,5,8] 是最低能量构象 (-4.891)
- 能量分布跨度: ~3.3 kcal/mol

---

## 6. 性能比较

| 指标 | K* | MARK* | 比率 |
|------|-------|---------|------|
| 时间 (ms) | 194 | 9 | 21.6x |
| 构象评估 | 24 | 0 | N/A |
| 节点展开 | N/A | 18 | N/A |
| 内存 | 643 MiB | 643 MiB | 1.0x |

**注意**: MARK*的速度优势是假象，因为它根本没完成计算。

---

## 7. 结论与建议

### 7.1 主要发现

1. **K*工作正常**：
   - 正确计算partition function
   - 能量边界关系一致
   - Bounds正常收敛

2. **MARK*存在严重问题**：
   - 过早终止（只完成初始化）
   - Partition function值不可信
   - 没有实际的最小化计算

3. **根本原因**：
   - `epsilonBound`可能初始就满足条件
   - 或者`runUntilNonZero()`后epsilon已达标
   - 导致主循环不执行

### 7.2 需要修复的问题

**优先级1（Critical）**：
```java
// 修复MARK*的qprime赋值
values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
// 而不是
values.qprime = rootNode.getUpperBound();
```

**优先级2（High）**：
- 调查为什么`epsilonBound ≤ targetEpsilon`在初始阶段就为真
- 确保MARK*真正进入`tightenBoundInPhases()`循环
- 添加调试输出显示epsilon的值

**优先级3（Medium）**：
- 在`runUntilNonZero()`中实际计算一些构象能量
- 验证节点bounds的正确性
- 比较rigid和minimizing能量矩阵的使用

### 7.3 测试改进建议

1. **添加调试模式**：
```java
markstar.debug = true;  // 启用详细输出
```

2. **强制进入主循环**：
```java
// 测试时使用更严格的epsilon
markstar.init(0.01);  // 而不是0.10
```

3. **捕获中间结果**：
```java
// 在runUntilNonZero()后检查bounds
log("After initial phase: epsilon=%f, lower=%e, upper=%e",
    markstar.epsilonBound,
    markstar.rootNode.getLowerBound(),
    markstar.rootNode.getUpperBound());
```

### 7.4 下一步行动

1. ✅ **测试框架已完成** - 可以用于比较
2. ❌ **MARK*需要修复** - 当前结果不可信
3. ⏳ **待验证** - 修复后重新运行测试
4. ⏳ **待扩展** - 在4和5个flexible residues上测试

### 7.5 预期正确结果

修复后，MARK*应该：
- **Partition function值** 与K*接近（同一数量级）
- **Gap更小** 比K*的gap小（更紧的bounds）
- **评估更多构象** 实际进行能量最小化
- **Free energy接近** 在K*的范围内或更精确

---

## 8. 测试代码质量评估

### 8.1 优点

✅ **结构清晰**：
- 分离的测试方法（3、4、5 flex）
- 清晰的结果数据结构
- 完整的比较输出

✅ **错误处理**：
- 正确使用TaskExecutor.ContextGroup
- 适当的try-with-resources
- 编译通过，运行稳定

✅ **可扩展性**：
- 易于添加新的测试系统
- ConformationBounds类可用于未来验证
- 输出格式易于解析

### 8.2 改进空间

1. **添加断言**：当前测试只输出，没有断言来自动验证结果
2. **保存结果**：可以将结果写入CSV用于后续分析
3. **绘图功能**：可视化bounds收敛过程
4. **错误检测**：自动检测异常的partition function值

---

## 附录A：Score计算的详细解析

### Score是下界吗？

**答案：是的！Score确实是能量的下界。**

### Score的计算方式

从代码分析，我们发现Score是通过A*算法的f-score计算的：

```java
f(node) = g(node) + h(node)
```

其中：

#### 1. **g-score（已定义部分的能量）**

定义在`PairwiseGScorer.java`:

```java
public double calc(ConfIndex confIndex, RCs rcs) {
    double gscore = emat.getConstTerm();  // 常数项

    // 单体能量（one-body）
    for (int i=0; i<confIndex.numDefined; i++) {
        int pos1 = confIndex.definedPos[i];
        int rc1 = confIndex.definedRCs[i];
        gscore += emat.getOneBody(pos1, rc1);
    }

    // 成对能量（pairwise）
    for (int i=0; i<confIndex.numDefined; i++) {
        int pos1 = confIndex.definedPos[i];
        int rc1 = confIndex.definedRCs[i];
        for (int j=0; j<i; j++) {
            int pos2 = confIndex.definedPos[j];
            int rc2 = confIndex.definedRCs[j];
            gscore += emat.getPairwise(pos1, rc1, pos2, rc2);
        }
    }

    return gscore;
}
```

**关键点**：g-score是已经确定的residue conformations之间的精确能量（来自能量矩阵）。

#### 2. **h-score（未定义部分的能量下界）**

定义在`TraditionalPairwiseHScorer.java`:

```java
public double calc(ConfIndex confIndex, RCs rcs) {
    double hscore = 0;

    // 对每个未定义的位置...
    for (int i=0; i<confIndex.numUndefined; i++) {
        int pos = confIndex.undefinedPos[i];

        // 在这个位置优化所有可能的RC
        double optRCEnergy = optimizer.initDouble();  // 初始化为+∞（minimize）或-∞（maximize）
        for (int j=0; j<rcs.get(pos).length; j++) {
            // 对于每个RC，计算它与：
            // 1. 已定义位置的最优pair能量
            // 2. 其他未定义位置的最优pair能量
            optRCEnergy = optimizer.opt(optRCEnergy, cachedEnergies[pos][j]);
        }

        hscore += optRCEnergy;
    }

    return hscore;
}
```

**关键点**：h-score是对未定义部分的**乐观估计**（optimistic estimate）：
- 对于minimize：选择每个位置的**最低**可能能量
- 这忽略了不同位置之间未来选择的相互影响
- 因此h-score ≤ 真实剩余能量

#### 3. **为什么f = g + h 是下界？**

**数学证明**：

设完整构象的能量为 E_total = E_defined + E_undefined + E_interactions

其中：
- E_defined：已定义residues之间的能量 = g-score（精确）
- E_undefined：未定义residues的内部能量
- E_interactions：未定义residues之间的相互作用

h-score的计算：
```
h-score = Σ min{E_i + min_j{pair(i,j)}} for all undefined positions i
```

由于h-score独立优化每个位置，忽略了未定义位置之间的耦合：
```
h-score ≤ E_undefined + E_interactions
```

因此：
```
f-score = g-score + h-score
        = E_defined + h-score
        ≤ E_defined + E_undefined + E_interactions
        = E_total
```

**∴ f-score ≤ E_total （真实能量）✅**

### Score在GradientDescent Pfunc中的使用

从`GradientDescentPfunc.java`，我们看到：

```java
// 计算energy的权重
result.scoreWeight = ctx.bcalc(result.econf.getScore());      // e^(-score/RT)
result.energyWeight = ctx.bcalc(result.econf.getEnergy());    // e^(-energy/RT)
```

**关键洞察**：

1. **Lower bound使用scoreWeight**:
```java
state.lowerScoreWeightSum += scoreWeight;  // Σ e^(-score/RT)
```

2. **Upper bound使用energyWeight**:
```java
state.energyWeightSum += energyWeight;     // Σ e^(-energy/RT)
```

3. **为什么这样做？**

由于 `score ≤ energy`：
```
e^(-score/RT) ≥ e^(-energy/RT)
```

因此：
```
Z_lower = Σ e^(-score/RT) ≥ Σ e^(-energy/RT) ≥ True Z
```

等等，这似乎反了！让我重新分析...

实际上，K*使用的是：
- **能量最小化后的conformations**用于计算true partition function
- **Score（A* f-score）**用于指导搜索顺序（最promising的先）

### 实际测试数据验证

从测试输出：

| Conf# | Score | Energy | Relationship | Boltzmann Factor Relationship |
|-------|---------|----------|--------------|------------------------------|
| 1 | -5.085 | -4.773 | Score < Energy ✅ | e^(-(-5.085)) < e^(-(-4.773)) ❌ |
| 2 | -5.261 | -4.891 | Score < Energy ✅ | e^(5.261) < e^(4.891) ❌ |
| 20 | -2.647 | -2.163 | Score < Energy ✅ | e^(2.647) < e^(2.163) ❌ |

**重要发现**：

**能量是负数**（结合能），所以：
- Score = -5.261 （更负 = 更稳定）
- Energy = -4.891 （较不负 = 较不稳定）
- Score < Energy ⟹ |Score| > |Energy| ⟹ Score更稳定

**Boltzmann因子**：
```
e^(-E/RT) where RT ≈ 0.6 kcal/mol at 300K

e^(-(-5.261)/0.6) = e^(8.768) ≈ 6.4×10^3
e^(-(-4.891)/0.6) = e^(8.152) ≈ 3.5×10^3

scoreWeight > energyWeight ✅
```

### 正确的Partition Function Bounds

从代码逻辑：

```java
// Lower bound: 只用已最小化的构象的score
lowerScoreWeightSum = Σ e^(-score_i/RT) for energied confs

// Upper bound: 已最小化的energy + 未最小化的score
energyWeightSum = Σ e^(-energy_i/RT) for energied confs
upperScoreWeightSum = Σ e^(-score_j/RT) for scored but not energied confs
```

**为什么这是正确的？**

1. **已最小化的构象**：
   - Score是下界：score ≤ true energy
   - 所以 e^(-score/RT) ≥ e^(-true_energy/RT)
   - 用score给出了乐观（较大）的Boltzmann weight

2. **未最小化的构象**：
   - 只有A*的score，还没做能量最小化
   - Score仍然是下界
   - 用score估计它们的最大可能贡献

3. **整体逻辑**：
   ```
   Q_lower = Σ_{energied} e^(-score/RT) + 0 (不算未energied的)
   P_upper = Σ_{energied} e^(-score/RT) + Σ_{scored} e^(-score/RT)
   ```

### 结论

✅ **Score确实是能量的下界**

✅ **Score来自A*的f-score = g-score + h-score**

✅ **所有测试数据都满足 Score ≤ Energy**

✅ **使用score计算的Boltzmann weight会高估partition function，因此适合作为upper bound的一部分**

✅ **GradientDescent通过逐步最小化conformations，将score替换为真实energy，从而收紧bounds**

## 附录B：完整测试输出

### K* 详细日志
```
[8  3  8] scores:0, confs:1, score:-5.085246, energy:-4.772786,
          bounds:[3.495281, Infinity] (log10p1), delta:1.000000
[8  5  8] scores:0, confs:2, score:-5.261103, energy:-4.891207,
          bounds:[3.841756, Infinity] (log10p1), delta:1.000000
[0  4  8] scores:0, confs:3, score:-3.597335, energy:-3.131962,
          bounds:[3.853876, Infinity] (log10p1), delta:1.000000
[0  3  8] scores:0, confs:4, score:-4.875070, energy:-4.396592,
          bounds:[3.944541, Infinity] (log10p1), delta:1.000000
[8  4  8] scores:0, confs:5, score:-3.819221, energy:-3.528431,
          bounds:[3.963070, Infinity] (log10p1), delta:1.000000
[0  5  8] scores:0, confs:6, score:-5.044793, energy:-4.507073,
          bounds:[4.048544, Infinity] (log10p1), delta:1.000000
[8  5  0] scores:0, confs:7, score:-3.399019, energy:-2.825649,
          bounds:[4.053075, Infinity] (log10p1), delta:1.000000
[2  5  8] scores:10, confs:8, score:-3.194442, energy:-2.932978,
          bounds:[4.058444, 5.057824] (log10p1), delta:0.899865
[4  5  8] scores:280, confs:9, score:-2.861508, energy:-2.619082,
          bounds:[4.061575, 4.123187] (log10p1), delta:0.132274
[8  5  6] scores:280, confs:10, score:-3.004802, energy:-2.330290,
          bounds:[4.063488, 4.119648] (log10p1), delta:0.121311
...
[0  2  8] scores:280, confs:20, score:-2.647334, energy:-2.162909,
          bounds:[4.083185, 4.097184] (log10p1), delta:0.031723
```

### MARK* 详细日志
```
Running until leaf is found...
MARK* expanded: 18, queued: 0, scored/sec: 27011,
      partial minimizations: 0, time: 2.81 ms,
      delta: 1.00000000000, heapMem: 21.2% of 2.0 GiB
Found a leaf!
```

---

**文档创建时间**: 2025年
**最后更新**: 本次测试执行后
**状态**: MARK*需要修复后重新测试
