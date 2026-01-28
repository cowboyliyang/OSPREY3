# MARK* vs K* 测试结果

## 测试配置

- 系统：Small System (3 flexible residues)
- Epsilon: 0.1
- Total conformations: 24

## K* GradientDescent 结果（正确的baseline）

```
Q* (lower):           1.208488e+04
P* (upper):           1.250764e+04
Free Energy Lower:   -5.594892 kcal/mol
Free Energy Upper:   -5.574500 kcal/mol
Confs Evaluated:      18
Time (ms):            90
```

**分析**：
- Partition function数量级：10⁴
- Free energy：约-5.6 kcal/mol
- 评估了18个conformations才收敛

## MARK* 结果（当前有bug）

```
Q* (lower):           2.617500e+26
P* (upper):           2.617500e+26
Free Energy Lower:   -36.074905 kcal/mol
Free Energy Upper:   -36.074905 kcal/mol
Confs Evaluated:      0
Time (ms):            34
```

**Debug输出**：
```
[MARK* DEBUG] negatedHScore=35.522459443876556, confUpperBound=-35.522459443876556
```

**分析**：
- Partition function：2.6×10²⁶ (**比K*大10²²倍！**)
- Free energy：-36.07 kcal/mol (**比K*低30 kcal/mol！**)
- **没有评估任何conformation**（0个）
- 时间更快（34 ms vs 90 ms），但结果完全错误

## 问题对比

| 指标 | K* | MARK* | 差异 |
|------|-----|--------|------|
| Z (partition function) | 1.2×10⁴ | 2.6×10²⁶ | **10²² 倍** |
| Free Energy | -5.6 kcal/mol | -36.1 kcal/mol | **-30.5 kcal/mol** |
| Confs Evaluated | 18 | 0 | **MARK*没有评估** |
| Time | 90 ms | 34 ms | MARK*更快（但错误） |

## 根本问题

### 1. confUpperBound计算错误

```
confUpperBound = rigidScore - negatedHScore
               = 0 - 35.52
               = -35.52 kcal/mol
```

这个值太低（负数太大），导致：

```
Boltzmann weight = e^(-confUpperBound/RT)
                 = e^(-(-35.52)/0.592)
                 = e^(60)
                 ≈ 10²⁶
```

### 2. 没有评估conformations

MARK*在`runUntilNonZero()`后就停止了，认为epsilon已经满足，因为：
- 初始bounds已经非常tight（upper ≈ lower）
- 但这些bounds是**完全错误**的

### 3. NegatedEnergyMatrix的作用

当前代码（修改后）：
```java
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs, MathTools.Optimizer.Maximize)
```

- `NegatedEnergyMatrix`把所有能量取反
- `Optimizer.Maximize`在negated值上求最大
- 结果：negatedHScore = +35.52

问题：这个+35.52是什么？为什么会导致confUpperBound = -35.52？

## 物理意义分析

Free energy = -RT ln(Z)

- K*: Z = 1.2×10⁴ → ΔG = -0.592 × ln(1.2×10⁴) = -5.59 kcal/mol ✓
- MARK*: Z = 2.6×10²⁶ → ΔG = -0.592 × ln(2.6×10²⁶) = -36.07 kcal/mol ✗

**物理上不可能**：
- 30 kcal/mol的差异相当于能量差了10²¹倍
- 这比宇宙中所有原子的数量还多
- 明显是代码bug，不是物理现象

## 结论

当前的"修复"（添加`Optimizer.Maximize`）**让问题变得更糟**，需要重新分析并正确修复。

原始代码（没有Optimizer参数）可能有其自己的问题，但当前的修复方向是错误的。

需要深入理解：
1. `confUpperBound`的物理意义
2. `NegatedEnergyMatrix`的设计意图
3. 为什么成功的MARK*测试能工作，而我们的测试不行
