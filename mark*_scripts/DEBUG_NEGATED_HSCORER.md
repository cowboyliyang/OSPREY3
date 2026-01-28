# Debug: NegatedHScorer Analysis

## 问题

MARKStarNode.java:405:
```java
double confUpperBound = rigidScore - negatedHScore;
```

当前：
- rigidScore = 0 (root node没有defined positions)
- negatedHScore = 35.5
- confUpperBound = -35.5 kcal/mol

这导致Boltzmann weight = e^(35.5/0.592) ≈ 10^26，完全错误。

## NegatedHScorer的实现

Line 368:
```java
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)
```

### TraditionalPairwiseHScorer的行为

使用默认构造函数 → `Optimizer.Minimize`

对于root node（所有positions都undefined）：
```java
hscore = 0;
for (each undefined position) {
    // 在这个position的所有RCs中选最优的
    optRCEnergy = optimizer.initDouble();  // Minimize → +∞
    for (each RC at this position) {
        energy = one-body + pairwise with other undefined positions
        optRCEnergy = optimizer.opt(optRCEnergy, energy);  // min(current, energy)
    }
    hscore += optRCEnergy;  // 加上这个position的最小能量
}
```

### NegatedEnergyMatrix的效果

```
Original rigidEmat (negative = favorable):
  pos0: one-body energies = [-50, -48, -49, -10, -9.5]
  pos1: one-body energies = [-30, -28, -29]
  pos2: one-body energies = [-20, -19, -18]

After NegatedEnergyMatrix (positive):
  pos0: one-body energies = [+50, +48, +49, +10, +9.5]
  pos1: one-body energies = [+30, +28, +29]
  pos2: one-body energies = [+20, +19, +18]

TraditionalPairwiseHScorer with Minimize:
  pos0: min(50, 48, 49, 10, 9.5) = 9.5  ← 选最小的正数
  pos1: min(30, 28, 29) = 28
  pos2: min(20, 19, 18) = 18

  negatedHScore = 9.5 + 28 + 18 = 55.5
```

## 问题分析

### 如果我们想要confUpperBound是能量的upper bound：

```
confUpperBound = rigidScore - negatedHScore
              = 0 - 55.5
              = -55.5 kcal/mol (非常低的能量)
```

这不对！Upper bound应该是**高能量**（不稳定），不是低能量。

### MARK*的意图是什么？

让我们考虑两种可能的设计意图：

#### 意图A：confUpperBound应该是能量的最大值（最不稳定）

那么我们需要：
```
confUpperBound = rigidScore + (max energies of undefined positions)
```

对于negatedHScore，如果我们想得到max energies：
```
Option 1: NegatedEnergyMatrix + Optimizer.Maximize
  - Negated: [+50, +48, +49, +10, +9.5]
  - Maximize: max(50, 48, 49, 10, 9.5) = 50
  - This gives the MOST POSITIVE value after negation
  - Which corresponds to MOST NEGATIVE original value (-50)
  - But we want MOST POSITIVE original value!

Option 2: Original EnergyMatrix + Optimizer.Maximize
  - Original: [-50, -48, -49, -10, -9.5]
  - Maximize: max(-50, -48, -49, -10, -9.5) = -9.5
  - This is LEAST NEGATIVE = HIGHEST energy
  - Then: confUpperBound = 0 - (-9.5) = +9.5
  - NO! Still wrong formula
```

#### 意图B：confUpperBound实际上是某种bound，不是直接的能量

也许公式中的减号有特殊含义？

让我查看K*如何计算conf upper bound...

## K*的做法（参考）

K* GradientDescentPfunc不使用这种公式。它用：
```
upperBound = (unscored confs) × (worst possible weight) + (scored upper bounds)
```

这是在Boltzmann weight空间计算，不是在能量空间。

## 可能的修复方案

### 方案1：修改optimizer为Maximize

```java
new TraditionalPairwiseHScorer(
    new NegatedEnergyMatrix(confSpace, rigidEmat),
    rcs,
    MathTools.Optimizer.Maximize  // ← 添加这个
)
```

但这会给出什么？
- Negated: [+50, +48, +49, +10, +9.5]
- Maximize: max(50, 48, 49, 10, 9.5) = 50
- negatedHScore = 50 + 30 + 20 = 100
- confUpperBound = 0 - 100 = -100

还是负数！

### 方案2：不使用NegatedEnergyMatrix

```java
new TraditionalPairwiseHScorer(
    rigidEmat,  // 直接使用原始矩阵
    rcs,
    MathTools.Optimizer.Maximize  // 选最高能量
)
```

这会给出：
- Original: [-50, -48, -49, -10, -9.5]
- Maximize: max(-50, -48, -49, -10, -9.5) = -9.5
- hScore = -9.5 + (-28) + (-18) = -55.5
- confUpperBound = 0 - (-55.5) = +55.5 kcal/mol

这是正数！而且物理意义：
- 选择每个position**最不稳定的RC**（能量最高）
- confUpperBound = +55.5是最差情况的能量

### 方案3：修改公式

也许公式本身就是错的？应该是：
```java
confUpperBound = rigidScore + negatedHScore
```

那么用原始代码：
- negatedHScore = 35.5
- confUpperBound = 0 + 35.5 = +35.5 kcal/mol

这也是正数，而且数值合理！

## 下一步

需要确定MARK*论文或设计文档中的公式是什么。

或者：测试方案2和方案3，看哪个给出合理的partition function结果。
