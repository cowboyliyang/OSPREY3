# OSPREY能量矩阵完整说明

## 概览

OSPREY使用多种能量矩阵来存储和计算蛋白质构象的能量。理解这些矩阵是理解MARK* bug的关键。

---

## 1. 基本概念

### 构象空间 (Conformation Space)
假设我们有一个简单的系统：
- **3个位置 (positions)**: pos=0, 1, 2
- **每个位置有多个rotamers (RCs - Rotational Conformers)**:
  - pos 0: 3个rotamers (RC 0, 1, 2)
  - pos 1: 2个rotamers (RC 0, 1)
  - pos 2: 4个rotamers (RC 0, 1, 2, 3)

**一个完整的conformation**：选择每个位置的一个rotamer
- 例如：`[0, 1, 2]` = pos0选RC0, pos1选RC1, pos2选RC2
- 总共有 3×2×4 = 24 种可能的conformations

### 能量计算
一个conformation的总能量包括：
1. **Constant term**: 固定的template能量（不变部分的能量）
2. **One-body terms**: 单个residue的内部能量
3. **Pairwise terms**: 两个residues之间的相互作用能量
4. **Higher-order terms**: 三个或更多residues的相互作用（通常忽略或近似）
  
**能量公式**：
```
E(conf) = E_const + Σ E_one(pos_i, RC_i) + Σ E_pair(pos_i, RC_i, pos_j, RC_j) + ...
```

---

## 2. EnergyMatrix (基础能量矩阵)

### 作用
存储预计算的pairwise能量（one-body和two-body terms）。

### 数据结构
```
constTerm: double                           # 常数项
oneBody[pos][rc]: double                    # 一体能量
pairwise[pos1][rc1][pos2][rc2]: double      # 两体能量
higherOrder[tuple]: double                  # 高阶能量（可选）
```

### 数字例子：Minimizing Energy Matrix

假设经过coordinate descent minimization后得到的能量：

```
constTerm = -10.0 kcal/mol

oneBody (one-body energies after minimization):
  pos 0, RC 0: -2.5 kcal/mol
  pos 0, RC 1: -1.8 kcal/mol
  pos 0, RC 2: -3.2 kcal/mol
  pos 1, RC 0: -4.1 kcal/mol
  pos 1, RC 1: -3.7 kcal/mol
  pos 2, RC 0: -1.2 kcal/mol
  pos 2, RC 1: -2.8 kcal/mol
  pos 2, RC 2: -3.5 kcal/mol
  pos 2, RC 3: -2.1 kcal/mol

pairwise (pairwise energies after minimization):
  (pos0,RC0) - (pos1,RC0): -0.5 kcal/mol
  (pos0,RC0) - (pos1,RC1): -0.3 kcal/mol
  (pos0,RC0) - (pos2,RC0): -0.8 kcal/mol
  (pos0,RC0) - (pos2,RC1): -0.4 kcal/mol
  (pos0,RC0) - (pos2,RC2): -0.6 kcal/mol
  ... (more pairs)
  (pos1,RC0) - (pos2,RC2): -1.2 kcal/mol
  ... (etc)
```

### 计算conformation能量
对于conformation `[0, 1, 2]` (pos0=RC0, pos1=RC1, pos2=RC2):

```
E = constTerm
  + oneBody[0][0] + oneBody[1][1] + oneBody[2][2]
  + pairwise[0][0][1][1] + pairwise[0][0][2][2] + pairwise[1][1][2][2]

E = -10.0
  + (-2.5) + (-3.7) + (-3.5)
  + (-0.3) + (-0.6) + (-0.9)

E = -10.0 - 9.7 - 1.8 = -21.5 kcal/mol
```

---

## 3. Rigid Energy Matrix vs Minimizing Energy Matrix

OSPREY计算**两套**能量矩阵：

### 3.1 Rigid Energy Matrix (rigidEmat)

**特点**：
- **不做minimization**
- 所有DOFs (degrees of freedom) 固定在初始位置
- 能量通常**更高**（不够稳定）
- 计算**快速**

**数字例子**：
```
constTerm = -10.0 kcal/mol

oneBody (rigid, no minimization):
  pos 0, RC 0: -1.8 kcal/mol  (比minimized的-2.5高)
  pos 0, RC 1: -1.2 kcal/mol
  pos 0, RC 2: -2.5 kcal/mol
  pos 1, RC 0: -3.2 kcal/mol
  pos 1, RC 1: -2.9 kcal/mol
  pos 2, RC 0: -0.8 kcal/mol
  pos 2, RC 1: -2.1 kcal/mol
  pos 2, RC 2: -2.8 kcal/mol
  pos 2, RC 3: -1.5 kcal/mol

pairwise (rigid):
  (pos0,RC0) - (pos1,RC0): -0.2 kcal/mol  (比minimized的-0.5高)
  (pos0,RC0) - (pos1,RC1): -0.1 kcal/mol
  (pos0,RC0) - (pos2,RC2): -0.3 kcal/mol
  ... (all values are higher/less favorable than minimized)
```

**用途**：
- 提供能量的**upper bound**（上界）
- 快速预估，用于A*搜索的heuristic

### 3.2 Minimizing Energy Matrix (minimizingEmat)

**特点**：
- **经过coordinate descent minimization**
- DOFs调整到局部能量最小值
- 能量通常**更低**（更稳定）
- 计算**慢**

**数字例子**：见上面Section 2的例子

**用途**：
- 提供能量的**lower bound**（下界）
- 更准确的能量估计
- 用于最终的partition function计算

### 比较

对于同一个conformation `[0, 1, 2]`:

```
Rigid Energy:      E_rigid     = -18.2 kcal/mol  (upper bound, less stable)
Minimized Energy:  E_minimized = -21.5 kcal/mol  (lower bound, more stable)

差距: 3.3 kcal/mol
```

**物理意义**：
- Rigid能量：假设所有原子固定，sterically strained
- Minimized能量：允许小幅度调整，释放strain，更接近真实能量

---

## 4. NegatedEnergyMatrix

### 作用
**将所有能量值取反号**的wrapper。

### 实现
```java
public class NegatedEnergyMatrix extends ProxyEnergyMatrix {
    @Override
    public Double getOneBody(int pos, int rc) {
        return -super.getOneBody(pos, rc);  // 取反！
    }

    @Override
    public Double getPairwise(int pos1, int rc1, int pos2, int rc2) {
        return -super.getPairwise(pos1, rc1, pos2, rc2);  // 取反！
    }
}
```

### 数字例子

**原始 Rigid Energy Matrix**:
```
oneBody[0][0] = -1.8 kcal/mol
oneBody[0][1] = -1.2 kcal/mol
oneBody[2][2] = -2.8 kcal/mol
pairwise[0][0][1][0] = -0.2 kcal/mol
```

**NegatedEnergyMatrix(rigidEmat)**:
```
oneBody[0][0] = +1.8 kcal/mol  (negated!)
oneBody[0][1] = +1.2 kcal/mol
oneBody[2][2] = +2.8 kcal/mol
pairwise[0][0][1][0] = +0.2 kcal/mol
```

### 为什么要取反？

**目的**：配合optimizer方向来计算upper/lower bounds

**Heuristic Scorer的逻辑**：
```java
// TraditionalPairwiseHScorer.calc()
for (int pos : undefinedPositions) {
    double optEnergy = optimizer.initDouble();  // Minimize: +∞, Maximize: -∞
    for (int rc : rcs[pos]) {
        double rcEnergy = getOneBody(pos, rc) + getPairwise(...);
        optEnergy = optimizer.opt(optEnergy, rcEnergy);  // min or max
    }
    hScore += optEnergy;
}
```

#### 场景1：计算Lower Bound (最小可能能量)
```java
hScorer = new TraditionalPairwiseHScorer(minimizingEmat, rcs, Optimizer.Minimize);

// pos2还未赋值，在RC 0,1,2,3中选最小的：
oneBody[2][0] = -0.8 kcal/mol
oneBody[2][1] = -2.1 kcal/mol
oneBody[2][2] = -2.8 kcal/mol ← 最小
oneBody[2][3] = -1.5 kcal/mol

hScore_pos2 = min(-0.8, -2.1, -2.8, -1.5) = -2.8 kcal/mol
```

#### 场景2：计算Upper Bound (最大可能能量)

**方法1（直接）**：
```java
hScorer = new TraditionalPairwiseHScorer(rigidEmat, rcs, Optimizer.Maximize);

// 在rigid能量上求最大：
oneBody[2][0] = -0.8 kcal/mol
oneBody[2][1] = -2.1 kcal/mol ← 最大（最不稳定）
oneBody[2][2] = -2.8 kcal/mol
oneBody[2][3] = -1.5 kcal/mol

hScore_pos2 = max(-0.8, -2.1, -2.8, -1.5) = -0.8 kcal/mol
```

**方法2（使用NegatedEnergyMatrix）**：
```java
negatedEmat = new NegatedEnergyMatrix(confSpace, rigidEmat);
hScorer = new TraditionalPairwiseHScorer(negatedEmat, rcs, Optimizer.Minimize);

// Negated energies:
oneBody[2][0] = +0.8 kcal/mol
oneBody[2][1] = +2.1 kcal/mol ← 在negated上最小
oneBody[2][2] = +2.8 kcal/mol
oneBody[2][3] = +1.5 kcal/mol

hScore_pos2 = min(+0.8, +2.1, +2.8, +1.5) = +0.8 kcal/mol

// 再取反回来：
actual_hScore = -hScore_pos2 = -0.8 kcal/mol  ✓ 与方法1一致！
```

**关键点**：
- `Minimize on negated values` = `Maximize on original values`
- 这允许代码使用**统一的Minimize optimizer**，通过取反来实现Maximize

---

## 5. UpdatingEnergyMatrix (Correction Matrix)

### 作用
在pairwise能量矩阵基础上，动态添加**higher-order corrections**（三体或更高阶的能量修正）。

### 数据结构
```
target: EnergyMatrix              # 底层的pairwise矩阵
corrections: Map<RCTuple, double> # 高阶修正项
```

### 为什么需要Correction？

**问题**：Pairwise近似不完美

```
真实的三体能量 ≠ 从pairwise加和估计的能量
```

**数字例子**：

Conformation `[0, 1, 2]`:

**从pairwise matrix估计**：
```
E_pairwise = constTerm + oneBody[0][0] + oneBody[1][1] + oneBody[2][2]
           + pairwise[0][0][1][1] + pairwise[0][0][2][2] + pairwise[1][1][2][2]
           = -10.0 + (-2.5) + (-3.7) + (-3.5) + (-0.3) + (-0.6) + (-0.9)
           = -21.5 kcal/mol
```

**实际triple minimization的结果**：
```
E_actual = -22.8 kcal/mol  (更低！)
```

**原因**：
- Pairwise minimization时，每对residues独立优化
- Triple minimization时，三个residues**同时**优化，可能找到更好的配置
- 三体cooperative effects未被pairwise捕捉

**Correction**：
```
correction = E_actual - E_pairwise
          = -22.8 - (-21.5)
          = -1.3 kcal/mol  (negative correction!)
```

### UpdatingEnergyMatrix的使用

```java
UpdatingEnergyMatrix correctionMatrix = new UpdatingEnergyMatrix(confSpace, minimizingEmat);

// 在MARK*运行时，当检测到pairwise估计不准确时：
correctionMatrix.setHigherOrder(tuple, correction);

// 之后计算能量时：
E_corrected = minimizingEmat.getInternalEnergy(tuple) + correctionMatrix.getHigherOrder(tuple);
            = -21.5 + (-1.3)
            = -22.8 kcal/mol  ✓ 正确！
```

### Negative Correction的重要性

**这就是MARK* bug的核心！**

当`correction < 0`时：
- 说明pairwise matrix的lower bound**过高**
- Triple minimization找到了**更低的能量**
- **必须**应用negative correction，否则bounds不准确

**Bug的表现**：
```java
// 原始错误代码
if (correction > 0) {
    correctionMatrix.setHigherOrder(tuple, correction);
}
else {
    System.err.println("Negative correction for "+tuple);
    // ❌ 没有应用correction！
}
```

**后果**：
- Lower bound保持过高
- Partition function被高估
- Boltzmann weights错误
- 整个算法失效

---

## 6. MARK*中的矩阵使用

### 在MARKStarBound.java中

```java
// Line 356-359
rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
    gscorerFactory.make(minimizingEmat),    // g-scorer (已确定部分的能量)
    hscorerFactory.make(minimizingEmat),    // h-scorer (lower bound)
    gscorerFactory.make(rigidEmat),         // rigid g-scorer
    new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)  // negated h-scorer
);
```

### 各个scorer的作用

#### 1. gScorer (minimizing)
**作用**：计算已确定部分的能量
```
conf = [0, 1, ?]  (pos2未定)
gScore = minimizingEmat.getInternalEnergy([0,1])
       = oneBody[0][0] + oneBody[1][1] + pairwise[0][0][1][1]
       = -2.5 + (-3.7) + (-0.3)
       = -6.5 kcal/mol
```

#### 2. hScorer (minimizing, minimize optimizer)
**作用**：计算未确定部分的lower bound
```
conf = [0, 1, ?]  (pos2未定)

对pos2，选择最小的one-body + pairwise：
pos2=RC0: -0.8 + pairwise[0][0][2][0] + pairwise[1][1][2][0] = -0.8 + (-0.8) + (-0.5) = -2.1
pos2=RC1: -2.1 + pairwise[0][0][2][1] + pairwise[1][1][2][1] = -2.1 + (-0.4) + (-0.7) = -3.2
pos2=RC2: -2.8 + pairwise[0][0][2][2] + pairwise[1][1][2][2] = -2.8 + (-0.6) + (-0.9) = -4.3 ← 最小
pos2=RC3: -1.5 + pairwise[0][0][2][3] + pairwise[1][1][2][3] = -1.5 + (-0.2) + (-0.6) = -2.3

hScore = min(-2.1, -3.2, -4.3, -2.3) = -4.3 kcal/mol

confLowerBound = gScore + hScore = -6.5 + (-4.3) = -10.8 kcal/mol
```

#### 3. rigidgScorer
**作用**：计算rigid能量（用于upper bound）
```
rigidScore = rigidEmat.getInternalEnergy([0,1])
           = oneBody_rigid[0][0] + oneBody_rigid[1][1] + pairwise_rigid[0][0][1][1]
           = -1.8 + (-2.9) + (-0.1)
           = -4.8 kcal/mol  (比minimized的-6.5高)
```

#### 4. negatedHScorer (NegatedEnergyMatrix + ?)
**作用**：计算upper bound的heuristic

**这就是bug所在！**

**原始设计意图（推测）**：
```java
negatedHScorer = new TraditionalPairwiseHScorer(
    new NegatedEnergyMatrix(confSpace, rigidEmat),
    rcs,
    Optimizer.Maximize  // ← 应该用Maximize！
);
```

这样：
- Maximize on negated rigid energies
- = Minimize on original rigid energies
- = 找到最不稳定的配置
- = Upper bound

**然后计算confUpperBound**：
```java
confUpperBound = rigidScore - negatedHScore
```

---

## 7. Upper Bound计算的数学

### 目标
找到conformation能量的**上界**（最不稳定的可能值）。

### 公式推导

**符号定义**：
- `E_rigid[i]`: position i的rigid能量
- `E_min[i]`: position i选择最差RC时的能量（upper bound）

**对于部分conformation `[0, 1, ?]`**:

**Lower bound** (已经确定的部分用minimized，未确定的选最好):
```
confLowerBound = E_minimized[0] + E_minimized[1] + min{E_minimized[2][rc]}
```

**Upper bound** (已经确定的部分用rigid，未确定的选最差):
```
confUpperBound = E_rigid[0] + E_rigid[1] + max{E_rigid[2][rc]}
```

### MARK*的实现

```java
// MARKStarNode.makeRoot Line 405
double confUpperBound = rigidScore - negatedHScore;
```

**这个公式意味着什么？**

假设：
- `rigidScore = E_rigid[0] + E_rigid[1]` = -4.8 kcal/mol
- `negatedHScore = ?`

**情况A：如果negatedHScorer使用Maximize optimizer**
```
在negated rigid energy上求最大 = 在原始rigid energy上求最小

negatedHScore (on negated values) = max{+0.8, +2.1, +2.8, +1.5} = +2.8
对应原始能量 = -2.8 kcal/mol  (最稳定的RC)

confUpperBound = rigidScore - negatedHScore
              = -4.8 - 2.8
              = -7.6 kcal/mol

这是错误的！我们想要upper bound（最不稳定），但得到了lower bound（最稳定）！
```

**情况B：如果negatedHScorer使用Minimize optimizer（当前代码）**
```
在negated rigid energy上求最小 = 在原始rigid energy上求最大

negatedHScore (on negated values) = min{+0.8, +2.1, +2.8, +1.5} = +0.8
对应原始能量 = -0.8 kcal/mol  (最不稳定的RC)

confUpperBound = rigidScore - negatedHScore
              = -4.8 - 0.8
              = -5.6 kcal/mol

这看起来对了！但等等...
```

**实际测试中的问题**：
```
测试输出显示：
negatedHScore = -23940 kcal/mol  (巨大的负数！)
confUpperBound = 0 - (-23940) = +23940 kcal/mol  (巨大的正数！)
```

**说明什么地方还是错了！**

---

## 8. MARK* Bug的真正原因（推测）

### 问题1：公式本身可能有问题

`confUpperBound = rigidScore - negatedHScore` 这个公式可疑。

**更自然的公式应该是**：
```java
confUpperBound = rigidScore + hScore_upper
```

其中`hScore_upper`应该是在rigid能量上求maximum的结果。

### 问题2：NegatedEnergyMatrix的使用混淆

**可能的正确用法**：

#### 选项A：不使用NegatedEnergyMatrix
```java
hScorer_upper = new TraditionalPairwiseHScorer(rigidEmat, rcs, Optimizer.Maximize);
confUpperBound = rigidScore + hScorer_upper.calc(...);
```

#### 选项B：使用NegatedEnergyMatrix但公式不同
```java
negatedHScorer = new TraditionalPairwiseHScorer(
    new NegatedEnergyMatrix(confSpace, rigidEmat),
    rcs,
    Optimizer.Minimize  // 在negated上minimize = 在original上maximize
);
double negatedHScore = negatedHScorer.calc(...);  // 返回positive值
double hScore_upper = -negatedHScore;  // 转回负值
confUpperBound = rigidScore + hScore_upper;
```

### 问题3：测试结果中的-23940来自哪里？

这个值太大了，不可能是单个position的能量。很可能是：
1. **累积了所有positions的能量** (不应该，h-score应该只算未定义的部分)
2. **包含了所有RCs的总和** (不应该，应该是min或max，不是sum)
3. **重复计算或其他bug**

---

## 9. 建议的调试步骤

### Step 1: 添加详细debug输出

在`MARKStarNode.makeRoot()`中：

```java
System.out.println("=== Debug negatedHScorer ===");
for (int pos = 0; pos < confSpace.positions.size(); pos++) {
    for (int rc = 0; rc < rcs.get(pos).length; rc++) {
        double rigidEnergy = rigidEmat.getOneBody(pos, rc);
        double negatedEnergy = new NegatedEnergyMatrix(confSpace, rigidEmat).getOneBody(pos, rc);
        System.out.println(String.format("  pos %d, rc %d: rigid=%f, negated=%f",
            pos, rc, rigidEnergy, negatedEnergy));
    }
}
double hScore = hScorer.calc(confIndex, rcs);
double negHScore = negatedHScorer.calc(confIndex, rcs);
System.out.println("hScore (minimize on minimizing): " + hScore);
System.out.println("negatedHScore (on negated rigid): " + negHScore);
```

### Step 2: 测试不同的配置

1. **移除NegatedEnergyMatrix**:
```java
new TraditionalPairwiseHScorer(rigidEmat, rcs, Optimizer.Maximize)
```

2. **修改confUpperBound公式**:
```java
confUpperBound = rigidScore + hScore_upper;  // 不是减法
```

3. **对比结果**

---

## 10. 总结

### 矩阵层次结构

```
EnergyMatrix (base class)
  ├─ 直接实例
  │    ├─ minimizingEmat (经过CCD minimization)
  │    └─ rigidEmat (没有minimization)
  │
  ├─ ProxyEnergyMatrix (wrapper base)
  │    ├─ NegatedEnergyMatrix (取反所有值)
  │    └─ UpdatingEnergyMatrix (添加higher-order corrections)
  │
  └─ 其他variants
       ├─ LazyEnergyMatrix (按需计算)
       └─ ReducedEnergyMatrix (K* specific)
```

### 关键数值关系

```
对于同一个conformation:
  E_rigid > E_minimized  (rigid能量总是更高/更不稳定)
  E_pairwise ≈ E_actual  (pairwise估计接近但不完美)
  correction = E_actual - E_pairwise  (可正可负！)

对于partition function bounds:
  Z_lower = Σ exp(-E_minimized / RT)  (lower bound)
  Z_upper = Σ exp(-E_rigid / RT)      (upper bound)
  Z_lower < Z_actual < Z_upper
```

### MARK* Bug核心

1. **Negative correction未被应用** → lower bounds过高
2. **confUpperBound计算可能有误** → 导致初始bounds错误
3. **两个bug共同作用** → partition function错误10²²倍

---

## 下一步

需要通过debug输出确认：
1. `negatedHScore = -23940` 是如何产生的？
2. 正确的`confUpperBound`公式应该是什么？
3. NegatedEnergyMatrix的设计意图是什么？（可能需要查看原始论文或提交历史）

