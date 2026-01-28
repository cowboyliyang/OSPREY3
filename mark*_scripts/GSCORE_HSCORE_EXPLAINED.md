# gScore和hScore详解：从矩阵到A*搜索

## 1. 核心概念：A*搜索在构象空间

### A*公式

```
f(node) = g(node) + h(node)

其中：
- g(node): 从起点到当前节点的**实际代价**（已知部分）
- h(node): 从当前节点到目标的**估计代价**（未知部分的启发式）
- f(node): 总估计代价
```

### 在OSPREY中的应用

**构象空间搜索**：
- **起点**: 空conformation `[]` (没有position被赋值)
- **目标**: 完整conformation `[RC0, RC1, RC2, ...]` (所有positions都被赋值)
- **中间节点**: 部分conformation `[RC0, RC1, ?]` (部分positions已赋值)

**能量作为代价**：
- **g(node)**: 已赋值部分的能量（从能量矩阵计算）
- **h(node)**: 未赋值部分的能量估计（启发式）
- **f(node)**: 完整conformation的能量估计

---

## 2. gScore：已确定部分的能量

### 定义

**gScore** = 已经赋值的positions的**实际能量**

### 从哪个矩阵计算？

**取决于scorer类型**：

#### MARK*中的gScore配置

```java
// MARKStarBound.java Line 356-358
rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
    gscorerFactory.make(minimizingEmat),      // ← minimizing g-scorer
    hscorerFactory.make(minimizingEmat),      // ← minimizing h-scorer
    gscorerFactory.make(rigidEmat),           // ← rigid g-scorer
    negatedHScorer                            // ← negated h-scorer
);
```

**两种gScorer**：
1. **Minimizing gScorer**: 使用`minimizingEmat`（更准确）
2. **Rigid gScorer**: 使用`rigidEmat`（用于upper bound）

### 计算过程（PairwiseGScorer.java）

```java
public double calc(ConfIndex confIndex, RCs rcs) {
    // Step 1: Constant term
    double gscore = emat.getConstTerm();

    // Step 2: One-body energies (已赋值的positions)
    for (int i=0; i<confIndex.numDefined; i++) {
        int pos = confIndex.definedPos[i];
        int rc = confIndex.definedRCs[i];
        gscore += emat.getOneBody(pos, rc);
    }

    // Step 3: Pairwise energies (已赋值positions之间的)
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

### 数字例子

**系统配置**：
- pos0有3个RCs, pos1有2个RCs, pos2有4个RCs

**能量矩阵**（minimizingEmat）:
```
constTerm = -10.0 kcal/mol

oneBody:
  pos0, RC0: -2.5
  pos0, RC1: -1.8
  pos0, RC2: -3.2
  pos1, RC0: -4.1
  pos1, RC1: -3.7
  pos2, RC0: -1.2
  pos2, RC1: -2.8
  pos2, RC2: -3.5
  pos2, RC3: -2.1

pairwise:
  (pos0=RC0, pos1=RC1): -0.3
  (pos0=RC0, pos2=RC2): -0.6
  (pos1=RC1, pos2=RC2): -0.9
  ... (more pairs)
```

#### 场景A：部分conformation `[RC0, RC1, ?]`

**已定义**: pos0=RC0, pos1=RC1
**未定义**: pos2

**gScore计算**：
```
gscore = constTerm
       + oneBody[pos0][RC0]
       + oneBody[pos1][RC1]
       + pairwise[pos0=RC0][pos1=RC1]

gscore = -10.0
       + (-2.5)
       + (-3.7)
       + (-0.3)

gscore = -16.5 kcal/mol
```

**物理意义**：pos0和pos1已经固定，它们的能量贡献是确定的。

#### 场景B：完整conformation `[RC0, RC1, RC2]`

**已定义**: pos0=RC0, pos1=RC1, pos2=RC2
**未定义**: 无

**gScore计算**：
```
gscore = constTerm
       + oneBody[pos0][RC0]
       + oneBody[pos1][RC1]
       + oneBody[pos2][RC2]
       + pairwise[pos0=RC0][pos1=RC1]
       + pairwise[pos0=RC0][pos2=RC2]
       + pairwise[pos1=RC1][pos2=RC2]

gscore = -10.0
       + (-2.5) + (-3.7) + (-3.5)
       + (-0.3) + (-0.6) + (-0.9)

gscore = -21.5 kcal/mol
```

**物理意义**：完整conformation的总能量（从pairwise矩阵）。

---

## 3. hScore：未确定部分的能量估计

### 定义

**hScore** = 未赋值positions的**估计能量**（启发式）

### 从哪个矩阵计算？

**取决于scorer类型和用途**：

#### MARK*中的hScore配置

1. **Normal hScorer** (minimizing, for lower bound):
   ```java
   hscorerFactory.make(minimizingEmat)  // 使用minimizing能量
   ```

2. **Negated hScorer** (for upper bound):
   ```java
   new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)
   ```

### 计算过程（TraditionalPairwiseHScorer.java）

这个计算比gScore复杂得多！

#### Step 1: 预计算（构造函数）

```java
// Line 62-85
// 预计算所有"未定义positions之间"的最优pairwise能量
undefinedEnergies = new double[numPos][][];

for (int pos1=0; pos1<numPos; pos1++) {
    for (int rc1 : rcs.get(pos1)) {
        for (int pos2=0; pos2<pos1; pos2++) {
            // 对于pos1的每个RC，找到与pos2的最优pairwise能量
            double optEnergy = optimizer.initDouble();  // Minimize: +∞, Maximize: -∞
            for (int rc2 : rcs.get(pos2)) {
                optEnergy = optimizer.opt(optEnergy, emat.getPairwise(pos1, rc1, pos2, rc2));
            }
            undefinedEnergies[pos1][rc1][pos2] = optEnergy;
        }
    }
}
```

**物理意义**：
- 对于每个(pos1, RC)组合，预先计算它与其他positions的最优pairwise能量
- 如果Minimize: 选择最有利的interaction
- 如果Maximize: 选择最不利的interaction

#### Step 2: 缓存计算（calcCachedEnergies）

```java
// Line 176-211
// 对每个未定义position的每个RC，计算它与已定义+未定义positions的总能量
for (int pos1 : undefinedPositions) {
    for (int rc1 : rcs.get(pos1)) {

        // 2.1: One-body energy
        double energy = emat.getOneBody(pos1, rc1);

        // 2.2: Pairwise with DEFINED positions（确定的）
        for (int pos2 : definedPositions) {
            int rc2 = definedRCs[pos2];
            energy += emat.getPairwise(pos1, rc1, pos2, rc2);
        }

        // 2.3: Optimal pairwise with UNDEFINED positions（预计算的最优值）
        for (int pos2 : undefinedPositions) {
            if (pos2 < pos1) {
                energy += undefinedEnergies[pos1][rc1][pos2];
            }
        }

        cachedEnergies[pos1][rc1] = energy;
    }
}
```

#### Step 3: 最终hScore计算（calc）

```java
// Line 100-121
double hscore = 0;

// 对每个未定义position
for (int pos : undefinedPositions) {

    // 在所有RCs中选择最优的（minimize或maximize）
    double optRCEnergy = optimizer.initDouble();
    for (int rc : rcs.get(pos)) {
        optRCEnergy = optimizer.opt(optRCEnergy, cachedEnergies[pos][rc]);
    }

    hscore += optRCEnergy;
}

return hscore;
```

### 数字例子

#### 场景：部分conformation `[RC0, RC1, ?]`

**已定义**: pos0=RC0, pos1=RC1
**未定义**: pos2
**使用**: minimizingEmat, Optimizer.Minimize

##### Step 1: 预计算（在构造时完成）

对于pos2未定义时，不需要预计算的"未定义-未定义"interactions（因为pos2是唯一未定义的）。

##### Step 2: 缓存计算

对pos2的每个RC，计算它与已定义positions的能量：

```
pos2=RC0:
  energy = oneBody[pos2][RC0]
         + pairwise[pos2=RC0][pos0=RC0]  (与已定义的pos0)
         + pairwise[pos2=RC0][pos1=RC1]  (与已定义的pos1)
  energy = -1.2 + (-0.8) + (-0.7) = -2.7 kcal/mol

pos2=RC1:
  energy = oneBody[pos2][RC1]
         + pairwise[pos2=RC1][pos0=RC0]
         + pairwise[pos2=RC1][pos1=RC1]
  energy = -2.8 + (-0.4) + (-0.9) = -4.1 kcal/mol

pos2=RC2:
  energy = oneBody[pos2][RC2]
         + pairwise[pos2=RC2][pos0=RC0]
         + pairwise[pos2=RC2][pos1=RC1]
  energy = -3.5 + (-0.6) + (-0.9) = -5.0 kcal/mol ← 最小

pos2=RC3:
  energy = oneBody[pos2][RC3]
         + pairwise[pos2=RC3][pos0=RC0]
         + pairwise[pos2=RC3][pos1=RC1]
  energy = -2.1 + (-0.3) + (-0.4) = -2.8 kcal/mol
```

##### Step 3: 选择最优RC

```
hScore = min(-2.7, -4.1, -5.0, -2.8)
       = -5.0 kcal/mol
```

**物理意义**：如果我们为pos2选择最优的RC（RC2），它将贡献-5.0 kcal/mol的能量。

#### 总的f-score

```
f-score = g-score + h-score
        = -16.5 + (-5.0)
        = -21.5 kcal/mol
```

这是对完整conformation `[RC0, RC1, RC2]`能量的估计！

#### 验证

实际上如果我们选择pos2=RC2，真实能量是：
```
真实能量 = -21.5 kcal/mol  (从前面的gScore例子B)
```

**完全一致！** 这是因为：
1. 只有一个未定义position
2. 我们的启发式恰好选中了最优的RC
3. Pairwise近似是精确的（没有三体效应）

---

## 4. 多个未定义Positions的复杂例子

### 场景：部分conformation `[RC0, ?, ?]`

**已定义**: pos0=RC0
**未定义**: pos1, pos2
**使用**: minimizingEmat, Optimizer.Minimize

#### Step 1: 预计算

```
对于pos1的每个RC，预计算它与pos2的最优pairwise能量：

pos1=RC0 与 pos2:
  undefinedEnergies[pos1=RC0][pos2] = min over pos2_RCs of pairwise[pos1=RC0][pos2=RC?]
  = min(pairwise[1,0][2,0], pairwise[1,0][2,1], pairwise[1,0][2,2], pairwise[1,0][2,3])
  = min(-0.7, -0.9, -1.2, -0.6)
  = -1.2 kcal/mol  (选择pos2=RC2最优)

pos1=RC1 与 pos2:
  undefinedEnergies[pos1=RC1][pos2] = min over pos2_RCs of pairwise[pos1=RC1][pos2=RC?]
  = min(pairwise[1,1][2,0], pairwise[1,1][2,1], pairwise[1,1][2,2], pairwise[1,1][2,3])
  = min(-0.5, -0.7, -0.9, -0.4)
  = -0.9 kcal/mol  (选择pos2=RC2最优)

pos2的RC与pos1（对称）：
  undefinedEnergies[pos2=RC0][pos1] = min(-0.7, -0.5) = -0.7
  undefinedEnergies[pos2=RC1][pos1] = min(-0.9, -0.7) = -0.9
  undefinedEnergies[pos2=RC2][pos1] = min(-1.2, -0.9) = -1.2
  undefinedEnergies[pos2=RC3][pos1] = min(-0.6, -0.4) = -0.6
```

#### Step 2: 缓存计算

**对pos1的每个RC**：

```
pos1=RC0:
  energy = oneBody[pos1][RC0]
         + pairwise[pos1=RC0][pos0=RC0]              (与已定义)
         + undefinedEnergies[pos1=RC0][pos2]         (与未定义pos2的最优)
  energy = -4.1 + (-0.5) + (-1.2) = -5.8 kcal/mol

pos1=RC1:
  energy = oneBody[pos1][RC1]
         + pairwise[pos1=RC1][pos0=RC0]
         + undefinedEnergies[pos1=RC1][pos2]
  energy = -3.7 + (-0.2) + (-0.9) = -4.8 kcal/mol
```

**对pos2的每个RC**：

```
pos2=RC0:
  energy = oneBody[pos2][RC0]
         + pairwise[pos2=RC0][pos0=RC0]
         + undefinedEnergies[pos2=RC0][pos1]
  energy = -1.2 + (-0.8) + (-0.7) = -2.7 kcal/mol

pos2=RC1:
  energy = -2.8 + (-0.4) + (-0.9) = -4.1 kcal/mol

pos2=RC2:
  energy = -3.5 + (-0.6) + (-1.2) = -5.3 kcal/mol

pos2=RC3:
  energy = -2.1 + (-0.3) + (-0.6) = -3.0 kcal/mol
```

#### Step 3: 选择每个position的最优RC

```
For pos1:
  optRCEnergy = min(-5.8, -4.8) = -5.8 kcal/mol  (选RC0)

For pos2:
  optRCEnergy = min(-2.7, -4.1, -5.3, -3.0) = -5.3 kcal/mol  (选RC2)

hScore = -5.8 + (-5.3) = -11.1 kcal/mol
```

#### gScore

```
gScore = constTerm + oneBody[pos0][RC0]
       = -10.0 + (-2.5)
       = -12.5 kcal/mol
```

#### f-score

```
f-score = gScore + hScore
        = -12.5 + (-11.1)
        = -23.6 kcal/mol
```

#### 验证

如果我们实际选择 `[pos0=RC0, pos1=RC0, pos2=RC2]`:

```
真实能量 = constTerm
         + oneBody[0,0] + oneBody[1,0] + oneBody[2,2]
         + pairwise[0,0][1,0] + pairwise[0,0][2,2] + pairwise[1,0][2,2]
         = -10.0 + (-2.5) + (-4.1) + (-3.5) + (-0.5) + (-0.6) + (-1.2)
         = -22.4 kcal/mol
```

**比较**：
- f-score估计: -23.6 kcal/mol
- 实际能量: -22.4 kcal/mol
- **差异**: 1.2 kcal/mol（估计更低，即更optimistic）

**为什么有差异？**

hScore中使用了"独立的最优选择"：
- 对pos1选RC0，假设pos2=RC2最优（pairwise = -1.2）
- 对pos2选RC2，假设pos1=RC0最优（pairwise = -1.2）

但这导致了**双重计数**pairwise能量！实际的pairwise[1,0][2,2] = -1.2只应该算一次，但在hScore中算了两次。

**这是A*启发式的本质**：它是**乐观的估计**（对minimize问题），保证不会高估代价，从而保证找到最优解。

---

## 5. MARK*中的多种Scorer组合

### MARKStarBound的Scorer配置

```java
// Line 356-359
rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
    gscorerFactory.make(minimizingEmat),     // gScorer (minimizing)
    hscorerFactory.make(minimizingEmat),     // hScorer (minimizing)
    gscorerFactory.make(rigidEmat),          // rigidgScorer
    new TraditionalPairwiseHScorer(
        new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)  // negatedHScorer
);
```

### 四种Scorers的用途

#### 1. gScorer (minimizing)
- **矩阵**: minimizingEmat
- **用途**: 计算已确定部分的**真实能量**（lower bound）
- **公式**: `E_defined = Σ E_1 + Σ E_2` (已定义positions)

#### 2. hScorer (minimizing, Optimizer.Minimize)
- **矩阵**: minimizingEmat
- **用途**: 估计未确定部分的**最低可能能量**（lower bound heuristic）
- **公式**: 对每个未定义pos，选择与已定义+未定义positions的最优RC

#### 3. rigidgScorer
- **矩阵**: rigidEmat
- **用途**: 计算已确定部分的**rigid能量**（upper bound）
- **公式**: 与gScorer相同，但使用rigid能量

#### 4. negatedHScorer
- **矩阵**: NegatedEnergyMatrix(rigidEmat)
- **用途**: 估计未确定部分的**最高可能能量**（upper bound heuristic）
- **Optimizer**: 应该是Maximize（在negated值上minimize = 在original值上maximize）

### Bounds计算

```java
// MARKStarNode.makeRoot() Line 402-406
double confLowerBound = gscore + hScore;           // minimizing + minimizing
double confUpperBound = rigidScore - negatedHScore;  // rigid - negated

rootNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound);
```

**这里就是bug所在！** `confUpperBound`的公式可疑。

---

## 6. 总结对比表

| Score类型 | 矩阵来源 | Optimizer | 包含的Positions | 用途 |
|-----------|----------|-----------|----------------|------|
| **gScore** | minimizingEmat | N/A | 已定义 | Lower bound的已知部分 |
| **hScore** | minimizingEmat | Minimize | 未定义 | Lower bound的估计部分 |
| **rigidgScore** | rigidEmat | N/A | 已定义 | Upper bound的已知部分 |
| **negatedHScore** | NegatedEnergyMatrix(rigidEmat) | Minimize (on negated) = Maximize (on original) | 未定义 | Upper bound的估计部分 |

### 能量矩阵取值层次

```
oneBody[pos][rc]        ← 直接从矩阵读取
pairwise[pos1][rc1][pos2][rc2]  ← 直接从矩阵读取

gScore = Σ (oneBody + pairwise) for defined positions  ← 简单求和

hScore = Σ min/max over RCs of (oneBody + pairwise + optimal_undefined_pairs)  ← 优化选择
```

### 数字范围示例

```
典型蛋白质conformation (3个positions):

oneBody: -5 ~ -1 kcal/mol per position
pairwise: -1 ~ -0.1 kcal/mol per pair

gScore (2 positions defined):
  ≈ -10 + (-2.5) + (-3.7) + (-0.3) = -16.5 kcal/mol

hScore (1 position undefined):
  ≈ min(-2.7, -4.1, -5.0, -2.8) = -5.0 kcal/mol

f-score = -16.5 + (-5.0) = -21.5 kcal/mol
```

---

## 7. 关键洞察

### gScore和hScore从哪个矩阵计算？

**答案取决于用途**：

1. **Lower Bound计算**:
   - gScore: minimizingEmat
   - hScore: minimizingEmat, Optimizer.Minimize

2. **Upper Bound计算**:
   - rigidgScore: rigidEmat
   - negatedHScore: NegatedEnergyMatrix(rigidEmat), Optimizer应该Maximize

### 为什么需要两套Scorers？

**Partition Function需要bounds**：

```
Z_lower ≤ Z_actual ≤ Z_upper

其中：
Z_lower = Σ exp(-(confLowerBound) / RT)
Z_upper = Σ exp(-(confUpperBound) / RT)
```

**confLowerBound**: 越低越好，使用minimizing能量
**confUpperBound**: 越高越好（越pessimistic），使用rigid能量

### hScore的精妙之处

hScore不是简单的矩阵查询，而是一个**优化过程**：
1. 预计算未定义positions之间的最优pairwise
2. 对每个未定义position的每个RC，计算与已定义和未定义positions的总能量
3. 对每个position选择最优RC
4. 求和

这保证了：
- **Admissible**: hScore不会高估真实代价（对minimize问题）
- **Consistent**: 满足A*的要求
- **Informative**: 提供了比简单和更好的指导

### 三体效应的影响

hScore基于pairwise矩阵，所以**不包括三体效应**。这就是为什么需要triple corrections：

```
真实能量 = gScore + hScore + triple_corrections
```

但在A*搜索时，我们用的是pairwise估计，之后再用corrections修正。
