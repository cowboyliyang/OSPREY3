# MARK* Triple Energy详解

## 1. 什么是Triple Energy？

Triple energy（三体能量）是**三个residues同时存在**时的能量，它捕捉了pairwise能量矩阵无法描述的**三体协同效应**。

---

## 2. 能量分解层次

### 完整的能量公式

```
E_total(conf) = E_const                           (constant term)
              + Σ E_1(pos_i, RC_i)               (one-body terms)
              + Σ E_2(pos_i, RC_i, pos_j, RC_j)  (pairwise terms)
              + Σ E_3(pos_i, RC_i, pos_j, RC_j, pos_k, RC_k)  (triple terms)
              + ...                               (higher-order)
```

### Pairwise近似

大多数算法（包括K*）只计算到pairwise level：

```
E_pairwise(conf) = E_const + Σ E_1 + Σ E_2
```

**假设**：忽略三体及更高阶项，或假设它们为0。

### MARK*的改进

MARK*动态计算**selected triples**的真实能量，并添加corrections。

---

## 3. 为什么Pairwise不够？

### 数字例子

考虑conformation `[pos0=RC0, pos1=RC1, pos2=RC2]`。

#### 从Pairwise矩阵估计

```
E_pairwise = E_const + E_1(0,0) + E_1(1,1) + E_1(2,2)
           + E_2(0,0,1,1) + E_2(0,0,2,2) + E_2(1,1,2,2)

代入数值（前面例子中的）：
E_pairwise = -10.0 + (-2.5) + (-3.7) + (-3.5)
           + (-0.3) + (-0.6) + (-0.9)
         = -10.0 - 9.7 - 1.8
         = -21.5 kcal/mol
```

#### 实际Triple Minimization

当我们**同时优化三个residues**时：

```
E_actual = -22.8 kcal/mol  (比pairwise估计低1.3 kcal/mol！)
```

#### 为什么不同？

**物理原因**：

1. **Pairwise minimization**:
   - 每对residues **独立**优化
   - pos0与pos1优化时，pos2不参与
   - pos0与pos2优化时，pos1不参与
   - pos1与pos2优化时，pos0不参与

2. **Triple minimization**:
   - 三个residues **同时**优化
   - 可能找到一个three-way configuration，其中：
     - pos0略微调整以accommodate pos1和pos2
     - pos1略微调整以accommodate pos0和pos2
     - pos2略微调整以accommodate pos0和pos1
   - 这种mutual accommodation在pairwise level无法捕捉

**类比**：
- Pairwise: 三个人两两握手，但彼此不看对方
- Triple: 三个人围成圈，同时握手并调整位置

---

## 4. MARK*中Triple Energy的计算

### 触发条件

在`MARKStarBound.computeDifferences()`中（Line 1028-1056）：

```java
double triplethreshold = 0.3;  // 只计算"可疑"的triples

for (each scored conformation conf) {
    for (each triple (pos1, pos2, pos3) in conf) {
        // 计算pairwise估计的误差范围
        double tupleBounds = rigidEmat.getInternalEnergy(tuple)
                           - minimizingEmat.getInternalEnergy(tuple);

        if (tupleBounds >= triplethreshold) {
            // 这个triple的rigid和minimized能量差异较大
            // 可能存在significant三体效应
            computeDifference(tuple, minimizingEcalc);
        }
    }
}
```

**逻辑**：
- 如果`rigid - minimized`差距大（≥0.3 kcal/mol）
- 说明minimization对这个triple有显著影响
- 值得做full triple minimization来获得更准确的能量

### Triple Energy计算过程

在`computeDifference()`中（Line 1061-1086）：

```java
private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
    // tuple = {pos: [pos1, pos2, pos3], RCs: [RC1, RC2, RC3]}

    // Step 1: 实际minimization（coordinate descent on all 3 residues simultaneously）
    minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
        double tripleEnergy = minimizedTuple.energy;  // 真实的三体最小化能量

        // Step 2: 从pairwise矩阵估计
        double lowerbound = minimizingEmat.getInternalEnergy(tuple);
        // 这包括：E_1(pos1,RC1) + E_1(pos2,RC2) + E_1(pos3,RC3)
        //        + E_2(pos1,RC1,pos2,RC2) + E_2(pos1,RC1,pos3,RC3) + E_2(pos2,RC2,pos3,RC3)

        // Step 3: 计算correction
        double correction = tripleEnergy - lowerbound;

        // Step 4: 存储correction（这就是我们修复的bug！）
        correctionMatrix.setHigherOrder(tuple, correction);
    });
}
```

### getInternalEnergy()做了什么？

在`EnergyMatrix.java`中：

```java
public double getInternalEnergy(RCTuple tup) {
    double energy = 0;

    // 加上所有one-body terms
    for (int i = 0; i < tup.pos.size(); i++) {
        int pos = tup.pos.get(i);
        int rc = tup.RCs.get(i);
        energy += getOneBody(pos, rc);
    }

    // 加上所有pairwise terms
    for (int i = 0; i < tup.pos.size(); i++) {
        for (int j = 0; j < i; j++) {
            int pos1 = tup.pos.get(i), rc1 = tup.RCs.get(i);
            int pos2 = tup.pos.get(j), rc2 = tup.RCs.get(j);
            energy += getPairwise(pos1, rc1, pos2, rc2);
        }
    }

    return energy;
}
```

**注意**：这**不包括**三体项！它只是把已有的one-body和pairwise能量加起来。

---

## 5. Triple Energy数字例子

### 场景

Conformation: `[pos0=RC0, pos1=RC1, pos2=RC2]`

### Step 1: Pairwise估计

从minimizingEmat：

```
One-body energies:
  E_1(0,0) = -2.5 kcal/mol
  E_1(1,1) = -3.7 kcal/mol
  E_1(2,2) = -3.5 kcal/mol

Pairwise energies:
  E_2(0,0,1,1) = -0.3 kcal/mol
  E_2(0,0,2,2) = -0.6 kcal/mol
  E_2(1,1,2,2) = -0.9 kcal/mol

lowerbound = minimizingEmat.getInternalEnergy([0,1,2])
          = (-2.5) + (-3.7) + (-3.5) + (-0.3) + (-0.6) + (-0.9)
          = -11.5 kcal/mol  (不包括E_const)
```

### Step 2: Triple Minimization

```java
minimizingEcalc.calcEnergy([pos0=RC0, pos1=RC1, pos2=RC2])
```

这会：
1. 创建一个molecule，只包含这三个residues的指定rotamers
2. 运行coordinate descent，**同时**优化所有三个residues的DOFs
3. 返回最小化后的能量

**结果**：
```
tripleEnergy = -12.8 kcal/mol  (不包括E_const)
```

### Step 3: Correction

```
correction = tripleEnergy - lowerbound
          = -12.8 - (-11.5)
          = -1.3 kcal/mol  (negative!)
```

**解释**：
- Pairwise估计：-11.5 kcal/mol
- 实际triple能量：-12.8 kcal/mol
- Triple minimization找到了**更低**的能量（多了1.3 kcal/mol的稳定化）
- 这1.3 kcal/mol就是**三体协同效应**的贡献

### Step 4: 存储Correction

```java
correctionMatrix.setHigherOrder(
    tuple = {pos:[0,1,2], RCs:[0,1,2]},
    correction = -1.3
);
```

### Step 5: 后续使用

当MARK*后续需要这个conformation的能量时：

```java
double energyFromMatrix = minimizingEmat.getInternalEnergy([0,1,2]);
double correctedEnergy = energyFromMatrix + correctionMatrix.getHigherOrder([0,1,2]);
                       = -11.5 + (-1.3)
                       = -12.8 kcal/mol  ✓ 正确！
```

---

## 6. Triple Energy在Partition Function中的作用

### Boltzmann Weight计算

```
w(conf) = exp(-E(conf) / RT)
```

**没有correction**：
```
E = -11.5 kcal/mol
w = exp(-(-11.5) / 0.592) = exp(19.4) ≈ 2.7 × 10^8
```

**有correction**：
```
E = -12.8 kcal/mol
w = exp(-(-12.8) / 0.592) = exp(21.6) ≈ 2.4 × 10^9
```

**差异**：
```
Ratio = 2.4×10^9 / 2.7×10^8 = 8.9倍
```

只是1.3 kcal/mol的差异，Boltzmann weight就差了**将近9倍**！

### 对Partition Function的影响

```
Z = Σ exp(-E_i / RT)
```

如果很多conformations都有类似的-1.3 kcal/mol的correction：

```
Z_corrected / Z_uncorrected ≈ 9  (如果所有confs都有类似correction)
```

在实际中，只有部分conformations有significant corrections，但累积效应仍然很大。

---

## 7. Negative vs Positive Corrections

### Negative Correction (更常见)

```
correction < 0
```

**含义**：
- Triple minimization找到了**更低**的能量
- Pairwise估计**过高**（too optimistic as upper bound）
- 三体cooperative stabilization

**例子**：
- 三个residues形成favorable三角hydrogen bond network
- 这种配置在pairwise level无法被完全捕捉

### Positive Correction (较少见)

```
correction > 0
```

**含义**：
- Triple minimization的能量**更高**
- Pairwise估计**过低**（too pessimistic as lower bound）
- 三体steric clash或unfavorable interactions

**例子**：
- 三个large side chains挤在一起
- Pairwise时每对都OK，但三个同时存在时空间不够

**为什么较少见？**
- Pairwise minimization通常已经能识别明显的clashes
- Pruning算法会提前排除这些不利的RCs
- 所以进入MARK*计算的conformations多数是reasonable的

---

## 8. MARK* Bug与Triple Energy的关系

### Bug的表现

**原始代码**（Line 1071-1077）：

```java
double correction = tripleEnergy - lowerbound;
if (correction > 0) {
    correctionMatrix.setHigherOrder(tuple, correction);
}
else {
    System.err.println("Negative correction for "+tuple.stringListing());
    // ❌ 没有设置correction！
}
```

### 后果

当`correction < 0`时（triple energy更低）：
1. **Correction没有被应用**
2. `minimizingEmat.getInternalEnergy(tuple)`返回pairwise估计（过高）
3. Lower bound过高
4. Boltzmann weight被**低估**
5. Partition function Z被**低估**

等等，这不对！测试显示Z被**高估**了10²²倍！

### 矛盾的解释

**关键洞察**：Bug不是在correction阶段，而是在**更早的阶段**！

让我重新分析测试输出：

```
MARK*:
  Q* (lower):           1.406097e+26  ← 这是lower bound!
  Confs Evaluated:      0              ← 没有评估任何conf!
```

**问题**：
- MARK*没有评估任何conformation
- 说明在`runUntilNonZero()`阶段就已经出错了
- Corrections根本没有机会被计算！

**真正的问题在root node initialization**：

在`MARKStarNode.makeRoot()`中计算的initial bounds就已经错了：
```
confUpperBound = rigidScore - negatedHScore
```

这个公式和`negatedHScore = -23940`导致了初始的巨大错误。

Triple energy corrections是**后期的修正机制**，但如果initial bounds就错得离谱，后期的corrections根本来不及应用。

---

## 9. Triple Energy的"Row/Column"结构

Triple energy**不存储在matrix**中！它存储在：

### UpdatingEnergyMatrix的corrections

```java
private TupleTrie corrections;  // Sparse storage
```

**数据结构**：

```
Map<RCTuple, Double>

RCTuple = {
    pos: [pos1, pos2, pos3],  // 3个positions
    RCs: [RC1, RC2, RC3]      // 对应的RCs
}

Double = correction value
```

### 为什么不用flat array？

Triple combinations太多了：

```
对于N个positions，每个有R个RCs：
Total triples = C(N,3) × R³
```

例如：
- 10个positions，平均5个RCs per position
- Triples = C(10,3) × 5³ = 120 × 125 = 15,000个triples

但实际上只有**很少的triples**需要correction（那些tupleBounds > threshold的）。

所以使用**sparse storage**（HashMap或Trie）更高效。

### 访问

```java
// 获取correction
RCTuple tuple = new RCTuple(pos=[0,1,2], RCs=[0,1,2]);
Double correction = correctionMatrix.getHigherOrder(tuple);

if (correction != null) {
    energy += correction;
}
```

---

## 10. 总结

### Triple Energy的本质

```
Triple Energy = "真实的三体最小化能量"

Triple Correction = Triple Energy - Pairwise估计
```

### 为什么重要？

1. **捕捉协同效应**：pairwise近似的不足
2. **提高精度**：partition function计算更准确
3. **选择性计算**：只计算"可疑"的triples（threshold筛选）

### 在MARK*中的作用

```
E_corrected(conf) = E_pairwise(conf) + Σ Correction_triple

Z = Σ exp(-E_corrected(conf) / RT)
```

### 数据结构

- **Pairwise**: 存储在flat 1D arrays（dense）
- **Triple corrections**: 存储在TupleTrie/HashMap（sparse）

### 与Bug的关系

- Negative corrections未被应用 → lower bounds过高
- 但更严重的bug在root node初始化 → 整个算法无法正常运行
- Triple corrections是第二道防线，但第一道防线（initial bounds）已经崩溃

### 物理意义

```
Negative correction (common):
  三体cooperative stabilization
  例如：三角hydrogen bond network

Positive correction (rare):
  三体steric clash
  例如：三个large side chains挤在一起
```

1.3 kcal/mol的correction → Boltzmann weight差9倍！能量微小差异对partition function影响巨大。
