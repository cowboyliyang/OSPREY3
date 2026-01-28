# Reference Energies详解

## 1. 什么是Reference Energy？

**Reference Energy (E_ref)** 是一个**能量归一化**技巧，用于：
1. 使不同amino acid types之间的能量**可比较**
2. 减少能量矩阵中的**数值范围**
3. 提高A*搜索的**效率**

---

## 2. 核心概念

### 问题：不同AA类型的能量baseline不同

考虑position 0可以是PHE或ALA：

```
PHE的所有rotamers:
  RC0 (PHE, rotamer 1): E_intra = -50.0 kcal/mol
  RC1 (PHE, rotamer 2): E_intra = -48.5 kcal/mol
  RC2 (PHE, rotamer 3): E_intra = -49.2 kcal/mol

ALA的所有rotamers:
  RC3 (ALA, rotamer 1): E_intra = -10.0 kcal/mol
  RC4 (ALA, rotamer 2): E_intra = -9.5 kcal/mol
```

**观察**：
- PHE的能量都在-50左右（大side chain，很多atoms）
- ALA的能量都在-10左右（小side chain，少atoms）

**问题**：
- 如果直接使用这些能量，PHE会**总是**被A*优先选择（因为能量更低）
- 但这个差异**不是因为PHE更好**，而是因为PHE simply有更多atoms！

### 解决方案：Reference Energy Normalization

**思路**：对每个AA type，找到它的"baseline"能量，然后减去。

```
E_corrected(RC) = E_original(RC) - E_ref(AA_type)
```

---

## 3. Reference Energy的计算

### 算法（SimplerEnergyMatrixCalculator.java Line 535-569）

```java
public SimpleReferenceEnergies calcReferenceEnergies() {
    SimpleReferenceEnergies eref = new SimpleReferenceEnergies();

    // 对每个position的每个conformation
    for (int pos : positions) {
        for (int conf : conformations[pos]) {

            // 计算intra-residue energy (不包括与其他residues的interaction)
            double intraE = confEcalc.calcIntraEnergy(pos, conf);

            String resType = getResidueType(pos, conf);  // e.g., "PHE", "ALA"

            // 对每个(pos, resType)组合，保留最小的intra energy
            Double currentMin = eref.get(pos, resType);
            if (currentMin == null || intraE < currentMin) {
                eref.set(pos, resType, intraE);
            }
        }
    }

    return eref;
}
```

### 具体例子

**Position 0的RCs**:
```
RC0 (PHE, rot1): intraE = -50.0 kcal/mol
RC1 (PHE, rot2): intraE = -48.5 kcal/mol
RC2 (PHE, rot3): intraE = -49.2 kcal/mol
RC3 (ALA, rot1): intraE = -10.0 kcal/mol
RC4 (ALA, rot2): intraE = -9.5 kcal/mol
```

**计算Reference Energies**:
```
E_ref(pos=0, AA=PHE) = min(-50.0, -48.5, -49.2) = -50.0 kcal/mol
E_ref(pos=0, AA=ALA) = min(-10.0, -9.5) = -10.0 kcal/mol
```

**物理意义**：
- E_ref是这个AA type在这个position上**最稳定的构象**的内部能量
- 作为这个AA type的"baseline"

---

## 4. Energy Matrix Correction

### 应用Reference Energies（ReferenceEnergies.java Line 96-117）

```java
void correctEnergyMatrix(EnergyMatrix emat) {
    for (int pos : positions) {
        for (int rc : RCs[pos]) {
            String AAType = getAAType(pos, rc);
            double eRef = eRefMatrix.get(pos).get(AAType);

            // 从one-body能量中减去reference
            double E_uncorrected = emat.getOneBody(pos, rc);
            double E_corrected = E_uncorrected - eRef;
            emat.setOneBody(pos, rc, E_corrected);
        }
    }
}
```

### 数字例子

**原始one-body energies**（包括intra + shell interactions）:
```
RC0 (PHE, rot1): E_one = -52.5 kcal/mol  (intra=-50.0, shell=-2.5)
RC1 (PHE, rot2): E_one = -50.3 kcal/mol  (intra=-48.5, shell=-1.8)
RC2 (PHE, rot3): E_one = -51.8 kcal/mol  (intra=-49.2, shell=-2.6)
RC3 (ALA, rot1): E_one = -11.2 kcal/mol  (intra=-10.0, shell=-1.2)
RC4 (ALA, rot2): E_one = -10.8 kcal/mol  (intra=-9.5, shell=-1.3)
```

**Correction**:
```
RC0: E_corrected = -52.5 - (-50.0) = -2.5 kcal/mol
RC1: E_corrected = -50.3 - (-50.0) = -0.3 kcal/mol
RC2: E_corrected = -51.8 - (-50.0) = -1.8 kcal/mol
RC3: E_corrected = -11.2 - (-10.0) = -1.2 kcal/mol
RC4: E_corrected = -10.8 - (-10.0) = -0.8 kcal/mol
```

**结果**：
- 所有corrected energies都在**相似的范围**（-3到0 kcal/mol）
- PHE和ALA现在**可比较**了
- 最好的PHE rotamer (RC0: -2.5) vs 最好的ALA (RC3: -1.2)
- 现在差异只反映了**relative stability**，不是absolute atom count

---

## 5. 为什么需要Reference Energies？

### 好处1：公平比较不同AA types

**没有E_ref**:
```
A*搜索会看到：
  PHE: -50 kcal/mol  ← 总是选这个！
  ALA: -10 kcal/mol
```

**有E_ref**:
```
A*搜索会看到：
  PHE: -2.5 kcal/mol  (相对于PHE的baseline)
  ALA: -1.2 kcal/mol  (相对于ALA的baseline)

现在可以根据实际的local environment来选择！
```

### 好处2：减小数值范围

**没有E_ref**:
```
oneBody energies: -100 到 -5 kcal/mol  (range: 95)
```

**有E_ref**:
```
corrected energies: -5 到 +2 kcal/mol  (range: 7)
```

更小的数值范围 → 更好的numerical stability

### 好处3：提高A*效率

A*的heuristic会更准确：
- 不会因为某个AA type inherently有更低的能量而偏向它
- 只考虑**相对优势**，not absolute baseline

---

## 6. Reference Energies在完整能量计算中的角色

### Energy Matrix中存储的是Corrected Energies

```
EnergyMatrix.oneBody[pos][rc] = E_intra + E_shell - E_ref
```

### 计算真实conformation能量时需要加回E_ref

```java
// ReferenceEnergies.java Line 120-128
public double confERef(int[] conf) {
    double totERef = 0;
    for (int pos=0; pos<numPos; pos++) {
        String AAType = getAAType(pos, conf[pos]);
        totERef += eRefMatrix.get(pos).get(AAType);
    }
    return totERef;
}
```

**完整公式**:
```
E_total(conf) = E_matrix(conf) + E_ref(conf)

其中:
  E_matrix(conf) = Σ (E_one_corrected + E_pair)  // 从矩阵读取
  E_ref(conf) = Σ E_ref(pos, AA_type)            // 加回reference
```

### 数字例子

**Conformation**: `[RC0 (PHE), RC3 (ALA)]`

**从能量矩阵计算**:
```
E_matrix = oneBody_corrected[0][RC0] + oneBody_corrected[1][RC3] + pairwise[0,0][1,3]
         = (-2.5) + (-1.2) + (-0.5)
         = -4.2 kcal/mol
```

**加回reference energies**:
```
E_ref = E_ref(pos=0, PHE) + E_ref(pos=1, ALA)
      = (-50.0) + (-10.0)
      = -60.0 kcal/mol
```

**真实总能量**:
```
E_total = E_matrix + E_ref
        = -4.2 + (-60.0)
        = -64.2 kcal/mol
```

---

## 7. 测试输出解读

### 从您的输出：

```
Calculating reference energies for 253 residue confs...
Finished in 75.8 ms
```

**这在做什么？**

1. **253 residue confs**: Complex state有253个RCs
2. 对每个RC计算**intra-residue energy**（不包括pairwise）
3. 对每个(position, AA_type)组合，找到最小的intra energy
4. 存储为E_ref

**时间很短（75.8 ms）** 因为：
- 只计算intra energy（单个residue内部）
- 不涉及pairwise interactions
- 相对简单的计算

### 为什么三个states都计算？

```
Complex:  253 residue confs → 75.8 ms
Protein:  110 residue confs → 13.5 ms
Ligand:   143 residue confs → 28.4 ms
```

因为：
- 每个state有不同的conf space
- 每个state需要自己的reference energies
- Complex包含protein+ligand的所有flexible positions
- Protein/Ligand各自只包含自己的flexible positions

---

## 8. Reference Energies vs Energy Matrix

### 计算顺序

```
1. Calculating reference energies...     ← 先计算E_ref
   (fast, only intra energies)

2. Calculating energy matrix...          ← 然后计算完整矩阵
   (slow, includes pairwise)

   在这个过程中：
   oneBody[pos][rc] -= E_ref[pos][AA_type]  ← 自动correction
```

### 数据结构

```
ReferenceEnergies:
  Map<(position, AA_type), Double>

  Example:
    (pos=0, "PHE") → -50.0
    (pos=0, "ALA") → -10.0
    (pos=1, "ASP") → -30.2
    ...

EnergyMatrix (after correction):
  oneBody[pos][rc] = corrected energies
  pairwise[pos1][rc1][pos2][rc2] = pairwise (unchanged)
```

### Pairwise不需要correction

**为什么？**

Pairwise能量是**两个residues之间的interaction**，不涉及单个residue的"baseline"：

```
E_pair(PHE-ASP) = -0.5 kcal/mol

这个值已经是"relative"的：
- 相对于两个residues不interact的情况
- 不需要进一步normalization
```

---

## 9. 关键总结

### Reference Energy的本质

```
E_ref(pos, AA_type) = min(intra energies of all rotamers of that AA type at that position)
```

**作用**：
1. ✅ **Normalize** 不同AA types的能量
2. ✅ **Fair comparison** 在protein design中
3. ✅ **Numerical stability** 减小数值范围
4. ✅ **A* efficiency** 更准确的heuristic

### 在能量计算中的使用

```
矩阵中存储：E_one_corrected = E_one_original - E_ref
A*搜索时：使用E_one_corrected（公平比较）
最终结果时：E_total = E_matrix + E_ref（真实能量）
```

### 与其他概念的区别

```
Reference Energy:  每个(pos, AA_type)的baseline
One-Body Energy:   单个RC的能量（intra + shell）
Pairwise Energy:   两个RCs之间的interaction
Constant Term:     整个系统的template能量
```

### 输出解读

```
"Calculating reference energies for 253 residue confs..."
= 计算253个RCs的intra energies
= 对每个(pos, AA_type)保留最小值
= 作为normalization baseline
= 快速（只有intra，不需要pairwise）
```

---

## 10. 为什么测试输出中先计算Reference Energies？

### 工作流程

```
1. Calculating reference energies
   ↓
2. Calculating minimizing energy matrix
   ├─ 计算时自动应用E_ref correction
   └─ 存储corrected one-body energies
   ↓
3. Calculating rigid energy matrix
   ├─ 同样应用E_ref correction
   └─ 两个矩阵使用相同的E_ref
```

### 重要性

**Reference energies是foundation**：
- 必须先计算，才能正确normalize能量矩阵
- 影响后续所有的A*搜索和partition function计算
- 确保不同AA types的能量可比较

没有reference energies → protein design算法会偏向large residues（因为它们有更多favorable interactions的机会）

有reference energies → 算法公平评估每个residue type的**相对优势**
