# OSPREY能量矩阵的Row和Column结构详解

## 核心概念

OSPREY的能量矩阵**不是**传统的2D矩阵！它使用**压缩的1D数组**来存储不同阶数的能量项。

---

## 1. 数据结构概览

```java
// TupleMatrixDouble.java (Line 48-49)
private double[] oneBody;   // 1D array for one-body energies
private double[] pairwise;  // 1D array for pairwise energies
```

**关键点**：
- 没有2D或4D数组！
- 所有能量存储在**flat 1D arrays**中
- 通过**索引计算**来访问特定的(position, RC)或(pos1, RC1, pos2, RC2)

---

## 2. One-Body矩阵结构

### 概念上的"Row"和"Column"

虽然物理上是1D数组，概念上可以看作：

```
          RC0    RC1    RC2    RC3    ...
pos0    [e0,0] [e0,1] [e0,2] [e0,3]  ...
pos1    [e1,0] [e1,1] [e1,2]         ...
pos2    [e2,0] [e2,1] [e2,2] [e2,3]  ...
...
```

- **Row**: Position (残基位置)
- **Column**: RC (Rotational Conformer，该位置的rotamer编号)

### 实际存储：连续的1D数组

```java
oneBody[] = [
    // pos0的所有RCs (假设pos0有3个RCs)
    e0,0, e0,1, e0,2,
    // pos1的所有RCs (假设pos1有2个RCs)
    e1,0, e1,1,
    // pos2的所有RCs (假设pos2有4个RCs)
    e2,0, e2,1, e2,2, e2,3,
    ...
]
```

### 索引计算

```java
// AbstractTupleMatrix.java Line 187-189
public int getOneBodyIndex(int res, int conf) {
    return oneBodyOffsets[res] + conf;
}
```

**oneBodyOffsets数组**记录每个position在1D数组中的起始位置：

```
oneBodyOffsets[0] = 0     // pos0从index 0开始
oneBodyOffsets[1] = 3     // pos1从index 3开始 (pos0有3个RCs)
oneBodyOffsets[2] = 5     // pos2从index 5开始 (pos1有2个RCs)
```

### 具体数字例子

**系统配置**：
- pos0: 3个RCs (RC 0, 1, 2)
- pos1: 2个RCs (RC 0, 1)
- pos2: 4个RCs (RC 0, 1, 2, 3)

**能量值（kcal/mol）**：

```
Conceptual Matrix (for visualization):
          RC0     RC1     RC2     RC3
pos0    [-2.5]  [-1.8]  [-3.2]   n/a
pos1    [-4.1]  [-3.7]   n/a     n/a
pos2    [-1.2]  [-2.8]  [-3.5]  [-2.1]
```

**实际存储（1D array）**：

```java
oneBody[] = {
    -2.5,  -1.8,  -3.2,    // pos0, RCs 0,1,2 (indices 0-2)
    -4.1,  -3.7,            // pos1, RCs 0,1   (indices 3-4)
    -1.2,  -2.8,  -3.5, -2.1  // pos2, RCs 0,1,2,3 (indices 5-8)
}

oneBodyOffsets[] = {0, 3, 5}
```

**访问示例**：
```java
// 获取 pos1, RC1 的能量
int index = oneBodyOffsets[1] + 1 = 3 + 1 = 4
double energy = oneBody[4] = -3.7 kcal/mol ✓

// 获取 pos2, RC2 的能量
int index = oneBodyOffsets[2] + 2 = 5 + 2 = 7
double energy = oneBody[7] = -3.5 kcal/mol ✓
```

---

## 3. Pairwise矩阵结构

### 概念上的"维度"

Pairwise能量有**4个索引**：`(pos1, RC1, pos2, RC2)`

概念上这是一个4D张量，但实际存储为1D数组。

### 对称性和存储优化

**重要约定**（AbstractTupleMatrix.java Line 47）：
```java
// tuples are sets not ordered pairs
// E(i_r, j_s) = E(j_s, i_r)
```

**存储规则**：只存储`pos1 > pos2`的pairs，因为能量是对称的。

### 维度示例

假设3个positions：
- pos0 有 3个RCs
- pos1 有 2个RCs
- pos2 有 4个RCs

**需要存储的pairs**：
1. `(pos1, pos0)` - pos1的每个RC × pos0的每个RC = 2×3 = 6个能量值
2. `(pos2, pos0)` - pos2的每个RC × pos0的每个RC = 4×3 = 12个能量值
3. `(pos2, pos1)` - pos2的每个RC × pos1的每个RC = 4×2 = 8个能量值

总共：6 + 12 + 8 = 26个pairwise能量值

### 实际存储

```java
pairwise[] = {
    // (pos1, pos0): pos1的2个RCs × pos0的3个RCs
    e[1,0][0,0], e[1,0][0,1], e[1,0][0,2],  // pos1=RC0 与 pos0=RC0,1,2
    e[1,1][0,0], e[1,1][0,1], e[1,1][0,2],  // pos1=RC1 与 pos0=RC0,1,2

    // (pos2, pos0): pos2的4个RCs × pos0的3个RCs
    e[2,0][0,0], e[2,0][0,1], e[2,0][0,2],  // pos2=RC0 与 pos0=RC0,1,2
    e[2,1][0,0], e[2,1][0,1], e[2,1][0,2],  // pos2=RC1 与 pos0=RC0,1,2
    e[2,2][0,0], e[2,2][0,1], e[2,2][0,2],  // pos2=RC2 与 pos0=RC0,1,2
    e[2,3][0,0], e[2,3][0,1], e[2,3][0,2],  // pos2=RC3 与 pos0=RC0,1,2

    // (pos2, pos1): pos2的4个RCs × pos1的2个RCs
    e[2,0][1,0], e[2,0][1,1],  // pos2=RC0 与 pos1=RC0,1
    e[2,1][1,0], e[2,1][1,1],  // pos2=RC1 与 pos1=RC0,1
    e[2,2][1,0], e[2,2][1,1],  // pos2=RC2 与 pos1=RC0,1
    e[2,3][1,0], e[2,3][1,1],  // pos2=RC3 与 pos1=RC0,1
}
```

### pairwiseOffsets计算

```java
// AbstractTupleMatrix.java Line 119-130
pairwiseOffsets = new int[numPos*(numPos - 1)/2];  // 对于3个pos：3*2/2 = 3

int pairwiseOffset = 0;
int pairwiseIndex = 0;
for (int res1=0; res1<numPos; res1++) {
    for (int res2=0; res2<res1; res2++) {
        pairwiseOffsets[pairwiseIndex++] = pairwiseOffset;
        pairwiseOffset += numConfAtPos[res1]*numConfAtPos[res2];
    }
}
```

**结果**：
```
pairwiseOffsets[0] = 0   // (pos1, pos0) 从 index 0 开始
pairwiseOffsets[1] = 6   // (pos2, pos0) 从 index 6 开始 (pos1有2个RCs×pos0有3个RCs=6)
pairwiseOffsets[2] = 18  // (pos2, pos1) 从 index 18 开始 (pos2有4个RCs×pos0有3个RCs=12)
```

### 索引计算

```java
// AbstractTupleMatrix.java Line 209-224
public int getPairwiseIndex(int res1, int conf1, int res2, int conf2) {
    // Ensure res1 > res2 (swap if necessary)
    if (res2 > res1) {
        swap(res1, res2);
        swap(conf1, conf2);
    }

    // Calculate index
    int pairIndex = res1*(res1 - 1)/2 + res2;  // 找到这对position的offset
    return pairwiseOffsets[pairIndex] + numConfAtPos[res2]*conf1 + conf2;
}
```

**公式解释**：
1. `res1*(res1-1)/2 + res2`: 计算position pair的索引（三角矩阵）
2. `pairwiseOffsets[pairIndex]`: 获取这对position的起始offset
3. `numConfAtPos[res2]*conf1 + conf2`: 在这对position的block内计算具体的RC pair位置

### 具体数字例子

**能量值**（kcal/mol）：

```
(pos1, pos0) pairs:
  (pos1=RC0, pos0=RC0): -0.5
  (pos1=RC0, pos0=RC1): -0.3
  (pos1=RC0, pos0=RC2): -0.4
  (pos1=RC1, pos0=RC0): -0.2
  (pos1=RC1, pos0=RC1): -0.3
  (pos1=RC1, pos0=RC2): -0.1

(pos2, pos0) pairs:
  (pos2=RC0, pos0=RC0): -0.8
  (pos2=RC0, pos0=RC1): -0.6
  (pos2=RC0, pos0=RC2): -0.7
  (pos2=RC1, pos0=RC0): -0.4
  (pos2=RC1, pos0=RC1): -0.5
  (pos2=RC1, pos0=RC2): -0.3
  (pos2=RC2, pos0=RC0): -0.6
  (pos2=RC2, pos0=RC1): -0.7
  (pos2=RC2, pos0=RC2): -0.5
  (pos2=RC3, pos0=RC0): -0.3
  (pos2=RC3, pos0=RC1): -0.4
  (pos2=RC3, pos0=RC2): -0.2

(pos2, pos1) pairs:
  (pos2=RC0, pos1=RC0): -0.7
  (pos2=RC0, pos1=RC1): -0.5
  (pos2=RC1, pos1=RC0): -0.9
  (pos2=RC1, pos1=RC1): -0.7
  (pos2=RC2, pos1=RC0): -1.2
  (pos2=RC2, pos1=RC1): -0.9
  (pos2=RC3, pos1=RC0): -0.6
  (pos2=RC3, pos1=RC1): -0.4
```

**实际1D array**：

```java
pairwise[] = {
    // (pos1, pos0): indices 0-5
    -0.5, -0.3, -0.4,  // pos1=RC0 与 pos0=RC0,1,2
    -0.2, -0.3, -0.1,  // pos1=RC1 与 pos0=RC0,1,2

    // (pos2, pos0): indices 6-17
    -0.8, -0.6, -0.7,  // pos2=RC0 与 pos0=RC0,1,2
    -0.4, -0.5, -0.3,  // pos2=RC1 与 pos0=RC0,1,2
    -0.6, -0.7, -0.5,  // pos2=RC2 与 pos0=RC0,1,2
    -0.3, -0.4, -0.2,  // pos2=RC3 与 pos0=RC0,1,2

    // (pos2, pos1): indices 18-25
    -0.7, -0.5,  // pos2=RC0 与 pos1=RC0,1
    -0.9, -0.7,  // pos2=RC1 与 pos1=RC0,1
    -1.2, -0.9,  // pos2=RC2 与 pos1=RC0,1
    -0.6, -0.4,  // pos2=RC3 与 pos1=RC0,1
}

pairwiseOffsets[] = {0, 6, 18}
```

**访问示例**：

```java
// 获取 (pos2=RC2, pos1=RC1) 的pairwise能量

// Step 1: Ensure pos2 > pos1 (already satisfied: 2 > 1)

// Step 2: 计算position pair的索引
int pairIndex = 2*(2-1)/2 + 1 = 2*1/2 + 1 = 1 + 1 = 2

// Step 3: 获取起始offset
int offset = pairwiseOffsets[2] = 18

// Step 4: 在block内计算RC pair的位置
// numConfAtPos[res2=1] = 2 (pos1有2个RCs)
// conf1=2 (pos2=RC2), conf2=1 (pos1=RC1)
int localIndex = numConfAtPos[1]*2 + 1 = 2*2 + 1 = 5

// Step 5: 最终索引
int finalIndex = offset + localIndex = 18 + 5 = 23

// 获取能量
double energy = pairwise[23] = -0.9 kcal/mol ✓
```

让我们验证：
- (pos2, pos1) block从index 18开始
- pos2=RC2的sub-block: indices 22-23 (每个sub-block有2个元素，因为pos1有2个RCs)
- pos2=RC2, pos1=RC1: index 23 ✓

---

## 4. 可视化总结

### One-Body Matrix (概念上)

```
Rows = Positions (残基位置)
Columns = RCs (rotamers)

          RC0    RC1    RC2    RC3
       ┌──────┬──────┬──────┬──────┐
pos0   │ -2.5 │ -1.8 │ -3.2 │  --  │
       ├──────┼──────┼──────┼──────┤
pos1   │ -4.1 │ -3.7 │  --  │  --  │
       ├──────┼──────┼──────┼──────┤
pos2   │ -1.2 │ -2.8 │ -3.5 │ -2.1 │
       └──────┴──────┴──────┴──────┘
```

**物理存储**：`[-2.5, -1.8, -3.2, -4.1, -3.7, -1.2, -2.8, -3.5, -2.1]`

### Pairwise Matrix (概念上)

这是一个4D tensor: `[pos1][RC1][pos2][RC2]`

但由于对称性，只存储`pos1 > pos2`的部分。

对于每个position pair (pos_i, pos_j) where i > j：
- 这是一个2D block
- Rows: pos_i的RCs
- Columns: pos_j的RCs

```
(pos1, pos0) block:  pos1有2个RCs, pos0有3个RCs
         pos0_RC0  pos0_RC1  pos0_RC2
       ┌─────────┬─────────┬─────────┐
pos1_RC0│  -0.5   │  -0.3   │  -0.4   │
       ├─────────┼─────────┼─────────┤
pos1_RC1│  -0.2   │  -0.3   │  -0.1   │
       └─────────┴─────────┴─────────┘

(pos2, pos0) block:  pos2有4个RCs, pos0有3个RCs
         pos0_RC0  pos0_RC1  pos0_RC2
       ┌─────────┬─────────┬─────────┐
pos2_RC0│  -0.8   │  -0.6   │  -0.7   │
       ├─────────┼─────────┼─────────┤
pos2_RC1│  -0.4   │  -0.5   │  -0.3   │
       ├─────────┼─────────┼─────────┤
pos2_RC2│  -0.6   │  -0.7   │  -0.5   │
       ├─────────┼─────────┼─────────┤
pos2_RC3│  -0.3   │  -0.4   │  -0.2   │
       └─────────┴─────────┴─────────┘

(pos2, pos1) block:  pos2有4个RCs, pos1有2个RCs
         pos1_RC0  pos1_RC1
       ┌─────────┬─────────┐
pos2_RC0│  -0.7   │  -0.5   │
       ├─────────┼─────────┤
pos2_RC1│  -0.9   │  -0.7   │
       ├─────────┼─────────┤
pos2_RC2│  -1.2   │  -0.9   │
       ├─────────┼─────────┤
pos2_RC3│  -0.6   │  -0.4   │
       └─────────┴─────────┘
```

**物理存储**：所有这些blocks的元素按row-major order连接成一个1D数组。

---

## 5. 为什么使用这种结构？

### 优点

1. **内存效率**：
   - 不存储对称的部分 (E[i,j] = E[j,i])
   - 不存储不可能的组合
   - 可变长度的RCs per position

2. **缓存友好**：
   - 连续内存访问
   - 对于iteration over all energies很高效

3. **灵活性**：
   - 每个position可以有不同数量的RCs
   - 容易扩展到higher-order terms

### 缺点

1. **索引复杂**：
   - 需要offset计算
   - 不直观

2. **调试困难**：
   - 不能直接"看到"矩阵
   - 需要理解索引映射

---

## 6. Higher-Order Terms

对于3-body或更高阶的能量项（如`UpdatingEnergyMatrix`的corrections），OSPREY使用：

1. **TupleTrie** 或 **HigherTupleFinder**: sparse数据结构
2. **RCTuple作为key**: `{pos: [0,1,2], RCs: [2,1,3]}`
3. **Double作为value**: correction能量

这些不存储在flat array中，而是使用HashMap或Trie结构，因为它们非常稀疏。

---

## 7. 总结

### One-Body矩阵
- **"Rows"**: Positions (残基位置)
- **"Columns"**: RCs (rotamers)
- **物理存储**: 按position顺序连接，每个position的所有RCs连续存储
- **访问**: `oneBody[oneBodyOffsets[pos] + rc]`

### Pairwise矩阵
- **"Dimensions"**: (pos1, RC1, pos2, RC2) - 4D tensor
- **只存储**: pos1 > pos2 的pairs（利用对称性）
- **每个pair block**: pos1的RCs × pos2的RCs 的2D矩阵
- **物理存储**: 所有blocks按position pair顺序连接，每个block内按row-major order
- **访问**: 复杂的offset计算，考虑position pair index和RC pair index

### 关键记忆点
```
oneBody[position][RC] → 实际是 oneBody[offset_pos + RC]
pairwise[pos1][RC1][pos2][RC2] → 实际是 pairwise[offset_pair + RC1*numRC2 + RC2]
```

这种压缩表示法在大规模蛋白质设计中节省大量内存，但代价是索引计算的复杂性。
