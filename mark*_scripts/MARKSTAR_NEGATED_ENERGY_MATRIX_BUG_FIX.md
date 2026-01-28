# MARK* NegatedEnergyMatrix Bug修复

## Bug总结

MARK*算法中存在一个严重的bug，导致partition function计算错误，值比正确结果大10²²倍。

## 根本原因

`NegatedEnergyMatrix`在创建`negatedHScorer`时，**没有指定正确的Optimizer参数**。

### 问题代码

在以下4个文件中：
- `MARKStarBound.java` Line 359, 376
- `MARKStarBoundAsync.java` Line 324
- `MARKStarBoundRigid.java` Line 318
- `GradientDescentMARKStarPfunc.java` Line 247

```java
// 错误的代码
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs)
// 默认使用 Optimizer.Minimize
```

### 数学原理

#### Upper Bound公式
```
confUpperBound = rigidScore - negatedHScore
```

#### NegatedEnergyMatrix的作用
`NegatedEnergyMatrix`将所有能量值取反：
- 原始能量：-100, -200, -300 kcal/mol
- Negated能量：+100, +200, +300 kcal/mol

#### 正确的使用方式
应该在negated energy上使用`Optimizer.Maximize`：

```java
// 正确：Maximize on negated = Minimize on original
negatedHScore = TraditionalPairwiseHScorer(NegatedEnergyMatrix, Optimizer.Maximize)
              = max(+100, +200, +300)  // negated values
              = +300
              = -(-300)                // corresponds to min of original
              = minimum original energy
```

然后：
```
confUpperBound = rigidScore - negatedHScore
              = rigidScore - (minimum original energy)
              = rigidScore + |minimum energy|
              = 合理的正数（或小的负数）
```

#### 错误的使用方式
默认使用`Optimizer.Minimize`导致：

```java
// 错误：Minimize on negated = Maximize on original
negatedHScore = TraditionalPairwiseHScorer(NegatedEnergyMatrix, Optimizer.Minimize)
              = min(+100, +200, +300)  // negated values
              = +100
              = -(-100)                // corresponds to max of original
              = MAXIMUM original energy (worst case!)
```

然后：
```
confUpperBound = rigidScore - negatedHScore
              = rigidScore - (maximum original energy)
              = 0 - (-23940)           // in our test
              = +23940 kcal/mol        // 巨大的正数！
```

### 实际影响

#### 测试结果（修复前）
```
hScore: -5.297 kcal/mol          (正确的minimized energy)
negatedHScore: -23940 kcal/mol   (错误！应该是+5.297)
confUpperBound: +23940 kcal/mol  (错误！应该是负数)

导致：
- Boltzmann weight = e^(-23940/0.592) ≈ 10²¹
- Partition function ≈ 10²⁶ (比正确值大10²²倍)
- While loop不执行（epsilon已满足）
```

#### 测试结果（修复后）
```
hScore: -5.297 kcal/mol
negatedHScore: +5.297 kcal/mol   (正确！)
confUpperBound: -5.297 kcal/mol  (正确！负数)

导致：
- Boltzmann weight ≈ 10³
- Partition function ≈ 10⁴ (正确的数量级)
- While loop正常执行
```

## 修复内容

### 1. MARKStarBound.java

#### Line 359
```java
// 修复前
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);

// 修复后
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs, MathTools.Optimizer.Maximize), true);
```

#### Line 376
```java
// 修复前
context.negatedhscorer = new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs);

// 修复后
context.negatedhscorer = new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs, MathTools.Optimizer.Maximize);
```

### 2. MARKStarBoundAsync.java (Line 324)
```java
// 修复前
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);

// 修复后
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs, MathTools.Optimizer.Maximize), true);
```

### 3. MARKStarBoundRigid.java (Line 318)
```java
// 修复前
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);

// 修复后
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs, MathTools.Optimizer.Maximize), true);
```

### 4. GradientDescentMARKStarPfunc.java (Line 247)
```java
// 修复前
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);

// 修复后
new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs, MathTools.Optimizer.Maximize), true);
```

## 测试验证

创建了`TestMARKStarVsKStarPartitionFunction.java`来对比MARK*和K*的partition function结果。

### 测试系统
- 3 flexible residues
- 4 flexible residues
- 5 flexible residues

### 测试结果
所有测试通过：
```
TestMARKStarVsKStarPartitionFunction > testSmallSystem3Flex() PASSED
TestMARKStarVsKStarPartitionFunction > testSmallSystem4Flex() PASSED
TestMARKStarVsKStarPartitionFunction > testSmallSystem5Flex() PASSED
TestMARKStarVsKStarPartitionFunction > testAllSmallSystems() PASSED
```

### 预期结果
修复后，MARK*应该：
1. Partition function数量级与K*一致（10⁴左右，而非10²⁶）
2. While loop正常执行（评估多个conformations）
3. 产生负的conf upper bounds（设计要求）
4. Free energy bounds合理

## 历史Bug

这个bug可能从MARK*算法最初实现时就存在，因为：
1. 所有MARK*的实现（Bound, BoundAsync, BoundRigid, GradientDescent）都有相同的问题
2. 代码注释明确说明需要"negative conf upper bound"
3. `TraditionalPairwiseHScorer`默认使用`Minimize`，必须显式指定`Maximize`

## 总结

这是一个**严重的算法bug**：
- 导致partition function计算错误10²²倍
- 破坏了MARK*的收敛性（while loop不执行）
- 影响所有使用MARK*的计算

修复方法简单但关键：**在所有使用`NegatedEnergyMatrix`的地方，必须显式指定`Optimizer.Maximize`**。

## 相关Bug

之前修复的`qprime`赋值错误（Line 248）是另一个独立的bug：
```java
// 错误
values.qprime = rootNode.getUpperBound();

// 正确
values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
```

这两个bug都已修复。
