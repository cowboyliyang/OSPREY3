# 如何发现 qprime 赋值错误

## 🔍 发现过程

### 第1步：查看测试代码

在测试代码 `TestMARKStarVsKStarPartitionFunction.java` 中，我看到：

```java
// K*的结果处理 (Line 258)
result.pstar = pfuncResult.values.qprime.add(pfuncResult.values.qstar);

// MARK*的结果处理 (Line 298)
result.pstar = pfuncResult.values.pstar;
```

**疑问1**：为什么K*需要 `qprime + qstar`，而MARK*直接用`pstar`？

### 第2步：查看 PartitionFunction.Values 的定义

**文件**: `PartitionFunction.java` Line 74-78

```java
public static class Values {
    public BigDecimal qstar;   // pfunc value of all evaluated confs
    public BigDecimal qprime;  // pfunc value of all unpruned, but unevaluated confs
    public BigDecimal pstar;   // pfunc value of all pruned confs
```

**关键发现**：
- `qstar` = 已评估构象的partition function值
- `qprime` = **未评估但未pruned构象的partition function值**（gap!）
- `pstar` = 已pruned构象的partition function值

**calcUpperBound()的实现** (Line 143-148):
```java
public BigDecimal calcUpperBound() {
    return new BigMath(decimalPrecision)
        .set(qstar)
        .add(qprime)     // ← qprime是gap
        .add(pstar)
        .get();
}
```

所以：**upper bound = qstar + qprime + pstar**

### 第3步：查看 K* (GradientDescentPfunc) 如何赋值

**文件**: `GradientDescentPfunc.java` Line 485-489

```java
values.qstar = state.getLowerBound();
values.qprime = bigMath()
    .set(state.getUpperBound())
    .sub(state.getLowerBound())    // ← qprime = upper - lower (gap!)
    .get();
```

✅ **K*正确**: qprime = gap

### 第4步：查看 MARK* (MARKStarBound) 如何赋值

**文件**: `MARKStarBound.java` Line 246-248

```java
values.pstar = rootNode.getUpperBound();
values.qstar = rootNode.getLowerBound();
values.qprime = rootNode.getUpperBound();  // ❌ 错误！
```

❌ **MARK*错误**: qprime = upper bound（应该是gap）

### 第5步：其他实现验证

**SimplePartitionFunction.java** Line 215:
```java
values.qprime = upperBound.totalBound.subtract(lowerBound.weightedScoreSum);
```
✅ qprime = gap

**ParallelConfPartitionFunction.java** Line 316-317:
```java
qprimeUnevaluated = qprimeUnevaluated.subtract(boltzmann.calc(econf.getScore()));
return qprimeUnevaluated.add(qprimeUnscored);
```
✅ qprime = 未评估的权重和

---

## 📊 数据验证

### MARK*当前的错误赋值

```java
qstar = 1.41×10²⁶  (lower bound)
qprime = 1.41×10²⁶ (upper bound) ❌ 应该是gap
pstar = 1.41×10²⁶  (upper bound)
```

### 测试代码计算pstar

```java
// K*:
pstar = qprime + qstar
      = gap + lower
      = (upper - lower) + lower
      = upper ✅

// MARK*:
pstar = pstar字段的值
      = upper ✅
```

所以在我们的测试中，**pstar的值其实是对的**，但qprime的值是错的。

### 影响分析

#### 1. **calcUpperBound() 会算错**

```java
// 当前错误的计算：
upperBound = qstar + qprime + pstar
           = lower + upper + upper  ❌
           = 1.41×10²⁶ + 1.41×10²⁶ + 1.41×10²⁶
           = 4.23×10²⁶

// 正确的计算：
upperBound = qstar + qprime + pstar
           = lower + gap + pruned
           = lower + (upper - lower) + 0
           = upper ✅
```

#### 2. **getEffectiveEpsilon() 会算错**

从`PartitionFunction.java` Line 129-131:
```java
BigDecimal s = MathTools.bigAdd(qprime, pstar, decimalPrecision);
BigDecimal qu = MathTools.bigAdd(s, qstar, decimalPrecision);
double delta = MathTools.bigDivide(s, qu, decimalPrecision).doubleValue();
```

当前错误：
```
s = qprime + pstar
  = upper + upper = 2 × upper
qu = s + qstar
   = 2×upper + lower
delta = s / qu
      = 2×upper / (2×upper + lower)
      ≈ 2/(2+0) = 1.0
```

正确应该：
```
s = gap + pruned
  = (upper - lower) + 0
  = gap
qu = gap + lower = upper
delta = gap / upper
      = epsilon value ✅
```

---

## 🐛 Bug的后果

### 1. **Epsilon计算错误**

如果qprime被错误赋值为upper bound：
- delta会接近1.0（而不是实际的epsilon）
- 这可能导致算法认为bounds已经足够紧

### 2. **可能影响while循环判断**

虽然while循环直接用`epsilonBound`而不是`getEffectiveEpsilon()`，但如果其他地方使用了epsilon值，可能会影响算法行为。

### 3. **Free energy计算**

Free energy使用的是`calcFreeEnergyBounds()`，它调用：
```java
public MathTools.DoubleBounds calcFreeEnergyBounds() {
    return new MathTools.DoubleBounds(
        calcLowerBoundFreeEnergy(),  // 用qstar
        calcUpperBoundFreeEnergy()   // 用calcUpperBound()
    );
}
```

如果`calcUpperBound()`错了，upper bound的free energy也会错。

---

## 🔍 为什么我怀疑这个？

### 线索1：命名不一致

K*的代码风格：
```java
values.qstar = lower;
values.qprime = upper - lower;  // 明确是gap
```

MARK*的代码风格：
```java
values.qstar = lower;
values.qprime = upper;  // 看起来像复制粘贴错误
values.pstar = upper;
```

**明显的复制粘贴错误！** qprime和pstar赋了同样的值。

### 线索2：API文档

`PartitionFunction.java`中明确写道：
```java
public BigDecimal qprime;  // pfunc value of all unpruned, but unevaluated confs
```

"unevaluated confs" 应该是一个**增量值**（gap），不是总的upper bound。

### 线索3：其他实现都是gap

所有其他partition function实现（SimplePartitionFunction, GradientDescentPfunc, ParallelConfPartitionFunction）都计算qprime为gap或未评估部分的和，**没有一个**直接赋值为upper bound。

---

## ✅ 修复方案

### 修复代码

**文件**: `MARKStarBound.java` Line 248

```java
// 从：
values.qprime = rootNode.getUpperBound();

// 改为：
values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
```

### 验证修复

修复后，检查：

1. **qprime的语义**：
```java
qprime = gap = upper - lower ✅
```

2. **calcUpperBound()**：
```java
upperBound = qstar + qprime + pstar
           = lower + gap + 0
           = lower + (upper - lower)
           = upper ✅
```

3. **getEffectiveEpsilon()**：
```java
s = qprime + pstar = gap + 0 = gap
qu = s + qstar = gap + lower = upper
delta = gap / upper ✅
```

---

## 🤔 这能解释10²²倍差异吗？

### 回答：不能完全解释

虽然qprime赋值错误，但：

1. **测试代码中我们用的是`pstar`字段**，不是`qprime`：
```java
result.pstar = pfuncResult.values.pstar;  // 这个值是对的
```

2. **主要问题仍然是while循环不执行**，导致：
   - 没有实际计算构象
   - bounds是初始值
   - 这些初始值才是10²²倍差异的根源

3. **但这个bug仍然需要修复**，因为：
   - 影响epsilon计算
   - 影响`calcUpperBound()`
   - 可能在其他使用场景中导致问题

---

## 📝 总结

### 我是如何发现的

1. ✅ **查看API文档** - qprime被定义为"unevaluated confs"的值
2. ✅ **对比K*实现** - K*明确计算qprime = upper - lower
3. ✅ **检查其他实现** - 所有实现都把qprime当作gap
4. ✅ **看到明显的复制粘贴** - qprime和pstar赋了同样的值
5. ✅ **理解partition function语义** - qprime应该是增量，不是总和

### Bug的优先级

🟡 **Medium** - 需要修复但不是主要问题

**原因**：
- ✅ 这是一个真实的bug
- ✅ 影响epsilon和upper bound计算
- ❌ 但不能解释10²²倍差异（主因是while循环不执行）
- ✅ 修复很简单（一行代码）

### 修复顺序

1. 🔴 **Critical**: 先修复while循环问题（主要bug）
2. 🟡 **Medium**: 然后修复qprime赋值（这个bug）
3. ✅ **Verify**: 运行测试确认两个修复都生效

---

## 附录：完整的Values语义

根据`PartitionFunction.java`的定义：

```
Partition Function Z 的分解：

Z = Z_evaluated + Z_unevaluated + Z_pruned
    ↓              ↓                ↓
  qstar          qprime           pstar

其中：
- qstar: 已经实际计算能量的构象的Boltzmann权重和
- qprime: 还没计算但在搜索空间中的构象的权重估计
- pstar: 被pruning算法剪枝的构象的权重估计

Bounds：
- Lower bound = qstar (确定的部分)
- Upper bound = qstar + qprime + pstar (所有可能的)

Epsilon：
- epsilon = (upper - lower) / upper
         = (qprime + pstar) / (qstar + qprime + pstar)
```

MARK*当前的错误：把qprime设为upper而不是gap，破坏了这个分解语义。
