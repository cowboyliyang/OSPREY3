# MARK* 修复计划

**备份位置**: `backups/markstar_before_fix_20251107_154200/`

---

## 问题总结

当前MARK*的主要问题：
1. ❌ **过早终止** - 只运行`runUntilNonZero()`，没有进入主循环
2. ❌ **Partition function异常大** - 1.41×10²⁶（应该≈10⁴）
3. ❌ **Free energy异常低** - -35.7 kcal/mol（应该≈-5.6）
4. ❌ **0个构象被评估** - 没有实际能量计算

---

## 需要修复的位置

### 🔴 **修复1: qprime赋值错误**

**文件**: `src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java`

**位置**: Line 248

**当前代码**:
```java
values.pstar = rootNode.getUpperBound();
values.qstar = rootNode.getLowerBound();
values.qprime = rootNode.getUpperBound();  // ❌ 错误！应该是gap
```

**问题**:
- `qprime`应该是gap（upper - lower），但当前代码赋值为upper bound
- 这导致后续使用`qprime`的地方都出错

**修复方案**:
```java
values.pstar = rootNode.getUpperBound();
values.qstar = rootNode.getLowerBound();
values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());  // ✅ 正确
```

**影响范围**:
- 影响partition function的bounds计算
- 影响free energy计算
- 但这不能解释10²²倍的差异

**优先级**: 🟡 Medium（需要修复但不是主要问题）

---

### 🔴 **修复2: While循环条件问题**

**文件**: `src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java`

**位置**: Line 222-235

**当前代码**:
```java
if(!nonZeroLower) {
    runUntilNonZero();
    updateBound();
}
while (epsilonBound > targetEpsilon &&           // ← 这个条件可能立即为false
       workDone()-previousConfCount < maxNumConfs
       && isStable(stabilityThreshold)) {
    debugPrint("Tightening from epsilon of "+epsilonBound);
    if(debug)
        debugHeap(queue);
    tightenBoundInPhases();
    debugPrint("Errorbound is now "+epsilonBound);
    ...
}
```

**问题分析**:

从测试输出，我们看到：
```
Running until leaf is found...
Found a leaf!
[立即停止]
```

这说明`runUntilNonZero()`执行后，while循环没有进入。可能原因：

1. **epsilonBound初始值问题**
   - `updateBound()`计算的epsilon可能已经≤0.10
   - 需要检查`rootNode.computeEpsilonErrorBounds()`的实现

2. **runUntilNonZero()的问题**
   ```java
   private void runUntilNonZero() {
       System.out.println("Running until leaf is found...");
       // ...
       boundLowestBoundConfUnderNode(rootNode,newNodes);
       queue.addAll(newNodes);
       newNodes.clear();
       System.out.println("Found a leaf!");
       nonZeroLower = true;
   }
   ```

   这个方法只是找到一个leaf就返回了！没有实际计算任何能量。

**需要调查**:
1. `epsilonBound`在`updateBound()`后的值
2. `rootNode.computeEpsilonErrorBounds()`的实现
3. 为什么找到一个leaf就认为完成了？

**调试方案**:
```java
// 在updateBound()后添加
System.out.println("After runUntilNonZero:");
System.out.println("  epsilonBound = " + epsilonBound);
System.out.println("  targetEpsilon = " + targetEpsilon);
System.out.println("  lowerBound = " + rootNode.getLowerBound());
System.out.println("  upperBound = " + rootNode.getUpperBound());
```

**优先级**: 🔴 Critical（这是主要问题！）

---

### 🔴 **修复3: runUntilNonZero()逻辑问题**

**文件**: `src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java`

**位置**: Line 460-476

**当前代码**:
```java
private void runUntilNonZero() {
    System.out.println("Running until leaf is found...");
    double bestConfUpper = Double.POSITIVE_INFINITY;

    List<MARKStarNode> newNodes = new ArrayList<>();
    List<MARKStarNode> leafNodes = new ArrayList<>();
    int numNodes = 0;
    Stopwatch leafLoop = new Stopwatch().start();
    Stopwatch overallLoop = new Stopwatch().start();
    boundLowestBoundConfUnderNode(rootNode,newNodes);
    queue.addAll(newNodes);

    newNodes.clear();
    System.out.println("Found a leaf!");
    nonZeroLower = true;
}
```

**问题**:
1. **没有实际计算任何能量** - 只是展开树节点
2. **没有循环** - 找到一个leaf就结束
3. **变量未使用** - `bestConfUpper`, `leafNodes`, `numNodes`, `leafLoop`, `overallLoop`都没用

**这看起来像未完成的代码！**

**可能的原始意图**:
```java
private void runUntilNonZero() {
    // 应该循环直到找到至少一个有非零lower bound的构象
    while (rootNode.getLowerBound().equals(BigDecimal.ZERO)) {
        // 展开最promising的节点
        // 计算leaf节点的能量
        // 更新bounds
    }
}
```

**需要检查**:
- 这个方法在其他测试中是否正常工作？
- 是否有其他版本的MARK*实现？

**优先级**: 🔴 Critical

---

### 🔴 **修复4: Bounds初始化检查**

**文件**: `src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java`

**需要检查的位置**:

1. **rootNode的初始化** - Line 195附近
   ```java
   // 需要确认rootNode的初始bounds是什么
   ```

2. **computeEpsilonErrorBounds()** - 在MARKStarNode中
   ```java
   // 需要看这个方法如何计算epsilon
   ```

3. **BigDecimal vs double转换**
   - MARK*使用BigDecimal存储bounds
   - 可能在转换过程中出错

**调查方案**:
```java
// 在init()后添加
System.out.println("Initial rootNode bounds:");
System.out.println("  lower: " + rootNode.getLowerBound());
System.out.println("  upper: " + rootNode.getUpperBound());
```

**优先级**: 🟡 High

---

### 🔴 **修复5: 能量单位和温度检查**

**需要检查的地方**:

1. **温度参数**
   ```java
   // 确认温度是否正确设置
   // 应该是298K，不是2.98K或其他
   ```

2. **能量单位**
   ```java
   // 确认是kcal/mol，不是kJ/mol或其他
   ```

3. **Boltzmann常数**
   ```java
   // RT应该 ≈ 0.592 kcal/mol at 298K
   ```

**优先级**: 🟡 Medium

---

## 修复顺序建议

### 阶段1: 调试和诊断 🔍

**目标**: 理解为什么while循环没有执行

1. **添加调试输出** (最优先)
   ```java
   // 在MARKStarBound.java的compute()方法中添加
   System.out.println("=== MARK* Compute Debug ===");
   System.out.println("targetEpsilon: " + targetEpsilon);

   if(!nonZeroLower) {
       runUntilNonZero();
       updateBound();
       System.out.println("After runUntilNonZero:");
       System.out.println("  epsilonBound: " + epsilonBound);
       System.out.println("  rootNode.lower: " + rootNode.getLowerBound());
       System.out.println("  rootNode.upper: " + rootNode.getUpperBound());
   }

   System.out.println("While loop condition check:");
   System.out.println("  epsilonBound > targetEpsilon? " + (epsilonBound > targetEpsilon));
   System.out.println("  workDone()-previousConfCount < maxNumConfs? "
       + ((workDone()-previousConfCount) + " < " + maxNumConfs));
   System.out.println("  isStable? " + isStable(stabilityThreshold));
   ```

2. **运行测试收集信息**
   ```bash
   ./gradlew test --tests "TestMARKStarVsKStarPartitionFunction.testSmallSystem3Flex"
   ```

3. **分析输出判断问题**

### 阶段2: 修复qprime赋值 ✏️

**文件**: `MARKStarBound.java` Line 248

```java
// 从：
values.qprime = rootNode.getUpperBound();

// 改为：
values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
```

### 阶段3: 修复主循环条件 🔧

根据阶段1的调试结果，可能需要：

**选项A**: 如果epsilonBound初始就很小
```java
// 确保初始epsilon足够大，强制进入循环
if(!nonZeroLower) {
    runUntilNonZero();
    updateBound();
    // 如果epsilon已经满足，强制设置一个大值
    if(epsilonBound <= targetEpsilon) {
        epsilonBound = 1.0;  // 强制进入循环
    }
}
```

**选项B**: 如果runUntilNonZero()不完整
```java
// 重新实现runUntilNonZero()，确保它真正计算一些构象
private void runUntilNonZero() {
    // 循环直到找到至少一个实际的leaf构象
    while (rootNode.getLowerBound().equals(BigDecimal.ZERO)) {
        // 展开节点
        // 计算能量
        // 更新bounds
    }
}
```

### 阶段4: 验证修复 ✅

1. **重新运行测试**
   ```bash
   ./gradlew test --tests "TestMARKStarVsKStarPartitionFunction.testSmallSystem3Flex"
   ```

2. **检查输出**:
   - MARK*应该进入主循环
   - 应该评估>0个构象
   - Partition function应该在10³-10⁵范围
   - Free energy应该在-7到-4 kcal/mol范围

3. **对比K*结果**:
   - Z_MARK* / Z_K* 应该在0.5-2.0之间
   - |G_MARK* - G_K*| 应该<1.5 kcal/mol

---

## 需要回答的问题

在开始修复之前，我需要你确认：

1. **是否要添加阶段1的调试输出？**
   - 这会帮助我们准确定位问题
   - 需要修改MARKStarBound.java

2. **是否先运行一次带调试输出的测试？**
   - 这样我们可以看到确切的epsilon值
   - 然后决定如何修复

3. **runUntilNonZero()看起来像未完成的代码，是否要重新实现它？**
   - 或者这个方法在其他地方有正确的实现？

4. **修复后是否要在所有3个测试系统上运行？**
   - 3-flex, 4-flex, 5-flex

---

## 预期结果

修复后，MARK*应该：

✅ **进入主计算循环**
```
Running until leaf is found...
Found a leaf!
Tightening from epsilon of 1.0
[进入循环，多次迭代]
Errorbound is now 0.09
```

✅ **评估多个构象**
```
MARK* evaluated: 15 conformations
```

✅ **产生合理的partition function**
```
Z_MARK* ≈ 1.0×10⁴ (与K*同数量级)
```

✅ **产生合理的free energy**
```
G_MARK* ≈ -5.5 kcal/mol (接近K*的-5.59)
```

✅ **展示更紧的bounds**
```
Gap_MARK* < Gap_K*
```

---

## 需要的文件

修复过程中可能需要修改：

1. ✏️ **MARKStarBound.java** - 主要修复目标
2. 🔍 **MARKStarNode.java** - 可能需要检查bounds计算
3. ✅ **TestMARKStarVsKStarPartitionFunction.java** - 已经准备好

已备份的文件：
- ✅ `backups/markstar_before_fix_20251107_154200/markstar/`
- ✅ `backups/markstar_before_fix_20251107_154200/TestMARKStarVsKStarPartitionFunction.java`

---

## 下一步

**请告诉我：你想先做什么？**

1. 📊 添加调试输出，运行测试看看确切的epsilon值？
2. ✏️ 直接修复qprime赋值错误？
3. 🔍 先检查computeEpsilonErrorBounds()的实现？
4. 🔧 其他方案？

我会根据你的选择，一步一步地进行修复。
