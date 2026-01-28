# MARK* Bounds Investigation

## 发现的问题

### 1. qprime赋值错误（已修复）
- **位置**: MARKStarBound.java Line 248
- **错误**: `values.qprime = rootNode.getUpperBound()`
- **修复**: `values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound())`
- **状态**: ✅ 已修复，但没有解决主要问题

### 2. While循环不执行（仍然存在）
- MARK*只运行`runUntilNonZero()`后就停止
- 没有进入主计算循环
- 结果：0个构象被评估

### 3. Bounds异常值（根本原因未知）
- rootNode.getLowerBound() = 1.41×10²⁶
- rootNode.getUpperBound() = 1.41×10²⁶
- 预期值应该在10³-10⁵范围

## 调查发现

### MARKStarBoundFastQueues vs MARKStarBound
- 测试了两种实现
- **结果**: 问题在两者中都存在
- **结论**: 问题不在于子类选择，而在于初始化或bounds计算

### Bounds计算逻辑

从`MARKStarNode.java`:

```java
// Line 388-409: makeRoot方法
public static MARKStarNode makeRoot(...) {
    Node rootNode = new Node(confSpace.positions.size());
    rootNode.index(confIndex);
    rootNode.gscore = gScorer.calc(confIndex, rcs);
    rootNode.rigidScore = rigidgScorer.calc(confIndex,rcs);
    double confUpperBound = rigidgScorer.calc(confIndex,rcs) - negatedHScorer.calc(confIndex, rcs);
    double confLowerBound = rootNode.gscore + hScorer.calc(confIndex, rcs);
    rootNode.computeNumConformations(rcs);
    rootNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound);
    return new MARKStarNode(rootNode, null);
}

// Line 465-477: setBoundsFromConfLowerAndUpper
public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
    if (lowerBound - upperBound > 1e-5){
        // 如果lower > upper，交换它们
        double temp = lowerBound;
        lowerBound = upperBound;
        upperBound = temp;
        lowerBound = Math.min(0, lowerBound);
        upperBound = Math.max(lowerBound, upperBound);
    }
    updateConfLowerBound(lowerBound);
    updateConfUpperBound(upperBound);
}

// Line 480-488: updateConfLowerBound
private void updateConfLowerBound(double tighterLower) {
    if(tighterLower > confLowerBound) {
        confLowerBound = tighterLower;
        updateSubtreeUpperBound(computeBoundsFromEnergy(confLowerBound));
    }
}

// Line 490-500: updateConfUpperBound
private void updateConfUpperBound(double tighterUpper) {
    if(tighterUpper == Double.POSITIVE_INFINITY)
        updateSubtreeLowerBound(BigDecimal.ZERO);
    if(tighterUpper < confUpperBound) {
        confUpperBound = tighterUpper;
        updateSubtreeLowerBound(computeBoundsFromEnergy(confUpperBound));
    }
}

// Line 502-504: computeBoundsFromEnergy
private BigDecimal computeBoundsFromEnergy(double energy) {
    return bc.calc(energy).multiply(new BigDecimal(getNumConformations()));
}
```

### 关键发现

**computeBoundsFromEnergy的逻辑**:
```
Partition Function Contribution = Boltzmann(energy) × NumConformations
                                 = e^(-E/RT) × N
```

对于root node（空赋值）：
- `energy` = g-score + h-score（A*的f-score）
- `NumConformations` = 24（对于3-flex系统）

**预期计算**:
```
假设 confLowerBound ≈ -5 kcal/mol
Boltzmann weight = e^(-(-5)/0.592) = e^8.45 ≈ 4700
SubtreeBound = 4700 × 24 ≈ 1.1×10⁵
```

**实际结果**:
```
SubtreeBound = 1.41×10²⁶
```

**差异**: 实际 / 预期 = 1.41×10²⁶ / 1.1×10⁵ ≈ **1.3×10²¹**

这意味着：
1. **能量值可能错误**: 如果能量是-60 kcal/mol而不是-5 kcal/mol
   - e^(60/0.592) ≈ 10⁴⁴（太大）
2. **温度参数可能错误**: 如果RT被设为错误的值
3. **NumConformations可能错误**: 如果计算出的构象数远大于24
4. **其他bug**: 在边界计算的某个地方

## 需要调查的方向

### 优先级1: 检查confLowerBound和confUpperBound的实际值
- 在`setBoundsFromConfLowerAndUpper`中添加打印
- 看看传入的能量值是多少
- 验证g-score和h-score的计算

### 优先级2: 检查BoltzmannCalculator
- 验证温度参数（应该是298K）
- 验证RT的值（应该≈0.592 kcal/mol）
- 检查是否有单位转换错误

### 优先级3: 检查NumConformations
- 验证getNumConformations()返回的值
- 应该是24，而不是10²¹

### 优先级4: 检查bounds更新逻辑
- `updateSubtreeLowerBound`和`updateSubtreeUpperBound`
- 看看是否有累积错误

## 下一步行动

1. ✅ 在MARKStarNode的makeRoot中添加调试输出
2. ✅ 在setBoundsFromConfLowerAndUpper中添加调试输出
3. ✅ 在computeBoundsFromEnergy中添加调试输出
4. ⬜ 运行测试，查看实际的能量值和NumConformations
5. ⬜ 根据输出定位问题
