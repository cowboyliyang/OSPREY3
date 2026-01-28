# MARKStar Subtree Cache 与原版差异分析

## 概述

本文档分析 MARKStar 中 **Subtree Cache（子树缓存）** 优化与原版实现的差异。

## 核心差异总结

### 1. 原版实现（无 Subtree Cache）

**最小化流程：**
```
EnergyCalculator.calcEnergy()
  ↓
创建 Minimizer
  ↓
minimizer.minimizeFrom(initialDOFs)  // 完整最小化整个构象
  ↓
返回结果
```

**特点：**
- 每次最小化都是**从头开始**，不利用历史结果
- 即使两个构象共享大部分相同的残基组合（RC assignments），也要重新最小化
- 简单直接，但存在大量重复计算

**代码位置：**
```721:726:src/main/java/edu/duke/cs/osprey/energy/EnergyCalculator.java
try (Minimizer minimizer = context.minimizers.make(f)) {
    // 原版：直接使用 minimizer
    Minimizer.Result result = minimizer.minimizeFrom(x);
    // ...
}
```

### 2. 带 Subtree Cache 的实现

**最小化流程：**
```
EnergyCalculator.calcEnergy()
  ↓
创建 Minimizer
  ↓
CachedMinimizer.wrapMinimizerIfNeeded()  // 包装为 CachedMinimizer
  ↓
CachedMinimizer.minimizeFrom()
  ↓
SubtreeDOFCache.minimizeWithCache()
  ├─ BranchDecomposition 分解构象为子树
  ├─ 检查每个子树的缓存
  ├─ 对未缓存的子树：ConstrainedMinimizer 只最小化该子树
  ├─ 组合所有子树的 DOF 值
  └─ 边界细化（refine boundaries）
  ↓
返回结果
```

**特点：**
- 使用 **BranchDecomposition** 将构象分解为子树
- 缓存每个子树的**最小化后的 DOF 值**
- 当新构象包含已缓存的子树时，**直接重用** DOF 值
- 只对**未缓存的子树**进行最小化
- 最后进行**边界细化**以优化子树之间的接口

**代码位置：**
```767:772:src/main/java/edu/duke/cs/osprey/energy/EnergyCalculator.java
private Minimizer wrapMinimizerIfNeeded(Minimizer minimizer, RCTuple conf, ObjectiveFunction objFunc) {
    if (edu.duke.cs.osprey.ematrix.CachedMinimizer.ENABLE_SUBTREE_CACHE && conf != null) {
        return new edu.duke.cs.osprey.ematrix.CachedMinimizer(minimizer, conf, objFunc);
    }
    return minimizer;
}
```

## 关键组件

### 1. CachedMinimizer

**作用：** Minimizer 的包装器，在最小化时启用子树缓存

**关键代码：**
```118:160:src/main/java/edu/duke/cs/osprey/ematrix/CachedMinimizer.java
@Override
public Result minimizeFrom(DoubleMatrix1D x) {
    if (!enableCache || conf == null || dofCache == null || objectiveFunction == null) {
        // 缓存未启用，使用标准最小化
        return delegate.minimizeFrom(x);
    }

    // 只对完整构象（size >= 3）使用缓存
    if (conf.size() < 3) {
        return delegate.minimizeFrom(x);
    }

    // 使用 SubtreeDOFCache 进行子树缓存最小化
    SubtreeDOFCache.MinimizationResult result =
        dofCache.minimizeWithCache(conf, delegate, x, objectiveFunction);

    return new Result(result.dofs, result.energy);
}
```

**关键特性：**
- 全局开关：`ENABLE_SUBTREE_CACHE`（默认 false）
- 需要 `ObjectiveFunction` 来创建 `ConstrainedMinimizer`
- 只对完整构象（size >= 3）启用缓存，片段（fragments）跳过

### 2. SubtreeDOFCache

**作用：** 实现真正的子树 DOF 缓存逻辑

**核心算法：**
```86:185:src/main/java/edu/duke/cs/osprey/ematrix/SubtreeDOFCache.java
public MinimizationResult minimizeWithCache(
        RCTuple conf,
        Minimizer minimizer,
        DoubleMatrix1D initialDOFs,
        ObjectiveFunction objectiveFunction) {

    // 1. 从 BranchDecomposition 获取子树
    List<Subtree> subtrees = getSubtrees(conf);

    // 2. 检查每个子树的缓存
    for (Subtree subtree : subtrees) {
        SubtreeKey key = new SubtreeKey(subtree, conf);
        MinimizedSubtree cached = cache.get(key);

        if (cached != null) {
            // 缓存命中：直接应用缓存的 DOF 值
            applySubtreeDOFs(combinedDOFs, cached.dofs, localDOFIndices);
        } else {
            // 缓存未命中：需要最小化这个子树
            uncachedSubtrees.add(subtree);
        }
    }

    // 3. 最小化未缓存的子树
    for (Subtree subtree : uncachedSubtrees) {
        // 使用 ConstrainedMinimizer 只优化这个子树的 DOF
        ConstrainedMinimizer constrainedMin = new ConstrainedMinimizer(
            minimizer, objectiveFunction, freeDOFIndices, combinedDOFs);
        
        Minimizer.Result subtreeResult = constrainedMin.minimizeFrom(combinedDOFs);
        
        // 缓存结果
        cache.put(key, new MinimizedSubtree(subtreeDOFs, subtreeResult.energy));
    }

    // 4. 边界细化
    if (ENABLE_BOUNDARY_REFINEMENT && subtrees.size() > 1) {
        finalEnergy = refineBoundaries(combinedDOFs, subtrees, minimizer, objectiveFunction, pmol);
    }

    return new MinimizationResult(combinedDOFs, finalEnergy, fullyCached);
}
```

**关键特性：**
- 使用 `BranchDecomposition` 分解构象
- 动态计算 **LOCAL DOF indices**（基于当前构象的 ParametricMolecule）
- 支持**部分缓存命中**（部分子树缓存，部分未缓存）
- **边界细化**：优化子树之间的接口 DOF

### 3. BranchDecomposition

**作用：** 将构象空间分解为子树结构

**算法：**
- 构建交互图（interaction graph）
- 使用贪心平衡分割将位置分成大致相等的两半
- 递归构建二叉树
- 每个节点代表一个位置子集

**代码位置：**
```75:93:src/main/java/edu/duke/cs/osprey/ematrix/BranchDecomposition.java
public BranchDecomposition(SimpleConfSpace confSpace) {
    this.confSpace = confSpace;

    // 获取所有位置
    Set<Integer> allPositions = new HashSet<>();
    for (int i = 0; i < confSpace.positions.size(); i++) {
        allPositions.add(i);
    }

    // 构建交互图
    Map<Integer, Set<Integer>> graph = buildInteractionGraph(allPositions);

    // 构建树
    this.root = buildTree(allPositions, graph);

    // 计算分支宽度
    this.branchWidth = computeBranchWidth(root);
}
```

### 4. ConstrainedMinimizer

**作用：** 只优化指定 DOF 的最小化器

**关键特性：**
- 只优化 `freeDOFIndices` 指定的 DOF
- 其他 DOF 保持固定
- 用于最小化单个子树，而不影响其他部分

## 数据流对比

### 原版数据流

```
构象 A = {pos0:RC2, pos1:RC5, pos2:RC1, pos3:RC3}
  ↓
完整最小化（所有 DOF）
  ↓
结果 DOF_A, Energy_A

构象 B = {pos0:RC2, pos1:RC5, pos2:RC1, pos3:RC7}  // 只有 pos3 不同
  ↓
完整最小化（所有 DOF）  // 重复计算 pos0,1,2
  ↓
结果 DOF_B, Energy_B
```

### Subtree Cache 数据流

```
初始化：
  BranchDecomposition 分解为：
    Subtree1: {pos0, pos1, pos2}
    Subtree2: {pos3}

构象 A = {pos0:RC2, pos1:RC5, pos2:RC1, pos3:RC3}
  ↓
检查缓存：
    Subtree1 {pos0:RC2, pos1:RC5, pos2:RC1} → 未缓存
    Subtree2 {pos3:RC3} → 未缓存
  ↓
最小化 Subtree1 → 缓存 DOF_subtree1
最小化 Subtree2 → 缓存 DOF_subtree2
  ↓
组合 + 边界细化
  ↓
结果 DOF_A, Energy_A

构象 B = {pos0:RC2, pos1:RC5, pos2:RC1, pos3:RC7}  // 只有 pos3 不同
  ↓
检查缓存：
    Subtree1 {pos0:RC2, pos1:RC5, pos2:RC1} → **缓存命中！** 重用 DOF_subtree1
    Subtree2 {pos3:RC7} → 未缓存
  ↓
最小化 Subtree2 → 缓存 DOF_subtree2_new
  ↓
组合（重用 Subtree1 + 新的 Subtree2）+ 边界细化
  ↓
结果 DOF_B, Energy_B
```

## 性能影响

### 预期加速

- **最佳情况**：30-50% 加速（当许多构象共享子树时）
- **最坏情况**：5-10% 开销（子树重用率低，缓存管理成本）

### 缓存统计

SubtreeDOFCache 提供详细的统计信息：
- **Full cache hits**：所有子树都缓存
- **Partial cache hits**：部分子树缓存
- **Full cache misses**：所有子树都未缓存
- **Subtree hit rate**：子树级别的缓存命中率

## 启用方式

### 1. 全局开关

```java
CachedMinimizer.ENABLE_SUBTREE_CACHE = true;
```

### 2. 初始化缓存

```java
CachedMinimizer.initializeGlobalCache(confSpace);
```

### 3. 传递 conf 参数

在 `EnergyCalculator.calcEnergy()` 中传递 `RCTuple conf`：
```java
ecalc.calcEnergy(pmol, inters, approximator, conf);  // conf 参数启用缓存
```

### 4. 查看统计

```java
CachedMinimizer.printGlobalStats();
```

## 关键设计决策

### 1. 为什么需要 ObjectiveFunction？

`ConstrainedMinimizer` 需要 `ObjectiveFunction` 来创建约束最小化器，只优化指定 DOF。

### 2. 为什么使用 LOCAL DOF indices？

不同氨基酸的 DOF 数量不同（如 GLY vs TRP），且构象的 DOF 顺序可能因突变而变化。使用 LOCAL indices 确保正确映射。

### 3. 为什么需要边界细化？

组合缓存的子树时，子树之间的接口 DOF 可能不是最优的。边界细化只优化接口 DOF，提高准确性。

### 4. 为什么跳过小片段（size < 3）？

小片段（单残基或残基对）通常用于能量矩阵计算，不适合子树分解。跳过它们避免不必要的开销。

## 代码修改位置总结

### 新增文件
1. `SubtreeDOFCache.java` - 子树缓存核心实现
2. `CachedMinimizer.java` - Minimizer 包装器
3. `BranchDecomposition.java` - 分支分解算法
4. `ConstrainedMinimizer.java` - 约束最小化器

### 修改的文件
1. `EnergyCalculator.java` - 添加 `wrapMinimizerIfNeeded()` 和 `calcEnergy(..., RCTuple conf)` 重载
2. `ConfEnergyCalculator.java` - 传递 `conf` 参数到 `EnergyCalculator`

## 总结

**Subtree Cache 的核心思想：**
- 识别构象之间的**共同子树**
- **缓存**这些子树的**最小化结果**
- **重用**缓存结果，只最小化**新的部分**
- 通过**边界细化**保证准确性

**与原版的主要区别：**
- 原版：每次都完整最小化
- Subtree Cache：智能重用，只最小化新部分

**适用场景：**
- 当许多构象共享相同的残基组合时效果最好
- 在 MARKStar 的 A* 搜索中，相邻节点经常共享子树，因此缓存效果显著
