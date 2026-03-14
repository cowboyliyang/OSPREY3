# BranchMARK* 实现指南

## 1. 目标

将 BWM* 的 branch decomposition 与 MARK* 结合，加速 MARK* 的 partition function 计算。

**核心思路**：MARK* 的 pairwise 能量打分（占 70-80% 计算量）沿 branch decomposition 完全可分解；
连续最小化（continuous minimization）不可分解，仍需全局处理，但仅占 10-20% 计算量。

**复杂度提升**：从 O(q^n) 降至 O(n·q^k)，其中 k = branchwidth << n。

## 2. 数学基础

### 2.1 Partition Function 沿 Branch Decomposition 的分解

对于 pairwise 能量矩阵：
```
Z = Σ_conf boltz(E(conf))

其中 E(conf) = Σ_i E_self(i, rc_i) + Σ_{i<j} E_pair(i, rc_i, j, rc_j)
```

给定 branch decomposition 树的某条边 e，其 M-set 为 M(e)，左右子树覆盖的位置集合为 L(e) 和 R(e)：
```
Z = Σ_{m ∈ assignments(M)} boltz(E_M(m)) · Z_left(m) · Z_right(m)
```

其中：
- `E_M(m)` = M-set 位置之间的 pairwise 能量（仅 M-set 内部的交互）
- `Z_left(m)` = 左子树位置在 M-set 赋值为 m 条件下的 partition function
- `Z_right(m)` = 右子树位置在 M-set 赋值为 m 条件下的 partition function

**这个分解对 pairwise 能量是精确的**。

### 2.2 连续最小化的处理

连续最小化（CCD）耦合所有位置，无法分解。处理策略：
1. 用分解后的 pairwise bounds 快速缩小搜索空间
2. 当所有子问题的 pairwise bounds 都紧但全局 ε 仍然大时，
   从各子问题的叶节点组装完整构象，做全局连续最小化
3. triple corrections 大部分停留在子问题内部（因为大多数三元组的位置都在同一子问题中）

## 3. 架构概览

```
                    BranchMARKStarBound (全局调度器)
                    implements PartitionFunction
                    extends MARKStarBound
                           |
              ┌────────────┼────────────┐
              |            |            |
        SubproblemMARKStar  ...   SubproblemMARKStar
        (左子树,M-set赋值m1)      (右子树,M-set赋值m2)
              |                         |
         本地 PQ                    本地 PQ
    (pairwise A* 展开)          (pairwise A* 展开)
```

**两级调度**：
- 全局调度器：选择 error 贡献最大的 (子问题, M-set赋值) 对
- 本地 MARK* 队列：在子问题内做 pairwise 展开

## 4. 新建文件清单

所有新文件位于：
```
src/main/java/edu/duke/cs/osprey/markstar/framework/branch/
```

### 4.1 已完成的文件 ✅

| 文件 | 状态 | 说明 |
|------|------|------|
| `BranchNode.java` | ✅ 已完成 | 分支分解树节点，改用 int 位置索引 |
| `BranchEdge.java` | ✅ 已完成 | 分支分解树边，M-set 用 `LinkedHashSet<Integer>` |
| `BranchTree.java` | ✅ 已完成 | 分支分解树数据结构，含 deepCopy 和 getBranchwidth |
| `BreadthFirstSearch.java` | ✅ 已完成 | BFS 求增广路径，用于最大流 |
| `MinVertexCut.java` | ✅ 已完成 | Ford-Fulkerson 最大流求最小顶点割 |
| `InteractionGraph.java` | ✅ 已完成 | 从 EnergyMatrix 构建交互图 |

### 4.2 待完成的文件

| 文件 | 说明 | 复杂度 | 依赖 |
|------|------|--------|------|
| `BranchDecomposition.java` | 核心：从交互图计算分支分解 | 高（~600行） | BranchTree, MinVertexCut |
| `MSetAssignment.java` | M-set 的一个赋值（位置→RC映射） | 低（~100行） | 无 |
| `SubproblemEnergyMatrix.java` | 将全局能量矩阵投影到子问题 | 中（~150行） | EnergyMatrix |
| `SubproblemMARKStar.java` | 子问题内的 MARK*-like 计算 | 高（~400行） | MSetAssignment, SubproblemEnergyMatrix, MARKStarNode |
| `BranchMARKStarNode.java` | 全局调度队列节点 | 低（~80行） | BranchEdge, MSetAssignment |
| `BranchMARKStarBound.java` | 主算法类，全局调度器 | 高（~500行） | 上述所有 |

### 4.3 需修改的现有文件

| 文件 | 修改内容 |
|------|----------|
| `MARKStarBound.java` | 将若干 private 字段/方法改为 protected |
| `MARKStar.java` | 添加 `useBranchDecomposition` 设置，条件性创建 BranchMARKStarBound |

## 5. 各文件详细设计

---

### 5.1 BranchDecomposition.java（待实现，~600行）

**来源**：从 `/home/users/lz280/BWM/src/BranchDecomposition/BranchDecompositionH.java`（1330行）移植并简化。

**主要改动**：
- 去掉文件 I/O，直接接受 `InteractionGraph` 对象
- 顶点类型 `String` → `int`（位置索引）
- M-set 类型 `LinkedHashSet<String>` → `LinkedHashSet<Integer>`
- 去掉 `GraphVertices` 类，直接用 `Set<Integer>`
- 用 Apache Commons Math 3 替代 JAMA 做特征值分解
- 去掉 timing 代码和 System.out 输出（改用 OSPREY 的日志框架或静默处理）

**关键方法**（从 BWM 移植）：

```java
package edu.duke.cs.osprey.markstar.framework.branch;

import org.apache.commons.math3.linear.*;

public class BranchDecomposition {

    private BranchTree bt;
    private Set<Integer> graphVertices;   // 所有位置索引
    private int numPositions;

    // === 构造与入口 ===

    /** 从交互图构建分支分解 */
    public BranchDecomposition(InteractionGraph graph) {
        // 1. 从 graph.getEdgeList() 构建初始叶节点
        // 2. 记录 graphVertices 和出现在多条边中的 duplicate 位置
    }

    /** 执行分支分解算法 */
    public void compute() {
        // 移植自 BranchDecompositionH 构造函数的主循环：
        // constructStarGraph(gvDup)
        // while (!done) {
        //     bigNodes = getNodesToSplit()
        //     if (bigNodes.isEmpty()) done = true
        //     else: try push → 2-sep → 3-sep → eigen, then performSplit
        // }
    }

    public BranchTree getTree() { return bt; }
    public int getBranchwidth() { return bt.getBranchwidth(); }

    // === 从 BWM 移植的私有方法 ===

    private void constructStarGraph(Set<Integer> gvDup);
    // 创建星形图：一个中心内部节点，连接所有叶节点
    // M-set 初始化：出现在多条边中的位置加入 M-set

    private LinkedHashSet<BranchNode> getNodesToSplit();
    // 返回度 > 3 的所有内部节点

    private BranchNode findPushNode(LinkedHashSet<BranchNode> bigNodes);
    // 找有未检查边的可推送节点

    private void pushNode(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy);
    // 推送操作：检查边对是否满足推送不等式
    // |M(i) ∪ M(j) ∩ M(其他)| ≤ max(|M(i)|, |M(j)|)

    private boolean find2sep(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy);
    // 查找 2-分离

    private boolean find3sep(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy);
    // 查找 3-分离（含三角形检测）

    private boolean find23sepHelper(BranchTree H, ...);
    // 2/3-分离的辅助方法

    private boolean find2sepHelper(BranchTree H, BranchNode v, BranchNode w, ...);
    // 2-分离辅助：求 v-w 最小顶点割

    private boolean findKsepPairV(BranchTree H, BranchNode v, BranchNode w, int k, ...);
    // k-分离：含顶点收缩

    private void doEigen(BranchNode pn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy);
    // 特征向量启发式分割
    // *** 改用 Apache Commons Math 3 ***
    // RealMatrix Fm = new Array2DRowRealMatrix(F);
    // EigenDecomposition ed = new EigenDecomposition(Fm);
    // 找第二小特征值对应的特征向量，按其排序做分割

    private void performSplit(BranchNode bn, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy);
    // 执行分割：将节点 bn 拆为 bx 和 by

    private void compMLRsets(BranchEdge e, LinkedHashSet<BranchEdge> sx, LinkedHashSet<BranchEdge> sy);
    // 计算新边的 M-set = M(sx) ∩ M(sy)

    private BranchTree constructGraphH(Object[] Sa, BranchEdge[] D);
    // 构造辅助图 H

    // contractV, unifyEdges, completeXY, modifyGraphH 等辅助方法...
    // isIndependentSetGraph：检查位置集合在原始交互图中是否独立
}
```

**关键移植点**：
1. `readGraphFile` → 替换为构造函数从 `InteractionGraph.getEdgeList()` 读取
2. `gv.getVind(String)` → 直接用 int 位置索引，无需查找
3. `gv.getGraphVertices()` → 改为 `graphVertices`（`Set<Integer>`）
4. `gv.getNumV()` → 改为 `numPositions`
5. `String` 比较 (`equalsIgnoreCase`) → `int` 比较 (`==`)
6. `BranchNode(true, v1, v2)` 中 v1/v2 从 `String` 改为 `int`
7. 辅助图 H 中节点的 `getv1()` → `getPos1()`
8. JAMA 的 `Matrix` 和 `EigenvalueDecomposition` → Apache Commons Math 3 的 `RealMatrix` 和 `EigenDecomposition`

---

### 5.2 MSetAssignment.java（待实现，~100行）

```java
package edu.duke.cs.osprey.markstar.framework.branch;

/**
 * 表示 M-set 位置的一个 RC 赋值。
 * positions[i] 是第 i 个 M-set 位置的全局索引，
 * rcs[i] 是对应位置的 RC 赋值。
 */
public class MSetAssignment implements Comparable<MSetAssignment> {

    private final int[] positions;  // M-set 位置（全局索引，已排序）
    private final int[] rcs;        // 各位置的 RC 赋值

    public MSetAssignment(int[] positions, int[] rcs);

    public int[] getPositions();
    public int[] getRCs();
    public int getRC(int globalPos);  // 给定全局位置返回 RC

    /** 计算此赋值在给定能量矩阵下的 pairwise 能量（仅 M-set 内部交互） */
    public double computeMSetEnergy(EnergyMatrix emat);

    /** 计算此赋值与某个 lambda 位置之间的 pairwise 能量贡献 */
    public double computeCrossEnergy(int lambdaPos, int lambdaRC, EnergyMatrix emat);

    /** 两个赋值是否在共有位置上一致 */
    public boolean isCompatible(MSetAssignment other);

    /** 合并两个在不相交位置上的赋值 */
    public MSetAssignment merge(MSetAssignment other);

    @Override public int hashCode();    // 基于 Arrays.hashCode(positions) ^ Arrays.hashCode(rcs)
    @Override public boolean equals(Object o);
    @Override public int compareTo(MSetAssignment o);  // 按能量排序
}
```

---

### 5.3 SubproblemEnergyMatrix.java（待实现，~150行）

```java
package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

/**
 * 将全局 EnergyMatrix 投影到子问题的位置子集上。
 * 不复制数据，通过索引映射代理到全局矩阵。
 *
 * 子问题位置 = lambdaSet（子问题私有位置）∪ mSet（M-set 位置）
 * M-set 位置的 RC 已固定（来自 MSetAssignment），
 * 其能量贡献折算到 lambda 位置的 one-body 项中。
 */
public class SubproblemEnergyMatrix {

    private final EnergyMatrix globalEmat;
    private final int[] localToGlobal;  // localPos → globalPos
    private final int[] globalToLocal;  // globalPos → localPos (-1 if not in subproblem)
    private final MSetAssignment mSetAssignment;  // 固定的 M-set 赋值
    private final RCs globalRCs;

    public SubproblemEnergyMatrix(EnergyMatrix globalEmat, Set<Integer> lambdaPositions,
                                   MSetAssignment mSetAssignment, RCs globalRCs);

    /** lambda 位置的数量（子问题的"自由"位置数） */
    public int getNumPos();

    /** 本地位置 pos 在全局矩阵中的 RC 数量 */
    public int getNumConfAtPos(int localPos);

    /**
     * 有效 one-body 能量 = 原始 E_self(i, rc_i)
     *   + Σ_{m ∈ M-set} E_pair(i, rc_i, m, rc_m)
     *
     * 即将 M-set 的交互折算到 one-body 中。
     */
    public double getOneBody(int localPos, int localRC);

    /** lambda 位置之间的 pairwise 能量（直接代理到全局矩阵） */
    public double getPairwise(int localPos1, int localRC1, int localPos2, int localRC2);

    /** 计算完整子问题构象的 pairwise 能量 */
    public double confE(int[] localConf);

    // 索引映射工具
    public int globalPos(int localPos);
    public int localPos(int globalPos);
}
```

**关键设计**：M-set 位置的 RC 已固定，其与 lambda 位置的 pairwise 能量折算到 one-body 项中。
这样子问题内部的搜索只在 lambda 位置上进行。

---

### 5.4 SubproblemMARKStar.java（待实现，~400行）

```java
package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.markstar.framework.MARKStarNode;

/**
 * 子问题内的 MARK*-like partition function 计算。
 * 仅处理 pairwise 展开（不做连续最小化）。
 *
 * 每个实例对应 branch decomposition 树中的一侧（左或右子树），
 * 可以有多个 M-set 赋值，每个赋值有独立的搜索树和 bounds。
 */
public class SubproblemMARKStar {

    // === 配置 ===
    private final Set<Integer> lambdaPositions;   // 子问题私有位置（全局索引）
    private final Set<Integer> mSetPositions;     // M-set 位置（全局索引）
    private final EnergyMatrix globalRigidEmat;
    private final EnergyMatrix globalMinEmat;
    private final RCs globalRCs;

    // === 每个 M-set 赋值的 bounds ===
    private final Map<MSetAssignment, SubproblemBound> boundsByMSet;

    public static class SubproblemBound {
        MARKStarNode rootNode;           // MARK* 搜索树根
        PriorityQueue<MARKStarNode> queue;  // 搜索队列
        SubproblemEnergyMatrix rigidEmat;   // 投影的 rigid 能量矩阵
        SubproblemEnergyMatrix minEmat;     // 投影的 minimizing 能量矩阵
        RCs localRCs;                       // 本地 RCs
        BigDecimal lowerBound;              // Boltzmann 加权下界
        BigDecimal upperBound;              // Boltzmann 加权上界
        int numExpansions;                  // 已展开节点数
    }

    // === 构造 ===
    public SubproblemMARKStar(Set<Integer> lambdaPositions, Set<Integer> mSetPositions,
                               EnergyMatrix globalRigidEmat, EnergyMatrix globalMinEmat,
                               RCs globalRCs);

    // === 核心方法 ===

    /** 为一个新的 M-set 赋值初始化搜索树 */
    public void initializeForMSet(MSetAssignment mSetAssignment) {
        // 1. 创建 SubproblemEnergyMatrix（rigid 和 minimizing 版本）
        // 2. 创建本地 RCs
        // 3. 创建根 MARKStarNode（所有 lambda 位置未赋值）
        //    - 用 pairwise scoring 计算初始 bounds
        // 4. 初始化搜索队列
        // 5. 存入 boundsByMSet
    }

    /** 执行一步 pairwise A* 展开，收紧指定 M-set 赋值下的 bounds */
    public void tightenBound(MSetAssignment mSetAssignment) {
        // 类似 MARKStarBound.processPartialConfNode，但：
        // - 仅做 pairwise 展开（不做连续最小化）
        // - 在 SubproblemEnergyMatrix 上操作
        // - 当节点到达叶（所有 lambda 位置都赋值），
        //   计算 rigid/minimizing pairwise bounds 作为该构象的上下界
        //
        // 流程：
        // 1. 从队列取 error 最大的节点
        // 2. 选择一个未赋值的 lambda 位置进行展开
        // 3. 对该位置的每个 RC 创建子节点
        // 4. 计算子节点的 pairwise bounds（h-score + g-score）
        // 5. 将子节点加入队列
        // 6. 更新 SubproblemBound 的 lowerBound 和 upperBound
    }

    /** 获取指定 M-set 赋值下的 bounds */
    public SubproblemBound getBound(MSetAssignment mSetAssignment);

    public BigDecimal getLowerBound(MSetAssignment mSetAssignment);
    public BigDecimal getUpperBound(MSetAssignment mSetAssignment);
    public double getEpsilon(MSetAssignment mSetAssignment);

    /** 获取完全赋值的叶节点构象（用于全局最小化） */
    public List<int[]> getLeafConformations(MSetAssignment mSetAssignment);
    // 返回的 int[] 是全局索引的完整构象（lambda + M-set 位置）

    public boolean isInitialized(MSetAssignment mSetAssignment);
}
```

**关键设计决策**：
- 不继承 `MARKStarBound`（因为它紧耦合到全局 confSpace）
- 复用 `MARKStarNode` 数据结构（它本身只是一个树节点 + bounds）
- 需要自己实现简化版的 pairwise scoring 和 h-score 计算
- 参考 `MARKStarBound.processPartialConfNode()`（891行）的逻辑

---

### 5.5 BranchMARKStarNode.java（待实现，~80行）

```java
package edu.duke.cs.osprey.markstar.framework.branch;

import java.math.BigDecimal;

/**
 * 全局调度队列中的节点。
 * 代表 (branch edge, M-set assignment) 对及其 error 贡献。
 */
public class BranchMARKStarNode implements Comparable<BranchMARKStarNode> {

    private final BranchEdge branchEdge;
    private final MSetAssignment mSetAssignment;
    private BigDecimal subtreeLowerBound;  // 此 M-set 赋值下的 factored lower bound
    private BigDecimal subtreeUpperBound;  // 此 M-set 赋值下的 factored upper bound

    public BranchMARKStarNode(BranchEdge edge, MSetAssignment assignment);

    public BranchEdge getBranchEdge();
    public MSetAssignment getMSetAssignment();
    public BigDecimal getErrorBound();  // upper - lower

    public void updateBounds(BigDecimal lower, BigDecimal upper);

    @Override
    public int compareTo(BranchMARKStarNode other);
    // 按 error bound 降序（最大 error 优先）
}
```

---

### 5.6 BranchMARKStarBound.java（待实现，~500行）

```java
package edu.duke.cs.osprey.markstar.framework.branch;

import edu.duke.cs.osprey.markstar.framework.MARKStarBound;

/**
 * BranchMARK* 主算法类。
 * 继承 MARKStarBound，重写 tightenBoundInPhases()。
 *
 * 两级调度：
 * 1. 全局调度器选择 error 贡献最大的 (子问题, M-set 赋值) 对
 * 2. 本地子问题做 pairwise A* 展开
 * 3. 当 pairwise bounds 紧但 ε 仍大时，做全局连续最小化
 */
public class BranchMARKStarBound extends MARKStarBound {

    // === 配置 ===
    private static final double INTERACTION_THRESHOLD = 0.1;  // kcal/mol
    private static final double DENSITY_THRESHOLD = 0.8;       // 图密度阈值
    private boolean fallbackToStandard = false;

    // === Branch decomposition ===
    private InteractionGraph interactionGraph;
    private BranchDecomposition branchDecomposition;
    private BranchTree branchTree;
    private int branchwidth;

    // === 子问题 ===
    // 对 branch tree 的每条边，存储左右两侧的子问题
    // key = 边索引 * 2 + side (0=left, 1=right)
    private Map<Integer, SubproblemMARKStar> subproblems;

    // === 全局调度 ===
    private PriorityQueue<BranchMARKStarNode> globalQueue;

    // 已枚举的 M-set 赋值
    private Map<Integer, List<MSetAssignment>> enumeratedAssignments;  // edgeIndex → assignments

    // === 构造函数 ===
    public BranchMARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat,
                                EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc,
                                RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);

        // 1. 构建交互图
        interactionGraph = InteractionGraph.buildFromEnergyMatrix(
            rigidEmat, minimizingEmat, rcs, INTERACTION_THRESHOLD);

        // 2. 检查是否值得分解
        if (interactionGraph.getDensity() > DENSITY_THRESHOLD) {
            fallbackToStandard = true;
            return;
        }

        // 3. 计算分支分解
        branchDecomposition = new BranchDecomposition(interactionGraph);
        branchDecomposition.compute();
        branchTree = branchDecomposition.getTree();
        branchwidth = branchDecomposition.getBranchwidth();

        // 4. branchwidth 太大则回退
        if (branchwidth >= rcs.getNumPos() / 2) {
            fallbackToStandard = true;
            return;
        }

        // 5. 初始化子问题
        initializeSubproblems();

        // 6. 枚举初始 M-set 赋值，播种全局队列
        globalQueue = new PriorityQueue<>();
        seedGlobalQueue();
    }

    // === 核心重写 ===

    @Override
    protected void tightenBoundInPhases() {
        if (fallbackToStandard) {
            super.tightenBoundInPhases();
            return;
        }

        // PHASE 1: 从全局队列取 error 最大的节点
        BranchMARKStarNode bestNode = globalQueue.poll();
        if (bestNode == null) return;

        // PHASE 2: 在对应子问题中做 pairwise 展开
        SubproblemMARKStar sub = getSubproblem(bestNode);
        MSetAssignment mAssign = bestNode.getMSetAssignment();
        if (!sub.isInitialized(mAssign))
            sub.initializeForMSet(mAssign);
        sub.tightenBound(mAssign);

        // PHASE 3: 重算 factored bounds
        recomputeGlobalBounds();

        // PHASE 4: 检查是否需要全局最小化
        if (allSubproblemsTight() && epsilonBound > targetEpsilon) {
            doGlobalMinimization();
        }

        // PHASE 5: 将更新后的节点放回队列
        if (sub.getEpsilon(mAssign) > 0) {
            bestNode.updateBounds(sub.getLowerBound(mAssign), sub.getUpperBound(mAssign));
            globalQueue.add(bestNode);
        }

        // PHASE 6: 必要时枚举新的 M-set 赋值
        maybeEnumerateNewAssignments();

        updateBound();
    }

    // === 子问题管理 ===

    private void initializeSubproblems() {
        // 遍历 branch tree，对每条边确定左右两侧的 lambda 位置
        // 为有 lambda 位置的一侧创建 SubproblemMARKStar
    }

    /** 确定一条边两侧各有哪些位置 */
    private Set<Integer>[] getPositionsOnSides(BranchEdge edge) {
        // BFS/DFS 从边的一端出发，收集所有叶节点的位置
        // 叶节点存储的是交互图的一条边 (pos1, pos2)
        // 一侧的位置 = 该侧所有叶节点涉及的位置
    }

    // === Factored bounds 计算 ===

    /**
     * Z_edge(lower) = Σ_m boltz(E_M(m)) · Z_left_lower(m) · Z_right_lower(m)
     */
    private BigDecimal computeFactoredLowerBound(int edgeIndex) { ... }
    private BigDecimal computeFactoredUpperBound(int edgeIndex) { ... }

    private void recomputeGlobalBounds() {
        // 从根边开始递归计算 factored bounds
        // 更新 rootNode 的 subtreeLowerBound 和 subtreeUpperBound
    }

    // === M-set 赋值枚举 ===

    private void seedGlobalQueue() {
        // 对根边的 M-set，枚举初始 top-k 赋值（按 pairwise 能量排序）
        // 懒枚举避免 q^k 的指数爆炸
    }

    private void maybeEnumerateNewAssignments() {
        // 当现有赋值的 bounds 都紧了，枚举新的赋值
    }

    // === 全局最小化 ===

    private boolean allSubproblemsTight() {
        // 检查所有子问题的 epsilon 是否足够小
    }

    private void doGlobalMinimization() {
        // 1. 从各子问题的叶节点组装完整构象
        // 2. 调用 super.processFullConfNode() 做连续最小化
        // 3. 应用 triple corrections
    }

    /** 从子问题叶节点组装前 k 个最优完整构象 */
    private List<int[]> assembleTopConformations(int maxConfs) { ... }
}
```

---

## 6. MARKStarBound.java 需要的可见性修改

以下字段/方法需从 `private` 改为 `protected`：

```java
// 字段
protected MARKStarNode rootNode;              // 根节点
protected double epsilonBound;                // 当前 epsilon
protected double targetEpsilon;               // 目标 epsilon
protected UpdatingEnergyMatrix correctionMatrix;  // 修正矩阵
protected ConfEnergyCalculator minimizingEcalc;   // 最小化能量计算器
protected Parallelism parallelism;
protected MathTools.Optimizer bc;             // Boltzmann 计算器
protected int numConfsEnergied;
protected int maxMinimizations;
protected TaskExecutor loopTasks;

// 方法
protected void processFullConfNode(List<MARKStarNode> newNodes, MARKStarNode node, Node confNode);
protected void processPartialConfNode(List<MARKStarNode> newNodes, MARKStarNode node, Node confNode);
protected void computeEnergyCorrection(ConfSearch.EnergiedConf eConf, ...);
protected void updateBound();
protected void loopCleanup(List<MARKStarNode> newNodes, Stopwatch watch, int numNodes);
```

需要通过阅读 `MARKStarBound.java` 确认具体行号和确切的字段名/方法签名。

---

## 7. MARKStar.java 修改

在 `MARKStar.Settings` 中添加：
```java
public boolean useBranchDecomposition = false;
// Builder 中添加:
public Builder setUseBranchDecomposition(boolean val) {
    useBranchDecomposition = val;
    return this;
}
```

在 `ConfSpaceInfo.calcPfunc()` 中（~340行），添加条件分支：
```java
PartitionFunction pfunc;
if (settings.useBranchDecomposition) {
    pfunc = new BranchMARKStarBound(confSpace, rigidEmat, minimizingEmat,
        minimizingConfEcalc, sequence.makeRCs(confSpace), settings.parallelism);
} else {
    pfunc = new MARKStarBoundFastQueues(confSpace, rigidEmat, minimizingEmat,
        minimizingConfEcalc, sequence.makeRCs(confSpace), settings.parallelism);
}
```

---

## 8. 实现顺序

```
Phase 1 (已完成 7/7):
  ✅ BranchNode.java
  ✅ BranchEdge.java
  ✅ BranchTree.java
  ✅ BreadthFirstSearch.java
  ✅ MinVertexCut.java
  ✅ InteractionGraph.java
  ✅ BranchDecomposition.java

Phase 2 (已完成 3/3):
  ✅ MSetAssignment.java
  ✅ SubproblemEnergyMatrix.java
  ✅ SubproblemMARKStar.java

Phase 3 (已完成 2/2):
  ✅ BranchMARKStarNode.java
  ✅ BranchMARKStarBound.java (在 framework 包中，非 branch 子包)

Phase 4 (已完成):
  ✅ MARKStar.java — 添加 useBranchDecomposition 设置 + calcPfunc 中条件创建 BranchMARKStarBound
  (MARKStarBound.java 无需修改 — BranchMARKStarBound 在同一包中，可访问 package-private 字段)
```

---

## 9. 关键参考文件路径

| 文件 | 用途 |
|------|------|
| `/home/users/lz280/BWM/src/BranchDecomposition/BranchDecompositionH.java` | 分支分解算法源码（1330行） |
| `/home/users/lz280/BWM/src/BDAStar/BWMStarNode.java` | BWM* 节点（参考递归结构） |
| `/home/users/lz280/IdeaProjects/OSPREY3/src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java` | MARK* 核心（1578行） |
| `/home/users/lz280/IdeaProjects/OSPREY3/src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBoundFastQueues.java` | 扩展模式参考 |
| `/home/users/lz280/IdeaProjects/OSPREY3/src/main/java/edu/duke/cs/osprey/markstar/MARKStar.java` | 入口点 |
| `/home/users/lz280/IdeaProjects/OSPREY3/src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarNode.java` | 节点数据结构 |
| `/home/users/lz280/IdeaProjects/OSPREY3/src/main/java/edu/duke/cs/osprey/ematrix/EnergyMatrix.java` | 能量矩阵 API |
| `/home/users/lz280/IdeaProjects/OSPREY3/src/main/java/edu/duke/cs/osprey/kstar/pfunc/PartitionFunction.java` | 接口定义 |

---

## 10. 测试策略

1. **单元测试**：InteractionGraph 构建、BranchDecomposition 在已知小图上的正确性
2. **集成测试**：SubproblemMARKStar 在 4 位置问题上与标准 MARK* 对比
3. **正确性测试**：BranchMARK* 与标准 MARK* 在 4-8 位置问题上收敛到相同 ε-近似
4. **回退测试**：稠密图自动回退到标准 MARK*
5. **性能测试**：10-25 位置稀疏交互图上比较计算时间
