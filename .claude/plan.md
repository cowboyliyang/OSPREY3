# Grid DP vs CCD Benchmark Implementation Plan

## Goal
Create一个 Grid DP minimizer，用 sparse graph + branch decomposition 替代 CCD，对固定 RC tuple 的连续 DOF 做离散化优化。然后 benchmark 对比两者的速度和能量质量。

## 新建文件

### File 1: `GridDPMinimizer.java`
**Location:** `src/main/java/edu/duke/cs/osprey/markstar/framework/GridDPMinimizer.java`

核心数据结构：
```java
public class GridDPMinimizer {
    // 输入
    SimpleConfSpace confSpace;
    InteractionGraph interactionGraph;
    RootedTreeEdge rootEdge;     // branch decomposition 的根
    int gridSize;                // g = 3 或 5

    // Position → DOF 映射
    int[][] positionDOFIndices;  // positionDOFIndices[pos] = {dofIdx1, dofIdx2, ...}
    int[] statesPerPosition;    // g^{numDOFs} per position

    // 预计算能量表
    double[][] oneBodyEnergy;    // oneBodyEnergy[pos][gridIdx] = E_intra + E_shell
    Map<Long, double[][]> pairEnergy;  // pairEnergy[edge_key][gridIdx_i][gridIdx_j] = E_pair

    // Grid 点坐标
    double[][][] gridPoints;     // gridPoints[pos][dofLocalIdx][gridPointIdx] = DOF value
}
```

核心方法：

#### 1. `buildPositionDOFMap(ParametricMolecule pmol)`
- 遍历 pmol.dofs，用 dof.getName() 或位置信息映射 position → DOF indices
- 计算 statesPerPosition[pos] = g^{numDOFs[pos]}

#### 2. `generateGridPoints(ParametricMolecule pmol)`
- 对每个 position 的每个 DOF：在 [dofBounds.min, dofBounds.max] 内均匀取 g 个点
- gridPoints[pos][dofLocalIdx] = {min, (min+max)/2, max} for g=3

#### 3. `precomputeOneBody(ParametricMolecule pmol, ResPairCache cache)`
- 对每个 position i:
  - 用 ResInterGen.of(confSpace).addIntra(i).addShell(i).make() 创建 interactions
  - 创建 ResidueForcefieldEnergy
  - 对每个 grid 点 (0 to statesPerPosition[i]-1):
    - 解码为 DOF values
    - 调用 dof.set(value) 设置每个 DOF
    - oneBodyEnergy[i][gridIdx] = efunc.getEnergy()

#### 4. `precomputePairwise(ParametricMolecule pmol, ResPairCache cache)`
- 对每条 interactionGraph 的边 (i, j):
  - 用 ResInterGen.of(confSpace).addInter(i, j).make() 创建 interactions
  - 创建 ResidueForcefieldEnergy
  - 对每个 grid 组合 (gridIdx_i, gridIdx_j):
    - 设置 position i 和 j 的 DOFs
    - pairEnergy[key][gridIdx_i][gridIdx_j] = efunc.getEnergy()

#### 5. `runDP()` - Bottom-up DP
```
postOrder(edge):
  if edge is leaf:
    // 直接从预计算表构建 DP table
    for each mState in [0, Π statesPerPosition[m] for m in M):
      dpTable[edge][mState] = 0
      // 加上 M-positions 的 one-body
      for m in M: dpTable[edge][mState] += oneBodyEnergy[m][gridIdx(m)]
      // 加上 M-positions 之间的 pairwise (如果是 edge)
      for (m1,m2) in edges∩M: dpTable += pairEnergy[m1,m2][gridIdx(m1)][gridIdx(m2)]
      // 枚举 lambda states，取 min
      best = +inf
      for each lambdaState:
        e = 0
        // lambda one-body
        for l in lambda: e += oneBodyEnergy[l][gridIdx(l)]
        // lambda-lambda pairwise
        for (l1,l2) in edges∩lambda: e += pairEnergy[...]
        // M-lambda pairwise
        for (m,l) in edges: e += pairEnergy[...]
        // children DP tables
        e += left_child_dp[project(mState+lambdaState, leftM)]
        e += right_child_dp[project(mState+lambdaState, rightM)]
        best = min(best, e)
      dpTable[edge][mState] += best

  else (internal):
    // 类似，但从 children 的 DP 表合并

  return dpTable[edge]
```

#### 6. `minimize(RCTuple conf)` - 主入口
```java
public Result minimize(RCTuple conf) {
    ParametricMolecule pmol = confSpace.makeMolecule(conf);
    buildPositionDOFMap(pmol);
    generateGridPoints(pmol);

    long t0 = System.nanoTime();
    precomputeOneBody(pmol, cache);
    precomputePairwise(pmol, cache);
    long precompTime = System.nanoTime() - t0;

    long t1 = System.nanoTime();
    double energy = runDP();
    long dpTime = System.nanoTime() - t1;

    return new Result(energy, precompTime, dpTime);
}
```

### File 2: `BenchmarkGridDPvsCCD.java`
**Location:** `src/test/java/edu/duke/cs/osprey/minimization/BenchmarkGridDPvsCCD.java`

```java
public class BenchmarkGridDPvsCCD {

    @Test
    public void benchmark() {
        // 1. 构建 conf space (1CC8, 7 flexible residues)
        //    和 TestMinimization 相同的设置

        // 2. 计算 energy matrix (rigid + minimizing)

        // 3. 构建 InteractionGraph (dual cutoff)

        // 4. 计算 BranchDecomposition, root tree

        // 5. 用 A* 获取 N 个 conformations

        // 6. 对每个 conf:
        //    a. CCD: 用 ConfEnergyCalculator.calcEnergy(rcTuple)
        //       记录时间和能量
        //    b. Grid DP: 用 GridDPMinimizer.minimize(conf)
        //       记录时间(precomp + DP 分开)和能量

        // 7. 输出对比表:
        //    Conf | CCD_E | GridDP_E | Gap | CCD_time | GridDP_precomp | GridDP_dp | Speedup
    }
}
```

## 实现顺序

1. **GridDPMinimizer**: position→DOF 映射 + grid 生成
2. **GridDPMinimizer**: one-body 和 pairwise 预计算
3. **GridDPMinimizer**: DP 算法（先不用 branch tree，用暴力枚举验证正确性）
4. **GridDPMinimizer**: DP 改为 branch tree 版本
5. **BenchmarkGridDPvsCCD**: 集成测试 + benchmark

## 关键注意事项

- DOF 设置会修改共享的 molecule 坐标，precompute 时每次都要重新设置
- 用 `ResInterGen` 创建 partial interactions，避免手动构建
- Mixed-radix encoding: gridIdx 解码为各 DOF 的 grid 点 index
- 验证正确性: g=很大 时 Grid DP 结果应接近 CCD 结果
- 先跑 g=3, 再试 g=5 看 bound 收紧多少
