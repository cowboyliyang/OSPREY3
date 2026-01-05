# CometsZ + BBK* + MARK* Integration

## 概述

这个集成组合了三种算法来优化蛋白质设计中的序列搜索和分区函数计算：

1. **CometsZ**: 序列空间搜索算法
2. **BBK***: K*算法的批处理优化版本
3. **MARK***: 高效的分区函数计算算法

## 架构层次

```
CometsZ (Sequence Search)
  └── BBK* (K* with Batching) [可选]
      └── MARK* (Partition Function)
          └── MARKStarBound (Efficient Bounds)
```

## 已实现的组合

### 1. CometsZ + GradientDescent (Baseline)
- 传统的CometsZ实现
- 使用GradientDescentPfunc计算分区函数
- **文件**: `TestCometsZWithMARKStarPerformance.java`

### 2. CometsZ + MARK*
- CometsZ序列搜索
- MARK*替代GradientDescentPfunc
- **文件**: `TestCometsZWithMARKStarPerformance.java`
- **测试**: ✅ 正在运行（cometsz_markstar_10135927.out）

### 3. CometsZ + BBK* + MARK* (NEW!)
- CometsZ序列搜索
- BBK*的批处理和剪枝策略
- MARK*用于分区函数计算
- **文件**: `TestCometsZWithBBKStarAndMARKStar.java`
- **测试**: ✅ 结构完成，等待运行

## 文件结构

```
src/test/java/edu/duke/cs/osprey/kstar/
├── TestCometsZWithMARKStarPerformance.java        # CometsZ + MARK*
├── TestCometsZWithBBKStarAndMARKStar.java         # CometsZ + BBK* + MARK* (NEW)
└── TestKStar.java                                  # Test case definitions

submit_scripts/
├── submit_cometsz_markstar_tests.sh               # CometsZ + MARK*
└── submit_cometsz_bbkstar_markstar_tests.sh       # CometsZ + BBK* + MARK* (NEW)
```

## 性能对比目标

### 预期加速

| 组合 | 预期 vs Baseline |
|------|------------------|
| CometsZ + MARK* | 1.5x - 3x |
| CometsZ + BBK* + GradientDescent | 1.2x - 2x |
| CometsZ + BBK* + MARK* | **2x - 5x** |

加速来源：
- **MARK***: 更高效的conformations minimize顺序
- **BBK***: 批处理减少overhead，更好的序列剪枝
- **组合**: 两者的乘法效应

## 运行测试

### 本地测试（快速验证）

```bash
# 测试CometsZ + MARK*
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsZWithMARKStarPerformance.compare2RL0FullGradientDescentVsMARKStar"

# 测试CometsZ + BBK* + MARK*结构
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestCometsZWithBBKStarAndMARKStar.testBBKStarIntegrationStructure"
```

### SLURM完整测试（几小时）

```bash
# 提交CometsZ + MARK*测试
sbatch submit_cometsz_markstar_tests.sh

# 提交CometsZ + BBK* + MARK*测试
sbatch submit_cometsz_bbkstar_markstar_tests.sh
```

### 查看结果

```bash
# 查看正在运行的作业
squeue -u $USER

# 查看输出
tail -f cometsz_bbkstar_markstar_JOBID.out

# 提取性能指标
grep "MARK\* Speedup:" cometsz_bbkstar_markstar_JOBID.out
grep "Total Time:" cometsz_bbkstar_markstar_JOBID.out
```

## 实现细节

### CometsZ集成点

CometsZ使用`pfuncFactory`来创建partition function实例：

```java
// GradientDescent版本
state.pfuncFactory = (rcs) -> {
    return new GradientDescentPfunc(
        state.confEcalc,
        minimizingConfTree,
        rigidConfTree,
        rcs.getNumPos()
    );
};

// MARK*版本
state.pfuncFactory = (rcs) -> {
    return new MARKStarBound(
        state.confSpace,
        rigidEmat,
        minimizingEmat,
        state.confEcalc,
        rcs,
        parallelism
    );
};
```

### BBK*特性

CometsZ已经内置了一些BBK*特性：
- `setMinNumConfTrees()`: 批处理
- Objective-based filtering: 序列剪枝
- Memory management: 高效内存使用

显式的BBK*集成可以添加：
- K* score-based sequence tree pruning
- 更激进的bounds-based filtering
- Sequence space partitioning

## 当前状态

### ✅ 完成
1. CometsZ + GradientDescent (baseline)
2. CometsZ + MARK* (正在SLURM上运行)
3. CometsZ + BBK* + MARK* 测试结构

### ⏳ 进行中
- CometsZ + MARK*性能测试（Job 10135927，已运行16小时）

### 📋 待办
1. 运行CometsZ + BBK* + MARK*完整测试
2. 对比所有组合的性能
3. 分析加速来源
4. 优化参数（batch size, epsilon, etc.）

## 测试配置

### 2RL0测试系统
- **Flexible positions**: 8
- **Conformation space**: ~10^16
- **Epsilon**: 0.95
- **Num sequences**: 25
- **Runtime**: 数小时到1天

### 资源需求
- **Memory**: 50 GB
- **CPUs**: 8
- **Time limit**: 48 hours

## 结果分析

测试完成后，会生成如下指标：

```
========================================
CometsZ + GradientDescent (Baseline)
========================================
Setup Time:      XXXX ms
Execution Time:  XXXX ms
Total Time:      XXXX ms
Sequences Found: 25

========================================
CometsZ + MARK*
========================================
Setup Time:      XXXX ms
Execution Time:  XXXX ms
Total Time:      XXXX ms
Sequences Found: 25

========================================
MARK* Speedup: X.XXx (Setup: X.XXx, Execution: X.XXx)
========================================
```

## 参考文献

1. **MARK***: Matrix A* for K* - 高效分区函数计算
2. **BBK***: Branch and Bound over K* - 批处理K*优化
3. **CometsZ**: Combinatorial Multi-state Exhaustive Tree Search - 序列空间搜索

## 联系方式

- **开发者**: Yuxi Long (lz280@duke.edu)
- **实验室**: Bruce Donald Lab, Duke University
- **项目**: OSPREY 3.0
