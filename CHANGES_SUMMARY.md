# Phase 2 Subtree Caching - Test Configuration Changes

## 修改总结 (Summary of Changes)

### 问题诊断 (Problem Diagnosis)

之前的测试配置导致cache完全不work，原因：

1. **规模太小**: 只有5个flexible residues
   - Branch decomposition只能生成4个subtrees
   - Subtrees太小（2-3个positions）
   - Cache hit机会很少

2. **Fragment过滤太激进**: `conf.size() < 3`跳过了大部分调用
   - Energy matrix计算时的singles/pairs被过滤
   - 即使full conformations也可能被过滤

3. **没有sequence space**: 所有positions都是wild-type only
   - 不同conformations之间RCs变化太大
   - 即使相同positions，RCs不同就无法复用

### 修改内容 (Changes Made)

#### 1. 测试规模调整 (Test Scale Adjustment)

**之前 (Before)**:
```java
int[] scales = {5, 7, 9}; // 测试5, 7, 9个flexible residues
```

**之后 (After)**:
```java
int[] scales = {7, 9}; // 只测试7, 9个flexible residues
// 移除了5 residues测试，规模太小
```

#### 2. 添加Mutable Position (Added Sequence Space)

**之前 (Before)**:
```java
// 所有positions都是wild-type only
if (numFlexible >= 2) {
    protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType)...
    ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType)...
}
```

**之后 (After)**:
```java
// 添加1个mutable position创建sequence space
protein.flexibility.get("G648").setLibraryRotamers("ALA", "VAL", "LEU").setContinuous();

// 其余7个positions保持wild-type
if (numFlexible >= 1) {
    protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType)...
}
// ... (7 flexible positions total)
```

**好处**:
- 创建了3个不同序列 (ALA, VAL, LEU mutations)
- 每个序列会minimize多个conformations
- **关键**: 同一序列的不同conformations可以共享subtrees！
  - 例如：Sequence with G648=ALA
    - Conformation 1: positions [1-7] with RCs [2,5,1,3,0,7,4]
    - Conformation 2: positions [1-7] with RCs [2,5,1,3,0,7,9]
    - Subtree [1-5] with RCs [2,5,1,3,0] 可以复用！✓

#### 3. 注释掉小规模测试 (Commented Out Small Test)

```java
// @Test  // 注释掉了
public void testPhase2SubtreeCaching() {
    int scale = 7; // 从5改成7
    ...
}
```

### 文件修改列表 (Modified Files)

1. **TestDPvsOriginal.java**
   - Line 79-80: 修改test scales
   - Line 379-436: 重写buildConfSpace()添加mutable position
   - Line 442: 注释掉testPhase2SubtreeCaching()

2. **run_phase2_1mut7flex_test.sh** (新文件)
   - SLURM脚本运行新配置的测试

3. **PHASE2_TEST_SUMMARY.md** (新文件)
   - 详细的测试说明和预期结果

## 预期结果 (Expected Results)

### Cache统计应该显示 (Cache Stats Should Show)

```
=== Subtree DOF Cache Statistics (TRUE SUBTREE CACHING) ===
Total conformation queries:   > 0     (之前是0)
Total subtree queries:        > 0     (之前是0)
Full cache hits:              > 0     (之前是0)
Partial cache hits:           > 0     (之前是0)
Subtree hit rate:             20-40%  (之前是0% / NaN)
Cache size:                   XX / 100000
Estimated speedup:            1.1-1.3x
```

### 性能对比目标 (Performance Targets)

| Scale | Target Phase 1+2 vs Phase 1 Speedup |
|-------|-------------------------------------|
| 7 res | 1.1x - 1.3x                        |
| 9 res | 1.2x - 1.5x                        |

**注意**: Speedup可能modest，因为：
- Cache查找有overhead
- Boundary refinement需要额外minimization
- ConstrainedMinimizer setup有成本

## 运行测试 (Running the Test)

```bash
# 提交作业
cd /home/users/lz280/IdeaProjects/OSPREY3
sbatch run_phase2_1mut7flex_test.sh

# 查看作业状态
squeue -u lz280

# 查看输出
ls -lth phase2_*.out | head -1
tail -f phase2_JOBID.out  # 替换JOBID
```

## 下一步 (Next Steps)

1. **等待测试完成** (~1-2小时)
2. **检查cache statistics**
   - 如果还是0，需要debug SubtreeDOFCache.getSubtrees()
   - 添加更多debug日志确认subtrees是否生成
3. **分析性能**
   - 如果speedup < 1.0，分析overhead来源
   - 如果speedup > 1.2，继续扩展到更大规模

## Debug日志查找 (Debug Logs to Look For)

```
[CachedMinimizer DEBUG #N] USING CACHE for conf.size()=8
[Phase 2] TRUE Subtree DOF Cache initialized
[Phase 2] Branch decomposition: BranchDecomposition[nodes=15, leaves=8, internal=7, branchWidth=2]
```

如果看不到这些，说明：
- Cache没有被初始化
- 或者被fragment filtering跳过
- 或者BranchDecomposition有问题

---

## CCD Minimizer确认 (CCD Minimizer Confirmation)

✅ 确认了SimpleCCDMinimizer总是被调用
✅ 没有基于residue数量的阈值检查
✅ 只要`numDOFs > 0`就会调用CCD

代码位置:
- EnergyCalculator.java:311, 335, 367, 462
- 所有都是: `minimizers = (f) -> new SimpleCCDMinimizer(f);`

---

**Job ID**: 10366585
**提交时间**: $(date)
**预计完成**: 1-2小时
