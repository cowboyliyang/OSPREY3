# 运行更大规模的MARK*测试

## 目的
测试更大系统（更多flexible残基）来确保minimize足够多的conformations，以便分析：
- Minimize顺序 vs Partition Function贡献的相关性

## 测试配置

### 系统大小
创建了不同规模的测试：

| 测试 | Flexible残基数 | 预计Conf Space | 预计耗时 | 预计minimize数量 |
|------|---------------|----------------|----------|------------------|
| test4FlexibleResiduesTouchstone | 4+1 (ligand) = 5 | ~10^10 | 5-10分钟 | 20-50 |
| test5FlexibleResidues | 5+1 = 6 | ~10^12 | 10-30分钟 | 50-100 |
| test6FlexibleResidues | 6+1 = 7 | ~10^14 | 30-60分钟 | 100-200 |
| test7FlexibleResidues | 7+1 = 8 | ~10^16 | 1-2小时 | 200-500 |

### 固定参数
- **Epsilon**: 0.01 (1% 相对误差)
- **Parallelism**: 4 CPUs
- **BigDecimal精度**: 64位有效数字

## 如何运行

### 快速测试（推荐先运行）
```bash
# 4个flexible残基 - 预计5-10分钟
./gradlew cleanTest test --tests "edu.duke.cs.osprey.markstar.TestLargerSystem.test4FlexibleResiduesTouchstone"
```

### 中等规模测试
```bash
# 5个flexible残基 - 预计10-30分钟
./gradlew cleanTest test --tests "edu.duke.cs.osprey.markstar.TestLargerSystem.test5FlexibleResidues"
```

### 大规模测试
```bash
# 6个flexible残基 - 预计30-60分钟
./gradlew cleanTest test --tests "edu.duke.cs.osprey.markstar.TestLargerSystem.test6FlexibleResidues"

# 7个flexible残基 - 预计1-2小时
./gradlew cleanTest test --tests "edu.duke.cs.osprey.markstar.TestLargerSystem.test7FlexibleResidues"
```

## 输出文件

每个测试会生成CSV文件：
```
markstar_order_vs_contribution_4flex_eps0p01.csv
markstar_order_vs_contribution_5flex_eps0p01.csv
markstar_order_vs_contribution_6flex_eps0p01.csv
markstar_order_vs_contribution_7flex_eps0p01.csv
```

CSV格式：
```
MinimizeOrder,Conformation,Score,LowerBound,Corrected,FinalEnergy,BoltzmannWeight,PercentContribution,ContributionRank
1,[4,11,7,10,3],-24.535469,-24.658562,-24.658562,-27.246422,1.234568e+08,12.3456,1
2,[4,11,7,10,4],-24.123456,-24.234567,-24.234567,-26.834512,9.876543e+07,9.8765,2
...
```

## 控制台输出示例

```
================================================================================
MARK* 快速测试: 4个flexible残基 (预计5-10分钟)
================================================================================

配置:
  - Flexible positions: 5
  - Target epsilon: 0.01
  - Total positions: 5
  - Estimated conf space size: ~10^10 (假设每个位置~100个rotamers)

设置protein flexible残基:
  - Residue 21: wild-type rotamers + continuous
  - Residue 22: wild-type rotamers + continuous
  - Residue 23: wild-type rotamers + continuous
  - Residue 24: wild-type rotamers + continuous
设置ligand flexible残基:
  - Residue 209: wild-type rotamers + continuous

创建energy calculators...
计算energy matrices...
Energy matrix计算完成 (耗时: 45.3s)

运行MARK* with tracking...
这可能需要几分钟到几十分钟，取决于系统大小...

[MARK*输出...]

================================================================================
MARK* 计算完成!
================================================================================
Partition function计算耗时: 387.5s
Status: Estimated
Lower bound (q*): 1.234567e+10
Upper bound (p*): 1.246789e+10
Effective epsilon: 0.0098
Conformations minimized: 47

解析到 47 个minimized conformations

================================================================================
分析: Minimize顺序 vs Partition Function贡献
================================================================================

总共minimize了 47 个conformations
Total Z (from minimized confs): 1.234567e+10

================================================================================
前20个被minimize的conformations:
--------------------------------------------------------------------------------
Order  | Conformation         | Energy     | Boltzmann    | % of Z     | Rank
--------------------------------------------------------------------------------
1      | [4,11,7,10,3]       |   -27.2464 |  1.234568e+08|    12.3456%| #1
2      | [4,11,7,10,4]       |   -26.8321 |  9.876543e+07|    9.8765% | #2
3      | [4,11,7,11,3]       |   -26.1234 |  7.654321e+07|    7.6543% | #3
...

前20个conformations的累积贡献: 68.42%

================================================================================
累积贡献分析:
--------------------------------------------------------------------------------
前  10% (   5 confs) 贡献:  45.23% of total Z
前  25% (  12 confs) 贡献:  72.15% of total Z
前  50% (  24 confs) 贡献:  89.34% of total Z
前  75% (  35 confs) 贡献:  96.78% of total Z
前 100% (  47 confs) 贡献: 100.00% of total Z

================================================================================
相关性分析:
--------------------------------------------------------------------------------
Spearman秩相关系数: 0.8234
Top-5 overlap: 5/5 (100.0%)
Top-10 overlap: 9/10 (90.0%)
Top-20 overlap: 17/20 (85.0%)

详细数据已保存到: markstar_order_vs_contribution_4flex_eps0p01.csv

================================================================================
结论:
================================================================================
✓ 优先minimize的conformations确实对Z贡献最大!
  - 前10%的confs贡献了 45.2% 的Z
  - Spearman相关系数: 0.823
================================================================================
```

## 分析指标解释

### 1. Spearman秩相关系数 (ρ)
- 衡量两个排序的相关性
- **ρ > 0.7**: 强相关 ✓ - minimize顺序与贡献大小强相关
- **ρ > 0.3**: 中等相关 ○
- **ρ < 0.3**: 弱相关 ✗

### 2. Top-K Overlap
- 前K个被minimize的有多少也在贡献最大的前K名中？
- **Top-10 overlap: 9/10 (90%)**
  - 意味着前10个minimize的，有9个也是贡献最大的前10名
  - 90%的overlap说明相关性很强

### 3. 累积贡献
- **前10%的conformations贡献了45%的Z**
  - 说明MARK*确实优先处理了重要的conformations
- 如果是随机顺序，前10%应该只贡献约10%

### 4. 期望结果
对于一个好的MARK*实现：
- Spearman相关系数应该 > 0.6
- Top-10 overlap应该 > 70%
- 前10%的conformations应该贡献 > 30% 的Z

## 监控进度

测试运行时，可以在另一个终端查看输出：
```bash
tail -f test_4flex_output.txt
```

或者查看是否完成：
```bash
ls -lh markstar_order_vs_contribution_*.csv
```

## 如果测试太慢

如果4个flexible残基的测试超过15分钟还没完成，可以：
1. 停止测试 (Ctrl+C)
2. 查看partial结果（如果有）
3. 考虑使用更少的flexible残基（3个）

## 预期发现

如果MARK*算法设计正确，我们应该看到：
1. **强正相关** (Spearman > 0.6)
2. **高Top-K overlap** (> 70%)
3. **前10%贡献大** (> 30% of Z)

这将验证：**MARK*确实优先minimize了对partition function贡献最大的conformations！**

## 当前测试正在运行

测试 `test4FlexibleResiduesTouchstone` 正在后台运行...
输出保存到: `test_4flex_output.txt`

检查进度：
```bash
tail -f test_4flex_output.txt | grep -E "(解析到|Spearman|结论|✓)"
```
