# Partition Function 打印数值 - 详细数据

## 数据来源

这些数据来自修复后的 TestMARKStar.java 测试，通过添加打印功能获得。
数据位置：`build/reports/tests/test/classes/edu.duke.cs.osprey.markstar.TestMARKStar.html`

---

## 测试 1: test1GUASmall (MARK* only)

**系统**: 1GUA, 3-flexible positions, ε=0.99  
**方法**: MARK* only (保守下界)  
**测试时间**: 3.395s  

### Partition Function Bounds

| State | Lower Bound | Upper Bound | Log10 Lower | Log10 Upper |
|-------|-------------|-------------|-------------|-------------|
| **Protein** | 1.021e+35 | 1.537e+35 | 35.009 | 35.187 |
| **Ligand** | 6.443e+03 | 1.327e+04 | 3.809 | 4.123 |
| **Complex** | 2.122e+46 | 5.569e+47 | 46.327 | 47.746 |

**完整精确值**:
```
Protein:  [102117205388813970480700096775914813.988..., 153662725780773488949267184945209612.324...]
Ligand:   [6442.986..., 13270.060...]
Complex:  [21222569113267543080046536186562860935461224052.258..., 556866885088918568161669041237604998054315873082.284...]
```

### Free Energy Bounds (kcal/mol)

| State | Lower Bound | Upper Bound |
|-------|-------------|-------------|
| **Protein** | -48.049 | -47.807 |
| **Ligand** | -5.630 | -5.201 |
| **Complex** | -65.199 | -63.261 |

### K* Score

- K* score: 3.226e+07
- K* score bounds: [1.041e+07, 8.464e+08]

---

## 测试 2: testMARKStarVsKStar (K* vs MARK* 直接对比)

**系统**: 1GUA, 8-flexible positions, ε=0.68  
**方法**: 同时运行 K* (Traditional) 和 MARK*  

### 即将从 HTML 报告中提取...

