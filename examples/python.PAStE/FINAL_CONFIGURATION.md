# 最终推荐配置 (Final Recommended Configuration)

## 回答：哪个版本更好？

**最佳配置（已实现）**：混合方法

```
✓ K03_template.in / K04_template.in  → 手工从CIF创建（准确几何结构）
✓ K03_coords.in / K04_coords.in      → ANTECHAMBER AM1-BCC电荷（量子化学准确度）
✓ K03_rotlib.dat / K04_rotlib.dat    → 手工创建（统计学rotamer库）
```

## 为什么这样配置最好？

### 1. Template.in（内坐标）：手工从CIF创建
**优点：**
- ✓ 100%准确的几何参数（键长、键角、二面角）来自PDB理想坐标
- ✓ 原子命名与coords.in一致（NZ, CH, OH, CQ, CS, CU等）
- ✓ OSPREY能正确识别所有原子
- ✓ 保持与AMBER force field的兼容性

**为什么不用ANTECHAMBER生成的template？**
- ✗ ANTECHAMBER使用CIF命名（CX, CK, CP, OP1, OP2）
- ✗ 与我们的coords.in命名不匹配
- ✗ OSPREY报错：`didn't find atom named CQ`

### 2. Coords.in（坐标+电荷）：ANTECHAMBER电荷
**优点：**
- ✓ AM1-BCC量子化学电荷（比手工估计准确得多）
- ✓ ALY总电荷：-0.001（几乎完美的0）
- ✓ SLL总电荷：-0.999（几乎完美的-1）
- ✓ 已发表的标准方法

**改进效果：**
```
ALY关键原子电荷变化：
  CH (羰基C):  +0.5973 → +0.2980  (差异: -0.30)
  OH (羰基O):  -0.5679 → -0.3660  (差异: +0.20)

SLL关键原子电荷变化：
  CU (羧基C):  +0.8014 → +0.4190  (差异: -0.38)
  OU1(羧基O):  -0.8188 → -0.2110  (差异: +0.61)
  OU2(羧基O):  -0.8188 → -0.2110  (差异: +0.61)
```

### 3. Rotlib.dat（rotamer库）：手工创建
**原因：**
- ✓ 基于统计数据（Dunbrack库）
- ✓ 不需要量子化学计算
- ✓ 已根据LYS骨架正确设置

## 当前状态确认

### 文件来源：
```bash
K03_template.in  → 手工（从ALY.cif，2025-01-07 12:28备份）
K04_template.in  → 手工（从SLL.cif，去质子化-COO⁻，2025-01-07 12:28备份）
K03_coords.in    → 混合（CIF坐标 + ANTECHAMBER电荷，2025-01-07 12:44更新）
K04_coords.in    → 混合（CIF坐标 + ANTECHAMBER电荷，2025-01-07 12:44更新）
K03_rotlib.dat   → 手工（基于LYS）
K04_rotlib.dat   → 手工（基于LYS）
```

### 关键决策：
1. **SLL质子化状态**：-COO⁻（去质子化，pH 7.4，总电荷-1）✓
2. **电荷来源**：ANTECHAMBER AM1-BCC（量子化学）✓
3. **几何来源**：PDB CIF（理想坐标）✓
4. **命名方式**：自定义（与OSPREY兼容）✓

## 验证步骤

### 1. 检查template和coords是否匹配：
```bash
# 检查K03原子名称一致性
grep "^   [0-9]" K03_template.in | awk '{print $2}' | sort > /tmp/template_atoms.txt
head -27 K03_coords.in | grep -v "^ALY" | awk '{print $1}' | sort > /tmp/coords_atoms.txt
diff /tmp/template_atoms.txt /tmp/coords_atoms.txt
# 应该无差异

# 检查K04原子名称一致性
grep "^   [0-9]" K04_template.in | awk '{print $2}' | sort > /tmp/template_atoms.txt
head -32 K04_coords.in | grep -v "^SLL" | awk '{print $1}' | sort > /tmp/coords_atoms.txt
diff /tmp/template_atoms.txt /tmp/coords_atoms.txt
# 应该无差异
```

### 2. 验证电荷总和：
```bash
# ALY应该≈0
python3 -c "
charges = []
with open('K03_coords.in') as f:
    for line in f:
        if line.startswith('ENDRES'): break
        if line.startswith('ALY'): continue
        parts = line.split()
        if len(parts) >= 5:
            charges.append(float(parts[4].rstrip('f')))
print(f'ALY total charge: {sum(charges):.4f}')
"

# SLL应该≈-1
python3 -c "
charges = []
with open('K04_coords.in') as f:
    for line in f:
        if line.startswith('ENDRES'): break
        if line.startswith('SLL'): continue
        parts = line.split()
        if len(parts) >= 5:
            charges.append(float(parts[4].rstrip('f')))
print(f'SLL total charge: {sum(charges):.4f}')
"
```

### 3. 运行OSPREY测试：
```bash
sbatch run_paste.slurm
# 检查sample-*.out中是否有错误
# 应该看到：
# - 成功加载K03/K04 templates
# - 无"NoWindowOverlap"警告
# - 合理的ΔΔG值（SLL: +3-4 kcal/mol，不是+25）
```

## 总结

**你问的"哪个版本更好"？答案是：**

当前使用的**混合版本**最好：
- Template来自手工创建（CIF几何）
- Coords电荷来自ANTECHAMBER（量子化学）
- 命名统一（自定义，OSPREY兼容）

这是最佳实践，结合了：
1. 最准确的几何结构（CIF理想坐标）
2. 最准确的部分电荷（AM1-BCC量子化学）
3. 最好的兼容性（命名一致）

**不要使用**纯ANTECHAMBER版本的template.in（会导致命名冲突）。

## 相关文件位置

```
当前工作文件：
  /home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE/K03_*.in
  /home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE/K04_*.in

备份（各个版本）：
  backup_before_antechamber/  - 手工CIF版本（template和coords）
  K03_coords_before_antechamber.in  - 更新电荷前的coords
  K04_coords_before_antechamber.in  - 更新电荷前的coords
  K03_template_from_prepi.in  - ANTECHAMBER版本（不推荐，命名冲突）
  K04_template_from_prepi.in  - ANTECHAMBER版本（不推荐，命名冲突）

ANTECHAMBER输出：
  ALY.prepi, SLL.prepi  - AM1-BCC电荷来源
```
