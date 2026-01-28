# Bug修复总结

## 发现的Bug

**错误信息:**
```
java.lang.IndexOutOfBoundsException: Index 1 out of bounds for length 0
```

**原因:**
OSPREY的模板文件解析器对文件格式有严格要求。模板文件的标题行格式不正确导致解析失败。

## 已修复的问题

### 1. K03_template.in (Acetyl-Lysine)
- ❌ 之前: `ACETYL-LYSINE`
- ✅ 修复: `N6-ACETYL-LYSINE`

### 2. K04_template.in (Succinyl-Lysine)  
- ❌ 之前: `SUCCINYL-LYSINE`
- ✅ 修复: `N6-SUCCINYL-LYSINE`

### 3. PASTE_6dv2.py
已更新以支持两种修饰:
- 添加了K04 (succinyl-lysine) 模板文件
- B781位置现在可以是: LYS (野生型), ALY (acetyl-lysine), SLL (succinyl-lysine)

## 文件清单

✅ 所有必需的文件已在当前目录:

**Acetyl-Lysine (ALY):**
- K03_template.in
- K03_coords.in  
- K03_rotlib.dat

**Succinyl-Lysine (SLL):**
- K04_template.in
- K04_coords.in
- K04_rotlib.dat

**主脚本:**
- PASTE_6dv2.py (已更新)

## 如何运行

### 方法1: 直接测试 (推荐先运行)
```bash
conda activate AmberTools22
python test_templates.py
```

这会测试模板是否能正确加载。

### 方法2: 提交完整任务
```bash
sbatch submit_paste.py
```

## 关键代码变更

### PASTE_6dv2.py 第12-16行:
```python
customizedTemplateLib = osprey.TemplateLibrary(
    extraTemplates=['K03_template.in', 'K04_template.in'],
    extraTemplateCoords=['K03_coords.in', 'K04_coords.in'],
    extraRotamers=['K03_rotlib.dat', 'K04_rotlib.dat']
)
```

### PASTE_6dv2.py 第21行:
```python
protein.flexibility['B781'].setLibraryRotamers(osprey.WILD_TYPE, 'ALY', 'SLL').addWildTypeRotamers().setContinuous()
```

## 设计空间说明

现在B781位置有3种选择:
1. **LYS** - 野生型赖氨酸
2. **ALY** - N6-acetyl-lysine (乙酰化赖氨酸)
3. **SLL** - N6-succinyl-lysine (琥珀酰化赖氨酸)

每种类型都有对应的rotamer库，允许侧链构象优化。

## 预期输出

运行成功后，PAStE会:
1. 计算energy matrix
2. 搜索最优序列
3. 输出结果到 `paste.results.tsv`
4. 生成PDB结构文件

## 故障排查

如果仍有错误:

1. **文件未找到**: 确认所有K03和K04文件在同一目录
2. **模板解析错误**: 检查.in文件格式没有被意外修改
3. **内存不足**: 调整submit_paste.py中的--mem参数

## 技术细节

- **Force Field**: AMBER
- **坐标来源**: PDB Chemical Component Dictionary (ALY, SLL)
- **Rotamer库**: 基于Dunbrack库，扩展了修饰基团的dihedral angles

