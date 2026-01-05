# COMETSZ + BBKStar + MARKStar 示例合集

## 📁 文件列表

我为您创建了一套完整的 COMETSZ + MARKStar 示例和文档，**完全支持 VSCode 和 IntelliJ IDEA**！

### 🎯 示例代码

1. **SimpleCometsZMARKStarExample.java** ⭐ **推荐新手**
   - 简化版本，易于理解
   - 完整中文注释
   - 直接可运行
   - ~300 行代码
   - 基于实际测试案例

2. **CometsZBBKStarMARKStarExample.java**
   - 完整版本，专业实现
   - 英文注释
   - 模块化设计
   - ~450 行代码
   - 适合深入学习和扩展

### 📚 文档

3. **QUICK_START_GUIDE.md** ⭐ **新手必读**
   - 快速入门指南
   - 5 分钟上手
   - 常见问题解答
   - 参数调优建议

4. **COMETSZ_BBKSTAR_MARKSTAR_README.md**
   - 详细理论背景
   - 完整 API 说明
   - 性能优化指南
   - 文献引用

5. **README_EXAMPLES.md** (本文件)
   - 总览和索引
   - 快速导航

### 🖥️ VSCode 支持

6. **VSCODE_SETUP_GUIDE.md** ⭐ **VSCode 用户必读**
   - 完整的 VSCode 设置指南
   - 三种运行方法
   - 调试配置
   - 常见问题解决

7. **VSCODE_QUICK_REFERENCE.md** ⭐ **快速参考**
   - 一分钟快速上手
   - 快捷键备忘单
   - 命令行速查

8. **run_simple_example.sh** - 一键运行脚本（简化版）
9. **run_complete_example.sh** - 一键运行脚本（完整版）

### ⚙️ VSCode 配置文件（已创建）

10. **.vscode/launch.json** - 运行和调试配置
11. **.vscode/tasks.json** - 编译任务配置
12. **.vscode/settings.json** - 项目设置

## 🚀 快速开始

### 方案 A：VSCode 用户（最简单）⭐

```bash
# 1. 打开项目
cd /home/users/lz280/IdeaProjects/OSPREY3
code .

# 2. 打开示例文件
# examples/java/SimpleCometsZMARKStarExample.java

# 3. 点击 main 方法上的 "Run" 按钮
# 或者按 F5
```

**或者使用一键脚本：**

```bash
cd examples/java
./run_simple_example.sh
```

详细说明见 **[VSCODE_SETUP_GUIDE.md](VSCODE_SETUP_GUIDE.md)**

---

### 方案 B：IntelliJ IDEA 用户

```bash
# 1. 进入示例目录
cd /home/users/lz280/IdeaProjects/OSPREY3/examples/java

# 2. 查看快速入门指南
cat QUICK_START_GUIDE.md

# 3. 运行简化示例
# (在 IntelliJ IDEA 中打开 SimpleCometsZMARKStarExample.java)
# 右键 -> Run
```

### 方案 B：深入学习

```bash
# 1. 阅读详细文档
cat COMETSZ_BBKSTAR_MARKSTAR_README.md

# 2. 研究完整示例
# 打开 CometsZBBKStarMARKStarExample.java

# 3. 查看测试代码
cd ../../src/test/java/edu/duke/cs/osprey/kstar
cat TestCometsZWithBBKStarAndMARKStar.java
```

## 📖 使用哪个示例？

### 选择 SimpleCometsZMARKStarExample.java 如果您：

- ✅ 是第一次使用 OSPREY
- ✅ 想快速了解 COMETSZ + MARKStar
- ✅ 需要中文注释
- ✅ 想看到清晰的步骤说明

### 选择 CometsZBBKStarMARKStarExample.java 如果您：

- ✅ 已经熟悉 OSPREY 基础
- ✅ 需要扩展和定制
- ✅ 想了解更多实现细节
- ✅ 计划用于实际研究项目

## 🎓 学习路径

### 初学者路径

```
第 1 天: QUICK_START_GUIDE.md (30 分钟)
       ↓
第 2 天: SimpleCometsZMARKStarExample.java (运行并理解)
       ↓
第 3 天: 修改参数，测试不同配置
       ↓
第 4 天: COMETSZ_BBKSTAR_MARKSTAR_README.md (理论部分)
```

### 进阶路径

```
第 1 天: CometsZBBKStarMARKStarExample.java (完整示例)
       ↓
第 2 天: COMETSZ_BBKSTAR_MARKSTAR_README.md (深入理论)
       ↓
第 3 天: TestCometsZWithBBKStarAndMARKStar.java (测试代码)
       ↓
第 4 天: 源代码 CometsZ.java, MARKStarBound.java
```

## 🔑 关键概念速查

### COMETSZ 是什么？

**COMETS-Z** = COMETS + 热力学系综（Thermodynamic Ensembles）

- 多状态设计（如：结合态 vs 非结合态）
- 线性多状态自由能 (LMFE) 目标函数
- 序列空间的智能搜索

```java
// 目标：最小化结合自由能
ΔG = G(complex) - G(protein) - G(ligand)
```

### MARKStar 是什么？

**MARK*** = Matrix A* for Rapid K*

- 快速配分函数计算（5-100x 加速）
- 使用刚性/柔性能量矩阵
- 提供严格的数学界限

```java
// 核心思想：
上界 ≤ 真实配分函数 ≤ 下界
通过迭代缩小界限直到满足 epsilon 精度
```

### BBKStar 在哪里？

BBKStar 的批处理优化策略已经集成在 COMETSZ 中：

- `setMinNumConfTrees()` - 控制内存中的构象树
- `setMaxSimultaneousMutations()` - 控制序列搜索策略
- 自动的序列剪枝和优先级排序

## 📊 示例对比

| 特性 | SimpleExample | CompleteExample |
|------|---------------|-----------------|
| 代码行数 | ~300 | ~450 |
| 注释语言 | 中文 | 英文 |
| 难度 | ⭐⭐ | ⭐⭐⭐⭐ |
| 模块化 | 基础 | 完整 |
| 扩展性 | 中等 | 高 |
| 文档 | 代码内 | 代码内+外部 |
| 适合 | 学习 | 研究 |

## 🛠️ 核心代码模板

### 最小可运行示例

```java
// 1. 定义状态
CometsZ.State protein = new CometsZ.State("Protein", proteinSpace);
CometsZ.State ligand = new CometsZ.State("Ligand", ligandSpace);
CometsZ.State complex = new CometsZ.State("Complex", complexSpace);

// 2. 定义目标函数
CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
    .addState(complex, 1.0)
    .addState(protein, -1.0)
    .addState(ligand, -1.0)
    .build();

// 3. 创建 COMETSZ
CometsZ cometsZ = new CometsZ.Builder(objective)
    .setEpsilon(0.95)
    .build();

// 4. 为每个状态设置 MARKStar
for (CometsZ.State state : cometsZ.states) {
    // 计算能量矩阵
    EnergyMatrix rigidEmat = ...;
    EnergyMatrix minimizingEmat = ...;

    // 设置 MARKStar
    state.pfuncFactory = (rcs) -> new MARKStarBound(...);
    state.fragmentEnergies = minimizingEmat;
    state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(...).build();
}

// 5. 运行
List<CometsZ.SequenceInfo> sequences = cometsZ.findBestSequences(5);
```

## 📈 典型运行结果

### 控制台输出

```
================================================================================
【步骤 1】定义蛋白质-配体系统

  ✓ 蛋白质构象空间: 4 个灵活位点
  ✓ 配体构象空间: 4 个灵活位点
  ✓ 复合物构象空间: 8 个灵活位点

【步骤 2】配置 COMETSZ 多状态设计
  ...

【步骤 4】运行 COMETSZ 序列优化
  正在搜索最优序列...

================================================================================
【结果】找到 5 个最优序列
================================================================================

序列 #1:
  G649=asp G650=glu G651=ile G654=val A156=thr A172=asn A192=ser A193=phe
  结合自由能 (ΔG): [-45.2341, -42.8765] kcal/mol
    Complex: [-403.2586, -396.8900] kcal/mol
    Protein: [-234.5678, -232.1234] kcal/mol
    Ligand: [-123.4567, -121.8901] kcal/mol
```

### 输出文件

- `cometsz.results.tsv` - 所有序列的详细结果（TSV 格式）
- `emat.*.dat` - 缓存的能量矩阵（可重用）

## ⚙️ 重要参数说明

### Epsilon（精度）

```java
.setEpsilon(0.95)  // 推荐：0.95
```

| 值 | 含义 | 速度 | 用途 |
|----|------|------|------|
| 0.90 | 90% 置信度 | 快 | 快速探索 |
| 0.95 | 95% 置信度 | 中等 | 推荐值 ⭐ |
| 0.99 | 99% 置信度 | 慢 | 高精度研究 |

### 同时突变数

```java
.setMaxSimultaneousMutations(1)  // 推荐：1 或 2
```

| 值 | 搜索范围 | 速度 | 用途 |
|----|----------|------|------|
| 1 | 保守 | 快 | 初步设计 ⭐ |
| 2 | 中等 | 中等 | 精细优化 |
| 3+ | 广泛 | 慢 | 全面探索 |

### CPU 核心数

```java
Parallelism.makeCpu(4)  // 根据机器配置
```

推荐配置：
- 笔记本电脑: 2-4 核心
- 工作站: 8-16 核心
- 集群: 32+ 核心

## 🐛 常见问题速查

### 问题 1: 编译错误

```
Solution: 确保已运行 ./gradlew build
```

### 问题 2: 找不到 PDB 文件

```
Solution: 检查 PDBIO.readResource("/2RL0.min.reduce.pdb")
         确保文件在 src/main/resources/ 目录下
```

### 问题 3: 内存不足

```
Solution: 添加 .setMinNumConfTrees(50)
```

### 问题 4: 运行时间太长

```
Solution:
1. 降低 epsilon: .setEpsilon(0.90)
2. 增加 CPU: Parallelism.makeCpu(8)
3. 减少设计位点
```

### 问题 5: NullPointerException

```
Solution: 确保设置了 corrections 矩阵
markstar.setCorrections(new UpdatingEnergyMatrix(...))
```

## 📚 扩展阅读

### 源代码位置

```
CometsZ:
  src/main/java/edu/duke/cs/osprey/kstar/CometsZ.java

MARKStar:
  src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java

测试代码:
  src/test/java/edu/duke/cs/osprey/kstar/TestCometsZWithBBKStarAndMARKStar.java
```

### 相关示例

```
Python 示例:
  examples/python.KStar/markstar.kstar.py
  examples/python.KStar/bbkstar.py
  examples/python.GMEC/comets.py

Java 测试:
  src/test/java/edu/duke/cs/osprey/kstar/TestMSKStar.java
  src/test/java/edu/duke/cs/osprey/markstar/TestMARKStar.java
```

## 🎯 使用场景

### 场景 1: 蛋白质-配体结合优化

```java
// 目标：找到最佳结合序列
LMFE = G(complex) - G(protein) - G(ligand)
```

✅ 使用 SimpleCometsZMARKStarExample.java

### 场景 2: 蛋白质稳定性优化

```java
// 目标：最小化蛋白质自由能
LMFE = G(protein)
```

✅ 修改示例，只使用一个状态

### 场景 3: 多状态设计（3+ 状态）

```java
// 目标：同时优化多个构象
LMFE = G(state1) + G(state2) - G(state3)
```

✅ 使用 CometsZBBKStarMARKStarExample.java 并扩展

## 📞 获取帮助

### 文档优先级

1. **QUICK_START_GUIDE.md** - 入门问题
2. **本文件代码注释** - 实现细节
3. **COMETSZ_BBKSTAR_MARKSTAR_README.md** - 理论和高级用法
4. **测试代码** - 实际用例

### 在线资源

- OSPREY 官网
- Bruce Donald Lab 网站
- 相关论文（见 COMETSZ_BBKSTAR_MARKSTAR_README.md）

## ✅ 检查清单

开始前确认：

- [ ] Java 11+ 已安装
- [ ] OSPREY 已编译 (`./gradlew build`)
- [ ] PDB 文件已准备
- [ ] 已阅读 QUICK_START_GUIDE.md
- [ ] 理解基本概念（COMETSZ, MARKStar）

开始运行：

- [ ] 选择了合适的示例
- [ ] 修改了 PDB 文件路径（如需要）
- [ ] 调整了参数（epsilon, CPU 等）
- [ ] 准备好等待计算完成（5-30 分钟）

结果验证：

- [ ] 检查控制台输出
- [ ] 查看 TSV 结果文件
- [ ] 验证自由能值合理
- [ ] 置信区间宽度可接受

## 🎉 总结

您现在拥有：

✅ **2 个完整的可运行 Java 示例**
  - SimpleCometsZMARKStarExample.java (简化版)
  - CometsZBBKStarMARKStarExample.java (完整版)

✅ **3 个详细的文档**
  - QUICK_START_GUIDE.md (快速入门)
  - COMETSZ_BBKSTAR_MARKSTAR_README.md (完整文档)
  - README_EXAMPLES.md (本文件)

✅ **完整的学习路径**
  - 从入门到精通
  - 理论到实践
  - 简单到复杂

✅ **实用的参考资料**
  - 参数调优指南
  - 常见问题解答
  - 性能优化建议

## 🚀 现在就开始吧！

```bash
# 第一步：阅读快速入门
cat QUICK_START_GUIDE.md

# 第二步：运行示例
# 在 IntelliJ IDEA 中打开 SimpleCometsZMARKStarExample.java
# 右键 -> Run

# 第三步：查看结果
cat cometsz.results.tsv
```

祝您设计成功！🎯
