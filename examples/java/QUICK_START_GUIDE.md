# COMETSZ + MARKStar 快速入门指南

## 文件说明

我为您创建了两个 Java 示例：

### 1. `SimpleCometsZMARKStarExample.java` ⭐ 推荐新手使用
- **简化版本**，易于理解
- 完整的中文注释
- 可直接运行
- 基于 TestCometsZWithBBKStarAndMARKStar 测试代码

### 2. `CometsZBBKStarMARKStarExample.java`
- **完整版本**，包含更多细节
- 英文注释，更专业
- 模块化设计，便于扩展
- 适合深入学习

### 3. `COMETSZ_BBKSTAR_MARKSTAR_README.md`
- **详细文档**，包含理论背景
- 参数调优指南
- 常见问题解答
- 性能优化建议

## 快速开始

### 第一步：确保环境准备就绪

```bash
# 检查 Java 版本 (需要 Java 11+)
java -version

# 检查项目是否编译
cd /home/users/lz280/IdeaProjects/OSPREY3
./gradlew build
```

### 第二步：准备 PDB 文件

确保您有以下文件：
```
src/main/resources/2RL0.min.reduce.pdb
```

或者修改代码中的 PDB 文件路径为您自己的结构。

### 第三步：运行简化示例

```bash
# 方法 1: 使用 IntelliJ IDEA
# 1. 打开 SimpleCometsZMARKStarExample.java
# 2. 右键点击 main 方法
# 3. 选择 "Run 'SimpleCometsZMARKStarExample.main()'"

# 方法 2: 使用命令行
cd examples/java
javac -cp ../../build/libs/osprey.jar SimpleCometsZMARKStarExample.java
java -cp ../../build/libs/osprey.jar:. SimpleCometsZMARKStarExample
```

## 代码核心逻辑

### 1. 定义三个状态

```java
CometsZ.State protein = new CometsZ.State("Protein", proteinSpace);
CometsZ.State ligand = new CometsZ.State("Ligand", ligandSpace);
CometsZ.State complex = new CometsZ.State("Complex", complexSpace);
```

### 2. 设置目标函数（结合自由能）

```java
CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
    .addState(complex, 1.0)    // +G_complex
    .addState(protein, -1.0)   // -G_protein
    .addState(ligand, -1.0)    // -G_ligand
    .build();
```

### 3. 配置 MARKStar

```java
// 为每个状态创建两种能量矩阵
EnergyMatrix rigidEmat = ...;      // 刚性（快速，上界）
EnergyMatrix minimizingEmat = ...; // 柔性（准确，下界）

// 创建 MARKStar 配分函数
state.pfuncFactory = (rcs) -> {
    MARKStarBound markstar = new MARKStarBound(
        state.confSpace,
        rigidEmat,
        minimizingEmat,
        state.confEcalc,
        rcs,
        parallelism
    );
    markstar.setCorrections(correctionMatrix);
    return markstar;
};
```

### 4. 运行设计

```java
List<CometsZ.SequenceInfo> sequences = cometsZ.findBestSequences(5);
```

## 输出解读

### 控制台输出示例

```
序列 #1:
  G649=asp G650=glu G651=ile G654=val A156=thr A172=asn A192=ser A193=phe
  结合自由能 (ΔG): [-45.2341, -42.8765] kcal/mol
    Complex: [-403.2586, -396.8900] kcal/mol
    Protein: [-234.5678, -232.1234] kcal/mol
    Ligand: [-123.4567, -121.8901] kcal/mol
```

### 结果含义

- **结合自由能** (ΔG): 越负越好（表示更强的结合）
- **置信区间**: [下界, 上界] 由 epsilon 控制
- **各状态自由能**: 独立计算，用于验证

### TSV 文件

`cometsz.results.tsv` 包含所有序列的详细结果，格式：

```
G649    G650    G651    ...    Objective Min    Objective Max    ...
asp     glu     ile     ...    -45.2341         -42.8765         ...
```

## 关键参数说明

### Epsilon（精度控制）

```java
.setEpsilon(0.95)
```

- **0.90**: 快速探索，较低精度
- **0.95**: 推荐值，平衡速度和精度 ⭐
- **0.99**: 高精度，计算时间长

### 最大同时突变数

```java
.setMaxSimultaneousMutations(1)
```

- **1**: 每次只突变一个位点（保守）⭐
- **2**: 同时考虑两个位点的突变（中等）
- **3+**: 更全面但更慢

### 能量窗口

```java
.setObjectiveWindowSize(100.0)   // 相对窗口
.setObjectiveWindowMax(100.0)    // 绝对窗口
```

- 控制搜索范围
- 值越小，搜索越集中在最优解附近
- 值越大，探索更多次优解

### CPU 核心数

```java
Parallelism.makeCpu(4)
```

- 根据您的机器配置调整
- 推荐使用 4-8 核心
- 更多核心 = 更快计算

## 常见问题

### Q1: 如何修改设计位点？

修改这部分代码：

```java
protein.flexibility.get("G649")
    .setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
    .addWildTypeRotamers()
    .setContinuous();
```

### Q2: 如何使用自己的 PDB 文件？

```java
// 从文件读取
Molecule mol = PDBIO.readFile("path/to/your.pdb");

// 从资源读取
Molecule mol = PDBIO.readResource("/your.pdb");
```

### Q3: 内存不足怎么办？

```java
cometsZ.Builder(objective)
    .setMinNumConfTrees(50)  // 限制内存中的树数量
    .build();
```

### Q4: 如何加快计算？

1. 增加 CPU 核心: `Parallelism.makeCpu(8)`
2. 降低精度: `setEpsilon(0.90)`
3. 减少设计位点
4. 减少每个位点的旋转异构体数量

### Q5: 结果不收敛怎么办？

检查：
1. epsilon 是否太高（如 0.99+）
2. 能量窗口是否太大
3. 是否有能量计算问题

## 与其他方法的对比

| 方法 | 速度 | 精度 | 多状态 | 严格界限 |
|------|------|------|--------|----------|
| 传统 K* | 慢 | 高 | ❌ | ✅ |
| BBK* | 快 | 高 | ❌ | ✅ |
| COMETS | 中等 | 高 | ✅ | ✅ |
| COMETSZ + MARKStar | **很快** | **高** | **✅** | **✅** |

## 性能提示

### 首次运行

```
预计时间：5-30 分钟（取决于系统大小）
```

- 需要计算能量矩阵（会自动缓存）
- 初始化 MARKStar

### 后续运行

```
预计时间：1-5 分钟
```

- 使用缓存的能量矩阵
- 跳过初始化步骤

### 大型系统

对于 10+ 个灵活位点的系统：

```java
// 使用 ConfDB 缓存
state.confDBFile = new File("confdb." + state.name + ".db");

// 限制内存
.setMinNumConfTrees(100)

// 使用更多 CPU
Parallelism.makeCpu(16)
```

## 下一步

### 学习更多

1. 阅读 `COMETSZ_BBKSTAR_MARKSTAR_README.md` 了解理论背景
2. 查看 `TestCometsZWithBBKStarAndMARKStar.java` 了解测试用例
3. 研究 `CometsZ.java` 了解实现细节

### 扩展示例

1. 添加约束条件：
```java
CometsZ.LMFE constraint = new CometsZ.LMFE.Builder()
    .addState(protein, 1.0)
    .constrainLessThan(-100.0)  // 蛋白质自由能必须 < -100
    .build();

cometsZ.Builder(objective)
    .addConstraint(constraint)
    .build();
```

2. 多个蛋白质状态：
```java
CometsZ.State proteinState1 = new CometsZ.State("Protein1", space1);
CometsZ.State proteinState2 = new CometsZ.State("Protein2", space2);
// ... 定义更复杂的目标函数
```

## 获得帮助

如果遇到问题：

1. 检查日志输出
2. 查看 TSV 结果文件
3. 参考测试代码 `TestCometsZWithBBKStarAndMARKStar.java`
4. 查阅 OSPREY 文档

## 总结

您现在拥有：

✅ 两个完整的可运行示例
✅ 详细的中文文档
✅ 快速入门指南（本文件）
✅ 参数调优建议
✅ 常见问题解答

开始您的蛋白质设计之旅吧！🚀
