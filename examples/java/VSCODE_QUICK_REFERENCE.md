# VSCode 快速参考卡

## 🚀 三种运行方式（从易到难）

### 方式 1: 图形界面运行 ⭐ **最简单**

```
1. 在 VSCode 中打开 SimpleCometsZMARKStarExample.java
2. 找到 main 方法
3. 点击上方的 "Run" 链接
4. 完成！
```

**优点**: 零配置，一键运行
**缺点**: 首次可能需要等待 Java 扩展加载

---

### 方式 2: 使用运行脚本 ⭐ **推荐**

在 VSCode 终端中（按 `` Ctrl+` ``）：

```bash
cd examples/java
./run_simple_example.sh
```

**优点**: 自动检查环境，有详细输出
**缺点**: 需要有执行权限

---

### 方式 3: 手动命令行 ⭐ **最灵活**

```bash
cd examples/java

# 编译
javac -cp ../../build/libs/*:../../build/classes/java/main \
    SimpleCometsZMARKStarExample.java

# 运行
java -Xmx4g \
    -cp ../../build/libs/*:../../build/classes/java/main:. \
    SimpleCometsZMARKStarExample
```

**优点**: 完全控制，易于调试
**缺点**: 需要记住命令

---

## 🎯 VSCode 快捷键

```
打开终端:         Ctrl + `
运行程序:         F5 (调试) 或 Ctrl+F5 (运行)
停止程序:         Shift+F5
命令面板:         Ctrl+Shift+P
查找文件:         Ctrl+P
全局搜索:         Ctrl+Shift+F
```

---

## 📋 使用 VSCode 任务

按 `Ctrl+Shift+P`，输入 "Tasks: Run Task"，选择：

- **Build OSPREY** - 编译整个项目
- **Compile Simple Example** - 只编译示例
- **Run Simple Example (Script)** - 运行简化版
- **Run Complete Example (Script)** - 运行完整版

---

## 🐛 调试

### 设置断点

1. 在代码行号左侧点击（出现红点）
2. 按 `F5` 启动调试
3. 程序会在断点处暂停

### 调试控制

```
继续运行:    F5
单步跳过:    F10
单步进入:    F11
单步跳出:    Shift+F11
```

### 查看变量

- 鼠标悬停在变量上查看值
- 在左侧 "Variables" 面板查看所有变量
- 在 "Watch" 面板添加表达式

---

## ⚡ 快速启动（新项目）

```bash
# 1. 安装 VSCode Java 扩展
# 在 VSCode 中: Ctrl+Shift+X
# 搜索: "Extension Pack for Java"
# 点击 Install

# 2. 打开项目
cd /home/users/lz280/IdeaProjects/OSPREY3
code .

# 3. 编译项目
./gradlew build

# 4. 运行示例
cd examples/java
./run_simple_example.sh
```

---

## 📁 文件位置

```
.vscode/
├── launch.json      ← 运行配置（已创建）
├── tasks.json       ← 任务配置（已创建）
└── settings.json    ← 项目设置（已创建）

examples/java/
├── SimpleCometsZMARKStarExample.java       ← 简化版示例
├── CometsZBBKStarMARKStarExample.java      ← 完整版示例
├── run_simple_example.sh                   ← 运行脚本
├── run_complete_example.sh                 ← 运行脚本
├── VSCODE_SETUP_GUIDE.md                   ← 详细设置指南
└── QUICK_START_GUIDE.md                    ← 快速入门
```

---

## 🔍 常见问题 1 分钟解决

### Q: 看不到 "Run" 按钮？

```bash
# 解决方案：
1. Ctrl+Shift+P
2. 输入 "Java: Clean Java Language Server Workspace"
3. 重启 VSCode
```

### Q: 编译失败？

```bash
# 解决方案：
cd /home/users/lz280/IdeaProjects/OSPREY3
./gradlew clean build
```

### Q: 找不到类？

```bash
# 解决方案：检查 classpath
ls build/libs/
# 应该看到 .jar 文件

# 如果没有，运行：
./gradlew build
```

### Q: 内存不足？

编辑 `.vscode/launch.json`：
```json
"vmArgs": "-Xmx8g -Xms2g"  // 增加到 8GB
```

---

## 📊 两个示例的区别

| 特性 | Simple | Complete |
|------|--------|----------|
| 代码行数 | ~300 | ~450 |
| 注释语言 | 中文 | 英文 |
| 运行时间 | 5-20分钟 | 10-30分钟 |
| 适合 | 学习 | 研究 |
| 推荐 | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |

---

## 💡 最佳实践

### 第一次运行

```bash
# 1. 先编译 OSPREY
./gradlew build

# 2. 运行简化版示例
cd examples/java
./run_simple_example.sh

# 3. 查看结果
cat cometsz.results.tsv
```

### 修改参数后

```bash
# 1. 只需重新编译示例
javac -cp ../../build/libs/*:../../build/classes/java/main \
    SimpleCometsZMARKStarExample.java

# 2. 运行
java -Xmx4g -cp ../../build/libs/*:../../build/classes/java/main:. \
    SimpleCometsZMARKStarExample
```

### 清理缓存

```bash
# 删除能量矩阵缓存（重新计算）
rm emat.*.dat

# 删除结果文件
rm cometsz*.tsv
```

---

## 🎓 学习路径

### Day 1: 环境设置（10 分钟）
```
✓ 安装 Java 扩展
✓ 编译项目
✓ 运行示例
```

### Day 2: 理解代码（30 分钟）
```
✓ 阅读 SimpleCometsZMARKStarExample.java
✓ 理解 5 个步骤
✓ 查看输出结果
```

### Day 3: 修改参数（1 小时）
```
✓ 修改 epsilon 值
✓ 修改设计位点
✓ 比较结果
```

### Day 4: 深入学习（2 小时）
```
✓ 阅读 COMETSZ_BBKSTAR_MARKSTAR_README.md
✓ 研究 CometsZBBKStarMARKStarExample.java
✓ 理解理论背景
```

---

## 🔗 相关文档

- **[VSCODE_SETUP_GUIDE.md](VSCODE_SETUP_GUIDE.md)** - 详细设置指南
- **[QUICK_START_GUIDE.md](QUICK_START_GUIDE.md)** - 快速入门
- **[COMETSZ_BBKSTAR_MARKSTAR_README.md](COMETSZ_BBKSTAR_MARKSTAR_README.md)** - 理论文档
- **[README_EXAMPLES.md](README_EXAMPLES.md)** - 示例总览

---

## ⌨️ VSCode 命令备忘单

```bash
# Java 相关
Ctrl+Shift+P → "Java: Clean Java Language Server Workspace"
Ctrl+Shift+P → "Java: Force Java Compilation"
Ctrl+Shift+P → "Java: Open Java Language Server Log File"

# 任务相关
Ctrl+Shift+P → "Tasks: Run Task"
Ctrl+Shift+P → "Tasks: Run Build Task"

# Git 相关
Ctrl+Shift+G → 打开 Git 面板
Ctrl+Shift+P → "Git: Commit"
Ctrl+Shift+P → "Git: Push"

# 终端相关
Ctrl+`       → 打开/关闭终端
Ctrl+Shift+` → 新建终端
```

---

## 🎉 快速开始命令（复制粘贴）

### 初次运行

```bash
# 一键运行（在 VSCode 终端中）
cd /home/users/lz280/IdeaProjects/OSPREY3/examples/java && ./run_simple_example.sh
```

### 后续运行

```bash
# 快速运行（在 examples/java 目录下）
./run_simple_example.sh
```

### 调试运行

```bash
# 编译后在 VSCode 中按 F5，选择 "Run Simple COMETSZ Example"
```

---

**提示**: 将此文件保存为书签，方便随时查阅！
