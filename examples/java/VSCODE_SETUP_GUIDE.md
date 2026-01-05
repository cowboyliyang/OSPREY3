# VSCode 运行指南

## 在 VSCode 中运行 COMETSZ + MARKStar 示例

完全可以！VSCode 是一个很好的选择。这里是完整的设置和运行指南。

## 🔧 前置要求

### 1. 安装 VSCode 扩展

打开 VSCode，安装以下扩展：

```
必需扩展：
- Extension Pack for Java (Microsoft)
  包含：
  - Language Support for Java(TM) by Red Hat
  - Debugger for Java
  - Test Runner for Java
  - Maven for Java
  - Project Manager for Java
  - Visual Studio IntelliCode

可选但推荐：
- Gradle for Java (Microsoft)
```

### 2. 确认 Java 环境

```bash
# 检查 Java 版本（需要 11+）
java -version

# 检查 JAVA_HOME 环境变量
echo $JAVA_HOME
```

如果没有设置 JAVA_HOME：
```bash
# Linux/Mac
export JAVA_HOME=/path/to/your/jdk
export PATH=$JAVA_HOME/bin:$PATH

# 或者添加到 ~/.bashrc 或 ~/.zshrc
```

## 📁 VSCode 项目设置

### 方法 1: 使用 Gradle（推荐）

#### 步骤 1: 在 VSCode 中打开项目

```bash
cd /home/users/lz280/IdeaProjects/OSPREY3
code .
```

#### 步骤 2: VSCode 会自动识别 Gradle 项目

VSCode 会在右下角显示 "Importing Gradle Project..."，等待完成。

#### 步骤 3: 编译项目

打开终端（Terminal → New Terminal）：

```bash
# 首次编译
./gradlew build

# 或者点击 VSCode 侧边栏的 Gradle 图标
# 展开 Tasks → build → build
```

### 方法 2: 直接编译和运行（快速方法）

#### 创建运行脚本

我为您创建了一个一键运行脚本：

**run_simple_example.sh**
```bash
#!/bin/bash

echo "==================================="
echo "编译并运行 SimpleCometsZMARKStarExample"
echo "==================================="

# 设置路径
PROJECT_ROOT="/home/users/lz280/IdeaProjects/OSPREY3"
EXAMPLE_DIR="$PROJECT_ROOT/examples/java"
BUILD_DIR="$PROJECT_ROOT/build"

# 进入项目根目录
cd "$PROJECT_ROOT"

# 编译 OSPREY（如果还没编译）
if [ ! -d "$BUILD_DIR/libs" ]; then
    echo "正在编译 OSPREY..."
    ./gradlew build
fi

# 编译示例
echo "正在编译示例..."
cd "$EXAMPLE_DIR"

javac -cp "$BUILD_DIR/libs/*:$BUILD_DIR/classes/java/main" SimpleCometsZMARKStarExample.java

# 运行示例
echo "正在运行示例..."
java -cp "$BUILD_DIR/libs/*:$BUILD_DIR/classes/java/main:." SimpleCometsZMARKStarExample

echo "运行完成！"
```

**使用方法：**
```bash
chmod +x run_simple_example.sh
./run_simple_example.sh
```

## 🚀 在 VSCode 中运行示例

### 方法 A: 使用 VSCode 的 Run 按钮（最简单）⭐

1. 在 VSCode 中打开 `SimpleCometsZMARKStarExample.java`

2. 你会看到 `main` 方法上方有一个 "Run | Debug" 的链接

3. 点击 **Run**

4. 查看终端输出

**如果没有看到 Run 按钮：**
- 确保安装了 "Extension Pack for Java"
- 右键点击代码 → "Run Java"

### 方法 B: 使用 VSCode 终端

在 VSCode 中打开终端（`` Ctrl+` `` 或 Terminal → New Terminal）：

```bash
# 进入示例目录
cd examples/java

# 编译
javac -cp ../../build/libs/*:../../build/classes/java/main SimpleCometsZMARKStarExample.java

# 运行
java -cp ../../build/libs/*:../../build/classes/java/main:. SimpleCometsZMARKStarExample
```

### 方法 C: 配置 launch.json（专业方法）

#### 步骤 1: 创建 launch.json

在项目根目录创建或编辑 `.vscode/launch.json`：

```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "type": "java",
            "name": "Run Simple COMETSZ Example",
            "request": "launch",
            "mainClass": "SimpleCometsZMARKStarExample",
            "projectName": "OSPREY3",
            "cwd": "${workspaceFolder}/examples/java",
            "classPaths": [
                "${workspaceFolder}/build/libs/*",
                "${workspaceFolder}/build/classes/java/main",
                "${workspaceFolder}/examples/java"
            ],
            "vmArgs": "-Xmx4g"
        },
        {
            "type": "java",
            "name": "Run Complete COMETSZ Example",
            "request": "launch",
            "mainClass": "CometsZBBKStarMARKStarExample",
            "projectName": "OSPREY3",
            "cwd": "${workspaceFolder}/examples/java",
            "classPaths": [
                "${workspaceFolder}/build/libs/*",
                "${workspaceFolder}/build/classes/java/main",
                "${workspaceFolder}/examples/java"
            ],
            "vmArgs": "-Xmx4g"
        }
    ]
}
```

#### 步骤 2: 运行

1. 按 `F5` 或点击左侧的 "Run and Debug" 图标
2. 选择 "Run Simple COMETSZ Example"
3. 点击绿色播放按钮

## 📝 创建 VSCode Tasks（可选）

创建 `.vscode/tasks.json` 来快速编译：

```json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "Compile OSPREY",
            "type": "shell",
            "command": "./gradlew build",
            "group": "build",
            "problemMatcher": []
        },
        {
            "label": "Compile Simple Example",
            "type": "shell",
            "command": "javac -cp ../../build/libs/*:../../build/classes/java/main SimpleCometsZMARKStarExample.java",
            "options": {
                "cwd": "${workspaceFolder}/examples/java"
            },
            "group": "build",
            "problemMatcher": ["$javac"]
        },
        {
            "label": "Run Simple Example",
            "type": "shell",
            "command": "java -cp ../../build/libs/*:../../build/classes/java/main:. SimpleCometsZMARKStarExample",
            "options": {
                "cwd": "${workspaceFolder}/examples/java"
            },
            "group": "test",
            "dependsOn": ["Compile Simple Example"],
            "problemMatcher": []
        }
    ]
}
```

**使用：**
1. 按 `Ctrl+Shift+P`
2. 输入 "Tasks: Run Task"
3. 选择 "Run Simple Example"

## 🔍 调试（Debug）

### 设置断点

1. 在代码左侧点击，设置红色断点
2. 按 `F5` 启动调试
3. 使用调试工具栏：
   - Continue (F5)
   - Step Over (F10)
   - Step Into (F11)
   - Step Out (Shift+F11)

### 调试技巧

```java
// 在关键位置设置断点，查看变量值
for (CometsZ.State state : cometsZ.states) {
    // 在这里设置断点，查看每个状态的处理
    state.confEcalc = ...;
}

// 在结果输出前设置断点
List<CometsZ.SequenceInfo> sequences = cometsZ.findBestSequences(5);
// 断点在这里，检查 sequences 的内容
```

## 📂 推荐的 VSCode 工作区结构

```
OSPREY3/
├── .vscode/
│   ├── launch.json       # 运行配置
│   ├── tasks.json        # 任务配置
│   └── settings.json     # 项目设置
├── examples/
│   └── java/
│       ├── SimpleCometsZMARKStarExample.java
│       ├── CometsZBBKStarMARKStarExample.java
│       ├── run_simple_example.sh          # 一键运行脚本
│       └── *.md                            # 文档
└── src/
    └── ...
```

## 🎯 一键运行脚本（最快方法）

我为您创建完整的运行脚本，直接使用即可：

### 脚本 1: `run_simple_example.sh`

```bash
#!/bin/bash
cd "$(dirname "$0")"
PROJECT_ROOT="../.."

echo "🚀 运行 Simple COMETSZ Example..."

# 编译 OSPREY（如果需要）
if [ ! -d "$PROJECT_ROOT/build/libs" ]; then
    echo "📦 首次运行，正在编译 OSPREY..."
    cd "$PROJECT_ROOT" && ./gradlew build && cd -
fi

# 编译示例
echo "🔨 编译示例..."
javac -cp "$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main" \
    SimpleCometsZMARKStarExample.java

# 运行
echo "▶️  运行示例..."
java -Xmx4g \
    -cp "$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main:." \
    SimpleCometsZMARKStarExample

echo "✅ 完成！"
```

### 脚本 2: `run_complete_example.sh`

```bash
#!/bin/bash
cd "$(dirname "$0")"
PROJECT_ROOT="../.."

echo "🚀 运行 Complete COMETSZ Example..."

# 编译 OSPREY（如果需要）
if [ ! -d "$PROJECT_ROOT/build/libs" ]; then
    echo "📦 首次运行，正在编译 OSPREY..."
    cd "$PROJECT_ROOT" && ./gradlew build && cd -
fi

# 编译示例
echo "🔨 编译示例..."
javac -cp "$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main" \
    CometsZBBKStarMARKStarExample.java

# 运行
echo "▶️  运行示例..."
java -Xmx4g \
    -cp "$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main:." \
    CometsZBBKStarMARKStarExample

echo "✅ 完成！"
```

**使用方法：**

```bash
# 在 VSCode 终端中
cd examples/java
chmod +x run_simple_example.sh run_complete_example.sh

# 运行简化示例
./run_simple_example.sh

# 运行完整示例
./run_complete_example.sh
```

## 🛠️ VSCode 设置优化

创建 `.vscode/settings.json`：

```json
{
    "java.configuration.updateBuildConfiguration": "automatic",
    "java.compile.nullAnalysis.mode": "automatic",
    "java.debug.settings.hotCodeReplace": "auto",
    "java.saveActions.organizeImports": true,

    // 内存设置
    "java.jdt.ls.vmargs": "-Xmx2g",

    // 文件关联
    "files.associations": {
        "*.gradle": "groovy"
    },

    // 终端设置
    "terminal.integrated.defaultProfile.linux": "bash",

    // 排除不必要的文件
    "files.exclude": {
        "**/.gradle": true,
        "**/build": false,  // 保持可见，因为需要查看编译结果
        "**/*.class": true
    }
}
```

## 📊 VSCode vs IntelliJ IDEA

| 特性 | VSCode | IntelliJ IDEA |
|------|--------|---------------|
| 启动速度 | ⚡ 快 | 慢 |
| 内存占用 | 💾 低 | 高 |
| Java 支持 | ✅ 好 | ⭐ 优秀 |
| 调试功能 | ✅ 完整 | ⭐ 更强大 |
| 免费 | ✅ 完全免费 | Community 版免费 |
| 插件生态 | 🔌 丰富 | 🔌 丰富 |
| 适合 | 轻量级开发 | 专业 Java 开发 |

**结论**: 两者都很好，VSCode 更轻量，IntelliJ 对 Java 支持更好。

## ⚡ 快捷方式

### VSCode Java 快捷键

```
运行程序:        F5 (Debug) 或 Ctrl+F5 (Run)
停止程序:        Shift+F5
重启程序:        Ctrl+Shift+F5
设置断点:        F9
单步跳过:        F10
单步进入:        F11
单步跳出:        Shift+F11

查找文件:        Ctrl+P
全局搜索:        Ctrl+Shift+F
命令面板:        Ctrl+Shift+P
终端:            Ctrl+`
```

## 🐛 常见问题

### Q1: VSCode 不识别 Java 类

**解决方案：**
```bash
# 1. 重新编译
./gradlew clean build

# 2. 在 VSCode 中
Ctrl+Shift+P → "Java: Clean Java Language Server Workspace"

# 3. 重启 VSCode
```

### Q2: 找不到依赖库

**解决方案：**
```bash
# 确保 build 目录存在
ls -la build/libs/

# 如果没有，运行
./gradlew build
```

### Q3: 内存不足错误

**解决方案：**

在 launch.json 中增加内存：
```json
"vmArgs": "-Xmx8g"  // 从 4g 增加到 8g
```

### Q4: 无法找到主类

**解决方案：**

确保 classPaths 正确：
```json
"classPaths": [
    "${workspaceFolder}/build/libs/*",
    "${workspaceFolder}/build/classes/java/main",
    "${workspaceFolder}/examples/java"  // 包含示例目录
]
```

## 📦 完整的一键设置脚本

创建 `setup_vscode.sh`：

```bash
#!/bin/bash

echo "🔧 设置 VSCode 环境..."

# 创建 .vscode 目录
mkdir -p .vscode

# 创建 launch.json
cat > .vscode/launch.json << 'EOF'
{
    "version": "0.2.0",
    "configurations": [
        {
            "type": "java",
            "name": "Run Simple COMETSZ Example",
            "request": "launch",
            "mainClass": "SimpleCometsZMARKStarExample",
            "projectName": "OSPREY3",
            "cwd": "${workspaceFolder}/examples/java",
            "classPaths": [
                "${workspaceFolder}/build/libs/*",
                "${workspaceFolder}/build/classes/java/main",
                "${workspaceFolder}/examples/java"
            ],
            "vmArgs": "-Xmx4g"
        }
    ]
}
EOF

# 创建 settings.json
cat > .vscode/settings.json << 'EOF'
{
    "java.configuration.updateBuildConfiguration": "automatic",
    "java.jdt.ls.vmargs": "-Xmx2g"
}
EOF

echo "✅ VSCode 配置完成！"
echo ""
echo "下一步："
echo "1. 在 VSCode 中打开项目"
echo "2. 安装 'Extension Pack for Java'"
echo "3. 打开 examples/java/SimpleCometsZMARKStarExample.java"
echo "4. 点击 'Run' 或按 F5"
```

**运行：**
```bash
chmod +x setup_vscode.sh
./setup_vscode.sh
```

## 🎉 总结

在 VSCode 中运行示例的最简单方法：

### 方法 1: 图形界面（推荐新手）

1. 安装 "Extension Pack for Java"
2. 打开 `SimpleCometsZMARKStarExample.java`
3. 点击 `main` 方法上的 "Run"
4. 完成！

### 方法 2: 使用脚本（推荐经验用户）

```bash
cd examples/java
./run_simple_example.sh
```

### 方法 3: 手动命令行（最灵活）

```bash
cd examples/java
javac -cp ../../build/libs/*:../../build/classes/java/main SimpleCometsZMARKStarExample.java
java -cp ../../build/libs/*:../../build/classes/java/main:. SimpleCometsZMARKStarExample
```

**所有方法都完全可行！选择您最喜欢的即可。** 🚀
