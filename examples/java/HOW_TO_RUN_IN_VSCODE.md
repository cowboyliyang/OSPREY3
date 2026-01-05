# 如何在 VSCode 中运行 CometsZBBKStarMARKStarExample

## 问题说明

由于这个例子需要完整的 OSPREY 运行时依赖（包括所有第三方库），直接使用 `javac` 和 `java` 命令会遇到 classpath 问题。

## ✅ 解决方案：使用现有的测试代码

最简单可靠的方法是**直接运行现有的测试代码**，它们已经在项目中配置好了所有依赖。

### 方法 1：运行现有的 MARKStar 测试 ⭐ **推荐**

在 VSCode 终端中运行：

```bash
cd /home/users/lz280/IdeaProjects/OSPREY3

# 运行 MARKStar 与 KStar 对比测试（大约 5-10 分钟）
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMARKStarVsKStarPartitionFunction"
```

这个测试会：
- ✅ 展示 MARKStar 如何工作
- ✅ 对比 MARKStar 与传统 K* 的性能
- ✅ 输出详细的分析结果

### 方法 2：运行 MSKStar 测试（类似 COMETSZ）

```bash
cd /home/users/lz280/IdeaProjects/OSPREY3

# 运行 MSKStar 测试（多状态设计）
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMSKStar.test2RL0PPI"
```

这个测试演示：
- ✅ 多状态设计（类似 COMETSZ）
- ✅ 蛋白质-配体结合优化
- ✅ LMFE 目标函数

### 方法 3：查看测试结果

```bash
# 测试报告位置
cat build/reports/tests/test/index.html

# 或者查看控制台输出
```

## 📖 理解示例代码

虽然直接运行示例有困难，但您可以：

1. **阅读示例代码** - 了解如何使用 API
   ```bash
   cat examples/java/CometsZBBKStarMARKStarExample.java
   ```

2. **阅读测试代码** - 查看实际运行的例子
   ```bash
   cat src/test/java/edu/duke/cs/osprey/kstar/TestMARKStarVsKStarPartitionFunction.java
   cat src/test/java/edu/duke/cs/osprey/kstar/TestCometsZWithBBKStarAndMARKStar.java
   ```

3. **修改测试代码** - 基于测试创建自己的版本
   - 测试代码已经配置好所有依赖
   - 可以复制测试方法并修改参数

## 🎯 实际使用建议

### 对于学习目的

**推荐路径：**
1. 运行测试查看输出
2. 阅读测试代码理解实现
3. 阅读示例代码理解API
4. 基于测试创建自己的版本

### 对于研究/生产使用

**推荐做法：**
1. 在测试目录创建自己的测试类
   ```java
   // 在 src/test/java/你的包名/MyCometsZTest.java
   public class MyCometsZTest {
       @Test
       public void myDesignTest() {
           // 你的代码
       }
   }
   ```

2. 运行自己的测试
   ```bash
   ./gradlew test --tests "你的包名.MyCometsZTest"
   ```

这样可以：
- ✅ 自动处理所有依赖
- ✅ 使用 JUnit 的测试框架
- ✅ 集成到现有构建系统

## 🔧 如果真的想运行独立示例

如果您坚持要运行独立的示例文件，需要：

### 步骤 1：获取完整的 classpath

```bash
cd /home/users/lz280/IdeaProjects/OSPREY3

# 方法 1：使用 Gradle 获取
./gradlew dependencies --configuration runtimeClasspath > deps.txt

# 方法 2：列出所有 jar
find ~/.gradle/caches -name "*.jar" > all_jars.txt
```

### 步骤 2：手动构建 classpath

```bash
# 这会非常长...
CLASSPATH="build/libs/*:build/classes/java/main"
CLASSPATH="$CLASSPATH:~/.gradle/caches/modules-2/files-2.1/com/beust/jcommander/1.72/..."
# ... 需要添加几十个依赖
```

### 步骤 3：编译和运行

```bash
javac -cp "$CLASSPATH" examples/java/CometsZBBKStarMARKStarExample.java
java -Xmx4g -cp "$CLASSPATH:examples/java" CometsZBBKStarMARKStarExample
```

**但这非常繁琐且容易出错！**

## 💡 最佳实践

### 推荐的工作流程

```bash
# 1. 创建测试文件
cat > src/test/java/edu/duke/cs/osprey/examples/MyExample.java << 'EOF'
package edu.duke.cs.osprey.examples;

import org.junit.jupiter.api.Test;
import edu.duke.cs.osprey.kstar.*;
// ... 其他 import

public class MyExample {
    @Test
    public void runMyDesign() {
        // 复制 CometsZBBKStarMARKStarExample 的代码到这里
        // 或者复制测试代码并修改
    }
}
EOF

# 2. 运行测试
./gradlew test --tests "edu.duke.cs.osprey.examples.MyExample"

# 3. 查看结果
cat build/reports/tests/test/index.html
```

## 📚 相关资源

### 可运行的测试

- `TestMARKStar` - MARKStar 基础测试
- `TestMARKStarVsKStarPartitionFunction` - 性能对比
- `TestCometsZWithBBKStarAndMARKStar` - COMETSZ + MARKStar
- `TestMSKStar` - 多状态设计

### 文档

- `COMETSZ_BBKSTAR_MARKSTAR_README.md` - 完整理论文档
- `QUICK_START_GUIDE.md` - 快速入门
- `README_EXAMPLES.md` - 示例总览

## 🎓 总结

**最简单的方法：**
```bash
# 直接运行现有测试
cd /home/users/lz280/IdeaProjects/OSPREY3
./gradlew test --tests "*MARKStar*"
```

**最灵活的方法：**
1. 在 `src/test/java` 创建自己的测试类
2. 复制并修改示例代码
3. 使用 `./gradlew test` 运行

**不推荐的方法：**
- ❌ 尝试手动配置 classpath
- ❌ 直接运行 `examples/java/*.java`（依赖问题）

记住：**测试就是可运行的示例！** OSPREY 的设计理念就是通过测试来展示功能。
