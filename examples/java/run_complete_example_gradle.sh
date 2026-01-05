#!/bin/bash

# 完整的 COMETSZ + BBKStar + MARKStar 示例运行脚本
# 使用 Gradle 的内置 Java，不依赖系统 Java

cd "$(dirname "$0")"
PROJECT_ROOT="../.."

echo ""
echo "================================================================================"
echo "🚀 运行 Complete COMETSZ + BBKStar + MARKStar 示例 (Gradle 版本)"
echo "================================================================================"
echo ""

# 获取 Gradle 使用的 Java 路径
echo "📋 检查环境..."
cd "$PROJECT_ROOT"

# 通过 Gradle 获取 Java 路径
JAVA_HOME_FROM_GRADLE=$(./gradlew -q javaToolchains 2>/dev/null | grep -A 1 "CURRENT" | grep "Location" | awk '{print $2}')

if [ -z "$JAVA_HOME_FROM_GRADLE" ]; then
    # 备选方案：使用 Gradle 自己的 Java
    GRADLE_JAVA=$(./gradlew properties 2>/dev/null | grep "java.home" | cut -d: -f2 | tr -d ' ')
    if [ -n "$GRADLE_JAVA" ]; then
        JAVA_HOME_FROM_GRADLE="$GRADLE_JAVA"
    fi
fi

# 设置 Java 路径
if [ -n "$JAVA_HOME_FROM_GRADLE" ]; then
    export JAVA_HOME="$JAVA_HOME_FROM_GRADLE"
    export PATH="$JAVA_HOME/bin:$PATH"
    echo "✅ 使用 Gradle 的 Java: $JAVA_HOME"
else
    # 最后的备选方案：直接从 Gradle wrapper 中提取
    echo "⚠️  无法自动检测 Java，尝试使用 Gradle wrapper..."
fi

cd - > /dev/null

# 检查 OSPREY 是否已编译
if [ ! -d "$PROJECT_ROOT/build/libs" ] || [ ! "$(ls -A $PROJECT_ROOT/build/libs)" ]; then
    echo ""
    echo "📦 首次运行，正在编译 OSPREY（这可能需要几分钟）..."
    cd "$PROJECT_ROOT"
    if ! ./gradlew build --no-daemon; then
        echo "❌ OSPREY 编译失败"
        exit 1
    fi
    cd - > /dev/null
    echo "✅ OSPREY 编译完成"
else
    echo "✅ OSPREY 已编译"
fi

# 设置 classpath
CLASSPATH="$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main"

# 编译示例 - 使用 Gradle 的 Java
echo ""
echo "🔨 编译示例文件..."
cd "$PROJECT_ROOT/examples/java"

if [ -n "$JAVA_HOME" ] && [ -x "$JAVA_HOME/bin/javac" ]; then
    # 使用检测到的 Java
    if ! "$JAVA_HOME/bin/javac" -cp "$CLASSPATH" CometsZBBKStarMARKStarExample.java 2>&1; then
        echo "❌ 示例编译失败"
        exit 1
    fi
else
    # 使用 Gradle 来编译
    echo "使用 Gradle 编译..."
    cd "$PROJECT_ROOT"
    cat > build/tmp/CompileExample.java << 'JAVAEOF'
import javax.tools.*;
import java.io.File;
import java.util.*;

public class CompileExample {
    public static void main(String[] args) {
        JavaCompiler compiler = ToolProvider.getSystemJavaCompiler();
        StandardJavaFileManager fileManager = compiler.getStandardFileManager(null, null, null);

        List<String> options = Arrays.asList(
            "-cp", args[0],
            "-d", "."
        );

        Iterable<? extends JavaFileObject> compilationUnits =
            fileManager.getJavaFileObjectsFromStrings(Arrays.asList(args[1]));

        JavaCompiler.CompilationTask task = compiler.getTask(
            null, fileManager, null, options, null, compilationUnits
        );

        System.exit(task.call() ? 0 : 1);
    }
}
JAVAEOF

    ./gradlew -q javaexec \
        --main=CompileExample \
        --args="$CLASSPATH $PROJECT_ROOT/examples/java/CometsZBBKStarMARKStarExample.java" \
        --no-daemon 2>&1

    cd - > /dev/null
fi

echo "✅ 示例编译成功"

# 运行示例
echo ""
echo "================================================================================"
echo "▶️  运行示例（这可能需要 10-60 分钟，取决于系统大小）"
echo "================================================================================"
echo ""

cd "$PROJECT_ROOT/examples/java"

if [ -n "$JAVA_HOME" ] && [ -x "$JAVA_HOME/bin/java" ]; then
    # 使用检测到的 Java
    JVM_OPTS="-Xmx4g -Xms1g"
    if ! "$JAVA_HOME/bin/java" $JVM_OPTS -cp "$CLASSPATH:." CometsZBBKStarMARKStarExample; then
        echo ""
        echo "❌ 运行失败"
        exit 1
    fi
else
    # 使用 Gradle 来运行
    echo "使用 Gradle 运行..."
    cd "$PROJECT_ROOT"
    ./gradlew -q javaexec \
        --main=CometsZBBKStarMARKStarExample \
        --classpath="$CLASSPATH:examples/java" \
        --jvm-args="-Xmx4g -Xms1g" \
        --no-daemon 2>&1
    cd - > /dev/null
fi

echo ""
echo "================================================================================"
echo "✅ 运行完成！"
echo "================================================================================"
echo ""
echo "📄 结果文件："
echo "  - cometsz.bbkstar.markstar.results.tsv (序列结果)"
if [ -f "cometsz.bbkstar.markstar.results.tsv" ]; then
    echo "    ✅ 文件已生成"
else
    echo "    ⚠️  文件未生成（可能是程序提前退出）"
fi
echo ""
echo "💡 提示："
echo "  - 查看结果: cat cometsz.bbkstar.markstar.results.tsv"
echo "  - 修改参数: 编辑 CometsZBBKStarMARKStarExample.java"
echo "  - 查看文档: cat COMETSZ_BBKSTAR_MARKSTAR_README.md"
echo ""
