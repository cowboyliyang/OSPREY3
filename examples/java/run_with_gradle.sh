#!/bin/bash

# 使用 Gradle 运行 COMETSZ + BBKStar + MARKStar 示例
# 这是最可靠的方法，因为 Gradle 会处理所有依赖

cd "$(dirname "$0")"
PROJECT_ROOT="../.."

echo ""
echo "================================================================================"
echo "🚀 运行 COMETSZ + BBKStar + MARKStar 示例 (使用 Gradle)"
echo "================================================================================"
echo ""

echo "📋 准备环境..."
cd "$PROJECT_ROOT"

# 确保 OSPREY 已编译
if [ ! -d "build/libs" ] || [ ! "$(ls -A build/libs 2>/dev/null)" ]; then
    echo "📦 首次运行，正在编译 OSPREY..."
    if ! ./gradlew build --no-daemon; then
        echo "❌ OSPREY 编译失败"
        exit 1
    fi
    echo "✅ OSPREY 编译完成"
else
    echo "✅ OSPREY 已编译"
fi

echo ""
echo "🔨 编译示例..."
cd examples/java

# 使用 Gradle 的 classpath 编译示例
../../gradlew --quiet :compileJava || true  # 确保 OSPREY 已编译

# 获取完整的 runtime classpath
CLASSPATH=$(cd "$PROJECT_ROOT" && ./gradlew -q printClasspath --no-daemon 2>/dev/null || echo "")

# 如果上面的命令失败，使用备用方法
if [ -z "$CLASSPATH" ]; then
    echo "使用备用 classpath 方法..."
    CLASSPATH="$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main"
    # 添加 Gradle 缓存中的所有依赖
    if [ -d "$HOME/.gradle/caches/modules-2/files-2.1" ]; then
        for jar in $(find "$HOME/.gradle/caches/modules-2/files-2.1" -name "*.jar" 2>/dev/null | grep -E "(jcommander|guava|commons)" | head -20); do
            CLASSPATH="$CLASSPATH:$jar"
        done
    fi
fi

echo "✅ Classpath 已准备"

echo ""
echo "🔨 编译示例文件..."

# 检查是否可以直接使用 javac
if command -v javac &> /dev/null; then
    javac -cp "$CLASSPATH" CometsZBBKStarMARKStarExample.java
    COMPILE_SUCCESS=$?
else
    # 使用 Gradle wrapper 中的 Java
    echo "使用 Gradle 的 Java 编译..."
    cd "$PROJECT_ROOT"
    ./gradlew --no-daemon -q javaexec \
        --main-class=javax.tools.ToolProvider \
        --classpath="$CLASSPATH:examples/java" 2>&1
    COMPILE_SUCCESS=$?
    cd - > /dev/null
fi

if [ $COMPILE_SUCCESS -ne 0 ]; then
    echo "⚠️  标准编译失败，尝试使用 Gradle 任务..."
    cd "$PROJECT_ROOT"
    # 直接使用 Gradle 运行，跳过单独编译
else
    echo "✅ 示例编译成功"
fi

echo ""
echo "================================================================================"
echo "▶️  运行示例（这可能需要 10-60 分钟）"
echo "================================================================================"
echo ""

cd "$PROJECT_ROOT"

# 直接使用 Gradle 运行，这会处理所有依赖
./gradlew --no-daemon javaexec \
    -Pmain=CometsZBBKStarMARKStarExample \
    -PjvmArgs="-Xmx4g -Xms1g" \
    -PworkingDir="examples/java" \
    -PextraClasspath="examples/java"

RESULT=$?

echo ""
if [ $RESULT -eq 0 ]; then
    echo "================================================================================"
    echo "✅ 运行完成！"
    echo "================================================================================"
    echo ""
    echo "📄 结果文件："
    if [ -f "examples/java/cometsz.bbkstar.markstar.results.tsv" ]; then
        echo "  ✅ cometsz.bbkstar.markstar.results.tsv"
    else
        echo "  ⚠️  结果文件未生成"
    fi
else
    echo "================================================================================"
    echo "❌ 运行失败 (退出码: $RESULT)"
    echo "================================================================================"
fi

echo ""
echo "💡 提示："
echo "  - 查看结果: cat examples/java/cometsz.bbkstar.markstar.results.tsv"
echo "  - 查看文档: cat examples/java/COMETSZ_BBKSTAR_MARKSTAR_README.md"
echo ""
