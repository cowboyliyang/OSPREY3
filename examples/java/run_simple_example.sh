#!/bin/bash

# 简化的 COMETSZ + MARKStar 示例运行脚本
# 适用于 VSCode 或任何终端环境

cd "$(dirname "$0")"
PROJECT_ROOT="../.."

echo ""
echo "================================================================================"
echo "🚀 运行 Simple COMETSZ + MARKStar 示例"
echo "================================================================================"
echo ""

# 检查 Java 版本
echo "📋 检查环境..."
if ! command -v java &> /dev/null; then
    echo "❌ 错误: 未找到 Java。请安装 Java 11 或更高版本。"
    exit 1
fi

JAVA_VERSION=$(java -version 2>&1 | grep -oP 'version "?(1\.)?\K\d+' | head -1)
if [ "$JAVA_VERSION" -lt 11 ]; then
    echo "❌ 错误: Java 版本过低 (当前: $JAVA_VERSION, 需要: 11+)"
    exit 1
fi
echo "✅ Java 版本: $JAVA_VERSION"

# 编译 OSPREY（如果需要）
if [ ! -d "$PROJECT_ROOT/build/libs" ] || [ ! "$(ls -A $PROJECT_ROOT/build/libs)" ]; then
    echo ""
    echo "📦 首次运行，正在编译 OSPREY（这可能需要几分钟）..."
    cd "$PROJECT_ROOT"
    if ! ./gradlew build; then
        echo "❌ OSPREY 编译失败"
        exit 1
    fi
    cd - > /dev/null
    echo "✅ OSPREY 编译完成"
else
    echo "✅ OSPREY 已编译"
fi

# 检查是否有编译后的 jar 文件
CLASSPATH="$PROJECT_ROOT/build/libs/*:$PROJECT_ROOT/build/classes/java/main"
if [ ! -d "$PROJECT_ROOT/build/libs" ]; then
    echo "❌ 错误: 未找到编译后的库文件"
    exit 1
fi

# 编译示例
echo ""
echo "🔨 编译示例文件..."
if ! javac -cp "$CLASSPATH" SimpleCometsZMARKStarExample.java 2>&1; then
    echo "❌ 示例编译失败"
    exit 1
fi
echo "✅ 示例编译成功"

# 运行示例
echo ""
echo "================================================================================"
echo "▶️  运行示例（这可能需要 5-30 分钟，取决于系统大小）"
echo "================================================================================"
echo ""

# 设置 JVM 参数
JVM_OPTS="-Xmx4g -Xms1g"

# 运行
if ! java $JVM_OPTS -cp "$CLASSPATH:." SimpleCometsZMARKStarExample; then
    echo ""
    echo "❌ 运行失败"
    exit 1
fi

echo ""
echo "================================================================================"
echo "✅ 运行完成！"
echo "================================================================================"
echo ""
echo "📄 结果文件："
echo "  - cometsz.results.tsv (序列结果)"
if [ -f "cometsz.results.tsv" ]; then
    echo "    ✅ 文件已生成"
else
    echo "    ⚠️  文件未生成（可能是程序提前退出）"
fi
echo ""
echo "💡 提示："
echo "  - 查看结果: cat cometsz.results.tsv"
echo "  - 修改参数: 编辑 SimpleCometsZMARKStarExample.java"
echo "  - 查看文档: cat QUICK_START_GUIDE.md"
echo ""
