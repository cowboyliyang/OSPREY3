#!/bin/bash
# 检查5个系统测试的进度

echo "=========================================="
echo "测试进度检查"
echo "=========================================="
echo ""

# 检查哪些系统已经完成
for i in {1..5}; do
    tag="system${i}_7flex"
    csv_file="markstar_order_vs_contribution_${tag}_eps0p010.csv"
    analysis_file="${tag}_analysis.txt"

    if [ -f "$csv_file" ] && [ -f "$analysis_file" ]; then
        num_confs=$(tail -n +2 "$csv_file" | wc -l)
        echo "✓ System $i: 完成 (${num_confs} conformations)"

        # 提取Spearman系数
        subtree_rho=$(grep "SubtreeUpper vs Contribution的Spearman相关系数" "$analysis_file" | grep -oP 'ρ = \K[-\d.]+')
        error_rho=$(grep "ErrorBound vs Contribution的Spearman相关系数" "$analysis_file" | grep -oP 'ρ = \K[-\d.]+')

        if [ -n "$subtree_rho" ] && [ -n "$error_rho" ]; then
            diff=$(python3 -c "print(f'{float('$subtree_rho') - float('$error_rho'):+.4f}')")
            echo "    SubtreeUpper: ρ = $subtree_rho"
            echo "    ErrorBound:   ρ = $error_rho"
            echo "    差异:         $diff"
        fi
    elif [ -f "${tag}_test.log" ]; then
        # 检查测试日志
        if grep -q "BUILD SUCCESSFUL" "${tag}_test.log"; then
            echo "⏳ System $i: 测试完成，正在分析..."
        elif grep -q "BUILD FAILED" "${tag}_test.log"; then
            echo "✗ System $i: 测试失败"
        else
            echo "⏳ System $i: 正在运行..."
            # 尝试显示当前进度
            last_line=$(tail -1 "${tag}_test.log" 2>/dev/null)
            if [ -n "$last_line" ]; then
                echo "    最新: ${last_line:0:80}"
            fi
        fi
    else
        echo "⏳ System $i: 等待中..."
    fi
    echo ""
done

echo "=========================================="
echo "后台进程状态："
echo "=========================================="
ps aux | grep -E "gradlew|run_all_5_systems" | grep -v grep | awk '{print $2, $11, $12, $13, $14}'
