#!/bin/bash
# 分析所有5个系统并生成汇总报告

echo "=========================================="
echo "分析所有5个系统"
echo "=========================================="
echo ""

# 分析每个系统
for i in {1..5}; do
    csv_file="markstar_order_vs_contribution_system${i}_7flex_eps0p010.csv"
    output_file="system${i}_7flex_analysis.txt"

    if [ -f "$csv_file" ]; then
        echo "分析 System $i..."
        python3 analyze_subtree_vs_error.py "$csv_file" > "$output_file"
        echo "✓ System $i 分析完成"
    else
        echo "✗ System $i: CSV文件不存在 ($csv_file)"
    fi
done

echo ""
echo "=========================================="
echo "生成汇总报告"
echo "=========================================="
echo ""

# 生成汇总报告
python3 << 'EOFPY'
import re
import csv

print("\n" + "="*100)
print(" "*30 + "5个系统SubtreeUpper vs ErrorBound对比汇总")
print("="*100 + "\n")

results = []

for i in range(1, 6):
    tag = f"system{i}_7flex"
    analysis_file = f"{tag}_analysis.txt"
    csv_file = f"markstar_order_vs_contribution_{tag}_eps0p010.csv"

    try:
        # 读取分析结果
        with open(analysis_file, 'r', encoding='utf-8') as f:
            content = f.read()

            # 提取Spearman系数
            subtree_match = re.search(r'SubtreeUpper vs Contribution的Spearman相关系数: ρ = ([-\d.]+)', content)
            error_match = re.search(r'ErrorBound vs Contribution的Spearman相关系数: ρ = ([-\d.]+)', content)
            actual_match = re.search(r'实际顺序 vs Contribution的Spearman相关系数: ρ = ([-\d.]+)', content)

            if subtree_match and error_match:
                subtree_rho = float(subtree_match.group(1))
                error_rho = float(error_match.group(1))
                actual_rho = float(actual_match.group(1)) if actual_match else 0

                # 读取CSV获取conformation数量
                with open(csv_file, 'r') as csvf:
                    reader = csv.reader(csvf)
                    next(reader)  # skip header
                    num_confs = sum(1 for _ in reader)

                # 提取Top-K数据
                top3_match = re.search(r'Top-3分析：.*?SubtreeUpper策略:\s+(\d+)/3', content, re.DOTALL)
                top5_match = re.search(r'Top-5分析：.*?SubtreeUpper策略:\s+(\d+)/5', content, re.DOTALL)

                top3_subtree = int(top3_match.group(1)) if top3_match else 0
                top5_subtree = int(top5_match.group(1)) if top5_match else 0

                results.append({
                    'system': f"System {i}",
                    'num_confs': num_confs,
                    'subtree_rho': subtree_rho,
                    'error_rho': error_rho,
                    'actual_rho': actual_rho,
                    'diff': subtree_rho - error_rho,
                    'top3': top3_subtree,
                    'top5': top5_subtree
                })

                print(f"✓ System {i}: {num_confs} conformations")
            else:
                print(f"✗ System {i}: 无法提取Spearman系数")

    except FileNotFoundError:
        print(f"✗ System {i}: 文件不存在")
    except Exception as e:
        print(f"✗ System {i}: 错误 - {e}")

if results:
    print("\n" + "="*100)
    print("Spearman相关系数对比")
    print("="*100)
    print(f"{'System':<12} {'N_confs':<10} {'SubtreeUpper ρ':<18} {'ErrorBound ρ':<18} {'MARK* ρ':<18} {'Diff (S-E)':<12}")
    print("-"*100)

    for r in results:
        print(f"{r['system']:<12} {r['num_confs']:<10} {r['subtree_rho']:<18.4f} {r['error_rho']:<18.4f} {r['actual_rho']:<18.4f} {r['diff']:+.4f}")

    # 计算平均值
    avg_subtree = sum(r['subtree_rho'] for r in results) / len(results)
    avg_error = sum(r['error_rho'] for r in results) / len(results)
    avg_actual = sum(r['actual_rho'] for r in results) / len(results)
    avg_diff = avg_subtree - avg_error

    print("-"*100)
    print(f"{'Average':<12} {'':<10} {avg_subtree:<18.4f} {avg_error:<18.4f} {avg_actual:<18.4f} {avg_diff:+.4f}")

    print("\n" + "="*100)
    print("Top-K匹配分析")
    print("="*100)
    print(f"{'System':<12} {'Top-3 (SubtreeUpper)':<25} {'Top-5 (SubtreeUpper)':<25}")
    print("-"*100)

    for r in results:
        print(f"{r['system']:<12} {r['top3']}/3 ({r['top3']*100/3:.0f}%){' '*15} {r['top5']}/5 ({r['top5']*100/5:.0f}%)")

    print("\n" + "="*100)
    print("结论")
    print("="*100)
    print()

    if abs(avg_diff) < 0.01:
        print(f"✓ 两种策略表现几乎相同 (平均差异: {avg_diff:+.4f}, {abs(avg_diff)*100/avg_error:.2f}%)")
        print()
        print("在5个不同蛋白区域的7-flex系统上，SubtreeUpper和ErrorBound作为prioritization")
        print("指标的表现没有显著差异。MARK*当前的ErrorBound策略已经非常有效。")
    elif avg_diff > 0.01:
        print(f"✓ SubtreeUpper策略平均优于ErrorBound {avg_diff:.4f} ({avg_diff*100/avg_error:.1f}%)")
        print()
        print("建议考虑修改MARKStarNode.compareTo()使用SubtreeUpper作为优先级指标。")
    else:
        print(f"✓ ErrorBound策略平均优于SubtreeUpper {-avg_diff:.4f} ({-avg_diff*100/avg_subtree:.1f}%)")
        print()
        print("MARK*当前的ErrorBound策略表现更好，应该保持。")

    print()
    print("="*100)
    print("详细分析结果保存在:")
    for i in range(1, 6):
        print(f"  - system{i}_7flex_analysis.txt")
    print("="*100)
    print()

else:
    print("\n没有成功分析任何系统！")

EOFPY

echo ""
echo "完成！"
