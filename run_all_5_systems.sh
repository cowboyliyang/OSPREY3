#!/bin/bash
# 运行5个不同系统的7-flex测试，并分析结果

echo "=========================================="
echo "运行5个不同系统的7-flex测试"
echo "=========================================="
echo ""

# 定义5个测试
tests=(
    "testSystem1_Residues21to27"
    "testSystem2_Residues30to36"
    "testSystem3_Residues40to46"
    "testSystem4_Residues50to56"
    "testSystem5_Residues60to66"
)

tags=(
    "system1_7flex"
    "system2_7flex"
    "system3_7flex"
    "system4_7flex"
    "system5_7flex"
)

# 运行每个测试
for i in {0..4}; do
    test_name="${tests[$i]}"
    tag="${tags[$i]}"

    echo ""
    echo "=========================================="
    echo "运行测试 $((i+1))/5: $test_name"
    echo "=========================================="
    echo ""

    # 运行测试
    ./gradlew test --tests "edu.duke.cs.osprey.markstar.TestMultipleSystems7Flex.$test_name" 2>&1 | tee "${tag}_test.log"

    if [ $? -eq 0 ]; then
        echo "✓ 测试 $test_name 完成"

        # 提取数据并生成CSV
        echo "提取数据到CSV..."
        python3 << EOF
import xml.etree.ElementTree as ET
import re
import csv
import math

# 读取XML
try:
    tree = ET.parse('build/test-results/test/TEST-edu.duke.cs.osprey.markstar.TestMultipleSystems7Flex.xml')
    root = tree.getroot()
    system_out = root.find('.//system-out').text

    # 提取所有conformation行
    pattern = re.compile(
        r'\[([^\]]+)\]conf:\s*(\d+),\s*score:\s*([-\d.]+),\s*lower:\s*([-\d.]+),\s*corrected:\s*([-\d.]+)\s*energy:\s*([-\d.]+),\s*confBounds:\[([-\d.\s]+),([-\d.\s]+)\],\s*subtreeBounds:\[([^\]]+)\],\s*errorBound:([^,]+)'
    )

    conformations = []
    for match in pattern.finditer(system_out):
        conf_str = match.group(1).strip()
        order = int(match.group(2))
        score = float(match.group(3))
        lower_bound = float(match.group(4))
        corrected = float(match.group(5))
        energy = float(match.group(6))
        conf_lower = float(match.group(7))
        conf_upper = float(match.group(8))
        subtree_bounds = match.group(9)
        error_bound_str = match.group(10).strip()

        # Parse subtree bounds
        subtree_parts = subtree_bounds.split(',')
        subtree_lower = float(subtree_parts[0].strip())
        subtree_upper = float(subtree_parts[1].strip())
        error_bound = float(error_bound_str)

        conformations.append({
            'order': order,
            'conf': conf_str,
            'score': score,
            'lower_bound': lower_bound,
            'corrected': corrected,
            'energy': energy,
            'conf_lower': conf_lower,
            'conf_upper': conf_upper,
            'subtree_lower': subtree_lower,
            'subtree_upper': subtree_upper,
            'error_bound': error_bound
        })

    print(f"Found {len(conformations)} conformations for ${tag}")

    if len(conformations) == 0:
        print("Warning: No conformations found!")
        exit(1)

    # 计算Boltzmann weights和贡献
    RT = 0.5962  # kcal/mol at 298K
    for conf in conformations:
        conf['boltzmann_weight'] = math.exp(-conf['energy'] / RT)

    total_weight = sum(c['boltzmann_weight'] for c in conformations)

    for conf in conformations:
        conf['contribution'] = conf['boltzmann_weight'] / total_weight * 100

    # 按贡献排序
    sorted_by_contrib = sorted(conformations, key=lambda x: x['contribution'], reverse=True)
    for rank, conf in enumerate(sorted_by_contrib, 1):
        conf['rank'] = rank

    # 写CSV
    csv_file = 'markstar_order_vs_contribution_${tag}_eps0p010.csv'
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'MinimizeOrder', 'Conformation', 'Score', 'LowerBound', 'Corrected', 'FinalEnergy',
            'OldConfLower', 'OldConfUpper', 'SubtreeLower', 'SubtreeUpper', 'ErrorBound',
            'BoltzmannWeight', 'PercentContribution', 'ContributionRank'
        ])

        for conf in sorted(conformations, key=lambda x: x['order']):
            writer.writerow([
                conf['order'],
                conf['conf'],
                f"{conf['score']:.6f}",
                f"{conf['lower_bound']:.6f}",
                f"{conf['corrected']:.6f}",
                f"{conf['energy']:.6f}",
                f"{conf['conf_lower']:.6f}",
                f"{conf['conf_upper']:.6f}",
                f"{conf['subtree_lower']:.6e}",
                f"{conf['subtree_upper']:.6e}",
                f"{conf['error_bound']:.6e}",
                f"{conf['boltzmann_weight']:.6e}",
                f"{conf['contribution']:.6f}",
                conf['rank']
            ])

    print(f"CSV written to {csv_file}")

    # 显示top 5
    print(f"\\nTop 5 by contribution:")
    for i, conf in enumerate(sorted_by_contrib[:5], 1):
        print(f"  #{i}: Order {conf['order']}, Energy {conf['energy']:.2f}, Contribution {conf['contribution']:.2f}%")

except Exception as e:
    print(f"Error processing test results: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
EOF

        if [ $? -eq 0 ]; then
            echo "✓ CSV生成成功"

            # 运行分析
            echo ""
            echo "运行分析..."
            python3 analyze_subtree_vs_error.py "markstar_order_vs_contribution_${tag}_eps0p010.csv" > "${tag}_analysis.txt"
            echo "✓ 分析完成，结果保存到 ${tag}_analysis.txt"
        else
            echo "✗ CSV生成失败"
        fi
    else
        echo "✗ 测试 $test_name 失败"
    fi

    echo ""
done

echo ""
echo "=========================================="
echo "所有测试完成！"
echo "=========================================="
echo ""
echo "生成的文件："
ls -lh system*_7flex*.csv system*_7flex*.txt 2>/dev/null

echo ""
echo "现在生成汇总报告..."

# 生成汇总报告
python3 << 'EOFPY'
import re
import csv

print("\n" + "="*80)
print("5个系统的SubtreeUpper vs ErrorBound对比汇总")
print("="*80 + "\n")

results = []

for i in range(1, 6):
    tag = f"system{i}_7flex"
    analysis_file = f"{tag}_analysis.txt"
    csv_file = f"markstar_order_vs_contribution_{tag}_eps0p010.csv"

    try:
        # 读取分析结果
        with open(analysis_file, 'r') as f:
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

                results.append({
                    'system': f"System {i}",
                    'num_confs': num_confs,
                    'subtree_rho': subtree_rho,
                    'error_rho': error_rho,
                    'actual_rho': actual_rho,
                    'diff': subtree_rho - error_rho
                })

                print(f"✓ {tag}: 找到 {num_confs} conformations")
            else:
                print(f"✗ {tag}: 无法提取Spearman系数")

    except FileNotFoundError:
        print(f"✗ {tag}: 文件不存在")
    except Exception as e:
        print(f"✗ {tag}: 错误 - {e}")

if results:
    print("\n" + "="*80)
    print("汇总结果")
    print("="*80)
    print(f"{'System':<15} {'N_confs':<10} {'SubtreeUpper ρ':<18} {'ErrorBound ρ':<18} {'Diff':<10}")
    print("-"*80)

    for r in results:
        print(f"{r['system']:<15} {r['num_confs']:<10} {r['subtree_rho']:<18.4f} {r['error_rho']:<18.4f} {r['diff']:+.4f}")

    # 计算平均值
    avg_subtree = sum(r['subtree_rho'] for r in results) / len(results)
    avg_error = sum(r['error_rho'] for r in results) / len(results)
    avg_diff = avg_subtree - avg_error

    print("-"*80)
    print(f"{'Average':<15} {'':<10} {avg_subtree:<18.4f} {avg_error:<18.4f} {avg_diff:+.4f}")

    print("\n" + "="*80)
    print("结论")
    print("="*80)

    if abs(avg_diff) < 0.01:
        print(f"两种策略表现几乎相同 (平均差异: {avg_diff:+.4f})")
    elif avg_diff > 0.01:
        print(f"SubtreeUpper策略平均优于ErrorBound {avg_diff:.4f} ({avg_diff*100/avg_error:.1f}%)")
    else:
        print(f"ErrorBound策略平均优于SubtreeUpper {-avg_diff:.4f} ({-avg_diff*100/avg_subtree:.1f}%)")

    print()

else:
    print("\n没有成功分析任何系统！")

EOFPY

echo ""
echo "完成！"
