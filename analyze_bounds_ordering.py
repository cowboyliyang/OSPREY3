#!/usr/bin/env python3
"""
分析minimized energy与初始bounds的排序关系

对于每个minimized conformation，我们有：
- Score: 初始的rigid energy matrix给出的分数 (上界)
- LowerBound: 初始的lower bound
- FinalEnergy: 最终minimize后的实际能量
- ErrorBound: 可以从upper和lower计算得出

问题：这些值的大小排序是什么样的？
"""

import csv
import sys

def analyze_bounds_ordering(csv_file):
    conformations = []

    # 读取CSV
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            conf = {
                'order': int(row['MinimizeOrder']),
                'conf': row['Conformation'],
                'score': float(row['Score']),
                'lower': float(row['LowerBound']),
                'corrected': float(row['Corrected']),
                'final': float(row['FinalEnergy']),
                'contribution': float(row['PercentContribution'])
            }
            conformations.append(conf)

    print("=" * 80)
    print("分析: Minimized Energy与初始Bounds的排序关系")
    print("=" * 80)
    print()

    # 统计不同的排序模式
    patterns = {
        'lower < final < score': 0,      # 理想情况：final在bounds之间
        'final < lower < score': 0,      # final低于lower bound (说明lower bound不准)
        'lower < score < final': 0,      # final高于score (说明score不准)
        'final < score < lower': 0,      # 顺序完全反了
        'score < lower < final': 0,      # score是lower，lower是upper？
        'score < final < lower': 0,      # score是lower，lower是upper，final在中间
        'other': 0
    }

    # 统计error bound (应该是 upper - lower)
    error_bounds = []

    print("详细数据 (前20个):")
    print("-" * 120)
    print(f"{'Order':<6} {'Score':<12} {'Lower':<12} {'Final':<12} {'排序模式':<30} {'贡献%':<10}")
    print("-" * 120)

    for i, conf in enumerate(conformations[:20]):
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        # 判断排序模式
        if lower < final < score:
            pattern = 'lower < final < score'
        elif final < lower < score:
            pattern = 'final < lower < score'
        elif lower < score < final:
            pattern = 'lower < score < final'
        elif final < score < lower:
            pattern = 'final < score < lower'
        elif score < lower < final:
            pattern = 'score < lower < final'
        elif score < final < lower:
            pattern = 'score < final < lower'
        else:
            pattern = 'other'

        patterns[pattern] += 1

        # 计算error bound (假设是 max - min)
        error = max(score, lower) - min(score, lower)
        error_bounds.append(error)

        print(f"{conf['order']:<6} {score:<12.4f} {lower:<12.4f} {final:<12.4f} {pattern:<30} {conf['contribution']:<10.2f}")

    print()
    print("=" * 80)
    print("统计结果 (所有conformations):")
    print("=" * 80)

    # 重新统计所有conformations
    patterns = {
        'lower < final < score': 0,
        'final < lower < score': 0,
        'lower < score < final': 0,
        'final < score < lower': 0,
        'score < lower < final': 0,
        'score < final < lower': 0,
        'other': 0
    }

    # 统计final相对于score和lower的位置
    final_between_bounds = 0
    final_below_both = 0
    final_above_both = 0

    # 统计score和lower谁更大
    score_greater_than_lower = 0
    lower_greater_than_score = 0

    for conf in conformations:
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        # 判断排序模式
        if lower < final < score:
            pattern = 'lower < final < score'
            final_between_bounds += 1
        elif final < lower < score:
            pattern = 'final < lower < score'
            final_below_both += 1
        elif lower < score < final:
            pattern = 'lower < score < final'
            final_above_both += 1
        elif final < score < lower:
            pattern = 'final < score < lower'
            final_between_bounds += 1
        elif score < lower < final:
            pattern = 'score < lower < final'
            final_above_both += 1
        elif score < final < lower:
            pattern = 'score < final < lower'
            final_between_bounds += 1
        else:
            pattern = 'other'

        patterns[pattern] += 1

        # 统计score vs lower
        if score > lower:
            score_greater_than_lower += 1
        else:
            lower_greater_than_score += 1

    total = len(conformations)

    print()
    print("排序模式分布:")
    for pattern, count in sorted(patterns.items(), key=lambda x: -x[1]):
        if count > 0:
            pct = count * 100.0 / total
            print(f"  {pattern:<30}: {count:3d} ({pct:5.1f}%)")

    print()
    print("Score vs LowerBound关系:")
    print(f"  Score > LowerBound: {score_greater_than_lower} ({score_greater_than_lower*100.0/total:.1f}%)")
    print(f"  Score < LowerBound: {lower_greater_than_score} ({lower_greater_than_score*100.0/total:.1f}%)")

    print()
    print("Final Energy位置:")
    print(f"  Final在两个bounds之间: {final_between_bounds} ({final_between_bounds*100.0/total:.1f}%)")
    print(f"  Final低于两个bounds:   {final_below_both} ({final_below_both*100.0/total:.1f}%)")
    print(f"  Final高于两个bounds:   {final_above_both} ({final_above_both*100.0/total:.1f}%)")

    print()
    print("=" * 80)
    print("距离分析:")
    print("=" * 80)

    # 计算平均距离
    dist_final_to_score = []
    dist_final_to_lower = []

    for conf in conformations:
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        dist_final_to_score.append(abs(final - score))
        dist_final_to_lower.append(abs(final - lower))

    avg_dist_score = sum(dist_final_to_score) / len(dist_final_to_score)
    avg_dist_lower = sum(dist_final_to_lower) / len(dist_final_to_lower)

    print(f"Final Energy到Score的平均距离:      {avg_dist_score:.4f} kcal/mol")
    print(f"Final Energy到LowerBound的平均距离: {avg_dist_lower:.4f} kcal/mol")

    if avg_dist_score < avg_dist_lower:
        print(f"\n✓ Final Energy平均而言更接近Score (比LowerBound近{avg_dist_lower/avg_dist_score:.2f}x)")
    else:
        print(f"\n✓ Final Energy平均而言更接近LowerBound (比Score近{avg_dist_score/avg_dist_lower:.2f}x)")

    print()
    print("=" * 80)
    print("结论:")
    print("=" * 80)

    # 理解字段含义
    if score_greater_than_lower > total * 0.9:
        print("\n根据数据分析：")
        print("  • Score字段 = 初始的UPPER bound (rigid energy)")
        print("  • LowerBound字段 = 初始的LOWER bound")
        print("  • Final Energy = minimize后的真实能量")
        print()
        print("能量排序关系：")
        if final_between_bounds > total * 0.5:
            print("  ✓ 大多数情况下: LowerBound ≤ FinalEnergy ≤ Score")
            print("    这是预期的理想情况：真实能量在bounds之间")
        elif final_below_both > total * 0.5:
            print("  ⚠ 大多数情况下: FinalEnergy < LowerBound < Score")
            print("    说明LowerBound估计不够紧（太高了）")
        elif final_above_both > total * 0.5:
            print("  ⚠ 大多数情况下: LowerBound < Score < FinalEnergy")
            print("    说明Score (upper bound)估计不够紧（太低了）")
    elif lower_greater_than_score > total * 0.9:
        print("\n⚠ 可能的字段名称问题：")
        print("  • 'Score'字段实际是LOWER bound")
        print("  • 'LowerBound'字段实际是UPPER bound")
        print("  需要检查MARK*输出的实际含义")

    print()

if __name__ == '__main__':
    csv_file = '/home/users/lz280/IdeaProjects/OSPREY3/markstar_order_vs_contribution_4flex_eps0p01.csv'
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    analyze_bounds_ordering(csv_file)
