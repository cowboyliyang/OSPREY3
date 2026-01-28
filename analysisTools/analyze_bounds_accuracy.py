#!/usr/bin/env python3
"""
分析minimized energy与初始upper/lower bounds的接近程度
"""

import csv
import sys

def analyze_bounds_accuracy(filename):
    print("="*80)
    print("分析: Minimized Energy vs 初始Upper/Lower Bounds")
    print("="*80)
    print()

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        data = list(reader)

    print(f"总共 {len(data)} 个minimized conformations\n")

    print("="*80)
    print("详细分析:")
    print("-"*80)
    print(f"{'Order':<6} | {'Lower':<10} | {'Upper':<10} | {'Final':<10} | {'→Lower':<8} | {'→Upper':<8} | 更接近")
    print("-"*80)

    closer_to_lower = 0
    closer_to_upper = 0

    total_lower_error = 0
    total_upper_error = 0

    details = []

    for row in data:
        order = int(row['MinimizeOrder'])
        lower = float(row['LowerBound'])
        upper = float(row['Score'])  # Score是upper bound (rigid energy)
        final = float(row['FinalEnergy'])

        # 计算到upper和lower的距离
        dist_to_lower = abs(final - lower)
        dist_to_upper = abs(final - upper)

        # 计算相对误差
        energy_range = abs(upper - lower)
        if energy_range > 0:
            rel_error_lower = dist_to_lower / energy_range
            rel_error_upper = dist_to_upper / energy_range
        else:
            rel_error_lower = 0
            rel_error_upper = 0

        total_lower_error += dist_to_lower
        total_upper_error += dist_to_upper

        if dist_to_lower < dist_to_upper:
            closer = "Lower ✓"
            closer_to_lower += 1
        else:
            closer = "Upper ✓"
            closer_to_upper += 1

        details.append({
            'order': order,
            'lower': lower,
            'upper': upper,
            'final': final,
            'dist_lower': dist_to_lower,
            'dist_upper': dist_to_upper,
            'closer': closer,
            'range': energy_range
        })

        print(f"{order:<6} | {lower:10.4f} | {upper:10.4f} | {final:10.4f} | "
              f"{dist_to_lower:8.4f} | {dist_to_upper:8.4f} | {closer}")

    print()
    print("="*80)
    print("统计汇总:")
    print("-"*80)
    print(f"更接近Lower Bound的: {closer_to_lower}/{len(data)} ({closer_to_lower*100/len(data):.1f}%)")
    print(f"更接近Upper Bound的: {closer_to_upper}/{len(data)} ({closer_to_upper*100/len(data):.1f}%)")
    print()
    print(f"平均到Lower Bound的距离: {total_lower_error/len(data):.4f} kcal/mol")
    print(f"平均到Upper Bound的距离: {total_upper_error/len(data):.4f} kcal/mol")

    # 分析bounds的质量
    print()
    print("="*80)
    print("Bounds质量分析:")
    print("-"*80)

    avg_range = sum(d['range'] for d in details) / len(details)
    print(f"平均bounds区间宽度: {avg_range:.4f} kcal/mol")

    # 计算final energy在bounds区间中的相对位置
    print()
    print("Final Energy在[Lower, Upper]区间中的位置:")
    print("-"*80)
    print(f"{'Order':<6} | {'Position':<10} | 解释")
    print("-"*80)

    positions = []
    for d in details:
        if d['range'] > 0:
            # position = 0 means at lower bound, 1 means at upper bound
            position = (d['final'] - d['lower']) / d['range']
            positions.append(position)

            if position < 0:
                interpretation = "低于Lower! (异常)"
            elif position > 1:
                interpretation = "高于Upper! (异常)"
            elif position < 0.33:
                interpretation = "接近Lower"
            elif position < 0.67:
                interpretation = "居中"
            else:
                interpretation = "接近Upper"

            print(f"{d['order']:<6} | {position:10.2%} | {interpretation}")
        else:
            print(f"{d['order']:<6} | {'N/A':<10} | Lower=Upper (已精确)")

    print()
    print("="*80)
    print("结论:")
    print("="*80)

    if positions:
        avg_position = sum(positions) / len(positions)
        print(f"\n平均位置: {avg_position:.2%}")
        print("  (0% = Lower Bound, 50% = 中点, 100% = Upper Bound)")
        print()

        if avg_position < 0.33:
            print("✓✓✓ Minimized energy显著更接近Lower Bound!")
            print("    说明: Lower bound (minimizing emat)是很好的估计")
            print("    说明: MARK*的lower bound估计质量很高")
        elif avg_position < 0.5:
            print("✓✓ Minimized energy更接近Lower Bound")
            print("   说明: Lower bound估计较准确")
        elif avg_position < 0.67:
            print("○ Minimized energy居中")
            print("  说明: Bounds估计一般")
        else:
            print("✗ Minimized energy更接近Upper Bound")
            print("  说明: Lower bound可能过于保守")

        # 检查异常值
        anomalies = [p for p in positions if p < 0 or p > 1]
        if anomalies:
            print(f"\n⚠ 警告: 发现 {len(anomalies)} 个异常值 (final energy超出bounds)")

    print("="*80)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "markstar_order_vs_contribution_4flex_eps0p01.csv"

    analyze_bounds_accuracy(filename)
