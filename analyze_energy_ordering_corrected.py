#!/usr/bin/env python3
"""
分析minimized energy与初始bounds的排序关系（考虑负数能量）

关键：在分子模拟中，能量通常是负数
- 更负 = 更低能量 = 更稳定
- 例如：-76 < -75 (能量关系), 但 -76更低/更好

所以在比较"哪个是上界/下界"时：
- Lower Bound = 能量的下界 = 更负的值
- Upper Bound = 能量的上界 = 更不负(更高)的值
"""

import csv
import sys

def analyze_energy_ordering(csv_file):
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

    print("=" * 100)
    print("分析: Minimized Energy与决定minimize时的Bounds的排序关系")
    print("=" * 100)
    print()

    print("⚠ 重要：能量是负数，更负 = 更低能量 = 更稳定")
    print()
    print("字段说明（根据代码分析）：")
    print("  • Score        = oldgscore = 决定minimize时node.gscore的值")
    print("  • LowerBound   = minimizingEmat.confE() = minimizing energy matrix的完整conformation能量")
    print("  • FinalEnergy  = 实际minimize后测得的真实能量")
    print("  • confUpperBound (未在CSV中) = rigidScore ≈ rigid energy matrix (应该更不负)")
    print()

    # 分析前20个conformations的详细数据
    print("\n详细数据 (前14个):")
    print("-" * 120)
    print(f"{'Order':<6} {'Score':<12} {'Lower':<12} {'Final':<12} {'排序(从最负到最不负)':<35} {'Final更接近':<15}")
    print("-" * 120)

    for conf in conformations:
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        # 按从负到正排序（即从低能量到高能量）
        values = [
            ('Score', score),
            ('Lower', lower),
            ('Final', final)
        ]
        values.sort(key=lambda x: x[1])  # 从最负到最不负

        order_str = ' < '.join([v[0] for v in values])

        # 计算距离
        dist_to_score = abs(final - score)
        dist_to_lower = abs(final - lower)
        if dist_to_score < dist_to_lower:
            closer = f"Score ({dist_to_score:.3f})"
        elif dist_to_lower < dist_to_score:
            closer = f"Lower ({dist_to_lower:.3f})"
        else:
            closer = "Equal"

        print(f"{conf['order']:<6} {score:<12.3f} {lower:<12.3f} {final:<12.3f} {order_str:<35} {closer:<15}")

    print("\n" + "=" * 100)
    print("统计分析:")
    print("=" * 100)

    # 统计Score vs Lower vs Final的排序
    patterns = {}

    # 统计Final相对于Score和Lower的位置
    final_lowest = 0      # Final最负（最低能量）
    final_middle = 0      # Final在中间
    final_highest = 0     # Final最不负（最高能量）

    # 统计Final更接近哪个
    closer_to_score = 0
    closer_to_lower = 0
    equal_dist = 0

    for conf in conformations:
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        # 找出最负和最不负的bound
        most_negative = min(score, lower)      # 最低能量的bound
        least_negative = max(score, lower)     # 最高能量的bound

        # 判断Final的位置
        if final < most_negative:
            final_lowest += 1
            position = "Final最低(最负)"
        elif final > least_negative:
            final_highest += 1
            position = "Final最高(最不负)"
        else:
            final_middle += 1
            position = "Final在中间"

        # 完整排序
        if score < lower < final:
            pattern = "Score < Lower < Final"
        elif score < final < lower:
            pattern = "Score < Final < Lower"
        elif lower < score < final:
            pattern = "Lower < Score < Final"
        elif lower < final < score:
            pattern = "Lower < Final < Score"
        elif final < score < lower:
            pattern = "Final < Score < Lower"
        elif final < lower < score:
            pattern = "Final < Lower < Score"
        else:
            pattern = "其他"

        patterns[pattern] = patterns.get(pattern, 0) + 1

        # 距离统计
        dist_to_score = abs(final - score)
        dist_to_lower = abs(final - lower)
        if abs(dist_to_score - dist_to_lower) < 1e-6:
            equal_dist += 1
        elif dist_to_score < dist_to_lower:
            closer_to_score += 1
        else:
            closer_to_lower += 1

    total = len(conformations)

    print("\n1. Final Energy相对于Bounds的位置:")
    print(f"   Final最负(能量最低):    {final_lowest:3d} ({final_lowest*100.0/total:5.1f}%)")
    print(f"   Final在两个bounds之间:  {final_middle:3d} ({final_middle*100.0/total:5.1f}%)")
    print(f"   Final最不负(能量最高):  {final_highest:3d} ({final_highest*100.0/total:5.1f}%)")

    print("\n2. 完整排序模式 (从最负到最不负 = 从最低到最高能量):")
    for pattern, count in sorted(patterns.items(), key=lambda x: -x[1]):
        pct = count * 100.0 / total
        print(f"   {pattern:<30}: {count:3d} ({pct:5.1f}%)")

    print("\n3. Final更接近哪个bound:")
    print(f"   更接近Score:      {closer_to_score:3d} ({closer_to_score*100.0/total:5.1f}%)")
    print(f"   更接近LowerBound: {closer_to_lower:3d} ({closer_to_lower*100.0/total:5.1f}%)")
    print(f"   等距:             {equal_dist:3d} ({equal_dist*100.0/total:5.1f}%)")

    # 计算平均距离
    avg_dist_score = sum(abs(c['final'] - c['score']) for c in conformations) / total
    avg_dist_lower = sum(abs(c['final'] - c['lower']) for c in conformations) / total

    print(f"\n   平均距离到Score:      {avg_dist_score:.4f} kcal/mol")
    print(f"   平均距离到LowerBound: {avg_dist_lower:.4f} kcal/mol")

    print("\n" + "=" * 100)
    print("结论:")
    print("=" * 100)

    print("\n▸ 关键发现：")

    if final_lowest > total * 0.5:
        print(f"\n  1. Final energy是最负的（最低能量）: {final_lowest}/{total}")
        print("     → minimize后得到比所有预估bounds更低的能量")
        print("     → 这是正常的！说明minimize确实找到了更稳定的结构")
        print("     → Score和Lower都是未minimize前的估计，会偏高")
    elif final_middle > total * 0.5:
        print(f"\n  1. Final energy在bounds之间: {final_middle}/{total}")
        print("     → minimize后的能量在预估范围内")
        print("     → 说明bounds估计比较准确")
    else:
        print(f"\n  1. Final energy是最不负的（最高能量）: {final_highest}/{total}")
        print("     → ⚠ 这很奇怪，minimize应该降低能量")

    # 分析Score vs Lower的关系
    score_lower_count = 0  # Score < Lower (Score更负)
    score_greater_count = 0  # Score > Lower (Score更不负)
    score_equal_count = 0

    for conf in conformations:
        if abs(conf['score'] - conf['lower']) < 1e-6:
            score_equal_count += 1
        elif conf['score'] < conf['lower']:
            score_lower_count += 1
        else:
            score_greater_count += 1

    print(f"\n  2. Score vs LowerBound:")
    print(f"     Score < Lower (Score更负):  {score_lower_count:3d} ({score_lower_count*100.0/total:5.1f}%)")
    print(f"     Score = Lower:             {score_equal_count:3d} ({score_equal_count*100.0/total:5.1f}%)")
    print(f"     Score > Lower (Score更不负): {score_greater_count:3d} ({score_greater_count*100.0/total:5.1f}%)")

    if score_lower_count > total * 0.5:
        print("     → oldgscore给出了更紧的下界（更负）")
        print("     → 理论上应该是: Final ≤ oldgscore ≤ Lower ≤ rigidScore")
    elif score_greater_count > total * 0.5:
        print("     → LowerBound (minimizingEmat)给出了更紧的下界（更负）")
        print("     → 理论上应该是: Final ≤ Lower ≤ oldgscore ≤ rigidScore")
    else:
        print("     → 两者大致相同")

    if closer_to_score > closer_to_lower:
        print(f"\n  3. Final更接近Score: {closer_to_score}/{total}")
        print(f"     → oldgscore是更准确的估计")
        print(f"     → 平均距离比LowerBound近 {avg_dist_lower/avg_dist_score:.2f}x")
    else:
        print(f"\n  3. Final更接近LowerBound: {closer_to_lower}/{total}")
        print(f"     → minimizingEmat.confE()是更准确的估计")
        print(f"     → 平均距离比Score近 {avg_dist_score/avg_dist_lower:.2f}x")

    # 根据最常见的pattern给出排序
    most_common = max(patterns.items(), key=lambda x: x[1])
    print(f"\n  4. 最常见的排序（{most_common[1]}/{total} = {most_common[1]*100.0/total:.1f}%）：")
    print(f"     {most_common[0]}")
    print("     (从左到右 = 从最低能量到最高能量)")

    # 推断完整的能量排序
    print("\n▸ 推断的完整能量排序（从最低到最高）：")

    # 计算平均值
    avg_final = sum(c['final'] for c in conformations) / total
    avg_score = sum(c['score'] for c in conformations) / total
    avg_lower = sum(c['lower'] for c in conformations) / total

    energy_order = [
        ('Final (minimize后真实能量)', avg_final),
        ('Score (oldgscore)', avg_score),
        ('LowerBound (minimizingEmat)', avg_lower),
    ]
    energy_order.sort(key=lambda x: x[1])

    for i, (name, val) in enumerate(energy_order):
        print(f"     {i+1}. {name:<35} = {val:10.3f} kcal/mol")

    print(f"     {len(energy_order)+1}. {'confUpperBound (rigidScore, 未记录)':<35} ≈ [应该更高]")

    print("\n▸ 关于Error Bound和minimize决策：")
    avg_error = sum(abs(c['score'] - c['lower']) for c in conformations) / total
    print(f"     • Score和LowerBound的平均差异: {avg_error:.4f} kcal/mol")
    print(f"     • MARK*使用 errorBound = (subtreeUpper - subtreeLower) × minimizationRatio")
    print(f"     • 优先minimize error bound最大的nodes")
    print(f"     • 这些就是对partition function上下界差距贡献最大的nodes")

    print("\n" + "=" * 100)

if __name__ == '__main__':
    csv_file = '/home/users/lz280/IdeaProjects/OSPREY3/markstar_order_vs_contribution_4flex_eps0p01.csv'
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    analyze_energy_ordering(csv_file)
