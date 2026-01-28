#!/usr/bin/env python3
"""
分析minimized energy与初始bounds的排序关系

根据代码分析：
- Score (CSV) = oldgscore = node.gscore在minimize前 = 基于minimizing energy matrix的confLowerBound估计
- LowerBound (CSV) = minimizingEmat.confE() = minimizing energy matrix给出的完整conformation能量（另一个下界）
- FinalEnergy (CSV) = 实际minimize后的能量
- confUpperBound = rigidScore ≈ rigid energy matrix给出的能量（上界）

注意：confUpperBound没有在CSV中，但我们可以推断它应该比Score和LowerBound都高（因为rigid>minimizing）

问题：这些值的排序关系是什么？
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
                'score': float(row['Score']),          # oldgscore (confLowerBound的一个版本)
                'lower': float(row['LowerBound']),     # minimizingEmat.confE()
                'corrected': float(row['Corrected']),  # correctionMatrix.confE()
                'final': float(row['FinalEnergy']),    # 真实minimize后的能量
                'contribution': float(row['PercentContribution'])
            }
            conformations.append(conf)

    print("=" * 100)
    print("分析: Minimized Energy与决定minimize时的Bounds的排序关系")
    print("=" * 100)
    print()

    print("字段说明（根据代码分析）：")
    print("  • Score        = oldgscore = 决定minimize时node.gscore的值（一个confLowerBound的估计）")
    print("  • LowerBound   = minimizingEmat.confE() = minimizing energy matrix的完整conformation能量")
    print("  • Corrected    = correctionMatrix.confE() = 已minimize过的conformations的correction")
    print("  • FinalEnergy  = 实际minimize后测得的真实能量")
    print("  • confUpperBound (未在CSV中) ≈ rigidScore = rigid energy matrix给出的能量")
    print()
    print("理论上的排序关系：")
    print("  1. 能量越低 → Boltzmann weight越高 → 对Z贡献越大")
    print("  2. rigidScore (上界) >= minimizingEmat (下界) >= 真实能量")
    print("  3. oldgscore 应该是confLowerBound的一个版本，理论上也是下界")
    print()
    print("=" * 100)

    # 分析前20个conformations的详细数据
    print("\n详细数据 (前20个, 按minimize顺序):")
    print("-" * 100)
    print(f"{'Order':<6} {'Score':<10} {'Lower':<10} {'Final':<10} {'排序':<25} {'距离分析':<20} {'贡献%':<8}")
    print("-" * 100)

    patterns = {}
    score_vs_lower = {'score>lower': 0, 'score<lower': 0, 'score=lower': 0}
    final_position = {
        'final < both': 0,      # final最小
        'final between': 0,     # final在中间
        'final > both': 0       # final最大
    }

    for i, conf in enumerate(conformations[:20]):
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        # 判断score vs lower
        if abs(score - lower) < 1e-6:
            sv = 'score=lower'
            score_vs_lower['score=lower'] += 1
        elif score > lower:
            sv = 'score>lower'
            score_vs_lower['score>lower'] += 1
        else:
            sv = 'score<lower'
            score_vs_lower['score<lower'] += 1

        # 判断final的位置
        min_bound = min(score, lower)
        max_bound = max(score, lower)

        if final < min_bound:
            fp = 'final < both'
            final_position['final < both'] += 1
            order_str = f"final < {sv}"
        elif final > max_bound:
            fp = 'final > both'
            final_position['final > both'] += 1
            order_str = f"{sv} < final"
        else:
            fp = 'final between'
            final_position['final between'] += 1
            if score > lower:
                order_str = "lower ≤ final ≤ score"
            else:
                order_str = "score ≤ final ≤ lower"

        # 距离分析
        dist_to_score = abs(final - score)
        dist_to_lower = abs(final - lower)
        if dist_to_score < dist_to_lower:
            closer = f"→score ({dist_to_score:.3f})"
        elif dist_to_lower < dist_to_score:
            closer = f"→lower ({dist_to_lower:.3f})"
        else:
            closer = "equal dist"

        patterns[order_str] = patterns.get(order_str, 0) + 1

        print(f"{conf['order']:<6} {score:<10.3f} {lower:<10.3f} {final:<10.3f} {order_str:<25} {closer:<20} {conf['contribution']:<8.2f}")

    # 统计所有conformations
    print("\n" + "=" * 100)
    print("统计结果 (所有 {} 个conformations):".format(len(conformations)))
    print("=" * 100)

    # 重新统计所有
    patterns_all = {}
    score_vs_lower_all = {'score>lower': 0, 'score<lower': 0, 'score=lower': 0}
    final_position_all = {'final < both': 0, 'final between': 0, 'final > both': 0}

    dist_to_score_list = []
    dist_to_lower_list = []

    closer_to_score = 0
    closer_to_lower = 0

    for conf in conformations:
        score = conf['score']
        lower = conf['lower']
        final = conf['final']

        # Score vs Lower
        if abs(score - lower) < 1e-6:
            sv = 'score=lower'
            score_vs_lower_all['score=lower'] += 1
        elif score > lower:
            sv = 'score>lower'
            score_vs_lower_all['score>lower'] += 1
        else:
            sv = 'score<lower'
            score_vs_lower_all['score<lower'] += 1

        # Final position
        min_bound = min(score, lower)
        max_bound = max(score, lower)

        if final < min_bound:
            final_position_all['final < both'] += 1
            if score > lower:
                order_str = "final < lower < score"
            else:
                order_str = "final < score < lower"
        elif final > max_bound:
            final_position_all['final > both'] += 1
            if score > lower:
                order_str = "lower < score < final"
            else:
                order_str = "score < lower < final"
        else:
            final_position_all['final between'] += 1
            if score > lower:
                order_str = "lower ≤ final ≤ score"
            else:
                order_str = "score ≤ final ≤ lower"

        patterns_all[order_str] = patterns_all.get(order_str, 0) + 1

        # 距离统计
        dist_to_score = abs(final - score)
        dist_to_lower = abs(final - lower)
        dist_to_score_list.append(dist_to_score)
        dist_to_lower_list.append(dist_to_lower)

        if dist_to_score < dist_to_lower:
            closer_to_score += 1
        elif dist_to_lower < dist_to_score:
            closer_to_lower += 1

    total = len(conformations)

    print("\n1. Score vs LowerBound 关系:")
    for key in ['score>lower', 'score=lower', 'score<lower']:
        count = score_vs_lower_all[key]
        pct = count * 100.0 / total
        print(f"   {key:<15}: {count:3d} ({pct:5.1f}%)")

    print("\n2. Final Energy 位置:")
    for key in ['final < both', 'final between', 'final > both']:
        count = final_position_all[key]
        pct = count * 100.0 / total
        print(f"   {key:<15}: {count:3d} ({pct:5.1f}%)")

    print("\n3. 完整排序模式分布:")
    for pattern, count in sorted(patterns_all.items(), key=lambda x: -x[1]):
        pct = count * 100.0 / total
        print(f"   {pattern:<30}: {count:3d} ({pct:5.1f}%)")

    print("\n4. Final Energy 距离分析:")
    avg_dist_score = sum(dist_to_score_list) / len(dist_to_score_list)
    avg_dist_lower = sum(dist_to_lower_list) / len(dist_to_lower_list)

    print(f"   到Score的平均距离:      {avg_dist_score:.4f} kcal/mol")
    print(f"   到LowerBound的平均距离: {avg_dist_lower:.4f} kcal/mol")
    print(f"   更接近Score:      {closer_to_score} ({closer_to_score*100.0/total:.1f}%)")
    print(f"   更接近LowerBound: {closer_to_lower} ({closer_to_lower*100.0/total:.1f}%)")

    # 计算error bound (假设是max - min)
    error_bounds = []
    for conf in conformations:
        score = conf['score']
        lower = conf['lower']
        error = abs(score - lower)
        error_bounds.append(error)

    avg_error = sum(error_bounds) / len(error_bounds)
    max_error = max(error_bounds)
    min_error = min(error_bounds)

    print("\n5. Error Bound 分析 (|Score - LowerBound|):")
    print(f"   平均 error bound: {avg_error:.4f} kcal/mol")
    print(f"   最大 error bound: {max_error:.4f} kcal/mol")
    print(f"   最小 error bound: {min_error:.4f} kcal/mol")

    print("\n" + "=" * 100)
    print("结论：")
    print("=" * 100)

    print("\n▸ 决定minimize时的排序关系：")

    if score_vs_lower_all['score=lower'] > total * 0.9:
        print("  ✓ Score = LowerBound (大部分情况)")
        print("    这意味着oldgscore和minimizingEmat.confE()给出了相同的下界估计")
    elif score_vs_lower_all['score>lower'] > total * 0.5:
        print("  ✓ Score > LowerBound (大部分情况)")
        print("    这意味着oldgscore比minimizingEmat.confE()的估计更高（更宽松的下界）")
    else:
        print("  ✓ Score < LowerBound (大部分情况)")
        print("    这意味着oldgscore比minimizingEmat.confE()的估计更低（更紧的下界）")

    print("\n▸ Minimized energy的位置：")

    if final_position_all['final < both'] > total * 0.5:
        print("  ✓ Final < min(Score, LowerBound) (大部分情况)")
        print("    说明minimize后的真实能量低于所有的下界估计")
        print("    → 下界估计不够紧，真实能量更优")
    elif final_position_all['final between'] > total * 0.5:
        print("  ✓ min(Score, LowerBound) ≤ Final ≤ max(Score, LowerBound) (大部分情况)")
        print("    说明minimize后的真实能量在两个估计之间")
        print("    → 这种情况表明bounds估计合理")
    else:
        print("  ✓ Final > max(Score, LowerBound) (大部分情况)")
        print("    说明minimize后的真实能量高于所有的下界估计")
        print("    → ⚠ 这违反了下界的定义！可能字段含义需要重新理解")

    print("\n▸ Final更接近哪个bound：")
    if closer_to_score > closer_to_lower:
        ratio = closer_to_score / max(closer_to_lower, 1)
        print(f"  ✓ Final更接近Score ({closer_to_score}/{total} = {closer_to_score*100.0/total:.1f}%)")
        print(f"    平均距离：Score {avg_dist_score:.4f} vs LowerBound {avg_dist_lower:.4f}")
        print(f"    → Score的估计比LowerBound更准确 ({avg_dist_lower/avg_dist_score:.2f}x)")
    else:
        ratio = closer_to_lower / max(closer_to_score, 1)
        print(f"  ✓ Final更接近LowerBound ({closer_to_lower}/{total} = {closer_to_lower*100.0/total:.1f}%)")
        print(f"    平均距离：Score {avg_dist_score:.4f} vs LowerBound {avg_dist_lower:.4f}")
        print(f"    → LowerBound的估计比Score更准确 ({avg_dist_score/avg_dist_lower:.2f}x)")

    print("\n▸ 关于Error Bound:")
    print(f"  • 平均 |Score - LowerBound| = {avg_error:.4f} kcal/mol")
    print(f"  • 这个差异可能来自于correction matrix的更新")
    print(f"  • Error bound = (subtreeUpper - subtreeLower) × minimizationRatio")
    print(f"  • 决定minimize的优先级取决于这个error bound的大小")

    print("\n▸ 最终排序（从低能量到高能量）：")

    # 采样一些conformations来确定一般模式
    samples = conformations[:5]
    avg_final = sum(c['final'] for c in samples) / len(samples)
    avg_score = sum(c['score'] for c in samples) / len(samples)
    avg_lower = sum(c['lower'] for c in samples) / len(samples)

    values = [
        ('Final (真实能量)', avg_final),
        ('Score (oldgscore)', avg_score),
        ('LowerBound (minimizingEmat)', avg_lower),
    ]

    # 按值排序
    values.sort(key=lambda x: x[1])

    print("  基于前5个conformations的平均值：")
    for i, (name, val) in enumerate(values):
        print(f"    {i+1}. {name:<30} = {val:10.3f} kcal/mol")

    # 添加理论上的上界
    print(f"    4. {'confUpperBound (rigidScore, 未在CSV)':<30} ≈ [更高的值，未记录]")

    print("\n" + "=" * 100)

if __name__ == '__main__':
    csv_file = '/home/users/lz280/IdeaProjects/OSPREY3/markstar_order_vs_contribution_4flex_eps0p01.csv'
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    analyze_energy_ordering(csv_file)
