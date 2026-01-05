#!/usr/bin/env python3
"""
基于现有CSV数据分析Error Bound vs Lower Bound的效果

虽然CSV中没有真正的error bound，但我们可以：
1. 用lower bound作为能量下界的估计
2. 计算基于lower bound的 prioritization vs contribution
3. 对比MARK*实际的minimize顺序 vs 理论最优顺序

关键问题：
- Lower bound prioritization: 能量最低的先minimize
- Error bound prioritization: error最大的先minimize (MARK*当前策略)
- 哪个更好？
"""

import csv
import sys
import math

# Boltzmann常数和温度
R = 1.9872036e-3  # kcal/(mol·K)
T = 298.15  # K (室温)
RT = R * T  # ≈ 0.592 kcal/mol

def boltzmann_weight(energy):
    """计算Boltzmann weight = exp(-E/RT)"""
    return math.exp(-energy / RT)

def analyze_prioritization_strategies(csv_file):
    conformations = []

    # 读取CSV
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            conf = {
                'order': int(row['MinimizeOrder']),
                'conf': row['Conformation'],
                'score': float(row['Score']),
                'lower_bound': float(row['LowerBound']),
                'final_energy': float(row['FinalEnergy']),
                'contribution': float(row['PercentContribution']),
                'rank': int(row['ContributionRank'])
            }
            conformations.append(conf)

    total = len(conformations)

    print("=" * 100)
    print("Prioritization策略对比分析")
    print("=" * 100)
    print()

    print(f"总共 {total} 个minimized conformations")
    print()

    # ========================================
    # 策略1：按Lower Bound排序（能量越低越优先）
    # ========================================
    print("=" * 100)
    print("策略1: Lower Bound Prioritization（能量最低的先minimize）")
    print("=" * 100)
    print()

    sorted_by_lower = sorted(conformations, key=lambda x: x['lower_bound'])

    print("如果按Lower Bound优先级minimize，顺序应该是：")
    print(f"{'LB Rank':<8} {'Order':<6} {'LowerBound':<12} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 70)
    for rank, conf in enumerate(sorted_by_lower[:10], 1):
        print(f"{rank:<8} {conf['order']:<6} {conf['lower_bound']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    # 计算Spearman相关系数
    def calculate_spearman(list1, list2):
        """list1和list2是两个排序后的conformation列表"""
        if len(list1) < 2:
            return 0
        # 创建rank映射
        rank1 = {conf['order']: i+1 for i, conf in enumerate(list1)}
        rank2 = {conf['order']: i+1 for i, conf in enumerate(list2)}
        sum_d_squared = 0
        for conf in conformations:
            order = conf['order']
            d = rank1[order] - rank2[order]
            sum_d_squared += d * d
        n = len(conformations)
        return 1.0 - (6.0 * sum_d_squared) / (n * (n * n - 1))

    # 按contribution排序
    sorted_by_contrib = sorted(conformations, key=lambda x: x['contribution'], reverse=True)

    rho_lb_contrib = calculate_spearman(sorted_by_lower, sorted_by_contrib)
    print(f"\n📊 Lower Bound排序 vs Contribution排序的Spearman相关系数: ρ = {rho_lb_contrib:.4f}")

    if rho_lb_contrib > 0.9:
        print("   ✓ 强相关！Lower bound是优秀的prioritization指标")
    elif rho_lb_contrib > 0.7:
        print("   ✓ 相关性高，lower bound是不错的prioritization指标")
    else:
        print("   ○ 相关性中等")

    # ========================================
    # 策略2：MARK*实际的minimize顺序
    # ========================================
    print()
    print("=" * 100)
    print("策略2: MARK*实际的minimize顺序（基于Error Bound）")
    print("=" * 100)
    print()

    sorted_by_actual_order = sorted(conformations, key=lambda x: x['order'])

    print("MARK*实际minimize的顺序：")
    print(f"{'Order':<8} {'Conf':<25} {'LowerBound':<12} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 80)
    for conf in sorted_by_actual_order[:10]:
        print(f"{conf['order']:<8} {conf['conf']:<25} {conf['lower_bound']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    rho_actual_contrib = calculate_spearman(sorted_by_actual_order, sorted_by_contrib)
    print(f"\n📊 实际顺序 vs Contribution排序的Spearman相关系数: ρ = {rho_actual_contrib:.4f}")

    if rho_actual_contrib > 0.9:
        print("   ✓ 强相关！MARK*的error bound prioritization很好")
    elif rho_actual_contrib > 0.7:
        print("   ✓ 相关性高，error bound prioritization不错")
    else:
        print("   ○ 相关性中等")

    # ========================================
    # 对比分析
    # ========================================
    print()
    print("=" * 100)
    print("对比分析：Lower Bound vs Error Bound Prioritization")
    print("=" * 100)
    print()

    print(f"Lower Bound策略相关系数: ρ = {rho_lb_contrib:.4f}")
    print(f"MARK* Error Bound策略相关系数: ρ = {rho_actual_contrib:.4f}")
    print()

    if rho_lb_contrib > rho_actual_contrib:
        diff = rho_lb_contrib - rho_actual_contrib
        print(f"🎯 Lower Bound策略更好！（高 {diff:.4f}）")
        print()
        print("结论：")
        print("  • 能量lower bound比error bound更能预测conformation的重要性")
        print("  • 建议改进MARK*的prioritization策略，增加lower bound的权重")
        print(f"  • 可能的改进：priority = f(errorBound, lowerBound)")
    elif rho_actual_contrib > rho_lb_contrib:
        diff = rho_actual_contrib - rho_lb_contrib
        print(f"🎯 MARK* Error Bound策略更好！（高 {diff:.4f}）")
        print()
        print("结论：")
        print("  • MARK*的error bound prioritization确实有效")
        print("  • Error bound不仅考虑能量，还考虑bounds的不确定性")
        print("  • 当前策略已经比纯lower bound更优")
    else:
        print("📊 两种策略效果相近")
        print()
        print("结论：")
        print("  • Lower bound和error bound都是合理的prioritization指标")
        print("  • 可能需要更大规模的测试来区分优劣")

    # ========================================
    # Top-K对比
    # ========================================
    print()
    print("=" * 100)
    print("Top-K overlap分析")
    print("=" * 100)
    print()

    for k in [3, 5, 10]:
        if k > total:
            continue

        # Lower Bound策略的Top-K
        top_k_lb = set(conf['order'] for conf in sorted_by_lower[:k])
        # 实际minimize的Top-K
        top_k_actual = set(conf['order'] for conf in sorted_by_actual_order[:k])
        # 实际贡献的Top-K
        top_k_contrib = set(conf['order'] for conf in sorted_by_contrib[:k])

        overlap_lb = len(top_k_lb & top_k_contrib)
        overlap_actual = len(top_k_actual & top_k_contrib)

        print(f"Top-{k}分析：")
        print(f"  Lower Bound策略匹配贡献Top-{k}: {overlap_lb}/{k} ({overlap_lb*100/k:.1f}%)")
        print(f"  MARK*实际顺序匹配贡献Top-{k}: {overlap_actual}/{k} ({overlap_actual*100/k:.1f}%)")

        if overlap_lb > overlap_actual:
            print(f"  ✓ Lower Bound更好（多{overlap_lb - overlap_actual}个）")
        elif overlap_actual > overlap_lb:
            print(f"  ✓ Error Bound更好（多{overlap_actual - overlap_lb}个）")
        else:
            print(f"  ○ 两者相同")
        print()

    # ========================================
    # 累积贡献分析
    # ========================================
    print("=" * 100)
    print("累积贡献分析")
    print("=" * 100)
    print()

    def cumulative_contribution(sorted_list):
        """计算按给定顺序minimize的累积贡献"""
        cumulative = []
        total_contrib = 0
        for conf in sorted_list:
            total_contrib += conf['contribution']
            cumulative.append(total_contrib)
        return cumulative

    cum_lb = cumulative_contribution(sorted_by_lower)
    cum_actual = cumulative_contribution(sorted_by_actual_order)

    print(f"{'N confs':<10} {'Lower Bound策略':<20} {'MARK*实际顺序':<20} {'差异':<15}")
    print("-" * 70)
    for i in [2, 4, 6, 8, 10, total-1]:
        if i >= total:
            i = total - 1
        diff = cum_lb[i] - cum_actual[i]
        sign = "✓" if diff > 0 else ("✗" if diff < 0 else "=")
        print(f"{i+1:<10} {cum_lb[i]:<20.2f}% {cum_actual[i]:<20.2f}% {sign} {abs(diff):<13.2f}%")

    print()
    print("解释：")
    print("  • 如果Lower Bound策略的累积贡献更高 → 说明它更快找到重要conformations")
    print("  • 如果MARK*实际顺序的累积贡献更高 → 说明error bound策略更有效")

    # ========================================
    # 最终建议
    # ========================================
    print()
    print("=" * 100)
    print("最终建议")
    print("=" * 100)
    print()

    if rho_lb_contrib > rho_actual_contrib + 0.05:  # 明显更好
        print("🎯 强烈建议修改MARK*的prioritization策略")
        print()
        print("当前策略：")
        print("  priority = (subtreeUpper - subtreeLower) × minimizationRatio")
        print("  = error bound（partition function bounds的差距）")
        print()
        print("建议改进：")
        print("  priority = α × errorBound + β × lowerBound_penalty")
        print("  其中 lowerBound_penalty = exp(-lowerBound/RT)")
        print("  （能量越低，penalty越大，优先级越高）")
        print()
        print(f"预期改进：相关系数从 {rho_actual_contrib:.4f} 提升到接近 {rho_lb_contrib:.4f}")

    elif rho_actual_contrib > rho_lb_contrib + 0.05:
        print("✓ MARK*的当前策略已经很好")
        print()
        print("Error bound prioritization比纯lower bound更优，说明：")
        print("  • 考虑bounds的不确定性是重要的")
        print("  • 不仅要找低能量conformations，还要减少整体的bounds gap")
        print("  • 这平衡了exploration和exploitation")

    else:
        print("📊 两种策略差异不大")
        print()
        print("可能的原因：")
        print("  • 测试系统较小（只有14个conformations）")
        print("  • Lower bound与error bound高度相关")
        print()
        print("建议：")
        print("  • 在更大的系统上测试（如6-7个flexible residues）")
        print("  • 分析error bound的具体组成")

    print()
    print("=" * 100)

if __name__ == '__main__':
    csv_file = '/home/users/lz280/IdeaProjects/OSPREY3/markstar_order_vs_contribution_4flex_eps0p01.csv'
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    analyze_prioritization_strategies(csv_file)
