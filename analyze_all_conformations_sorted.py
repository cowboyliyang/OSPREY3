#!/usr/bin/env python3
"""
提取所有conformations的bounds信息，并按照不同的维度排序
"""

import csv
import sys

def analyze_conformations_sorted(csv_file):
    conformations = []

    # 读取CSV
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            conf = {
                'order': int(row['MinimizeOrder']),
                'conf': row['Conformation'],
                'score': float(row['Score']),          # oldgscore
                'lower_bound': float(row['LowerBound']),  # minimizingEmat.confE()
                'final_energy': float(row['FinalEnergy']),  # 实际minimize后的能量
                'contribution': float(row['PercentContribution']),
                'rank': int(row['ContributionRank'])
            }

            # 计算error bound (假设是score和lower_bound之间的差异)
            conf['error_bound'] = abs(conf['score'] - conf['lower_bound'])

            # 对于能量bounds，我们需要理解：
            # - lower_bound (minimizingEmat) 应该是能量的下界（最负）
            # - score (oldgscore) 可能被correction更新
            # - 根据之前的分析，实际上 lower_bound < score < final_energy
            # 所以我们设：
            conf['energy_lower'] = min(conf['score'], conf['lower_bound'])
            conf['energy_upper'] = max(conf['score'], conf['lower_bound'])

            conformations.append(conf)

    total = len(conformations)
    print("=" * 120)
    print(f"所有 {total} 个Conformations的详细数据和排序分析")
    print("=" * 120)
    print()

    # ========================================
    # 1. 按Lower Bound排序（能量下界，最负的排最前）
    # ========================================
    print("=" * 120)
    print("1️⃣  按 LOWER BOUND 排序（从最低能量到最高能量）")
    print("=" * 120)
    print(f"{'Rank':<5} {'Order':<6} {'Conformation':<25} {'LowerBound':<12} {'Score':<12} {'Final':<12} {'Contribution%':<14} {'ContribRank':<6}")
    print("-" * 120)

    sorted_by_lower = sorted(conformations, key=lambda x: x['lower_bound'])
    for rank, conf in enumerate(sorted_by_lower, 1):
        print(f"{rank:<5} {conf['order']:<6} {conf['conf']:<25} {conf['lower_bound']:<12.4f} "
              f"{conf['score']:<12.4f} {conf['final_energy']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<5}")

    print()
    print("📊 统计：")
    print(f"   最低 Lower Bound: {sorted_by_lower[0]['lower_bound']:.4f} kcal/mol (Order {sorted_by_lower[0]['order']})")
    print(f"   最高 Lower Bound: {sorted_by_lower[-1]['lower_bound']:.4f} kcal/mol (Order {sorted_by_lower[-1]['order']})")
    print(f"   范围: {sorted_by_lower[-1]['lower_bound'] - sorted_by_lower[0]['lower_bound']:.4f} kcal/mol")
    print()

    # ========================================
    # 2. 按Upper Bound排序（能量上界，这里用score或者max(score, lower)）
    # ========================================
    print("=" * 120)
    print("2️⃣  按 UPPER BOUND 排序（使用Score作为上界估计，从最低能量到最高能量）")
    print("=" * 120)
    print(f"{'Rank':<5} {'Order':<6} {'Conformation':<25} {'Score(Upper)':<12} {'LowerBound':<12} {'Final':<12} {'Contribution%':<14} {'ContribRank':<6}")
    print("-" * 120)

    sorted_by_upper = sorted(conformations, key=lambda x: x['score'])
    for rank, conf in enumerate(sorted_by_upper, 1):
        print(f"{rank:<5} {conf['order']:<6} {conf['conf']:<25} {conf['score']:<12.4f} "
              f"{conf['lower_bound']:<12.4f} {conf['final_energy']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<5}")

    print()
    print("📊 统计：")
    print(f"   最低 Score: {sorted_by_upper[0]['score']:.4f} kcal/mol (Order {sorted_by_upper[0]['order']})")
    print(f"   最高 Score: {sorted_by_upper[-1]['score']:.4f} kcal/mol (Order {sorted_by_upper[-1]['order']})")
    print(f"   范围: {sorted_by_upper[-1]['score'] - sorted_by_upper[0]['score']:.4f} kcal/mol")
    print()

    # ========================================
    # 3. 按Error Bound排序（从大到小）
    # ========================================
    print("=" * 120)
    print("3️⃣  按 ERROR BOUND 排序（|Score - LowerBound|，从大到小）")
    print("=" * 120)
    print(f"{'Rank':<5} {'Order':<6} {'Conformation':<25} {'ErrorBound':<12} {'Score':<12} {'LowerBound':<12} {'Contribution%':<14} {'ContribRank':<6}")
    print("-" * 120)

    sorted_by_error = sorted(conformations, key=lambda x: x['error_bound'], reverse=True)
    for rank, conf in enumerate(sorted_by_error, 1):
        print(f"{rank:<5} {conf['order']:<6} {conf['conf']:<25} {conf['error_bound']:<12.4f} "
              f"{conf['score']:<12.4f} {conf['lower_bound']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<5}")

    print()
    print("📊 统计：")
    print(f"   最大 Error Bound: {sorted_by_error[0]['error_bound']:.4f} kcal/mol (Order {sorted_by_error[0]['order']})")
    print(f"   最小 Error Bound: {sorted_by_error[-1]['error_bound']:.4f} kcal/mol (Order {sorted_by_error[-1]['order']})")
    print(f"   平均 Error Bound: {sum(c['error_bound'] for c in conformations)/total:.4f} kcal/mol")
    print()

    # ========================================
    # 4. 按Contribution排序（从大到小）
    # ========================================
    print("=" * 120)
    print("4️⃣  按 PARTITION FUNCTION 贡献排序（从大到小）")
    print("=" * 120)
    print(f"{'Rank':<5} {'Order':<6} {'Conformation':<25} {'Contribution%':<14} {'Final Energy':<12} {'LowerBound':<12} {'Score':<12} {'ErrorBound':<12}")
    print("-" * 120)

    sorted_by_contrib = sorted(conformations, key=lambda x: x['contribution'], reverse=True)
    cumulative = 0
    for rank, conf in enumerate(sorted_by_contrib, 1):
        cumulative += conf['contribution']
        print(f"{rank:<5} {conf['order']:<6} {conf['conf']:<25} {conf['contribution']:<14.4f} "
              f"{conf['final_energy']:<12.4f} {conf['lower_bound']:<12.4f} {conf['score']:<12.4f} {conf['error_bound']:<12.4f}")

    print()
    print("📊 统计：")
    print(f"   最大贡献: {sorted_by_contrib[0]['contribution']:.4f}% (Order {sorted_by_contrib[0]['order']})")
    print(f"   最小贡献: {sorted_by_contrib[-1]['contribution']:.4f}% (Order {sorted_by_contrib[-1]['order']})")
    top3_contrib = sum(c['contribution'] for c in sorted_by_contrib[:3])
    print(f"   Top-3贡献总和: {top3_contrib:.2f}%")
    top5_contrib = sum(c['contribution'] for c in sorted_by_contrib[:5])
    print(f"   Top-5贡献总和: {top5_contrib:.2f}%")
    print()

    # ========================================
    # 5. 相关性分析
    # ========================================
    print("=" * 120)
    print("5️⃣  相关性分析")
    print("=" * 120)

    # Error Bound vs Contribution
    print("\n🔍 Error Bound vs Contribution 的关系：")
    print(f"{'ErrorBound Rank':<16} {'Order':<6} {'Error Bound':<12} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 70)
    for rank, conf in enumerate(sorted_by_error[:10], 1):
        print(f"{rank:<16} {conf['order']:<6} {conf['error_bound']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    # Lower Bound vs Contribution
    print("\n🔍 Lower Bound vs Contribution 的关系：")
    print(f"{'LowerBound Rank':<16} {'Order':<6} {'Lower Bound':<12} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 70)
    for rank, conf in enumerate(sorted_by_lower[:10], 1):
        print(f"{rank:<16} {conf['order']:<6} {conf['lower_bound']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    # Final Energy vs Contribution
    print("\n🔍 Final Energy vs Contribution 的关系：")
    sorted_by_final = sorted(conformations, key=lambda x: x['final_energy'])
    print(f"{'Final Energy Rank':<18} {'Order':<6} {'Final Energy':<12} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 72)
    for rank, conf in enumerate(sorted_by_final[:10], 1):
        print(f"{rank:<18} {conf['order']:<6} {conf['final_energy']:<12.4f} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    print()

    # ========================================
    # 6. Spearman相关系数计算
    # ========================================
    print("=" * 120)
    print("6️⃣  Spearman 秩相关系数")
    print("=" * 120)
    print()

    # 创建rank映射
    lower_bound_ranks = {conf['order']: rank for rank, conf in enumerate(sorted_by_lower, 1)}
    error_bound_ranks = {conf['order']: rank for rank, conf in enumerate(sorted_by_error, 1)}
    final_energy_ranks = {conf['order']: rank for rank, conf in enumerate(sorted_by_final, 1)}
    contribution_ranks = {conf['order']: conf['rank'] for conf in conformations}

    def calculate_spearman(ranks1, ranks2, confs):
        if len(confs) < 2:
            return 0
        sum_d_squared = 0
        for conf in confs:
            order = conf['order']
            d = ranks1[order] - ranks2[order]
            sum_d_squared += d * d
        n = len(confs)
        return 1.0 - (6.0 * sum_d_squared) / (n * (n * n - 1))

    spearman_lower_contrib = calculate_spearman(lower_bound_ranks, contribution_ranks, conformations)
    spearman_error_contrib = calculate_spearman(error_bound_ranks, contribution_ranks, conformations)
    spearman_final_contrib = calculate_spearman(final_energy_ranks, contribution_ranks, conformations)

    print(f"Lower Bound (最低能量) vs Contribution (最大贡献):  ρ = {spearman_lower_contrib:7.4f}")
    print(f"Error Bound (最大) vs Contribution (最大贡献):     ρ = {spearman_error_contrib:7.4f}")
    print(f"Final Energy (最低能量) vs Contribution (最大贡献): ρ = {spearman_final_contrib:7.4f}")
    print()

    print("📊 解释：")
    print("   • ρ > 0.7: 强相关")
    print("   • ρ > 0.3: 中等相关")
    print("   • ρ < 0.3: 弱相关")
    print()

    if spearman_lower_contrib > 0.7:
        print("   ✓ Lower Bound与Contribution强相关：能量越低的conformations对Z贡献越大")
    if spearman_final_contrib > 0.7:
        print("   ✓ Final Energy与Contribution强相关：最终能量越低的conformations对Z贡献越大")
    if spearman_error_contrib > 0.3:
        print("   ✓ Error Bound与Contribution中等相关：但这是间接关系")

    print()

    # ========================================
    # 7. 总结
    # ========================================
    print("=" * 120)
    print("7️⃣  总结")
    print("=" * 120)
    print()

    print("🎯 关键发现：")
    print()

    # 检查最低lower bound的conformation是否贡献最大
    lowest_lower_order = sorted_by_lower[0]['order']
    lowest_lower_contrib_rank = sorted_by_lower[0]['rank']
    print(f"1. 最低 Lower Bound 的 conformation (Order {lowest_lower_order}):")
    print(f"   - Lower Bound: {sorted_by_lower[0]['lower_bound']:.4f} kcal/mol")
    print(f"   - Contribution: {sorted_by_lower[0]['contribution']:.4f}%")
    print(f"   - Contribution Rank: #{lowest_lower_contrib_rank}")

    # 检查最大error bound的conformation
    largest_error_order = sorted_by_error[0]['order']
    largest_error_contrib_rank = sorted_by_error[0]['rank']
    print(f"\n2. 最大 Error Bound 的 conformation (Order {largest_error_order}):")
    print(f"   - Error Bound: {sorted_by_error[0]['error_bound']:.4f} kcal/mol")
    print(f"   - Contribution: {sorted_by_error[0]['contribution']:.4f}%")
    print(f"   - Contribution Rank: #{largest_error_contrib_rank}")

    # 检查最大贡献的conformation
    largest_contrib_order = sorted_by_contrib[0]['order']
    print(f"\n3. 最大 Contribution 的 conformation (Order {largest_contrib_order}):")
    print(f"   - Contribution: {sorted_by_contrib[0]['contribution']:.4f}%")
    print(f"   - Lower Bound: {sorted_by_contrib[0]['lower_bound']:.4f} kcal/mol")
    print(f"   - Final Energy: {sorted_by_contrib[0]['final_energy']:.4f} kcal/mol")
    print(f"   - Error Bound: {sorted_by_contrib[0]['error_bound']:.4f} kcal/mol")

    print()
    print("🔑 核心洞察：")
    if spearman_final_contrib > 0.7:
        print("   ✓ 能量越低 → 对Z贡献越大（强相关）")
    print(f"   • Top-3 conformations 贡献了 {top3_contrib:.1f}% 的 partition function")
    print(f"   • Top-5 conformations 贡献了 {top5_contrib:.1f}% 的 partition function")
    print(f"   • MARK* minimize了 {total} 个conformations 达到 ε ≤ 0.01")

    print()
    print("=" * 120)

if __name__ == '__main__':
    csv_file = '/home/users/lz280/IdeaProjects/OSPREY3/markstar_order_vs_contribution_4flex_eps0p01.csv'
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    analyze_conformations_sorted(csv_file)
