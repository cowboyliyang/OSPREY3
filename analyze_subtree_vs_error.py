#!/usr/bin/env python3
"""
分析Subtree Upper Bound vs Error Bound作为prioritization指标的效果

现在CSV包含了真正的：
- SubtreeUpper: 对Z的贡献上界
- ErrorBound: MARK*用于prioritization的error bound
"""

import csv
import sys

def calculate_spearman(list1, list2, conformations):
    """计算Spearman秩相关系数"""
    if len(conformations) < 2:
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

def analyze_prioritization(csv_file):
    conformations = []

    # 读取CSV
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            conf = {
                'order': int(row['MinimizeOrder']),
                'conf': row['Conformation'],
                'subtree_lower': float(row['SubtreeLower']),
                'subtree_upper': float(row['SubtreeUpper']),
                'error_bound': float(row['ErrorBound']),
                'contribution': float(row['PercentContribution']),
                'rank': int(row['ContributionRank'])
            }
            conformations.append(conf)

    total = len(conformations)

    print("=" * 100)
    print(f"Subtree Upper vs Error Bound Prioritization分析 ({total} conformations)")
    print("=" * 100)
    print()

    # 按Contribution排序
    sorted_by_contrib = sorted(conformations, key=lambda x: x['contribution'], reverse=True)

    # 按SubtreeUpper排序（大的优先）
    sorted_by_subtree = sorted(conformations, key=lambda x: x['subtree_upper'], reverse=True)

    # 按ErrorBound排序（大的优先）
    sorted_by_error = sorted(conformations, key=lambda x: x['error_bound'], reverse=True)

    # 实际minimize顺序
    sorted_by_actual = sorted(conformations, key=lambda x: x['order'])

    # ========================================
    # 1. Subtree Upper策略
    # ========================================
    print("=" * 100)
    print("策略1: Subtree Upper Bound (对Z的贡献上界，大的优先)")
    print("=" * 100)
    print()

    print("如果按SubtreeUpper prioritize，顺序应该是：")
    print(f"{'Rank':<6} {'Order':<6} {'SubtreeUpper':<15} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 70)
    for rank, conf in enumerate(sorted_by_subtree[:15], 1):
        print(f"{rank:<6} {conf['order']:<6} {conf['subtree_upper']:<15.6e} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    rho_subtree = calculate_spearman(sorted_by_subtree, sorted_by_contrib, conformations)
    print(f"\n📊 SubtreeUpper vs Contribution的Spearman相关系数: ρ = {rho_subtree:.4f}")

    # ========================================
    # 2. Error Bound策略
    # ========================================
    print()
    print("=" * 100)
    print("策略2: Error Bound (MARK*当前策略，大的优先)")
    print("=" * 100)
    print()

    print("如果按ErrorBound prioritize，顺序应该是：")
    print(f"{'Rank':<6} {'Order':<6} {'ErrorBound':<15} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 70)
    for rank, conf in enumerate(sorted_by_error[:15], 1):
        print(f"{rank:<6} {conf['order']:<6} {conf['error_bound']:<15.6e} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    rho_error = calculate_spearman(sorted_by_error, sorted_by_contrib, conformations)
    print(f"\n📊 ErrorBound vs Contribution的Spearman相关系数: ρ = {rho_error:.4f}")

    # ========================================
    # 3. 实际顺序
    # ========================================
    print()
    print("=" * 100)
    print("策略3: MARK*实际minimize顺序")
    print("=" * 100)
    print()

    print("MARK*实际minimize的顺序：")
    print(f"{'Order':<6} {'SubtreeUpper':<15} {'ErrorBound':<15} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 80)
    for conf in sorted_by_actual[:15]:
        print(f"{conf['order']:<6} {conf['subtree_upper']:<15.6e} {conf['error_bound']:<15.6e} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    rho_actual = calculate_spearman(sorted_by_actual, sorted_by_contrib, conformations)
    print(f"\n📊 实际顺序 vs Contribution的Spearman相关系数: ρ = {rho_actual:.4f}")

    # ========================================
    # 对比分析
    # ========================================
    print()
    print("=" * 100)
    print("对比分析")
    print("=" * 100)
    print()

    print(f"SubtreeUpper策略:  ρ = {rho_subtree:.4f}")
    print(f"ErrorBound策略:    ρ = {rho_error:.4f}")
    print(f"MARK*实际顺序:     ρ = {rho_actual:.4f}")
    print()

    # 判断最好的策略
    best_rho = max(rho_subtree, rho_error, rho_actual)
    if rho_subtree == best_rho:
        print("🎯 SubtreeUpper策略最好！")
        if rho_subtree > rho_error:
            diff = rho_subtree - rho_error
            print(f"   比ErrorBound好 {diff:.4f} ({diff*100/rho_error:.1f}%)")
        if rho_subtree > rho_actual:
            diff = rho_subtree - rho_actual
            print(f"   比MARK*实际顺序好 {diff:.4f} ({diff*100/rho_actual:.1f}%)")
    elif rho_error == best_rho:
        print("🎯 ErrorBound策略最好！")
        if rho_error > rho_subtree:
            diff = rho_error - rho_subtree
            print(f"   比SubtreeUpper好 {diff:.4f} ({diff*100/rho_subtree:.1f}%)")
        if rho_error > rho_actual:
            diff = rho_error - rho_actual
            print(f"   比MARK*实际顺序好 {diff:.4f} ({diff*100/rho_actual:.1f}%)")
    else:
        print("🎯 MARK*实际顺序最好！")
        if rho_actual > rho_subtree:
            diff = rho_actual - rho_subtree
            print(f"   比SubtreeUpper好 {diff:.4f} ({diff*100/rho_subtree:.1f}%)")
        if rho_actual > rho_error:
            diff = rho_actual - rho_error
            print(f"   比ErrorBound好 {diff:.4f} ({diff*100/rho_error:.1f}%)")

    # ========================================
    # Top-K分析
    # ========================================
    print()
    print("=" * 100)
    print("Top-K Overlap分析")
    print("=" * 100)
    print()

    for k in [3, 5, 10, min(20, total)]:
        if k > total:
            continue

        top_k_subtree = set(conf['order'] for conf in sorted_by_subtree[:k])
        top_k_error = set(conf['order'] for conf in sorted_by_error[:k])
        top_k_actual = set(conf['order'] for conf in sorted_by_actual[:k])
        top_k_contrib = set(conf['order'] for conf in sorted_by_contrib[:k])

        overlap_subtree = len(top_k_subtree & top_k_contrib)
        overlap_error = len(top_k_error & top_k_contrib)
        overlap_actual = len(top_k_actual & top_k_contrib)

        print(f"Top-{k}分析：")
        print(f"  SubtreeUpper策略:  {overlap_subtree}/{k} ({overlap_subtree*100/k:.1f}%)")
        print(f"  ErrorBound策略:    {overlap_error}/{k} ({overlap_error*100/k:.1f}%)")
        print(f"  MARK*实际顺序:     {overlap_actual}/{k} ({overlap_actual*100/k:.1f}%)")

        best = max(overlap_subtree, overlap_error, overlap_actual)
        if overlap_subtree == best:
            print(f"  ✓ SubtreeUpper最好")
        elif overlap_error == best:
            print(f"  ✓ ErrorBound最好")
        else:
            print(f"  ✓ 实际顺序最好")
        print()

    # ========================================
    # 累积贡献分析
    # ========================================
    print("=" * 100)
    print("累积贡献分析")
    print("=" * 100)
    print()

    def cumulative_contribution(sorted_list):
        cumulative = []
        total_contrib = 0
        for conf in sorted_list:
            total_contrib += conf['contribution']
            cumulative.append(total_contrib)
        return cumulative

    cum_subtree = cumulative_contribution(sorted_by_subtree)
    cum_error = cumulative_contribution(sorted_by_error)
    cum_actual = cumulative_contribution(sorted_by_actual)

    print(f"{'N confs':<10} {'SubtreeUpper':<15} {'ErrorBound':<15} {'实际顺序':<15}")
    print("-" * 60)
    for i in [2, 4, 9, 19, min(49, total-1)]:
        if i >= total:
            i = total - 1
        print(f"{i+1:<10} {cum_subtree[i]:<15.2f}% {cum_error[i]:<15.2f}% {cum_actual[i]:<15.2f}%")

    print()

    # ========================================
    # 结论
    # ========================================
    print("=" * 100)
    print("最终结论")
    print("=" * 100)
    print()

    if rho_subtree > rho_error + 0.05:
        print("✅ 强烈推荐使用SubtreeUpper代替ErrorBound")
        print()
        print("建议修改MARKStarNode.compareTo():")
        print("  当前: return -getErrorBound().compareTo(other.getErrorBound())")
        print("  改为: return -getSubtreeUpperBound().compareTo(other.getSubtreeUpperBound())")
        print()
        print(f"预期改进: 相关系数从 {rho_error:.4f} 提升到 {rho_subtree:.4f}")
    elif rho_error > rho_subtree + 0.05:
        print("✅ 当前ErrorBound策略已经很好")
        print()
        print("Error bound考虑了bounds的gap，在这个系统上表现更优。")
    else:
        print("📊 两种策略差异不大")
        print()
        print("可能需要在更多系统上测试。")

    print()
    print("=" * 100)

if __name__ == '__main__':
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'markstar_order_vs_contribution_15flex_eps0p010.csv'

    try:
        analyze_prioritization(csv_file)
    except FileNotFoundError:
        print(f"错误：找不到文件 {csv_file}")
        print("请先运行测试生成CSV文件。")
    except KeyError as e:
        print(f"错误：CSV文件缺少必需的列 {e}")
        print("请确保CSV包含: SubtreeLower, SubtreeUpper, ErrorBound")
