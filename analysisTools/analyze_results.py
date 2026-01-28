#!/usr/bin/env python3
"""
分析MARK*的CSV结果
"""

import pandas as pd
import sys
from scipy.stats import spearmanr

def analyze_csv(filename):
    # 读取CSV
    df = pd.read_csv(filename)

    print("="*80)
    print(f"分析文件: {filename}")
    print("="*80)
    print(f"\n总共minimize了 {len(df)} 个conformations\n")

    # 基本统计
    total_contribution = df['PercentContribution'].sum()
    print(f"Total contribution: {total_contribution:.2f}%")
    print(f"Total Boltzmann weight: {df['BoltzmannWeight'].sum():.6e}\n")

    # 显示前10个被minimize的conformations
    print("="*80)
    print("前10个被minimize的conformations:")
    print("-"*80)
    print(f"{'Order':<6} | {'Energy':<10} | {'Boltzmann':<12} | {'% of Z':<10} | {'Rank':<6}")
    print("-"*80)

    for i, row in df.head(10).iterrows():
        print(f"{row['MinimizeOrder']:<6} | {row['FinalEnergy']:<10.4f} | {row['BoltzmannWeight']:<12.4e} | {row['PercentContribution']:<10.4f} | #{row['ContributionRank']:<5}")

    # 累积贡献分析
    print("\n" + "="*80)
    print("累积贡献分析:")
    print("-"*80)

    for pct in [10, 25, 50, 75, 100]:
        n_confs = max(1, len(df) * pct // 100)
        cumulative = df.head(n_confs)['PercentContribution'].sum()
        print(f"前 {pct:3d}% ({n_confs:3d} confs) 贡献: {cumulative:6.2f}% of total Z")

    # Spearman相关系数
    print("\n" + "="*80)
    print("相关性分析:")
    print("-"*80)

    # 计算Spearman
    minimize_order = df['MinimizeOrder'].values
    contribution_rank = df['ContributionRank'].values

    spearman_corr, p_value = spearmanr(minimize_order, contribution_rank)
    print(f"Spearman秩相关系数: {spearman_corr:.4f} (p-value: {p_value:.4e})")

    # Top-K overlap分析
    print()
    for k in [3, 5, 10, min(20, len(df))]:
        if k > len(df):
            continue

        top_k_by_order = set(df.head(k)['Conformation'])
        top_k_by_contribution = set(df.nsmallest(k, 'ContributionRank')['Conformation'])

        overlap = len(top_k_by_order.intersection(top_k_by_contribution))
        overlap_pct = (overlap * 100.0) / k

        print(f"Top-{k} overlap: {overlap}/{k} ({overlap_pct:.1f}%)")

    # 详细对比前5个
    print("\n" + "="*80)
    print("前5个minimize的 vs 前5个贡献最大的:")
    print("-"*80)

    print("\n按Minimize顺序:")
    for i, row in df.head(5).iterrows():
        print(f"  Order {row['MinimizeOrder']}: Energy={row['FinalEnergy']:.4f}, "
              f"Contribution={row['PercentContribution']:.2f}%, Rank=#{row['ContributionRank']}")

    print("\n按贡献大小:")
    top_5_contrib = df.nsmallest(5, 'ContributionRank')
    for i, row in top_5_contrib.iterrows():
        print(f"  Rank #{row['ContributionRank']}: Energy={row['FinalEnergy']:.4f}, "
              f"Contribution={row['PercentContribution']:.2f}%, Order={row['MinimizeOrder']}")

    # 结论
    print("\n" + "="*80)
    print("结论:")
    print("="*80)

    top_10_pct = df.head(max(1, len(df)//10))['PercentContribution'].sum()

    if spearman_corr > 0.7 or top_10_pct > 50:
        print("✓ 优先minimize的conformations确实对Z贡献最大!")
        print(f"  - 前10%的confs贡献了 {top_10_pct:.1f}% 的Z")
        print(f"  - Spearman相关系数: {spearman_corr:.3f}")
        if spearman_corr > 0.8:
            print("  - 相关性非常强！")
    elif spearman_corr > 0.3 or top_10_pct > 30:
        print("○ 有一定相关性，但不是非常强")
        print(f"  - 前10%的confs贡献了 {top_10_pct:.1f}% 的Z")
        print(f"  - Spearman相关系数: {spearman_corr:.3f}")
    else:
        print("✗ minimize顺序和贡献大小关系不大")
        print(f"  - 前10%的confs只贡献了 {top_10_pct:.1f}% 的Z")
        print(f"  - Spearman相关系数: {spearman_corr:.3f}")

    print("="*80)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "markstar_order_vs_contribution_4flex_eps0p01.csv"

    analyze_csv(filename)
