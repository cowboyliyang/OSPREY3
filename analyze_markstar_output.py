#!/usr/bin/env python3
"""
分析MARK*的输出来回答两个问题：
1. error bound在运行过程中是否会变大？
2. 优先被minimize的conformation是否对最终partition function贡献最大？
"""

import re
import sys
from collections import defaultdict

class MARKStarAnalyzer:
    def __init__(self):
        self.epsilon_history = []
        self.bound_history = []
        self.minimized_confs = []  # (order, conf, energy, initial_lower_bound)

    def parse_output(self, lines):
        """解析MARK*的输出"""
        for line in lines:
            # 查找epsilon bound的更新
            # Pattern: "Current overall error bound: 0.xxxx, spread of [lower, upper]"
            match = re.search(r'Current overall error bound:\s+([\d.]+),\s+spread of\s+\[([\d.e+-]+),\s+([\d.e+-]+)\]', line)
            if match:
                epsilon = float(match.group(1))
                lower = float(match.group(2))
                upper = float(match.group(3))
                self.epsilon_history.append(epsilon)
                self.bound_history.append((lower, upper))

            # 查找被minimize的conformation
            # Pattern: conf info with energy bounds
            match = re.search(r'\[([\d,\s]+)\]conf:\s*(\d+).*energy:\s*([\d.e+-]+).*bounds:\[([\d.e+-]+),\s*([\d.e+-]+)\].*delta:\s*([\d.]+)', line)
            if match:
                conf_str = match.group(1)
                conf_num = int(match.group(2))
                energy = float(match.group(3))
                lower_bound = float(match.group(4))
                upper_bound = float(match.group(5))
                delta = float(match.group(6))

                self.minimized_confs.append({
                    'order': conf_num,
                    'conf': conf_str,
                    'energy': energy,
                    'lower_bound': lower_bound,
                    'upper_bound': upper_bound,
                    'delta': delta
                })

    def analyze_epsilon_monotonicity(self):
        """分析epsilon bound是否单调递减"""
        print("\n" + "="*70)
        print("问题1: Error Bound是否会在运行过程中变大？")
        print("="*70)

        if len(self.epsilon_history) < 2:
            print("没有足够的epsilon history数据")
            return

        print(f"\n总共记录了 {len(self.epsilon_history)} 次epsilon更新\n")
        print("Iteration | Epsilon    | Change    | Status")
        print("-"*50)

        increases = []
        for i, eps in enumerate(self.epsilon_history):
            if i == 0:
                print(f"{i:9d} | {eps:10.8f} | {'N/A':9s} | Initial")
            else:
                prev_eps = self.epsilon_history[i-1]
                change = eps - prev_eps
                status = "INCREASE!" if change > 1e-6 else "decrease" if change < -1e-6 else "~same"

                print(f"{i:9d} | {eps:10.8f} | {change:+9.6f} | {status}")

                if change > 1e-6:
                    increases.append((i, prev_eps, eps, change))

        print("\n" + "="*70)
        if increases:
            print(f"✗ 发现 {len(increases)} 次epsilon INCREASE:")
            for iter_num, prev, curr, change in increases:
                print(f"  Iteration {iter_num}: {prev:.8f} -> {curr:.8f} (增加了 {change:.8f})")
            print("\n结论: Error bound在某些迭代中确实变大了!")
        else:
            print("✓ Error bound从未增大 - 完全单调递减!")
            print("\n结论: Error bound严格单调递减")
        print("="*70)

    def analyze_contribution_correlation(self):
        """分析早期minimize的conf是否贡献更大"""
        print("\n" + "="*70)
        print("问题2: 优先minimize的conformation是否对partition function贡献最大？")
        print("="*70)

        if len(self.minimized_confs) < 2:
            print("没有足够的minimized conformation数据")
            return

        # 计算每个conf的Boltzmann weight (exp(-E/RT))
        # 假设RT=0.001987*300 ≈ 0.6 kcal/mol
        RT = 0.5961621

        for conf in self.minimized_confs:
            conf['boltzmann'] = 2.71828 ** (-conf['energy'] / RT)  # e^(-E/RT)

        total_weight = sum(c['boltzmann'] for c in self.minimized_confs)

        # 按minimize顺序排序
        sorted_confs = sorted(self.minimized_confs, key=lambda x: x['order'])

        print(f"\n总共minimize了 {len(sorted_confs)} 个conformations")
        print(f"Total Boltzmann weight: {total_weight:.6e}\n")

        print("前20个被minimize的conformations:")
        print("Order | Energy    | Boltzmann Weight | % Contribution | Bounds Spread")
        print("-"*75)

        cumulative_contribution = 0
        for i, conf in enumerate(sorted_confs[:20]):
            contribution_pct = (conf['boltzmann'] / total_weight) * 100
            cumulative_contribution += contribution_pct
            bounds_spread = conf['upper_bound'] - conf['lower_bound']

            print(f"{conf['order']:5d} | {conf['energy']:9.3f} | {conf['boltzmann']:16.6e} | {contribution_pct:13.6f}% | {bounds_spread:13.6e}")

        # 分析前10%, 25%, 50%的贡献
        percentiles = [10, 25, 50]
        print(f"\n累积贡献分析:")
        print("-"*50)

        for pct in percentiles:
            n_confs = max(1, len(sorted_confs) * pct // 100)
            top_confs = sorted_confs[:n_confs]
            top_contribution = sum(c['boltzmann'] for c in top_confs)
            contribution_pct = (top_contribution / total_weight) * 100

            print(f"前 {pct:2d}% ({n_confs:3d} confs) 贡献: {contribution_pct:6.2f}% of total Z")

        # 结论
        print("\n" + "="*70)
        top_10_pct_contribution = (sum(c['boltzmann'] for c in sorted_confs[:max(1, len(sorted_confs)//10)]) / total_weight) * 100

        if top_10_pct_contribution > 50:
            print(f"✓ 前10%的conformations贡献了 {top_10_pct_contribution:.1f}% 的partition function")
            print("结论: 优先minimize的conformations确实对Z贡献最大!")
        elif top_10_pct_contribution > 30:
            print(f"○ 前10%的conformations贡献了 {top_10_pct_contribution:.1f}% 的partition function")
            print("结论: 优先minimize的conformations对Z有显著贡献，但不是压倒性的")
        else:
            print(f"✗ 前10%的conformations只贡献了 {top_10_pct_contribution:.1f}% 的partition function")
            print("结论: 优先minimize的conformations并非主要贡献者")
        print("="*70)

def main():
    if len(sys.argv) > 1:
        # 从文件读取
        with open(sys.argv[1], 'r') as f:
            lines = f.readlines()
    else:
        # 从stdin读取
        print("请提供MARK*的输出文件，或者通过stdin输入")
        print("用法: python analyze_markstar_output.py <output_file>")
        print("或者: ./gradlew test --tests TestMARKStar.test1GUASmall 2>&1 | python analyze_markstar_output.py")
        return

    analyzer = MARKStarAnalyzer()
    analyzer.parse_output(lines)
    analyzer.analyze_epsilon_monotonicity()
    analyzer.analyze_contribution_correlation()

if __name__ == "__main__":
    main()
