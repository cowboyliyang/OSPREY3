#!/usr/bin/env python3
"""
计算真正的Error Bound用于MARK*的prioritization

根据代码分析：
errorBound = (subtreeUpperBound - subtreeLowerBound) × minimizationRatio

其中：
- subtreeUpperBound = exp(-confLowerBound/RT) × numConfs
- subtreeLowerBound = exp(-confUpperBound/RT) × numConfs

注意：能量下界 → Z上界，能量上界 → Z下界（反转）

从CSV我们有：
- Score = oldgscore ≈ confLowerBound的一个估计
- LowerBound = minimizingEmat.confE() ≈ 另一个confLowerBound估计
- 但我们没有confUpperBound（应该是rigid score，未记录）

我们需要推断：
1. confLowerBound = 决定minimize时的能量下界
2. confUpperBound = 决定minimize时的能量上界（rigid score）
3. 从这两个计算errorBound
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

def analyze_true_error_bounds(csv_file):
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
                'final_energy': float(row['FinalEnergy']),
                'contribution': float(row['PercentContribution']),
                'rank': int(row['ContributionRank'])
            }
            conformations.append(conf)

    print("=" * 100)
    print("真正的Error Bound分析")
    print("=" * 100)
    print()

    print("⚠️ 问题：CSV中没有记录confUpperBound（rigid score）")
    print()
    print("根据代码分析：")
    print("  • subtreeUpperBound = exp(-confLowerBound/RT) × numConfs")
    print("  • subtreeLowerBound = exp(-confUpperBound/RT) × numConfs")
    print("  • errorBound = (subtreeUpper - subtreeLower) × minimizationRatio")
    print()
    print("从CSV我们有：")
    print("  • Score ≈ confLowerBound的一个估计")
    print("  • LowerBound ≈ minimizingEmat给出的confLowerBound")
    print("  • confUpperBound = ??? (rigid score, 未记录)")
    print()

    # 策略1：假设Score和LowerBound中较低的是confLowerBound
    # 策略2：估计confUpperBound比confLowerBound高约1-2 kcal/mol（基于经验）

    print("=" * 100)
    print("策略：基于能量bounds推断Error Bound")
    print("=" * 100)
    print()

    # 对于每个conformation，我们用已知的信息推断
    print("假设1：confLowerBound = min(Score, LowerBound)")
    print("假设2：confUpperBound ≈ confLowerBound + ΔE")
    print("         (ΔE = rigid与minimizing的差异，通常1-2 kcal/mol)")
    print()

    # 从final energy vs lower bounds的差异来估计ΔE
    delta_estimates = []
    for conf in conformations:
        conf_lower = min(conf['score'], conf['lower_bound'])
        # final energy高于lower bounds，说明有一定gap
        # rigid应该更高
        delta = abs(conf['final_energy'] - conf_lower)
        delta_estimates.append(delta)

    avg_delta = sum(delta_estimates) / len(delta_estimates)
    print(f"从数据估计的平均ΔE: {avg_delta:.4f} kcal/mol")
    print()

    # 实际上，让我们重新思考：
    # 根据之前的分析，lower_bound < score < final
    # 如果rigid score是上界，它应该在哪里？

    print("=" * 100)
    print("重新审视：什么是真正的error bound？")
    print("=" * 100)
    print()

    print("根据代码第422行：")
    print("  errorBound = (subtreeUpperBound - subtreeLowerBound) × minimizationRatio")
    print()
    print("其中subtreeUpperBound和subtreeLowerBound是partition function Z的bounds，")
    print("通过Boltzmann变换从能量bounds得到：")
    print("  subtreeUpper = exp(-energyLower/RT) × numConfs")
    print("  subtreeLower = exp(-energyUpper/RT) × numConfs")
    print()

    print("问题的关键：我们需要在minimize**之前**的energyLower和energyUpper")
    print()
    print("从CSV我们只有minimize**之后**的信息，此时：")
    print("  • Score = oldgscore = minimize前的某个估计")
    print("  • LowerBound = minimizingEmat.confE() = matrix预测")
    print("  • FinalEnergy = minimize后的真实能量")
    print()

    print("关键洞察：")
    print("  在决定是否minimize这个conformation时，MARK*使用的bounds是：")
    print("  • confLowerBound = 基于minimizing energy matrix的估计")
    print("  • confUpperBound = 基于rigid energy matrix的估计")
    print()
    print("但CSV中**没有记录rigid score**！")
    print()

    print("=" * 100)
    print("可能的解决方案")
    print("=" * 100)
    print()

    print("方案1：修改MARKStarBound.java，在minimize之前记录error bound")
    print("  需要在processFullConfNode中，minimize之前调用getErrorBound()")
    print()

    print("方案2：根据理论关系推断")
    print("  假设rigid score比minimizing score高固定百分比（如10-20%）")
    print()

    print("方案3：使用现有数据做近似分析")
    print("  用|Score - LowerBound|作为error bound的代理指标")
    print("  虽然不精确，但可以看趋势")
    print()

    print("=" * 100)
    print("方案3的分析：用|Score - LowerBound|作为代理")
    print("=" * 100)
    print()

    # 计算代理error bound
    for conf in conformations:
        conf['proxy_error'] = abs(conf['score'] - conf['lower_bound'])

    # 排序
    sorted_by_proxy_error = sorted(conformations, key=lambda x: x['proxy_error'], reverse=True)
    sorted_by_contrib = sorted(conformations, key=lambda x: x['contribution'], reverse=True)

    print(f"{'Proxy Error Rank':<18} {'Order':<6} {'|Score-Lower|':<14} {'Contribution%':<14} {'Contrib Rank':<12}")
    print("-" * 80)
    for rank, conf in enumerate(sorted_by_proxy_error, 1):
        print(f"{rank:<18} {conf['order']:<6} {conf['proxy_error']:<14.4f} {conf['contribution']:<14.4f} #{conf['rank']:<11}")

    print()
    print("⚠️ 这只是代理指标，不是真正的error bound！")
    print("真正的error bound需要subtreeUpper和subtreeLower，需要:")
    print("  1. confUpperBound (rigid score) - 未记录")
    print("  2. confLowerBound (minimizing score)")
    print("  3. numConfs (该node的conformations数量) - 未记录")
    print("  4. minimizationRatio - 未记录")
    print()

    print("=" * 100)
    print("结论")
    print("=" * 100)
    print()

    print("要正确分析error bound vs contribution的关系，我们需要：")
    print()
    print("✓ 修改代码以记录：")
    print("  1. 每个conformation在minimize**之前**的error bound")
    print("  2. confUpperBound (rigid score)")
    print("  3. confLowerBound (minimizing estimate)")
    print("  4. numConfs和minimizationRatio")
    print()
    print("当前CSV数据不足以计算真正的error bound。")
    print()
    print("建议：修改MARKStarBound.java的printMinimizationOutput()，添加：")
    print("  • node.confUpperBound")
    print("  • node.confLowerBound")
    print("  • curNode.getErrorBound() (在minimize前)")
    print("  • node.getNumConformations()")
    print()

if __name__ == '__main__':
    csv_file = '/home/users/lz280/IdeaProjects/OSPREY3/markstar_order_vs_contribution_4flex_eps0p01.csv'
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    analyze_true_error_bounds(csv_file)
