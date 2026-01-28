#!/usr/bin/env python3
"""
从PDB CIF文件提取真实的ALY几何参数
"""
import math

# 从ALY.cif提取的ideal坐标（第45-47列）
aly_atoms = {
    'NZ':  (-0.509, 0.231, 3.324),
    'HZ':  (-1.357, 0.683, 3.191),
    'CE':  (0.275, -0.208, 2.167),
    'CH':  (-0.053, 0.007, 4.572),  # acetyl carbonyl
    'OH':  (1.006, -0.557, 4.738),  # acetyl oxygen
    'CH3': (-0.861, 0.460, 5.761),  # methyl
}

def distance(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

def angle(a, b, c):
    """计算b为顶点的角度"""
    ba = tuple(a[i]-b[i] for i in range(3))
    bc = tuple(c[i]-b[i] for i in range(3))
    dot = sum(ba[i]*bc[i] for i in range(3))
    len_ba = math.sqrt(sum(x**2 for x in ba))
    len_bc = math.sqrt(sum(x**2 for x in bc))
    return math.degrees(math.acos(dot / (len_ba * len_bc)))

print("=== 真实的ALY几何参数（从PDB CIF） ===\n")

print("1. 关键键长:")
print(f"   CE-NZ:  {distance(aly_atoms['CE'], aly_atoms['NZ']):.3f} Å")
print(f"   NZ-HZ:  {distance(aly_atoms['NZ'], aly_atoms['HZ']):.3f} Å")
print(f"   NZ-CH:  {distance(aly_atoms['NZ'], aly_atoms['CH']):.3f} Å (酰胺C-N)")
print(f"   CH-OH:  {distance(aly_atoms['CH'], aly_atoms['OH']):.3f} Å (C=O)")
print(f"   CH-CH3: {distance(aly_atoms['CH'], aly_atoms['CH3']):.3f} Å (C-C)")

print("\n2. 关键键角:")
print(f"   CE-NZ-CH:  {angle(aly_atoms['CE'], aly_atoms['NZ'], aly_atoms['CH']):.1f}°")
print(f"   NZ-CH-OH:  {angle(aly_atoms['NZ'], aly_atoms['CH'], aly_atoms['OH']):.1f}°")
print(f"   NZ-CH-CH3: {angle(aly_atoms['NZ'], aly_atoms['CH'], aly_atoms['CH3']):.1f}°")
print(f"   CE-NZ-HZ:  {angle(aly_atoms['CE'], aly_atoms['NZ'], aly_atoms['HZ']):.1f}°")

print("\n=== 我当前template中的值（估计值）===")
print("键长:")
print("   NZ-CH: 1.522 Å  (我猜的，基于C-C)")
print("   CH-OH: 1.229 Å  (我猜的，基于C=O)")
print("\n键角:")
print("   CE-NZ-CH: 109.47° (我猜的，基于四面体)")
print("   NZ-CH-OH: 120.00° (我猜的)")

print("\n=== 结论 ===")
print("我的template使用的是【猜测值】，不是真实几何参数！")
print("这会导致:")
print("- 初始结构可能有应变")
print("- 能量计算不够准确")
print("- 但continuous optimization会部分修正")
