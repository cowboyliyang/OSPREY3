#!/usr/bin/env python3
"""
从PDB CIF文件中提取理想几何参数，生成准确的OSPREY模板文件
"""

import re
import math
from collections import defaultdict

def parse_cif_atoms(cif_file):
    """解析CIF文件中的原子坐标"""
    atoms = []
    with open(cif_file, 'r') as f:
        lines = f.readlines()

    # 找到原子数据部分
    start_idx = None
    for i, line in enumerate(lines):
        if '_chem_comp_atom.pdbx_ordinal' in line:
            start_idx = i + 1
            break

    if start_idx is None:
        return atoms

    # 解析原子
    for line in lines[start_idx:]:
        if line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
            break
        parts = line.split()
        if len(parts) >= 17:
            # CIF格式: comp_id atom_id alt_atom_id type_symbol charge ... model_x model_y model_z ideal_x ideal_y ideal_z ...
            # ALY OH OH O 0 1 N N N N N N 9.990 20.290 -7.423 1.006 -0.557 4.738 ...
            try:
                atom_name = parts[1]
                x_ideal = float(parts[15])
                y_ideal = float(parts[16])
                z_ideal = float(parts[17])
                atoms.append({
                    'name': atom_name,
                    'x': x_ideal,
                    'y': y_ideal,
                    'z': z_ideal
                })
            except (ValueError, IndexError):
                continue
    return atoms

def parse_cif_bonds(cif_file):
    """解析CIF文件中的键连接信息"""
    bonds = []
    with open(cif_file, 'r') as f:
        in_bond_section = False
        for line in f:
            if '_chem_comp_bond.pdbx_ordinal' in line:
                in_bond_section = True
                continue
            if in_bond_section:
                if line.startswith('#') or line.startswith('loop_') or line.startswith('_pdbx'):
                    break
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append((parts[1], parts[2], parts[3]))  # atom1, atom2, order
    return bonds

def distance(atom1, atom2):
    """计算两原子间距离"""
    dx = atom1['x'] - atom2['x']
    dy = atom1['y'] - atom2['y']
    dz = atom1['z'] - atom2['z']
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def angle(atom1, atom2, atom3):
    """计算三原子间的键角（度）"""
    # atom2是中心原子
    v1 = [atom1['x']-atom2['x'], atom1['y']-atom2['y'], atom1['z']-atom2['z']]
    v2 = [atom3['x']-atom2['x'], atom3['y']-atom2['y'], atom3['z']-atom2['z']]

    # 点积和模长
    dot = sum(a*b for a,b in zip(v1,v2))
    len1 = math.sqrt(sum(a*a for a in v1))
    len2 = math.sqrt(sum(a*a for a in v2))

    cos_angle = dot / (len1 * len2)
    cos_angle = max(-1, min(1, cos_angle))  # 防止浮点误差
    return math.degrees(math.acos(cos_angle))

def dihedral(atom1, atom2, atom3, atom4):
    """计算四原子间的二面角（度）"""
    # 基于向量计算
    b1 = [atom2['x']-atom1['x'], atom2['y']-atom1['y'], atom2['z']-atom1['z']]
    b2 = [atom3['x']-atom2['x'], atom3['y']-atom2['y'], atom3['z']-atom2['z']]
    b3 = [atom4['x']-atom3['x'], atom4['y']-atom3['y'], atom4['z']-atom3['z']]

    # 法向量
    n1 = [b1[1]*b2[2]-b1[2]*b2[1], b1[2]*b2[0]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]]
    n2 = [b2[1]*b3[2]-b2[2]*b3[1], b2[2]*b3[0]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]]

    # 归一化b2
    len_b2 = math.sqrt(sum(x*x for x in b2))
    m1 = [n1[1]*b2[2]/len_b2-n1[2]*b2[1]/len_b2,
          n1[2]*b2[0]/len_b2-n1[0]*b2[2]/len_b2,
          n1[0]*b2[1]/len_b2-n1[1]*b2[0]/len_b2]

    x = sum(a*b for a,b in zip(n1,n2))
    y = sum(a*b for a,b in zip(m1,n2))

    angle_rad = math.atan2(y, x)
    return math.degrees(angle_rad)

def analyze_ptm(cif_file, residue_name):
    """分析PTM的几何参数"""
    print(f"\n=== Analyzing {residue_name} from {cif_file} ===\n")

    atoms = parse_cif_atoms(cif_file)
    bonds = parse_cif_bonds(cif_file)

    # 创建原子字典
    atom_dict = {a['name']: a for a in atoms}

    # 打印关键键长
    print("Key bond lengths:")
    key_bonds = [
        ('NZ', 'CE'),
        ('NZ', 'CH' if residue_name == 'ALY' else 'CX'),
        ('CH' if residue_name == 'ALY' else 'CX', 'OH' if residue_name == 'ALY' else 'OX'),
        ('CH' if residue_name == 'ALY' else 'CX', 'CH3' if residue_name == 'ALY' else 'CK'),
    ]

    for atom1, atom2 in key_bonds:
        if atom1 in atom_dict and atom2 in atom_dict:
            d = distance(atom_dict[atom1], atom_dict[atom2])
            print(f"  {atom1}-{atom2}: {d:.3f} Å")

    if residue_name == 'SLL':
        # Succinyl链的额外键长
        extra_bonds = [
            ('CK', 'CL'),
            ('CL', 'CP'),
            ('CP', 'OP1'),
            ('CP', 'OP2'),
        ]
        for atom1, atom2 in extra_bonds:
            if atom1 in atom_dict and atom2 in atom_dict:
                d = distance(atom_dict[atom1], atom_dict[atom2])
                print(f"  {atom1}-{atom2}: {d:.3f} Å")

    # 打印关键键角
    print("\nKey bond angles:")
    key_angles = [
        ('CE', 'NZ', 'CH' if residue_name == 'ALY' else 'CX'),
        ('NZ', 'CH' if residue_name == 'ALY' else 'CX', 'OH' if residue_name == 'ALY' else 'OX'),
        ('NZ', 'CH' if residue_name == 'ALY' else 'CX', 'CH3' if residue_name == 'ALY' else 'CK'),
    ]

    for atom1, atom2, atom3 in key_angles:
        if atom1 in atom_dict and atom2 in atom_dict and atom3 in atom_dict:
            a = angle(atom_dict[atom1], atom_dict[atom2], atom_dict[atom3])
            print(f"  {atom1}-{atom2}-{atom3}: {a:.1f}°")

    # 打印关键二面角
    print("\nKey dihedral angles:")
    key_dihedrals = [
        ('CD', 'CE', 'NZ', 'CH' if residue_name == 'ALY' else 'CX'),
        ('CE', 'NZ', 'CH' if residue_name == 'ALY' else 'CX', 'CH3' if residue_name == 'ALY' else 'CK'),
    ]

    for atom1, atom2, atom3, atom4 in key_dihedrals:
        if all(a in atom_dict for a in [atom1, atom2, atom3, atom4]):
            d = dihedral(atom_dict[atom1], atom_dict[atom2], atom_dict[atom3], atom_dict[atom4])
            print(f"  {atom1}-{atom2}-{atom3}-{atom4}: {d:.1f}°")

# 分析ALY和SLL
analyze_ptm('ALY.cif', 'ALY')
analyze_ptm('SLL.cif', 'SLL')

print("\n=== Comparison with my templates ===")
print("\nMy approximations:")
print("  NZ-CH bond: 1.522 Å (should be ~1.347 Å for amide)")
print("  NZ-CH-O angle: ~109.5° (should be ~120° for sp2)")
print("  All bonds assumed tetrahedral geometry")
print("\nThe CIF files provide the correct geometry for accurate templates.")
