#!/usr/bin/env python3
"""
从CIF文件生成准确的coords.in文件
使用真实的笛卡尔坐标，但保留AMBER电荷参数
"""

# ALY的AMBER电荷（从之前的模板）
ALY_CHARGES = {
    'N': -0.3479, 'H': 0.2747, 'CA': -0.2400, 'HA': 0.1426,
    'CB': -0.0094, 'HB2': 0.0362, 'HB3': 0.0362,
    'CG': 0.0187, 'HG2': 0.0103, 'HG3': 0.0103,
    'CD': -0.0479, 'HD2': 0.0621, 'HD3': 0.0621,
    'CE': -0.0143, 'HE2': 0.1135, 'HE3': 0.1135,
    'NZ': -0.4157, 'HZ': 0.2719,
    'CH': 0.5973, 'OH': -0.5679,
    'CH3': -0.1825, 'HH31': 0.0603, 'HH32': 0.0603, 'HH33': 0.0603,
    'C': 0.7341, 'O': -0.5894,
}

# AMBER原子类型
ALY_TYPES = {
    'N': 'N', 'H': 'H', 'CA': 'CT', 'HA': 'H1',
    'CB': 'CT', 'HB2': 'HC', 'HB3': 'HC',
    'CG': 'CT', 'HG2': 'HC', 'HG3': 'HC',
    'CD': 'CT', 'HD2': 'HC', 'HD3': 'HC',
    'CE': 'CT', 'HE2': 'HP', 'HE3': 'HP',
    'NZ': 'N', 'HZ': 'H',
    'CH': 'C', 'OH': 'O',
    'CH3': 'CT', 'HH31': 'HC', 'HH32': 'HC', 'HH33': 'HC',
    'C': 'C', 'O': 'O',
}

# 从ALY.cif提取的理想坐标
ALY_COORDS = {
    'OH': (1.006, -0.557, 4.738),
    'CH': (-0.053, 0.007, 4.572),
    'CH3': (-0.861, 0.460, 5.761),
    'NZ': (-0.509, 0.231, 3.324),
    'CE': (0.275, -0.208, 2.167),
    'CD': (-0.461, 0.164, 0.879),
    'CG': (0.358, -0.294, -0.327),
    'CB': (-0.377, 0.078, -1.615),
    'CA': (0.442, -0.381, -2.823),
    'N': (1.755, 0.276, -2.802),
    'C': (-0.283, -0.013, -4.092),
    'O': (-0.079, 1.056, -4.616),
    'OXT': (-1.155, -0.873, -4.640),
    'HH31': (-0.336, 0.194, 6.679),
    'HH32': (-0.994, 1.541, 5.719),
    'HH33': (-1.835, -0.027, 5.746),
    'HZ': (-1.357, 0.683, 3.191),
    'HE3': (1.250, 0.278, 2.183),
    'HE2': (0.409, -1.289, 2.209),
    'HD3': (-1.435, -0.322, 0.864),
    'HD2': (-0.594, 1.245, 0.837),
    'HG3': (1.333, 0.192, -0.312),
    'HG2': (0.492, -1.375, -0.285),
    'HB3': (-1.352, -0.409, -1.631),
    'HB2': (-0.511, 1.159, -1.657),
    'HA': (0.575, -1.462, -2.781),
    'H': (1.583, 1.269, -2.842),
    'H2': (2.159, 0.088, -1.896),
    'HXT': (-1.620, -0.637, -5.454),
}

# SLL的AMBER电荷
SLL_CHARGES = {
    'N': -0.3479, 'H': 0.2747, 'CA': -0.2400, 'HA': 0.1426,
    'CB': -0.0094, 'HB2': 0.0362, 'HB3': 0.0362,
    'CG': 0.0187, 'HG2': 0.0103, 'HG3': 0.0103,
    'CD': -0.0479, 'HD2': 0.0621, 'HD3': 0.0621,
    'CE': -0.0143, 'HE2': 0.1135, 'HE3': 0.1135,
    'NZ': -0.4157, 'HZ': 0.2719,
    'CQ': 0.5973, 'OQ1': -0.5679,
    'CS': -0.0094, 'HS2': 0.0362, 'HS3': 0.0362,
    'CT': -0.0094, 'HT2': 0.0362, 'HT3': 0.0362,
    'CU': 0.7341, 'OU1': -0.5894, 'OU2': -0.5894,
    'C': 0.7341, 'O': -0.5894,
}

SLL_TYPES = {
    'N': 'N', 'H': 'H', 'CA': 'CT', 'HA': 'H1',
    'CB': 'CT', 'HB2': 'HC', 'HB3': 'HC',
    'CG': 'CT', 'HG2': 'HC', 'HG3': 'HC',
    'CD': 'CT', 'HD2': 'HC', 'HD3': 'HC',
    'CE': 'CT', 'HE2': 'HP', 'HE3': 'HP',
    'NZ': 'N', 'HZ': 'H',
    'CQ': 'C', 'OQ1': 'O',
    'CS': 'CT', 'HS2': 'HC', 'HS3': 'HC',
    'CT': 'CT', 'HT2': 'HC', 'HT3': 'HC',
    'CU': 'C', 'OU1': 'O2', 'OU2': 'O2',
    'C': 'C', 'O': 'O',
}

# 从SLL.cif提取的理想坐标 (CIF的命名 -> 我的命名映射)
SLL_COORDS_CIF = {
    'C': (-6.032, -0.453, -0.056),     # backbone C
    'N': (-4.914, 1.564, 0.710),       # backbone N
    'O': (-6.681, -0.404, 0.962),      # backbone O
    'CA': (-4.826, 0.435, -0.226),
    'CB': (-3.557, -0.369, 0.062),
    'CD': (-1.059, -0.308, 0.058),
    'CE': (0.169, 0.557, -0.235),
    'CG': (-2.329, 0.496, -0.230),
    'CK': (3.846, -0.441, 0.137),      # -> CS in my template
    'CL': (5.074, 0.423, -0.155),      # -> CT in my template
    'CP': (6.325, -0.368, 0.128),      # -> CU in my template
    'CX': (2.595, 0.350, -0.147),      # -> CQ in my template
    'OX': (2.678, 1.493, -0.544),      # -> OQ1 in my template
    'NZ': (1.385, -0.212, 0.041),
    'OP1': (6.242, -1.507, 0.523),     # -> OU1 in my template
    'OP2': (7.530, 0.192, -0.059),     # -> OU2 in my template
    'OXT': (-6.382, -1.304, -1.034),
    'H': (-4.946, 1.240, 1.665),
    'H2': (-4.150, 2.208, 0.572),
    'HA': (-4.794, 0.812, -1.248),
    'HB': (-3.547, -0.671, 1.109),
    'HBA': (-3.536, -1.255, -0.572),
    'HD': (-1.038, -1.194, -0.576),
    'HDA': (-1.050, -0.609, 1.105),
    'HE': (0.160, 0.859, -1.282),
    'HEA': (0.148, 1.444, 0.399),
    'HG': (-2.338, 0.798, -1.278),
    'HGA': (-2.349, 1.382, 0.404),
    'HK': (3.866, -1.328, -0.497),     # -> HS2 in my template
    'HKA': (3.855, -0.743, 1.184),     # -> HS3 in my template
    'HL': (5.064, 0.725, -1.203),      # -> HT2 in my template
    'HLA': (5.053, 1.310, 0.478),      # -> HT3 in my template
    'HNZ': (1.319, -1.127, 0.358),
    'HOP2': (8.303, -0.355, 0.137),    # OU2的氢
    'HXT': (-7.161, -1.856, -0.878),
}

# 创建映射
SLL_COORDS = {}
name_map = {
    'C': 'C', 'N': 'N', 'O': 'O', 'CA': 'CA', 'CB': 'CB',
    'CD': 'CD', 'CE': 'CE', 'CG': 'CG',
    'CK': 'CS', 'CL': 'CT', 'CP': 'CU',
    'CX': 'CQ', 'OX': 'OQ1', 'NZ': 'NZ',
    'OP1': 'OU1', 'OP2': 'OU2',
    'H': 'H', 'H2': 'H2', 'HA': 'HA',
    'HB': 'HB2', 'HBA': 'HB3',
    'HD': 'HD2', 'HDA': 'HD3',
    'HE': 'HE2', 'HEA': 'HE3',
    'HG': 'HG2', 'HGA': 'HG3',
    'HK': 'HS2', 'HKA': 'HS3',
    'HL': 'HT2', 'HLA': 'HT3',
    'HNZ': 'HZ',
}

for cif_name, my_name in name_map.items():
    if cif_name in SLL_COORDS_CIF:
        SLL_COORDS[my_name] = SLL_COORDS_CIF[cif_name]

# 生成ALY coords
aly_atom_order = ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3',
                  'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ',
                  'CH', 'OH', 'CH3', 'HH31', 'HH32', 'HH33', 'C', 'O']

print("Generating accurate ALY coords...")
with open('K03_coords_accurate.in', 'w') as f:
    f.write(f"ALY {len(aly_atom_order)}\n")
    for atom in aly_atom_order:
        if atom in ALY_COORDS:
            x, y, z = ALY_COORDS[atom]
            charge = ALY_CHARGES[atom]
            atype = ALY_TYPES[atom]
            f.write(f"{atom:4s}  {x:6.3f}f  {y:6.3f}f  {z:6.3f}f  {charge:7.4f}f  {atype}\n")
    f.write("ENDRES\n")
print("Generated K03_coords_accurate.in")

# 生成SLL coords
sll_atom_order = ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3',
                  'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ',
                  'CQ', 'OQ1', 'CS', 'HS2', 'HS3', 'CT', 'HT2', 'HT3',
                  'CU', 'OU1', 'OU2', 'C', 'O']

print("\nGenerating accurate SLL coords...")
with open('K04_coords_accurate.in', 'w') as f:
    f.write(f"SLL {len(sll_atom_order)}\n")
    for atom in sll_atom_order:
        if atom in SLL_COORDS:
            x, y, z = SLL_COORDS[atom]
            charge = SLL_CHARGES[atom]
            atype = SLL_TYPES[atom]
            f.write(f"{atom:4s}  {x:6.3f}f  {y:6.3f}f  {z:6.3f}f  {charge:7.4f}f  {atype}\n")
    f.write("ENDRES\n")
print("Generated K04_coords_accurate.in")

print("\nCoordinates are now from PDB CIF files with AMBER charges and atom types!")
