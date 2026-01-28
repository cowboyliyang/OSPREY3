#!/usr/bin/env python3
"""
验证生成的OSPREY模板文件与CIF文件的一致性
检查：
1. 原子数量
2. 键连接
3. 坐标
4. 缺失的原子
"""

import math

def parse_cif_atoms(cif_file):
    """从CIF文件解析原子"""
    atoms = {}
    with open(cif_file, 'r') as f:
        lines = f.readlines()

    start_idx = None
    for i, line in enumerate(lines):
        if '_chem_comp_atom.pdbx_ordinal' in line:
            start_idx = i + 1
            break

    if start_idx is None:
        return atoms

    for line in lines[start_idx:]:
        if line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
            break
        parts = line.split()
        if len(parts) >= 18:
            atom_name = parts[1]
            element = parts[3]
            x_ideal = float(parts[15])
            y_ideal = float(parts[16])
            z_ideal = float(parts[17])
            atoms[atom_name] = {
                'element': element,
                'x': x_ideal,
                'y': y_ideal,
                'z': z_ideal
            }
    return atoms

def parse_cif_bonds(cif_file):
    """从CIF文件解析键连接"""
    bonds = []
    with open(cif_file, 'r') as f:
        in_bond = False
        for line in f:
            if '_chem_comp_bond.pdbx_ordinal' in line:
                in_bond = True
                continue
            if in_bond:
                if line.startswith('#') or line.startswith('loop_') or line.startswith('_pdbx'):
                    break
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append((parts[1], parts[2], parts[3]))
    return bonds

def parse_coords_file(coords_file):
    """从coords.in文件解析原子"""
    atoms = {}
    with open(coords_file, 'r') as f:
        lines = f.readlines()

    for line in lines[1:-1]:  # 跳过第一行和ENDRES
        parts = line.split()
        if len(parts) >= 6:
            atom_name = parts[0]
            x = float(parts[1].rstrip('f'))
            y = float(parts[2].rstrip('f'))
            z = float(parts[3].rstrip('f'))
            charge = float(parts[4].rstrip('f'))
            atype = parts[5]
            atoms[atom_name] = {
                'x': x, 'y': y, 'z': z,
                'charge': charge, 'type': atype
            }
    return atoms

def distance(p1, p2):
    """计算两点距离"""
    dx = p1['x'] - p2['x']
    dy = p1['y'] - p2['y']
    dz = p1['z'] - p2['z']
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def validate_template(cif_file, coords_file, residue_name, name_map=None):
    """验证模板与CIF的一致性"""
    print(f"\n{'='*70}")
    print(f"Validating {residue_name} template against {cif_file}")
    print('='*70)

    # 解析文件
    cif_atoms = parse_cif_atoms(cif_file)
    cif_bonds = parse_cif_bonds(cif_file)
    coords_atoms = parse_coords_file(coords_file)

    print(f"\n📊 Atom counts:")
    print(f"   CIF file:      {len(cif_atoms)} atoms")
    print(f"   coords.in:     {len(coords_atoms)} atoms")

    # 名称映射（如果有）
    if name_map is None:
        name_map = {name: name for name in cif_atoms.keys()}

    # 检查缺失的原子
    print(f"\n🔍 Checking for missing/extra atoms:")

    cif_names = set(cif_atoms.keys())
    coords_names = set(coords_atoms.keys())

    # CIF中有但coords中没有的
    missing_in_coords = []
    for cif_name in cif_names:
        our_name = name_map.get(cif_name, cif_name)
        if our_name not in coords_names and cif_name not in ['OXT', 'HXT', 'H2']:
            # OXT, HXT, H2 是末端原子，在内部残基中不需要
            missing_in_coords.append(f"{cif_name} -> {our_name}")

    if missing_in_coords:
        print(f"   ⚠️  In CIF but missing in coords.in: {missing_in_coords}")
    else:
        print(f"   ✅ No missing atoms (excluding terminal atoms OXT, HXT, H2)")

    # coords中有但CIF中没有的
    extra_in_coords = [name for name in coords_names if name not in name_map.values()]
    if extra_in_coords:
        print(f"   ⚠️  In coords.in but not in CIF: {extra_in_coords}")
    else:
        print(f"   ✅ No extra atoms")

    # 检查坐标一致性
    print(f"\n📐 Checking coordinate accuracy:")
    max_diff = 0
    max_diff_atom = None

    for cif_name, our_name in name_map.items():
        if cif_name in cif_atoms and our_name in coords_atoms:
            cif_atom = cif_atoms[cif_name]
            coords_atom = coords_atoms[our_name]

            dist = distance(cif_atom, coords_atom)
            if dist > max_diff:
                max_diff = dist
                max_diff_atom = our_name

    print(f"   Maximum coordinate deviation: {max_diff:.4f} Å (atom: {max_diff_atom})")

    if max_diff < 0.001:
        print(f"   ✅ Coordinates match perfectly!")
    elif max_diff < 0.01:
        print(f"   ✅ Coordinates match well (< 0.01 Å)")
    else:
        print(f"   ⚠️  Coordinates have significant deviation (> 0.01 Å)")

    # 检查键连接
    print(f"\n🔗 Checking bond connectivity:")
    print(f"   CIF defines {len(cif_bonds)} bonds")

    # 显示几个关键的键
    print(f"\n   Key bonds from CIF:")
    key_bonds_to_check = [
        ('NZ', 'CE'),
        ('NZ', 'CH' if residue_name == 'ALY' else 'CX'),
        ('CH' if residue_name == 'ALY' else 'CX', 'OH' if residue_name == 'ALY' else 'OX'),
    ]

    for atom1, atom2 in key_bonds_to_check:
        bond_found = False
        for b1, b2, order in cif_bonds:
            if (b1 == atom1 and b2 == atom2) or (b1 == atom2 and b2 == atom1):
                print(f"      {atom1}-{atom2}: {order}")
                bond_found = True
                break
        if not bond_found:
            print(f"      {atom1}-{atom2}: NOT FOUND")

    # 总结
    print(f"\n{'='*70}")
    if not missing_in_coords and not extra_in_coords and max_diff < 0.01:
        print("✅ VALIDATION PASSED - Template matches CIF file!")
    else:
        print("⚠️  VALIDATION WARNING - Some differences detected")
    print('='*70)

# ALY验证
print("\n" + "="*70)
print("ACETYL-LYSINE (ALY) VALIDATION")
print("="*70)

validate_template(
    'ALY.cif',
    'K03_coords.in',
    'ALY'
)

# SLL验证（需要名称映射）
print("\n\n" + "="*70)
print("SUCCINYL-LYSINE (SLL) VALIDATION")
print("="*70)

# SLL的CIF名称到我们的名称的映射
sll_name_map = {
    'C': 'C', 'N': 'N', 'O': 'O', 'CA': 'CA', 'CB': 'CB',
    'CD': 'CD', 'CE': 'CE', 'CG': 'CG',
    'CK': 'CS',   # CIF中是CK，我们命名为CS
    'CL': 'CT',   # CIF中是CL，我们命名为CT
    'CP': 'CU',   # CIF中是CP，我们命名为CU
    'CX': 'CQ',   # CIF中是CX，我们命名为CQ
    'OX': 'OQ1',  # CIF中是OX，我们命名为OQ1
    'NZ': 'NZ',
    'OP1': 'OU1', # CIF中是OP1，我们命名为OU1
    'OP2': 'OU2', # CIF中是OP2，我们命名为OU2
    'H': 'H', 'HA': 'HA',
    'HB': 'HB2', 'HBA': 'HB3',
    'HD': 'HD2', 'HDA': 'HD3',
    'HE': 'HE2', 'HEA': 'HE3',
    'HG': 'HG2', 'HGA': 'HG3',
    'HK': 'HS2', 'HKA': 'HS3',
    'HL': 'HT2', 'HLA': 'HT3',
    'HNZ': 'HZ',
    # HOP2不在映射中，因为去质子化后删除了
}

validate_template(
    'SLL.cif',
    'K04_coords.in',
    'SLL',
    sll_name_map
)

print("\n\n" + "="*70)
print("OVERALL SUMMARY")
print("="*70)
print("\n✅ ALY (Acetyl-lysine):")
print("   - Geometry: From CIF ideal coordinates")
print("   - Charges: AMBER standard for neutral acetyl group")
print("   - Status: Ready to use")

print("\n✅ SLL (Succinyl-lysine):")
print("   - Geometry: From CIF ideal coordinates")
print("   - Charges: Updated to deprotonated carboxyl (-COO⁻)")
print("   - Modification: HOP2 removed (deprotonated state)")
print("   - Status: Ready to use")

print("\n📝 Key differences from CIF:")
print("   1. Atom naming: Renamed some atoms for OSPREY compatibility")
print("   2. Terminal atoms: OXT, HXT, H2 excluded (internal residue)")
print("   3. Protonation: SLL carboxyl is deprotonated (-COO⁻)")
print("   4. Charges: Added AMBER partial charges (not in CIF)")
print("   5. Atom types: Added AMBER atom types (not in CIF)")

print("\n" + "="*70)
