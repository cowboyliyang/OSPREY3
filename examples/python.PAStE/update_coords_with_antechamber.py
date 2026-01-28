#!/usr/bin/env python3
"""
用ANTECHAMBER计算的电荷更新coords.in文件
保持CIF的几何坐标，只更新电荷参数
"""

import os
import shutil

def parse_prepi_charges(prepi_file):
    """从prepi文件解析电荷"""
    charges = {}
    with open(prepi_file, 'r') as f:
        in_atom_section = False
        for line in f:
            if 'INT' in line and 'CORR' not in line:
                in_atom_section = True
                continue
            if in_atom_section:
                if line.strip().startswith('LOOP') or line.strip().startswith('IMPROPER'):
                    break
                parts = line.split()
                if len(parts) >= 11:
                    atom_name = parts[1]
                    try:
                        charge = float(parts[10])
                        charges[atom_name] = charge
                    except (ValueError, IndexError):
                        pass
    return charges

def update_coords_file(old_coords, new_coords, prepi_charges, name_map=None):
    """更新coords.in文件的电荷"""
    if name_map is None:
        name_map = {}

    with open(old_coords, 'r') as f:
        lines = f.readlines()

    # 备份原文件
    backup = old_coords.replace('.in', '_before_antechamber.in')
    shutil.copy(old_coords, backup)
    print(f"  Backed up to: {backup}")

    updated_lines = []
    updated_count = 0

    for line in lines:
        if line.startswith('ENDRES') or not line.strip():
            updated_lines.append(line)
            continue

        parts = line.split()
        if len(parts) >= 6:
            # 第一行是残基名和原子数
            if parts[0] in ['ALY', 'SLL'] and len(parts) == 2:
                updated_lines.append(line)
                continue

            atom_name = parts[0]
            x = parts[1]
            y = parts[2]
            z = parts[3]
            old_charge = parts[4]
            atom_type = parts[5]

            # 检查是否需要名称映射
            prepi_atom_name = name_map.get(atom_name, atom_name)

            # 如果ANTECHAMBER有这个原子的电荷，就更新
            if prepi_atom_name in prepi_charges:
                new_charge = prepi_charges[prepi_atom_name]
                new_line = f"{atom_name:4s}  {x:>7s}  {y:>7s}  {z:>7s}  {new_charge:7.4f}f  {atom_type}\n"
                updated_lines.append(new_line)
                updated_count += 1
            else:
                # 保持原电荷
                updated_lines.append(line)
        else:
            updated_lines.append(line)

    # 写入新文件
    with open(new_coords, 'w') as f:
        f.writelines(updated_lines)

    print(f"  Updated {updated_count} atom charges")
    print(f"  Saved to: {new_coords}")

print("="*70)
print("Updating coords.in files with ANTECHAMBER charges")
print("="*70)

os.chdir('/home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE')

# 更新ALY
print("\n1. Updating ALY (K03_coords.in)...")
aly_charges = parse_prepi_charges('ALY.prepi')
print(f"  Loaded {len(aly_charges)} charges from ALY.prepi")

update_coords_file('K03_coords.in', 'K03_coords_antechamber.in', aly_charges)

# 计算总电荷
total_aly = sum(aly_charges.values())
print(f"  Total charge: {total_aly:.4f} (expected: ~0)")

# 更新SLL - 需要名称映射
print("\n2. Updating SLL (K04_coords.in)...")
sll_charges = parse_prepi_charges('SLL.prepi')
print(f"  Loaded {len(sll_charges)} charges from SLL.prepi")

# SLL的名称映射：我们的名称 -> ANTECHAMBER的名称
sll_name_map = {
    'CQ': 'CX',   # 我们叫CQ，ANTECHAMBER叫CX
    'OQ1': 'OX',  # 我们叫OQ1，ANTECHAMBER叫OX
    'CS': 'CK',   # 我们叫CS，ANTECHAMBER叫CK
    'CT': 'CL',   # 我们叫CT，ANTECHAMBER叫CL
    'CU': 'CP',   # 我们叫CU，ANTECHAMBER叫CP
    'OU1': 'OP1', # 我们叫OU1，ANTECHAMBER叫OP1
    'OU2': 'OP2', # 我们叫OU2，ANTECHAMBER叫OP2
}

update_coords_file('K04_coords.in', 'K04_coords_antechamber.in', sll_charges, sll_name_map)

# 计算总电荷
total_sll = sum(sll_charges.values())
print(f"  Total charge: {total_sll:.4f} (expected: ~-1)")

print("\n" + "="*70)
print("✅ Updated coords files created!")
print("="*70)

print("\nGenerated files:")
print("  - K03_coords_antechamber.in (ALY with ANTECHAMBER charges)")
print("  - K04_coords_antechamber.in (SLL with ANTECHAMBER charges)")

print("\nBackup files:")
print("  - K03_coords_before_antechamber.in")
print("  - K04_coords_before_antechamber.in")

print("\nNext steps:")
print("  1. Review the new coords files")
print("  2. Replace the old ones:")
print("     cp K03_coords_antechamber.in K03_coords.in")
print("     cp K04_coords_antechamber.in K04_coords.in")
print("  3. Re-run OSPREY calculation")
print("  4. Compare new ΔΔG with old results")
