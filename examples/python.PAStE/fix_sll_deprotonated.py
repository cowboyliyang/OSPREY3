#!/usr/bin/env python3
"""
修正SLL模板为去质子化形式 (-COO⁻)
1. 删除HOP2氢原子
2. 修改CU、OU1、OU2的电荷为去质子化值
"""

# 去质子化羧酸根的AMBER电荷（参考ASP/GLU）
# 在-COO⁻中，两个氧原子是等价的（共振结构）
DEPROTONATED_CHARGES = {
    'CU': 0.8014,      # 羧基碳
    'OU1': -0.8188,    # 羧基氧1（共振）
    'OU2': -0.8188,    # 羧基氧2（共振）
}

def generate_deprotonated_template():
    """生成去质子化SLL的template.in"""
    # 读取当前的template
    with open('K04_template.in', 'r') as f:
        lines = f.readlines()

    # 删除HOP2相关的行，并调整后续原子编号
    new_lines = []
    atom_num_offset = 0

    for line in lines:
        # 跳过HOP2定义行
        if 'HOP2' in line or 'OU2 HOP2' in line:
            atom_num_offset = -1
            continue

        new_lines.append(line)

    # 写入新文件
    with open('K04_template_deprotonated.in', 'w') as f:
        f.writelines(new_lines)

    print("✓ Generated K04_template_deprotonated.in (removed HOP2)")

def generate_deprotonated_coords():
    """生成去质子化SLL的coords.in"""

    # SLL的AMBER电荷（更新羧基部分）
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
        # 去质子化羧基（-COO⁻）
        'CU': 0.8014, 'OU1': -0.8188, 'OU2': -0.8188,
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

    # 从当前coords文件读取坐标
    coords_dict = {}
    with open('K04_coords.in', 'r') as f:
        lines = f.readlines()
        for line in lines[1:-1]:  # 跳过第一行和ENDRES
            parts = line.split()
            if len(parts) >= 6:
                atom_name = parts[0]
                x = parts[1].rstrip('f')
                y = parts[2].rstrip('f')
                z = parts[3].rstrip('f')
                coords_dict[atom_name] = (float(x), float(y), float(z))

    # 原子顺序（删除HOP2，但保留OU2）
    sll_atom_order = ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3',
                      'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ',
                      'CQ', 'OQ1', 'CS', 'HS2', 'HS3', 'CT', 'HT2', 'HT3',
                      'CU', 'OU1', 'OU2', 'C', 'O']  # 注意：没有HOP2

    print("\nGenerating deprotonated SLL coords...")
    with open('K04_coords_deprotonated.in', 'w') as f:
        f.write(f"SLL {len(sll_atom_order)}\n")
        for atom in sll_atom_order:
            if atom in coords_dict:
                x, y, z = coords_dict[atom]
                charge = SLL_CHARGES[atom]
                atype = SLL_TYPES[atom]
                f.write(f"{atom:4s}  {x:6.3f}f  {y:6.3f}f  {z:6.3f}f  {charge:7.4f}f  {atype}\n")
        f.write("ENDRES\n")

    print("✓ Generated K04_coords_deprotonated.in")

    # 计算总电荷
    total_charge = sum(SLL_CHARGES[atom] for atom in sll_atom_order)
    print(f"\n📊 Total charge: {total_charge:.4f}")
    print(f"   (Expected: ~-1.0 for deprotonated carboxyl)")

    # 显示羧基电荷变化
    print("\n🔧 Carboxyl group charge changes:")
    print("   Before (protonated -COOH):")
    print("     CU:  +0.7341")
    print("     OU1: -0.5894")
    print("     OU2: -0.5894")
    print("     HOP2: (exists)")
    print("     Subtotal: -0.4447")
    print("\n   After (deprotonated -COO⁻):")
    print(f"     CU:  {SLL_CHARGES['CU']:+.4f}")
    print(f"     OU1: {SLL_CHARGES['OU1']:+.4f}")
    print(f"     OU2: {SLL_CHARGES['OU2']:+.4f}")
    print("     HOP2: (deleted)")
    print(f"     Subtotal: {SLL_CHARGES['CU'] + SLL_CHARGES['OU1'] + SLL_CHARGES['OU2']:+.4f}")

# 运行生成
print("=" * 60)
print("Fixing SLL template for deprotonated state (-COO⁻)")
print("=" * 60)

generate_deprotonated_template()
generate_deprotonated_coords()

print("\n" + "=" * 60)
print("✅ COMPLETED!")
print("=" * 60)
print("\nGenerated files:")
print("  - K04_template_deprotonated.in (removed HOP2)")
print("  - K04_coords_deprotonated.in (updated charges, no HOP2)")
print("\nNext steps:")
print("  1. Review the charge changes above")
print("  2. Replace K04_template.in and K04_coords.in with these files")
print("  3. Re-run OSPREY calculation")
print("  4. Compare new ΔΔG with +25.61 kcal/mol")
