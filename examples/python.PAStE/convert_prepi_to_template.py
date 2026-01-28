#!/usr/bin/env python3
"""
从ANTECHAMBER的prepi文件转换为OSPREY的template.in格式
"""

def convert_prepi_to_template(prepi_file, template_file, residue_name):
    """
    转换prepi到template.in

    prepi格式（AMBER）：
        序号 原子名 类型 ... ... ... ... 键长 键角 二面角 电荷

    template.in格式（OSPREY）：
        序号 原子名 类型 tree ... ... 键长 键角 二面角 电荷
    """

    with open(prepi_file, 'r') as f:
        lines = f.readlines()

    # 解析prepi文件
    atom_lines = []
    in_atom_section = False
    header_lines = []
    footer_lines = []

    for i, line in enumerate(lines):
        # 检测头部
        if 'INT' in line and 'CORR' not in line:
            in_atom_section = True
            # 保存头部信息
            header_lines = lines[:i+2]  # 包括INT和CORR行
            continue

        # 检测尾部
        if in_atom_section and (line.strip().startswith('LOOP') or line.strip().startswith('IMPROPER')):
            footer_lines = lines[i:]
            break

        # 收集原子行
        if in_atom_section and line.strip():
            parts = line.split()
            if len(parts) >= 11 and parts[0].isdigit():
                atom_lines.append(line)

    # 写入template文件
    with open(template_file, 'w') as f:
        # 写入头部注释
        f.write("\n")
        f.write(f"Template generated from ANTECHAMBER prepi file\n")
        f.write(f"Based on AM1-BCC quantum chemistry calculations\n")

        # 写入残基名称和格式标记
        f.write(f"{residue_name}\n")
        f.write(f" {residue_name}  INT     1\n")
        f.write(" CORR OMIT DU   BEG\n")
        f.write("   0.00000\n")

        # 写入原子行（prepi格式已经是正确的）
        for line in atom_lines:
            f.write(line)

        # 写入尾部（LOOP和IMPROPER）
        for line in footer_lines:
            f.write(line)

    print(f"✓ Converted {prepi_file} → {template_file}")
    print(f"  Total atoms: {len(atom_lines)}")

import os
os.chdir('/home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE')

print("="*70)
print("Converting ANTECHAMBER prepi to OSPREY template.in")
print("="*70)

# 备份当前template文件
import shutil
for old, backup in [('K03_template.in', 'K03_template_manual.in'),
                     ('K04_template.in', 'K04_template_manual.in')]:
    if os.path.exists(old):
        shutil.copy(old, backup)
        print(f"\n✓ Backed up {old} → {backup}")

print("\n" + "-"*70)
print("Converting ALY...")
print("-"*70)
convert_prepi_to_template('ALY.prepi', 'K03_template_from_prepi.in', 'ALY')

print("\n" + "-"*70)
print("Converting SLL...")
print("-"*70)
convert_prepi_to_template('SLL.prepi', 'K04_template_from_prepi.in', 'SLL')

print("\n" + "="*70)
print("✅ Conversion complete!")
print("="*70)

print("\nGenerated files:")
print("  - K03_template_from_prepi.in (ALY from ANTECHAMBER)")
print("  - K04_template_from_prepi.in (SLL from ANTECHAMBER)")

print("\nBackup files:")
print("  - K03_template_manual.in (手工创建的)")
print("  - K04_template_manual.in (手工创建的)")

print("\nNext steps:")
print("  1. Compare template files:")
print("     diff K03_template.in K03_template_from_prepi.in")
print("  2. If satisfied, replace:")
print("     cp K03_template_from_prepi.in K03_template.in")
print("     cp K04_template_from_prepi.in K04_template.in")
print("  3. Test with OSPREY")
