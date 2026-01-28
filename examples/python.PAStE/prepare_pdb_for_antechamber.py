#!/usr/bin/env python3
"""
从CIF文件生成ANTECHAMBER需要的PDB文件
- ALY: 中性形式（-COOH, 总电荷0）
- SLL: 去质子化形式（-COO⁻, 总电荷-1，无HOP2）
"""

def parse_cif_atoms(cif_file):
    """从CIF文件解析原子"""
    atoms = []
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
            atoms.append({
                'name': atom_name,
                'element': element,
                'x': x_ideal,
                'y': y_ideal,
                'z': z_ideal
            })
    return atoms

def write_pdb(atoms, filename, residue_name, exclude_atoms=None):
    """写PDB文件"""
    if exclude_atoms is None:
        exclude_atoms = []

    with open(filename, 'w') as f:
        f.write(f"REMARK   Generated from CIF for ANTECHAMBER\n")
        f.write(f"REMARK   Residue: {residue_name}\n")

        atom_num = 1
        for atom in atoms:
            # 跳过末端原子和排除的原子
            # 注意：保留H2（N端氢），只删除OXT, HXT
            if atom['name'] in ['OXT', 'HXT'] + exclude_atoms:
                continue

            # PDB格式：
            # ATOM序号 原子名 残基名 链 残基号 x y z 占有率 B因子 元素
            line = f"ATOM  {atom_num:5d}  {atom['name']:<4s}{residue_name} A   1    "
            line += f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
            line += f"  1.00  0.00          {atom['element']:>2s}\n"
            f.write(line)
            atom_num += 1

        f.write("END\n")

# 生成ALY的PDB（中性）
print("Generating ALY.pdb (neutral acetyl-lysine, charge=0)...")
aly_atoms = parse_cif_atoms('ALY.cif')
write_pdb(aly_atoms, 'ALY.pdb', 'ALY')
aly_count = len([a for a in aly_atoms if a['name'] not in ['OXT', 'HXT']])
print(f"✓ Generated ALY.pdb with {aly_count} atoms (including H2)")

# 生成SLL的PDB（去质子化，不含HOP2，保留H2）
print("\nGenerating SLL.pdb (deprotonated succinyl-lysine, charge=-1)...")
sll_atoms = parse_cif_atoms('SLL.cif')
write_pdb(sll_atoms, 'SLL.pdb', 'SLL', exclude_atoms=['HOP2'])  # 去质子化：删除HOP2
sll_count = len([a for a in sll_atoms if a['name'] not in ['OXT', 'HXT', 'HOP2']])
print(f"✓ Generated SLL.pdb with {sll_count} atoms (including H2, excluding HOP2)")
print("  (HOP2 excluded for deprotonated -COO⁻ state)")

print("\n" + "="*60)
print("PDB files ready for ANTECHAMBER!")
print("="*60)
print("\nNext step:")
print("  ./run_antechamber.sh")
