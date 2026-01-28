#!/usr/bin/env python3
"""
比较ANTECHAMBER生成的电荷与我们当前模板的电荷
"""

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
                    # prepi格式：序号 原子名 类型 ... ... ... ... ... ... 电荷
                    atom_name = parts[1]
                    try:
                        charge = float(parts[10])
                        charges[atom_name] = charge
                    except (ValueError, IndexError):
                        pass
    return charges

def parse_coords_charges(coords_file):
    """从coords.in文件解析电荷"""
    charges = {}
    with open(coords_file, 'r') as f:
        for line in f:
            if line.startswith('ENDRES'):
                break
            parts = line.split()
            if len(parts) >= 5 and not line.startswith('ALY') and not line.startswith('SLL'):
                atom_name = parts[0]
                try:
                    charge = float(parts[4].rstrip('f'))
                    charges[atom_name] = charge
                except (ValueError, IndexError):
                    pass
    return charges

import os

# 确定工作目录
work_dir = '/home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE'
os.chdir(work_dir)

print("="*70)
print("ANTECHAMBER vs Current Template Charge Comparison")
print("="*70)
print(f"Working directory: {os.getcwd()}")

# 检查文件是否存在
print("\nChecking files:")
files_to_check = [
    'ALY.prepi',
    'SLL.prepi',
    'K03_coords.in',
    'K04_coords.in'
]

for f in files_to_check:
    exists = "✓" if os.path.exists(f) else "✗"
    print(f"  {exists} {f}")

print()

# ALY比较
print("="*70)
print("ALY (Acetyl-lysine)")
print("="*70)

try:
    antechamber_aly = parse_prepi_charges('ALY.prepi')
    current_aly = parse_coords_charges('K03_coords.in')

    print(f"\nANTECHAMBER: {len(antechamber_aly)} atoms")
    print(f"Current:     {len(current_aly)} atoms")

    print("\nKey atoms comparison:")
    key_atoms = ['NZ', 'CH', 'OH', 'CH3', 'HH31', 'HH32', 'HH33']

    print(f"\n{'Atom':<6} {'ANTECHAMBER':>12} {'Current':>12} {'Difference':>12}")
    print("-" * 50)

    for atom in key_atoms:
        ante_charge = antechamber_aly.get(atom, None)
        curr_charge = current_aly.get(atom, None)
        if ante_charge is not None and curr_charge is not None:
            diff = ante_charge - curr_charge
            print(f"{atom:<6} {ante_charge:>12.4f} {curr_charge:>12.4f} {diff:>+12.4f}")
        elif ante_charge is not None:
            print(f"{atom:<6} {ante_charge:>12.4f} {'N/A':>12} {'N/A':>12}")
        elif curr_charge is not None:
            print(f"{atom:<6} {'N/A':>12} {curr_charge:>12.4f} {'N/A':>12}")

    # 总电荷
    total_ante = sum(antechamber_aly.values())
    total_curr = sum(current_aly.values())
    print("-" * 50)
    print(f"{'TOTAL':<6} {total_ante:>12.4f} {total_curr:>12.4f} {total_ante-total_curr:>+12.4f}")

except FileNotFoundError as e:
    print(f"⚠️  File not found: {e}")
    print("   Please run ANTECHAMBER first!")

# SLL比较
print("\n" + "="*70)
print("SLL (Succinyl-lysine, deprotonated)")
print("="*70)

try:
    antechamber_sll = parse_prepi_charges('SLL.prepi')
    current_sll = parse_coords_charges('K04_coords.in')

    print(f"\nANTECHAMBER: {len(antechamber_sll)} atoms")
    print(f"Current:     {len(current_sll)} atoms")

    print("\nKey atoms comparison (carboxyl group):")
    key_atoms = ['NZ', 'CU', 'OU1', 'OU2']

    # SLL的ANTECHAMBER使用CIF命名，需要映射
    name_map = {
        'CP': 'CU',
        'OP1': 'OU1',
        'OP2': 'OU2',
    }

    print(f"\n{'Atom':<6} {'ANTECHAMBER':>12} {'Current':>12} {'Difference':>12}")
    print("-" * 50)

    for atom in key_atoms:
        # 检查是否需要映射
        ante_atom = atom
        for cif_name, our_name in name_map.items():
            if atom == our_name:
                ante_atom = cif_name
                break

        ante_charge = antechamber_sll.get(ante_atom, None)
        curr_charge = current_sll.get(atom, None)

        if ante_charge is not None and curr_charge is not None:
            diff = ante_charge - curr_charge
            display_name = f"{atom}"
            if ante_atom != atom:
                display_name += f" ({ante_atom})"
            print(f"{display_name:<6} {ante_charge:>12.4f} {curr_charge:>12.4f} {diff:>+12.4f}")
        elif ante_charge is not None:
            print(f"{atom:<6} {ante_charge:>12.4f} {'N/A':>12} {'N/A':>12}")
        elif curr_charge is not None:
            print(f"{atom:<6} {'N/A':>12} {curr_charge:>12.4f} {'N/A':>12}")

    # 总电荷
    total_ante = sum(antechamber_sll.values())
    total_curr = sum(current_sll.values())
    print("-" * 50)
    print(f"{'TOTAL':<6} {total_ante:>12.4f} {total_curr:>12.4f} {total_ante-total_curr:>+12.4f}")
    print(f"\n(Expected total: ALY≈0, SLL≈-1)")

except FileNotFoundError as e:
    print(f"⚠️  File not found: {e}")
    print("   Please run ANTECHAMBER first!")

print("\n" + "="*70)
print("Next steps:")
print("  1. Review the charge differences above")
print("  2. If ANTECHAMBER charges are significantly different,")
print("     update K03_coords.in and K04_coords.in with new charges")
print("  3. Re-run OSPREY calculation to see if ΔΔG improves")
print("="*70)
