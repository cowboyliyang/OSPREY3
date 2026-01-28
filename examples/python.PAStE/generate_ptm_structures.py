#!/usr/bin/env python3
"""
生成B781位置替换为ALY和SLL的PDB结构
保持骨架原子（N, CA, C, O），替换侧链
"""

def parse_coords_file(coords_file):
    """从coords.in文件读取PTM的坐标"""
    atoms = {}
    with open(coords_file) as f:
        for line in f:
            if line.startswith('ENDRES') or not line.strip():
                break
            parts = line.split()
            if len(parts) >= 6 and not line.startswith('ALY') and not line.startswith('SLL'):
                atom_name = parts[0]
                x = float(parts[1].rstrip('f'))
                y = float(parts[2].rstrip('f'))
                z = float(parts[3].rstrip('f'))
                atoms[atom_name] = (x, y, z)
    return atoms

def get_backbone_atoms(pdb_file, chain, resnum):
    """从PDB文件提取指定残基的骨架原子坐标"""
    backbone = {}
    with open(pdb_file) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if f' {chain} {resnum:>3}' in line or f' {chain}{resnum:>4}' in line:
                atom_name = line[12:16].strip()
                if atom_name in ['N', 'CA', 'C', 'O', 'H', 'HA']:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    backbone[atom_name] = (x, y, z, line)
    return backbone

def calculate_transformation(ref_coords, ptm_coords):
    """
    计算从PTM模板坐标到实际蛋白坐标的转换
    基于N, CA, C三个骨架原子
    """
    import numpy as np

    # 提取骨架原子坐标
    ref_points = []
    ptm_points = []
    for atom in ['N', 'CA', 'C']:
        if atom in ref_coords and atom in ptm_coords:
            ref_points.append(ref_coords[atom][:3])
            ptm_points.append(ptm_coords[atom])

    if len(ref_points) < 3:
        print(f"Warning: Not enough backbone atoms for alignment")
        return None, None

    ref_points = np.array(ref_points)
    ptm_points = np.array(ptm_points)

    # 计算质心
    ref_center = ref_points.mean(axis=0)
    ptm_center = ptm_points.mean(axis=0)

    # 中心化
    ref_centered = ref_points - ref_center
    ptm_centered = ptm_points - ptm_center

    # 计算旋转矩阵 (Kabsch算法)
    H = ptm_centered.T @ ref_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # 确保右手坐标系
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # 平移向量
    t = ref_center - R @ ptm_center

    return R, t

def transform_point(point, R, t):
    """应用旋转和平移到一个点"""
    import numpy as np
    return R @ np.array(point) + t

def write_ptm_pdb(input_pdb, output_pdb, chain, resnum, ptm_name, ptm_coords, backbone):
    """
    生成替换了PTM的新PDB文件
    """
    import numpy as np

    # 计算坐标转换
    R, t = calculate_transformation(backbone, ptm_coords)

    if R is None:
        print("Error: Cannot calculate transformation")
        return

    # 读取原始PDB，准备写入新PDB
    with open(input_pdb) as f:
        pdb_lines = f.readlines()

    with open(output_pdb, 'w') as out:
        atom_num = 0
        skip_residue = False
        wrote_ptm = False

        for line in pdb_lines:
            # 记录原子编号
            if line.startswith('ATOM'):
                atom_num = int(line[6:11])

            # 检查是否是要替换的残基
            if line.startswith('ATOM') and f' {chain} {resnum:>3}' in line:
                # 只在第一次遇到时写入PTM
                if not wrote_ptm:
                    # 写入PTM的所有原子
                    atom_count = 0
                    for atom_name in ['N', 'H', 'CA', 'HA', 'CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3',
                                     'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ']:
                        if atom_name in ptm_coords:
                            atom_count += 1
                            coord = ptm_coords[atom_name]

                            # 如果是骨架原子，使用原始坐标
                            if atom_name in backbone:
                                x, y, z = backbone[atom_name][:3]
                            else:
                                # 侧链原子需要转换
                                transformed = transform_point(coord, R, t)
                                x, y, z = transformed

                            # 确定元素
                            element = atom_name[0]
                            if element == 'H':
                                element = 'H'
                            elif atom_name.startswith('HB') or atom_name.startswith('HG') or \
                                 atom_name.startswith('HD') or atom_name.startswith('HE') or \
                                 atom_name.startswith('HZ') or atom_name.startswith('HH'):
                                element = 'H'

                            # 写入ATOM行
                            out.write(f"ATOM  {atom_num:>5}  {atom_name:<4}{ptm_name} {chain}{resnum:>4}    ")
                            out.write(f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")
                            atom_num += 1

                    # 写入PTM特异原子
                    if ptm_name == 'ALY':
                        # 乙酰基: CH, OH, CH3, HH31-33
                        for atom_name in ['CH', 'OH', 'CH3', 'HH31', 'HH32', 'HH33']:
                            if atom_name in ptm_coords:
                                coord = ptm_coords[atom_name]
                                transformed = transform_point(coord, R, t)
                                x, y, z = transformed

                                if atom_name == 'CH':
                                    element = 'C'
                                elif atom_name == 'OH':
                                    element = 'O'
                                elif atom_name == 'CH3':
                                    element = 'C'
                                else:  # HH31-33
                                    element = 'H'

                                out.write(f"ATOM  {atom_num:>5}  {atom_name:<4}{ptm_name} {chain}{resnum:>4}    ")
                                out.write(f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")
                                atom_num += 1

                    elif ptm_name == 'SLL':
                        # 琥珀酰基: CQ, HQ2-3, CS, HS2-3, CT, HT2-3, CU, OU1, OU2
                        for atom_name in ['CQ', 'HQ2', 'HQ3', 'CS', 'HS2', 'HS3',
                                        'CT', 'HT2', 'HT3', 'CU', 'OU1', 'OU2']:
                            if atom_name in ptm_coords:
                                coord = ptm_coords[atom_name]
                                transformed = transform_point(coord, R, t)
                                x, y, z = transformed

                                if atom_name[0] == 'C':
                                    element = 'C'
                                elif atom_name[0] == 'O':
                                    element = 'O'
                                else:  # H atoms
                                    element = 'H'

                                out.write(f"ATOM  {atom_num:>5}  {atom_name:<4}{ptm_name} {chain}{resnum:>4}    ")
                                out.write(f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")
                                atom_num += 1

                    # 写入骨架C和O
                    for atom_name in ['C', 'O']:
                        if atom_name in backbone:
                            x, y, z = backbone[atom_name][:3]
                            element = atom_name
                            out.write(f"ATOM  {atom_num:>5}  {atom_name:<4}{ptm_name} {chain}{resnum:>4}    ")
                            out.write(f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")
                            atom_num += 1

                    wrote_ptm = True
                    skip_residue = True
                continue  # 跳过原始残基的这一行

            # 检查是否结束了要替换的残基
            if skip_residue and line.startswith('ATOM'):
                current_res = line[21:26].strip()
                target_res = f"{chain} {resnum}"
                if current_res != target_res:
                    skip_residue = False

            # 如果不是要跳过的行，写入原始内容
            if not skip_residue:
                # 更新ATOM编号
                if line.startswith('ATOM'):
                    atom_num += 1
                    out.write(f"ATOM  {atom_num:>5}{line[11:]}")
                else:
                    out.write(line)

# 主程序
print("="*80)
print("生成B781位置的PTM结构")
print("="*80)

import os
os.chdir('/home/users/lz280/IdeaProjects/OSPREY3/examples/python.PAStE')

# 读取PTM坐标
print("\n1. 读取PTM模板坐标...")
aly_coords = parse_coords_file('K03_coords.in')
sll_coords = parse_coords_file('K04_coords.in')
print(f"   ALY: {len(aly_coords)} 原子")
print(f"   SLL: {len(sll_coords)} 原子")

# 读取B781的骨架原子
print("\n2. 读取B781位置的骨架原子...")
backbone = get_backbone_atoms('6dv2_strip_reduce_prep_rc_add_rc.pdb', 'B', 781)
print(f"   找到 {len(backbone)} 个骨架原子")

# 生成ALY结构
print("\n3. 生成ALY结构...")
write_ptm_pdb(
    '6dv2_strip_reduce_prep_rc_add_rc.pdb',
    '6dv2_B781_ALY.pdb',
    'B', 781, 'ALY',
    aly_coords,
    backbone
)
print("   ✓ 生成: 6dv2_B781_ALY.pdb")

# 生成SLL结构
print("\n4. 生成SLL结构...")
write_ptm_pdb(
    '6dv2_strip_reduce_prep_rc_add_rc.pdb',
    '6dv2_B781_SLL.pdb',
    'B', 781, 'SLL',
    sll_coords,
    backbone
)
print("   ✓ 生成: 6dv2_B781_SLL.pdb")

print("\n" + "="*80)
print("✅ 完成！")
print("="*80)
print("\n生成的文件:")
print("  - 6dv2_B781_ALY.pdb (B781位置替换为acetyl-lysine)")
print("  - 6dv2_B781_SLL.pdb (B781位置替换为succinyl-lysine)")
print("\n注意:")
print("  侧链坐标通过叠合N-CA-C骨架原子从模板转换而来")
print("  可能需要在PyMOL或Chimera中检查几何合理性")
