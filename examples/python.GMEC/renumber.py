# from Bio.PDB import PDBParser, PDBIO
# import sys
#
# def renumber_pdb(input_pdb, output_pdb):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("pdb_structure", input_pdb)
#
#     new_resnum = 1
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.id[0] != ' ':
#                     continue
#                 old_id = residue.id
#                 new_id = (' ', new_resnum, ' ')
#                 residue.id = new_id
#                 new_resnum += 1
#
#     io = PDBIO()
#     io.set_structure(structure)
#     io.save(output_pdb)
#
# if __name__ == "__main__":
#     if len(sys.argv) != 3:
#         print("python renumber_pdb.py input.pdb output.pdb")
#         sys.exit(1)
#
#     input_pdb = sys.argv[1]
#     output_pdb = sys.argv[2]
#     renumber_pdb(input_pdb, output_pdb)
import sys

def renumber_pdb(input_file, output_file):
    current_chain = None
    current_residue_key = None
    new_res_num = 1

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')):
                # 解析关键字段
                chain = line[21]
                res_seq = line[22:26].strip()  # 残基编号部分（4字符去空格）
                i_code = line[26].strip()       # 插入码（单独处理）

                # 构建残基唯一标识
                residue_key = (chain, res_seq, i_code)
                

                # 链变化检测
                if chain != current_chain:
                    current_chain = chain
                    current_residue_key = residue_key
                    new_res_num = 1
                else:
                    # 残基变化检测
                    if residue_key != current_residue_key:
                        new_res_num += 1
                        current_residue_key = residue_key

                # 生成新残基编号（4字符右对齐）
                new_res_seq = f"{new_res_num:4}"
                # 重建PDB行
                new_line = line[:22] + new_res_seq + line[26] + line[27:]
                f_out.write(new_line)
            else:
                # 非ATOM/HETATM行直接写入
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("使用方法: python pdb_renumber.py 输入文件.pdb 输出文件.pdb")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    renumber_pdb(input_pdb, output_pdb)
