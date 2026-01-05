from Bio.PDB import PDBParser, PDBIO


def remove_duplicate_residues_ignore_icode(input_pdb, output_pdb):
    parser = PDBParser()
    structure = parser.get_structure('processed_structure', input_pdb)

    for model in structure:
        for chain in model:
            seen_residues = set()
            duplicates = []

            last_res = (0, 0, 0)
            for residue in chain:
                # 定义唯一性：链ID + 残基序号（忽略插入码）
                residue_key = (residue.id[1],)  # 仅用残基序号
                if residue_key in seen_residues:
                    duplicates.append(last_res)
                else:
                    seen_residues.add(residue_key)
                last_res = residue.id

            for res_id in duplicates:
                chain.detach_child(res_id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

remove_duplicate_residues_ignore_icode("3bte_clean.pdb", "3bte_clean_no_dup.pdb")



# def remove_icode_keep_plain(input_pdb, output_pdb):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure('processed_structure', input_pdb)
#
#     for model in structure:
#         for chain in model:
#             resseq_map = {}
#             for residue in chain:
#                 resseq = residue.id[1]
#                 icode = residue.id[2]
#                 resseq_map.setdefault(resseq, []).append(residue.id)
#
#             for resseq, residue_ids in resseq_map.items():
#                 # 查找是否有无icode的
#                 plain_id = None
#                 for res_id in residue_ids:
#                     if res_id[2] == '':
#                         plain_id = res_id
#                         break
#                 # 如果有无icode的，保留它，删掉其余
#                 if plain_id:
#                     for res_id in residue_ids:
#                         if res_id != plain_id:
#                             chain.detach_child(res_id)
#                 else:
#                     # 如果全是带icode的，保留第一个，删剩下的
#                     for res_id in residue_ids[1:]:
#                         chain.detach_child(res_id)
#
#     io = PDBIO()
#     io.set_structure(structure)
#     io.save(output_pdb)
#
# remove_icode_keep_plain("3bte.pdb", "output_fixed.pdb")
