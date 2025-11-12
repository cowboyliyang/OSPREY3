from Find_Doublets import SCOPE, rank_flex_overlap, rank_design_overlap
from MONTAGE import run_MONTAGE
from ARISE import maintain_gly_pro, run_ARISE

# -------------------------- SCOPE --------------------------

intrachain_pairs, interchain_pairs = SCOPE('6ov7.pdb', 'pdb_hulls', 'C',
                                           ['VAL', 'CYS', 'LEU', 'ILE', 'MET', 'TRP',
                                            'PHE', 'LYS', 'ARG', 'HID', 'HIE', 'HIP', 'SER', 'THR', 'TYR',
                                            'ASN', 'GLN', 'ASP', 'GLU', 'ALA', 'GLY', 'PRO'], True,
                                           'L', [])


# optional: order the flexible residues by volume overlap with design chain hulls
# this is useful for prioritizing flexible residues over a large search space
# returns: {design res, target res : cubic angstrom overlap}
flex_order = rank_flex_overlap('C', 'A', intrachain_pairs, interchain_pairs, 'pdb_hulls')
print(flex_order)

# optional: order the intra-chain volume overlaps
# returns: {design res, design res : cubic angstrom overlap}
design_order = rank_design_overlap('C', intrachain_pairs, "pdb_hulls")
print(design_order)


# -------------------------- MONTAGE --------------------------

run_MONTAGE("L_to_D", "L", "D", 20, 4)

# -------------------------- ARISE --------------------------

# maintain GLY and PRO residues
gly_pro_doublets = maintain_gly_pro('B')
print("Maintaining the following GLY/PRO residues: %s" % gly_pro_doublets)

# using the MONTAGE GMEC PDBs, construct full sequences
run_ARISE(1, 12, set(), gly_pro_doublets, 'B', 'A',
          'completed-matches', 0.2, 'L')
