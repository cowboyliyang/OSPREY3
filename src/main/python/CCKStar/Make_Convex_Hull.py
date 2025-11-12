import sys

from scipy.spatial import ConvexHull as ScipyConvexHull

from resources.RotamerConstructs import get_rotamers
from resources.lookup_dicts import VALp_L_origin, VALp_D_origin
from Lovell_to_PDB import translate_rotamer, rotate_rotamer, reflect_rotamer

def make_convex_hull(AA_types: list, outname, savePDB, chirality):

    all_rotamers = []

    # build AA rotamer list
    for a in AA_types:
        # ALA and GLY only have 1 rotamer, so ignore for CH
        if a in ['GLY', 'ALA']:
            continue
        amino = get_rotamers(a)
        all_rotamers.append(amino)

    # reflect L-lovell library to D space
    if chirality == 'D':
        for iden in all_rotamers:
            for rot in iden:
                reflect_rotamer(rot)

    # translate the rotamers to origin
    for iden in all_rotamers:
        for rot in iden:
            translate_rotamer(rot, [0, 0, 0])


    # we'll arbitrarily use transposed VAL rotamer p as the target alignment (at origin)
    if chirality == 'D':
        target_rotation = VALp_D_origin
    elif chirality == 'L':
        target_rotation = VALp_L_origin
    else:
        print("ERROR! Chirality not properly defined", file=sys.stderr)
        exit()

    # rotate the rotamers (CA must be at origin)
    for iden in all_rotamers:
        for rot in iden:
            rotate_rotamer(rot, target_rotation)

    # union list of all atoms in all rotamers
    coords_holder = []
    for iden in all_rotamers:
        for rot in iden:
            rot_atoms = [a for a in dir(rot) if not a.startswith('__') and a not in ["name", "print_pdb", "H", "HA",
                                                                                     "CA_back", "C_back", "N_back"]]
            coords = []
            for atom in rot_atoms:
                xyz = getattr(rot, atom)
                coords.append(xyz)
            for item in coords:
                if item not in coords_holder:
                    coords_holder.append(item)

    # create convex hull
    hull = ScipyConvexHull(coords_holder)
    hull_points = []

    for index in hull.vertices:
        coord = coords_holder[index]
        hull_points.append(coord)

    # format the data to visualize in PYMOL
    if savePDB:
        f = open(outname, "w")
        atom_num = 1
        for point in hull_points:
            pdb_info = [""] * 9
            pdb_info[0] = "ATOM".ljust(6)
            pdb_info[1] = str(atom_num).rjust(5)
            atom_name = ("C%s" % atom_num)
            pdb_info[2] = atom_name.center(4)
            pdb_info[3] = "VAL".ljust(3)
            pdb_info[4] = "A".rjust(1)
            pdb_info[5] = "1".rjust(4)
            pdb_info[6] = str('%8.3f' % (point[0])).rjust(8)
            pdb_info[7] = str('%8.3f' % (point[1])).rjust(8)
            pdb_info[8] = str('%8.3f' % (point[2])).rjust(8)
            f.write(("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                      pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                      pdb_info[8])))

            atom_num += 1
        f.close()
        print("Saved convex hull in PDB format to %s" % outname)

    return hull_points

# example = make_convex_hull(['VAL', 'CYS', 'LEU', 'ILE', 'MET', 'TRP', 'PHE', 'LYS', 'ARG', 'HID', 'HIE', 'HIP',
#                             'SER', 'THR', 'TYR', 'ASN', 'GLN', 'ASP', 'GLU'], "CH_points_test.pdb", True, "L")

