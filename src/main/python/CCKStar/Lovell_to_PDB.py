# this script serves to parse Lovell et al. v2 template coordinate libraries into PDB format

import re
from resources.lookup_dicts import *
from resources.RotamerConstructs import *
import os
import numpy as np
import math


# given an identity's rotamer, returns the backbone coords (CA, N, C) from lovell anchorCoords
def get_anchors(identity: str, rotamer: str):

    lovell_name = ClassToLovellRotamer[rotamer]
    lovell_start = ("[frag.%s.conf.%s]" % (identity, lovell_name))
    at_coords = False
    at_anchor = False

    ca_coords = []
    n_coords = []
    c_coords = []

    file = open("resources/lovell.conflib", "r")
    for line in file:
        if lovell_start in line:
            at_coords = True
        elif at_coords and "anchorCoords" in line:
            at_anchor = True
        elif at_coords and at_anchor:
            ca_anchor = line.split('=')[2].replace(" ", '')[1:-3]
            for item in ca_anchor.split(","):
                ca_coords.append(float(item))
            n_anchor = line.split('=')[3].replace(" ", '')[1:-3]
            for item in n_anchor.split(","):
                n_coords.append(float(item))
            if lovell_name in ["down", "up"]:
                c_anchor = line.split('=')[4].replace(" ", '')[1:-3]
            else:
                c_anchor = line.split('=')[4].replace(" ", '')[1:-11]
            for item in c_anchor.split(","):
                c_coords.append(float(item))
            break

    if not at_coords or not at_anchor:
        print("ERROR! Was unable to find anchors for %s %s" % (identity, rotamer))

    return ca_coords, n_coords, c_coords


def construct_lovell(identity: str, coords: list):

    rotname = coords[0].replace("-", "_")
    CA, N, C = get_anchors(identity, rotname)
    coords.append(CA)
    coords.append(N)
    coords.append(C)
    coords.pop(0)

    sconstruct = ("%s%s = %sRotamer(\'%s\',\n" % (identity, rotname, identity, rotname))

    counter = 1
    for xyz in coords:
        if counter % 2 == 0:
            sconstruct += ("%s,\n" % xyz)
        else:
            sconstruct += ("\t\t%s, " % xyz)
        counter += 1

    result = sconstruct[:-2]

    result += ")"

    print(result)

    return rotname

# parse a substring of lovell info (called by get_lovell)
def parse_lovell(ident: str, resinfo: str):

    at_coords = False
    first_name = True

    rot_info = []

    rot_names = []

    for line in resinfo.splitlines():
        start = re.search("^name.*", line)
        coordinates = re.search("^coords.*", line)
        end = re.search("^anchorCoords.*", line)
        if start:
            if first_name:
                first_name = False
                continue
            current_name = start.string[8:-1]
            at_coords = False
            rot_info.append(current_name)
        elif coordinates:
            at_coords = True
        elif end:
            at_coords = False
            rotn = construct_lovell(ident, rot_info)
            rot_names.append(rotn)
            rot_info = []
        elif at_coords:
            end_coord = re.search("^]", line)
            if end_coord:
                continue
            repls = {']': '', '}': '', " ": "", "#": ''}
            holder = line.split('[')[1]
            for k, v in repls.items():
                holder = holder.replace(k, v)
            xyz_str = holder.split(",")[0:3]
            xyz_float = []
            for element in xyz_str:
                xyz_float.append(float(element))
            rot_info.append(xyz_float)

    rot_name_list = ("%s_rotamers = [" % ident)
    for name in rot_names:
        rot_name_list += ("%s%s, " % (ident, name))
    rot_name_list = rot_name_list[:-2]
    rot_name_list += "]"
    print(rot_name_list)
    print('\n')

# given an AA identity, parse the atom info from the lovell.conflib library (ignores N-term conformations)
def get_lovell(identity: str):

    if identity == 'GLY' or identity == 'ALA':
        print("ERROR! Can't parse Gly or Ala; they only have 1 rotamer and aren't in the convex hull.")
        exit()

    # get all the conflib info
    info_start = "[frag." + identity + "]"

    file = open("resources/lovell.conflib", 'r')

    reached_info = False

    info = ""

    for line in file:
        if reached_info and line == '\n':
            break
        elif info_start in line:
            reached_info = True
        elif reached_info:
            info += line

    # parse to get coords for each rotamer
    parse_lovell(identity, info)


# create a lookup dictionary for rotamer atom names to general element names (e.g., HG11 = H)
# used to generate NameToElements in lookup_dicts
def make_element_dict():

    dict_holder = {}
    final_dict = {}

    file = open("resources/lovell.conflib", "r")

    at_atoms = False

    for line in file:
        astart = re.search("^atoms.=.*", line)
        aend = re.search("^]", line)
        if aend:
            at_atoms = False
        elif astart:
            at_atoms = True
        elif at_atoms:
            holder = line.split('\'')
            rot_name = holder[1]
            atom_name = holder[3]
            dict_holder[rot_name] = atom_name

    for key, value in dict_holder.items():
        if key not in final_dict.keys():
            final_dict[key] = value

    print(final_dict)

# generate rotamer classes for each amino acid identity
def make_acid_class(identity: str):

    atoms = ResAtoms[identity]

    cleanup = {']': '', "[": ''}
    for k, v in cleanup.items():
        atoms = atoms.replace(k, v)

    sclass = ("class %sRotamer: \n" % identity)
    sclass += ("\t def __init__(self, name, %s):\n" % atoms)
    sclass += "\t\tself.name = name\n"

    for item in atoms.split(','):
        cleaned = item.replace(" ", "")
        sclass += ("\t\tself.%s = %s\n" % (cleaned, cleaned))

    print(sclass)

# writes a PDB for all rotamers for the given identity
def print_ident_pdb(identity: str, rotamers: list):

    num_atoms = len(ResAtoms[identity].split(",")) - 1
    chain_tracker = 0
    resnum = 1
    filename = ("out_pdb/%s.pdb" % identity)
    if os.path.exists(filename):
        os.remove(filename)

    file = open(filename, "a")

    for rot in rotamers:
        curr_atom = resnum * num_atoms - (num_atoms - 1)
        rot.print_pdb(curr_atom, ChainAlphabet[chain_tracker], resnum, file)
        chain_tracker += 1
        resnum += 1

    file.write("END")
    file.close()

    print("saved %s" % filename)


# translate a rotamer to a specific xyz coordinate based on the alpha carbon
def translate_rotamer(rotamer, target: list):

    rot_atoms = [a for a in dir(rotamer) if not a.startswith('__') and a not in ["name", "print_pdb", "CA_back"]]

    CA_anchor = getattr(rotamer, "CA_back")

    x_diff = target[0] - CA_anchor[0]
    y_diff = target[1] - CA_anchor[1]
    z_diff = target[2] - CA_anchor[2]

    rotamer.CA_back = target
    for atom in rot_atoms:
        xyz = getattr(rotamer, atom)
        x_new = xyz[0] + x_diff
        y_new = xyz[1] + y_diff
        z_new = xyz[2] + z_diff
        rotamer.__setattr__(atom, [x_new, y_new, z_new])


# make a dictionary from no dash name (in class) -> lovell rotamer name
def make_lovell_name():

    name_dict = {}
    name_final = {}

    file = open("resources/lovell.conflib", "r")

    for line in file:
        name = re.search("^name.*", line)
        if name and line[8].islower():
            lovell_name = line.split("=")[1].strip().replace('\'', '')
            class_name = lovell_name.replace("-", '_')
            name_dict[class_name] = lovell_name

    for key, value in name_dict.items():
        if key not in name_final.values():
            name_final[key] = value

    print(name_final)


# find the Euclidean distance between two atoms in a rotamer
def euclidean_distance(rot, atom1, atom2):

    a1_xyz = np.array(getattr(rot, atom1))
    a2_xyz = np.array(getattr(rot, atom2))

    x_diff = (a1_xyz[0] - a2_xyz[0]) ** 2
    y_diff = (a1_xyz[1] - a2_xyz[1]) ** 2
    z_diff = (a1_xyz[2] - a2_xyz[2]) ** 2

    return math.sqrt(x_diff + y_diff + z_diff)


# calculates the length of a vector
def find_vec_length(vector: list):

    return math.sqrt((vector[0] ** 2) + (vector[1] ** 2) + (vector[2] ** 2))


# return the vector that points from root -> target (matrix subtraction)
def make_result_vector(root: list, target: list):

    return [(target[0] - root[0]), (target[1] - root[1]), (target[2] - root[2])]


# calculate scaled vector from vector root -> vector target
def make_norm_vector(root: list, target: list):

    new_vec = make_result_vector(root, target)
    vec_length = find_vec_length(new_vec)

    norm_vector = []

    for coord in new_vec:
        norm_coord = coord / vec_length
        norm_vector.append(norm_coord)

    return norm_vector


def calc_cross_product(vec1: list, vec2: list):

    first = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1])
    second = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2])
    third = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0])

    answer = [first, second, third]

    return answer


def calc_vec_angle(vec1: list, vec2: list):

    numerator = ((vec1[0] * vec2[0]) + (vec1[1] * vec2[1]) + (vec1[2] * vec2[2]))
    denominator = (find_vec_length(vec1)) * (find_vec_length(vec2))

    rad = math.acos(numerator / denominator)

    deg = math.degrees(rad)

    return deg


# move a rotamer to new CA, N, C anchor position
# we can find inverse, R, by solving R = V' * V^-1 where RV = V'
def rotate_rotamer(rotamer, target):

    # we ignore PRO
    if rotamer.name in ["down", "up"]:
        return

    # get the current backbone atoms (we use HA instead of CA so we can find inv matrix)
    old_N = rotamer.N_back
    old_C = rotamer.C_back
    old_HA = rotamer.HA

    # format the V matrix and find inverse
    V_matrix = np.array([old_N, old_C, old_HA])
    V_matrix_formatted = np.transpose(V_matrix)
    V_matrix_inv = np.linalg.inv(V_matrix_formatted)

    # format the V' matrix
    V_prime_matrix = np.array([target["N"], target["C"], target["HA"]])
    V_prime_formatted = np.transpose(V_prime_matrix)

    # calculate rotation matrix: target matrix x inverse original
    rotation_matrix = np.matmul(V_prime_formatted, V_matrix_inv)

    # calculate the new coords + update each class atom
    rot_atoms = [a for a in dir(rotamer) if not a.startswith('__') and a not in ["name", "print_pdb"]]
    for atom in rot_atoms:
        xyz = getattr(rotamer, atom)
        xyz_matrix = np.array(xyz)
        new_coords = np.matmul(rotation_matrix, xyz_matrix).tolist()
        rotamer.__setattr__(atom, new_coords)


# generate D-rotamers from L-rotamer library by reflection operation
def reflect_rotamer(rotamer):

    # we ignore PRO
    if rotamer.name in ["down", "up"]:
        return

    # for each atom, reflect over z axis
    rot_atoms = [a for a in dir(rotamer) if not a.startswith('__') and a not in ["name", "print_pdb"]]
    for atom in rot_atoms:
        xyz = getattr(rotamer, atom)
        reflected_coords = [xyz[0], xyz[1], (xyz[2] * -1)]
        rotamer.__setattr__(atom, reflected_coords)


# print all the xyz coordinates for each atom of a rotamer
# excludes distant H
def get_rot_coords(rotamer):

    coords = []

    rot_atoms = [a for a in dir(rotamer) if not a.startswith('__') and a not in ["name", "print_pdb", "H"]]
    for atom in rot_atoms:
        xyz = getattr(rotamer, atom)
        coords.append(xyz)

    return coords

# --------------------------------------------------------------------------------------------------------------
# # generate resources/IdentityClasses.py
# for i in AllIdentities:
#     make_acid_class(i)
#
# # make class constructors
# for i in AllIdentities:
#     get_lovell(i)

# # translate the rotamers to origin
# for iden in all_identity_rotamers:
#     for rot in iden:
#         translate_rotamer(rot, [0, 0, 0])
#
#
# # we'll arbitrarily use transposed VAL rotamer p as the target alignment
# target_rotation = {'CA': [0, 0, 0],
#                    'HA': [0.8059999999999974, -0.5640000000000001, -0.4930000000000003],
#                    'N': [0.605000000000004, 0.8560000000000016, 1.016],
#                    'C': [-0.7420000000000044, 0.8540000000000028, -1.0280000000000005]}
#
# # rotate the rotamers (assumes CA at origin)
# for iden in all_identity_rotamers:
#     for rot in iden:
#         rotate_rotamer(rot, target_rotation)
#
# # save the PDBs
# print_ident_pdb('VAL', VAL_rotamers)
# print_ident_pdb('CYS', CYS_rotamers)
# # print_ident_pdb('PRO', PRO_rotamers)
# print_ident_pdb('LEU', LEU_rotamers)
# print_ident_pdb('ILE', ILE_rotamers)
# print_ident_pdb('MET', MET_rotamers)
# print_ident_pdb('TRP', TRP_rotamers)
# print_ident_pdb('PHE', PHE_rotamers)
# print_ident_pdb('LYS', LYS_rotamers)
# print_ident_pdb('ARG', ARG_rotamers)
# print_ident_pdb('HID', HID_rotamers)
# print_ident_pdb('HIE', HIE_rotamers)
# print_ident_pdb('HIP', HIP_rotamers)
# print_ident_pdb('SER', SER_rotamers)
# print_ident_pdb('THR', THR_rotamers)
# print_ident_pdb('TYR', TYR_rotamers)
# print_ident_pdb('ASN', ASN_rotamers)
# print_ident_pdb('GLN', GLN_rotamers)
# print_ident_pdb('ASP', ASP_rotamers)
# print_ident_pdb('GLU', GLU_rotamers)
