import sys
import time

from Bio.PDB import PDBParser, PDBIO
from Find_Doublets import SCOPE
from KStarPrep import osprey_fileprep_kstar, ConfSpaceSpecs
import shutil
import glob
import subprocess

import math
import numpy as np
import os
import csv

from Bio.PDB import PDBParser
import warnings

warnings.simplefilter('ignore')


def euclidean_distance(coord1, coord2):

    x_diff = coord1[0] - coord2[0]
    y_diff = coord1[1] - coord2[1]
    z_diff = coord1[2] - coord2[2]

    return math.sqrt((x_diff ** 2) + (y_diff ** 2) + (z_diff ** 2))


def good_disulfide_vectors(res1, res2):

    res1_vector = res1['CB'].get_coord() - res1['CA'].get_coord()
    res2_vector = res2['CB'].get_coord() - res2['CA'].get_coord()

    dot_product = np.dot(res1_vector, res2_vector)

    if dot_product > 0:
        return True

    return False


def find_candidate_pairs(candidate_residues: list):

    candidate_pairs = []

    for res_idx1 in range(0, len(candidate_residues)-1):
        for res_idx2 in range(res_idx1+1, len(candidate_residues)):

            residue1 = candidate_residues[res_idx1]
            residue2 = candidate_residues[res_idx2]

            # check not neighbors by chain index
            if abs(residue1.get_id()[1] - residue2.get_id()[1]) <= 2:
                continue

            # CA distances range from 3.0A to 7.5A (https://www.nature.com/articles/s41598-020-67230-z)
            CA_distance = euclidean_distance(residue1['CA'].get_coord(), residue2['CA'].get_coord())
            if CA_distance < 3.0 or CA_distance > 7.5:
                continue

            # check residues face in same general direction (don't want disulfide formation to twist chain)
            if not good_disulfide_vectors(residue1, residue2):
                continue

            candidate_pairs.append([residue1.get_id()[1], residue2.get_id()[1]])

    return candidate_pairs


def find_cys_pairs(fixed_pdb: str, chain_id: str, ignored_positions: list):

    parser = PDBParser(PERMISSIVE=1)

    structure = parser.get_structure("complex", fixed_pdb)
    model = structure[0]
    chain = model[chain_id]

    candidate_residues = []

    for residue in chain:
        resname = residue.get_resname()
        resnum = residue.get_id()[1]

        # skip GLY and PRO -> mutation to CYS can twist backbone
        if resname == 'GLY' or resname == 'PRO' or resnum in ignored_positions:
            print("Ignoring residue %s %s for disulfide pair generation" % (resname, resnum))
            continue

        candidate_residues.append(residue)

    # determine all viable pairs based on geometry
    candidate_pairs = find_candidate_pairs(candidate_residues)

    print("Found %s possible disulfide bonds: %s" % (len(candidate_pairs), candidate_pairs))

    return candidate_pairs


def build_disulfides(input_directory: str, chain_id: str, ignore_posns: list, output_directory):

    # refactor the fxn good_disulfide_vector to account for diverse motifs
    print("WARNING: This script assumes a beta hairpin motif")

    os.mkdir(output_directory)
    print("Saving K* fileprep to %s" % output_directory)

    # find all pairs that have geometry for disulfide bond if mutated
    for complex_pdb in os.listdir(input_directory):

        print("Checking disulfide geometry for %s" % complex_pdb)
        entire_path = os.path.join(input_directory, complex_pdb)
        all_pairs = find_cys_pairs(entire_path, chain_id, ignore_posns)

        if len(all_pairs) == 0:
            print("WARNING: No candidates pairs exist for this structure. Moving to next structure.")
            continue

        # determine flexible sets for pairs
        # approximate CYS + vDw using nonpolar residues
        shutil.copy(entire_path, complex_pdb)
        doublets, interchains = SCOPE(complex_pdb, "", chain_id,
                                                        ["CYS", "LEU", "ILE", "PHE"], False,
                                                        "D", [])

        # for each disulfide pair, define flexibility
        disulfide_confspace = []
        for pair in all_pairs:
            total_flex = list(set(interchains[pair[0] - 1]) | set((interchains[pair[1] - 1])))
            newConfspace = ConfSpaceSpecs(pair, total_flex, [], [])
            disulfide_confspace.append(newConfspace)

        # make kstar storage folders
        match_storage = complex_pdb.split('-')[0] + '-kstar'
        os.mkdir(match_storage)

        # prep files for K*
        osprey_fileprep_kstar(complex_pdb, match_storage, ["CYS"], disulfide_confspace, chain_id,
                              "D", False, True, True,
                              False, "resources/K_bash.sh", [])

        # file cleanup
        old_PDB_template = "./" + match_storage.split('-')[0] + "*.pdb"
        for to_delete in glob.glob(old_PDB_template):
            os.remove(to_delete)
        shutil.move(match_storage, "%s/%s" % (output_directory, match_storage))

        print("\n\nCompleted fileprep for %s\n\n" % complex_pdb)

    print("\n\nFileprep for all candidate pairs for all matches completed!\n\n")


def cluster_runner(infolder: str):

    print("Submitting all disulfide candidates to cluster")

    runner = subprocess.run(["./cluster_runner.sh", infolder], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    errors = runner.stderr.decode('utf-8')
    if errors:
        print("ERROR! Unable to submit K* runs", file=sys.stderr)
        print(errors)
        exit()


def good_disulfide_bond(res1, res2):

    # C-S-S-C typical S-S bond length <= 2.05A
    res1_SG = res1['SG'].get_coord()
    res2_SG = res2['SG'].get_coord()
    SG_dist = euclidean_distance(res1_SG, res2_SG)
    if SG_dist <= 2.05:
        return True

    return False


def find_disulfides(kstar_directory, chain_id: str, get_best: bool, out_directory: str):

    parser = PDBParser(PERMISSIVE=1)

    match_directory = "%s/*-kstar" % kstar_directory

    if get_best:
        os.mkdir(out_directory)
        print("Saving best CYS bond for each match to %s" % out_directory)

    for match_path in glob.glob(match_directory):

        valid_disulfides = []

        print("\n\nNow checking %s" % match_path)

        pdb_names = "%s/kstar-*/ensembles/*pdb" % match_path

        for pdb in glob.glob(pdb_names):
            residue_names = pdb.split('/')[2].split('-')[1][1:-1].split('_')

            structure = parser.get_structure('complex', pdb)
            model = structure[0]
            chain = model[chain_id]

            # from PDB GMEC, get both CYS residue CA, CB
            residue_coords = []
            for residue in chain:
                resid = residue.get_id()[1]
                if resid == int(residue_names[0]) or resid == int(residue_names[1]):
                    residue_coords.append(residue)

            if len(residue_coords) != 2:
                print("ERROR! Something went wrong getting the CYS coordinates", file=sys.stderr)

            # check sulfur distances
            is_disulfide = good_disulfide_bond(residue_coords[0], residue_coords[1])

            if is_disulfide:
                valid_disulfides.append(pdb.split('/')[2])

        print("K* predicted valid disulfides for %s: %s" % (match_path.split('/')[1], valid_disulfides))

        # if only want bond for each match w/ best K* score
        if get_best:

            best_score = 0
            best_sequence = ""
            doublet_name = ""

            for screened in valid_disulfides:
                log_path = "%s/%s/submit.out" % (match_path, screened)

                with open(log_path, 'r', newline='\n') as logfile:
                    reader = csv.reader(logfile)
                    for row in reader:
                        if len(row) > 1 and row[0] == '1':
                            sequence = row[1]
                            score = float(row[2])
                            if score > best_score:
                                best_score = score
                                best_sequence = sequence
                                doublet_name = log_path.split('/')[2]

            print("Best bond for match %s is %s: %s [%s]" % (match_path.split('/')[1], doublet_name, best_score, best_sequence))

            pdb_name = "%s/%s/ensembles/seq.%s.pdb" % (match_path, doublet_name, best_sequence.replace(' ', '-'))
            target_location = "%s/%s.pdb" % (out_directory, match_path.split('/')[1])

            shutil.copy(pdb_name, target_location)


# find disulfide candidates + fileprep K*
build_disulfides("completed-matches", "B", [1, 2, 17, 18], "disulfide-builder-fileprep")

# run on cluster using cluster_runner.sh
# cluster_runner("disulfide-builder-fileprep")

# after running on cluster, analyze ensemble directory to find predicted disulfide bonds
# find_disulfides("disulfide-builder-fileprep", "B", True, "best-CYS")
