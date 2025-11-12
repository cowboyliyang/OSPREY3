import sys
import shutil
import glob
import os
import subprocess
import time
import random

from Bio.PDB import PDBParser

from KStarPrep import ConfSpaceSpecs, doublet_confspace_info, reduce_doublet, osprey_fileprep_kstar, save_doublets

def find_best_log(logfile, apo_tolerance):

    best_score = 0
    best_sequence = ""
    apo_lig_WT = 0

    on_header = False
    got_WT_apo_lig = False

    f = open(logfile, 'r')
    for line in f:
        if "error" in line or "Error" in line or "Java" in line or "java" in line:
            print("ERROR! Found logfile errors in %s" % logfile, file=sys.stderr)
            exit()
        if "seconds" in line:
            break
        if on_header:
            fields = line.split(',')

            # get the mutations
            mut_field = fields[1].replace(' ', '-')

            # skip if mutant destabilizes apo ligand beyond tolerance
            if got_WT_apo_lig:
                pfunc_mut = fields[9]
                pfunc_diff = apo_lig_WT - float(pfunc_mut)
                if pfunc_diff > apo_tolerance:
                    continue

            # get upper bound of WT apo ligand pfunc
            if not got_WT_apo_lig:
                pfunc_field = fields[10]
                if pfunc_field == 'none' or '-' in pfunc_field:
                    print("ERROR! No WT pfunc for %s" % logfile)
                    exit()
                apo_lig_WT = float(pfunc_field)
                got_WT_apo_lig = True

            # get the kstar score
            score = fields[2]
            if score == 'none' or '-' in score:
                continue
            score = float(score)
            if score > best_score:
                best_score = score
                best_sequence = mut_field

        if "Assignments" in line:
            on_header = True

    if best_sequence == "":
        print("ERROR! Unable to find best mutant for %s" % logfile, file=sys.stderr)
        exit()

    return best_sequence, best_score


def find_best_scans(infolder, outfolder, visited_dict, apo_tolerance):

    print("Storing the GMEC PDB from each match into %s" % outfolder)

    for matchfolder in os.listdir(infolder):
        log_path = os.path.join(infolder, matchfolder, 'kstar-*', 'submit.out')

        best_doublet = ""
        doublet_score = 0
        best_seq = ""

        for logfile in glob.glob(log_path):
            best_sequence, high_score = find_best_log(logfile, apo_tolerance)
            if high_score > doublet_score:
                best_doublet = logfile
                doublet_score = high_score
                best_seq = best_sequence

        if doublet_score == 0:
            print("ERROR! No doublet runs had a positive K* score for %s" % matchfolder, file=sys.stderr)
            print("If adequate flexibility was previously set, then this may indicate the need to prune this backbone")
            exit()
        print("Best doublet for %s was %s with a score of %s" % (matchfolder,
                                                                 best_doublet.split('/')[2], doublet_score))

        best_pdb_inname = 'seq.' + best_seq + '.pdb'
        best_pdb_outname = matchfolder.split('-')[0] + '-' + best_doublet.split('/')[2].split('-')[1] + '.pdb'

        pdb_inpath = os.path.join(infolder, matchfolder, best_doublet.split('/')[2], 'ensembles', best_pdb_inname)
        pdb_outpath = os.path.join(outfolder, best_pdb_outname)

        shutil.copy(pdb_inpath, pdb_outpath)

        # add this new doublet to our dictionary
        match_num = matchfolder.split('-')[0][5:]
        doublet_str = best_doublet.split('/')[2].split('-')[1].replace('[', '').replace(']', '')
        doublet_str_list = doublet_str.split('_')
        doublet_numbers = [int(x) for x in doublet_str_list]
        if match_num in visited_dict:
            old_doublets = visited_dict[match_num]
            old_doublets.append(doublet_numbers)
            visited_dict[match_num] = old_doublets
        else:
            visited_dict[match_num] = [doublet_numbers]

    print("Updated visited dictionary:")
    print(visited_dict)
    return visited_dict


def is_match_fully_assigned(visited_nodes, completed_matches, chain_length):

    for match, visited in visited_nodes.items():
        fixed_residues = set()
        for d in visited:
            for res in d:
                fixed_residues.add(res)
        if len(fixed_residues) == chain_length:
            print("Match %s is complete!" % match)
            completed_matches.add(match)
        else:
            print("Match %s is unfinished. The following res are fixed: %s" % (match, fixed_residues))

    return completed_matches


def move_finished_matches(completed_matches, visited_doublets, new_scan_PDB_folder, round_number,
                          final_designs_outfolder):

    moved_matches = []

    # figure out which matches we've already moved
    for finished_match in os.listdir(final_designs_outfolder):
        matchnum = finished_match.split('-')[0][5:]
        moved_matches.append(matchnum)

    for match_number in completed_matches:
        if match_number in moved_matches:
            continue

        # get the last doublet to find the PDB file
        last_doublet = visited_doublets[match_number][-1]
        formatted_doublet = "[" + str(last_doublet[0]) + '_' + str(last_doublet[1]) + '].pdb'
        old_filename = "match" + str(match_number) + '-' + formatted_doublet
        old_path = os.path.join(new_scan_PDB_folder,  old_filename)

        # specify new path to completed outfolder and store PDB with scan round it ended on
        new_filename = old_filename.split('.')[0] + "-Scan" + str(round_number) + '.pdb'
        new_path = os.path.join(final_designs_outfolder, new_filename)

        # copy over the PDB and updated log
        shutil.copyfile(old_path, new_path)
        print("Saved finished match %s to %s" % (match_number, new_path))
        print("This match will be skipped in future scans")


def setup_scan_round(round_number, previous_visits, apo_tolerance):

    new_scan_PDBs = "ARISE-Scan%s-PDB" % round_number
    os.mkdir(new_scan_PDBs)
    out_kstar_files = "ARISE-Scan%s-kstar" % round_number
    out_hull_files = "ARISE-Scan%s-hulls" % round_number
    os.mkdir(out_kstar_files)
    os.mkdir(out_hull_files)
    print("\nSaving kstar prep files to %s, hulls to %s" % (out_kstar_files, out_hull_files))

    if round_number == 1:
        for gmec in os.listdir("MONTAGE_GMEC"):
            full_pdb_path = os.path.join("MONTAGE_GMEC", gmec)
            shutil.copy(full_pdb_path, new_scan_PDBs)
        return new_scan_PDBs, previous_visits, out_kstar_files, out_hull_files

    else:
        old_scan_outfolder = "ARISE-Scan%s-kstar" % (round_number - 1)
        visited_doublets = find_best_scans(old_scan_outfolder, new_scan_PDBs, previous_visits, apo_tolerance)
        return new_scan_PDBs, visited_doublets, out_kstar_files, out_hull_files


def choose_next_doublet(all_visited_doublets, all_specs, fixed_residues, length_design):

    # first round: choose any doublet if none yet visited
    if len(fixed_residues) == 0:
        for spec in all_specs:
            if len(spec.doublet) == 2:
                return spec

    # general case: choose a doublet that involves a fixed residue with previous scan
    previous_doublet = all_visited_doublets[-1]
    print("Previous doublet for this match was %s" % previous_doublet)
    for spec in all_specs:
        shared_res = set(previous_doublet) & set(spec.doublet)
        if len(shared_res) == 1 and len(spec.doublet) == 2:
            if spec.doublet[0] in fixed_residues and spec.doublet[1] in fixed_residues:
                continue
            else:
                print("Selected doublet %s for next run" % spec.doublet)
                return spec

    # alternative case: choose a doublet that shares fixed res with any previous run
    print("No next doublet could be selected based on shared res with the previous doublet")
    for spec in all_specs:
        shared_res = fixed_residues & set(spec.doublet)
        if len(shared_res) == 1 and len(spec.doublet) == 2:
            print("Selected doublet %s for next run" % spec.doublet)
            return spec

    # rare case: no doublets have a fixed res, so choose a doublet with 0 fixed residues
    print("No next doublet could be selected based on shared res with any previous doublets")
    for spec in all_specs:
        shared_res = fixed_residues & set(spec.doublet)
        if len(shared_res) == 0 and len(spec.doublet) == 2:
            if spec.doublet[0] not in fixed_residues and spec.doublet[1] not in fixed_residues:
                print("Selected doublet %s for next run" % spec.doublet)
                return spec

    # special case: unfixed residues are not in any doublets (only island(s) remain)
    ordered_fixed_residues = sorted(list(fixed_residues))
    unfixed_residues = [x for x in range(1, length_design+1) if x not in fixed_residues]

    if len(unfixed_residues) >= 1:

        print("Only islands remain, so pairing with random fixed residue.")
        selected_island = unfixed_residues[0]
        other_residue = random.choice(ordered_fixed_residues)

        print("Prepping island residue %s with fixed residue %s" % (selected_island, other_residue))
        island_flexset = []
        for spec in all_specs:
            if len(spec.doublet) == 1 and spec.doublet[0] == selected_island:
                island_flexset = spec.flexset
                print("Island residue %s intersects with target residues %s" % (selected_island, island_flexset))

        new_doublet = ConfSpaceSpecs([selected_island, other_residue], island_flexset,
                                     [], [])
        return new_doublet

    # catch if some weird geometry slips past us
    print("ERROR! Was unable to select the next doublet for this match", file=sys.stderr)
    exit()


def fileprep_scan_round(new_scan_PDB, visited_doublets, design_ID, target_ID, out_kstar_files, out_hull_files,
                        finished_matches, chain_length, design_chirality):

    for file in os.listdir(new_scan_PDB):

        curr_matchnum = file.split('-')[0][5:]
        if curr_matchnum in finished_matches:
            print("Match %s is fully assigned, so not prepping" % curr_matchnum)
            continue

        print("\n--Now prepping K* files using %s--\n" % file)

        # set up the output directories for hulls and OSPREY kstar files
        full_filepath = os.path.join(new_scan_PDB, file)
        shutil.copy(full_filepath, file)
        hull_foldername = ("%s-hulls" % (file.split("-")[0]))
        kstar_foldername = ("%s-kstar" % (file.split("-")[0]))
        os.mkdir(hull_foldername)
        os.mkdir(kstar_foldername)

        # find which residues identities on the design need to be fixed
        all_visited_doublets = visited_doublets[curr_matchnum]
        fixed_residues = set()
        for d in all_visited_doublets:
            for res in d:
                fixed_residues.add(res)
        print("Fixed residues are %s" % fixed_residues)
        print("Remaining residues are %s" % [x for x in range(1, chain_length+1) if x not in fixed_residues])

        design_chain_muts = ['VAL', 'CYS', 'LEU', 'ILE', 'MET', 'TRP',
                             'PHE', 'LYS', 'ARG', 'HID', 'HIE', 'HIP', 'SER', 'THR', 'TYR',
                             'ASN', 'GLN', 'ASP', 'GLU', 'ALA', 'GLY', 'PRO']

        # find inter and intra-chain contacts
        all_specs = doublet_confspace_info(file, design_ID, design_chain_muts, design_chirality,
                                           hull_foldername, list(fixed_residues))

        # choose the next doublet based on geometry of previous doublet run
        selected_spec = choose_next_doublet(all_visited_doublets, all_specs, fixed_residues, chain_length)

        # reduce number of flexible target residues
        reduced_spec = reduce_doublet(design_ID, target_ID, selected_spec, 2, hull_foldername)

        # prepare compiled conformation space files for Kstar
        osprey_fileprep_kstar(file, kstar_foldername, design_chain_muts, reduced_spec, design_ID, design_chirality,
                              False, True, False, False,
                              "resources/K_bash.sh", list(fixed_residues))

        # save the confspecs
        doublet_filename = ("%s-doublet-info.txt" % (file.split("-")[0]))
        save_doublets(reduced_spec, doublet_filename)
        os.rename(doublet_filename, ("%s/%s" % (hull_foldername, doublet_filename)))
        print("Moved doublet info to %s" % hull_foldername)

        # organize files into directory tree and delete intermediate PDBs
        old_PDB_template = "./" + file.split('-')[0] + "*.pdb"
        for to_delete in glob.glob(old_PDB_template):
            os.remove(to_delete)
        shutil.move(hull_foldername, out_hull_files)
        print("Moved %s to %s" % (hull_foldername, out_hull_files))
        shutil.move(kstar_foldername, out_kstar_files)
        print("Moved %s to %s" % (kstar_foldername, out_kstar_files))


def start_new_run(new_scan_outfolder, round_number):

    print("\nSubmitting K* runs for Round %s" % round_number)

    # we have to use bash to submit via slurm
    runner = subprocess.run(["./cluster_runner.sh", new_scan_outfolder], stdout=subprocess.DEVNULL,
                            stderr=subprocess.PIPE)
    errors = runner.stderr.decode('utf-8')
    if errors:
        print("ERROR! Unable to submit K* runs", file=sys.stderr)
        print(errors)
        exit()

    print("All jobs submitted")


def are_kstar_finished(kstar_folder):

    match_dirs = os.path.join(kstar_folder, "match*-kstar")
    total_matches = 0
    for match in glob.glob(match_dirs):
        total_matches += 1

    print("Checking the %s matches in %s for completion..." % (total_matches, kstar_folder))

    log_files = os.path.join(match_dirs, 'kstar-*', 'submit.out')
    finished_matches = 0
    unfinished_matches = []
    for log in glob.glob(log_files):
        is_running = True
        for line in open(log, 'r'):
            if 'completed' in line:
                finished_matches += 1
                is_running = False
                break
            if 'error' in line or 'Error' in line or 'java' in line:
                print("ERROR! Something broke during the K* search for %s" % log, file=sys.stderr)
                exit()
        if is_running:
            unfinished_matches.append(log)

    if finished_matches != total_matches:
        print("%s out of %s matches are not finished. %s matches are running: %s" % ((total_matches-finished_matches),
                                                                                     total_matches,
                                                                                     len(unfinished_matches),
                                                                                     unfinished_matches))
        return False

    print("All %s matches completed for this round!" % total_matches)
    return True


def maintain_gly_pro(design_ID):

    parser = PDBParser(PERMISSIVE=1)

    gly_pro_doubls = {}

    for gmec in os.listdir("MONTAGE_GMEC"):
        full_pdb_path = os.path.join("MONTAGE_GMEC", gmec)

        curr_gly_pro = []

        structure = parser.get_structure("complex", full_pdb_path)
        model = structure[0]

        for chain in model:
            if chain.id == design_ID:
                for residue in chain:
                    resname = residue.resname
                    resid = residue.id[1]
                    if resname == 'GLY' or resname == 'PRO':
                        curr_gly_pro.append([resid, resid])

        match_num = gmec.split('-')[0].split('match')[1]
        gly_pro_doubls[match_num] = curr_gly_pro

    return gly_pro_doubls


def run_ARISE(round_number: int, length_chain: int, finished_matches: set, visited_doubs: dict, designID: str,
              targetID: str, final_designs_outfolder: str, apo_tolerance: float, design_chirality: str):

    # setup outfolder if it doesn't already exist
    if not os.path.isdir(final_designs_outfolder):
        print("Created directory %s to store completed matches" % final_designs_outfolder)
        os.mkdir(final_designs_outfolder)

    print("\n\n\n\n\n\n--Now starting ARISE round %s--\n\n\n\n\n\n" % round_number)

    # setup GMEC PDB directory and storage directories
    new_scan_PDB_folder, visited_doublets, out_kstar_files, out_hull_files = setup_scan_round(round_number,
                                                                                              visited_doubs,
                                                                                              apo_tolerance)

    # check if any matches are done, and if so update our tracker
    completed_matches = is_match_fully_assigned(visited_doublets, finished_matches, length_chain)
    print("%s of %s matches are fully assigned: %s" % (len(completed_matches), len(visited_doublets),
                                                       completed_matches))

    # if done, rename + move to completed outdirectory
    move_finished_matches(completed_matches, visited_doublets, new_scan_PDB_folder, round_number,
                          final_designs_outfolder)

    # exit if all fully assigned
    if len(completed_matches) == len(visited_doublets) and len(visited_doublets) != 0:
        print("\nSuccess: ARISE complete for all %s matches!" % len(finished_matches))
        print("All fully assigned sequences saved to %s" % final_designs_outfolder)
        exit()

    # prep all K* files
    fileprep_scan_round(new_scan_PDB_folder, visited_doublets, designID, targetID, out_kstar_files, out_hull_files,
                        completed_matches, length_chain, design_chirality)

    # start the runs on the cluster
    start_new_run(out_kstar_files, round_number)
    print("Jobs are running...")

    # periodically check our scan progress
    runs_done = False
    while not runs_done:
        time.sleep(1800)
        runs_done = are_kstar_finished(out_kstar_files)

    # if we finished this round, recursively start another
    run_ARISE(round_number+1, length_chain, completed_matches, visited_doublets, designID, targetID,
              final_designs_outfolder, apo_tolerance, design_chirality)

