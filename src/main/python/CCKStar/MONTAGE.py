import glob
import shutil
import time
import tempfile

from Bio.PDB import PDBParser, PDBIO
import sys
import warnings
import os
import fileinput
import subprocess

from Find_Doublets import SCOPE, rank_flex_overlap
from KStarPrep import find_gly_pro_index, osprey_fileprep_kstar, ConfSpaceSpecs

# start OSPREY (if not already running)
import jpype
if not jpype.isJVMStarted():
    import osprey
    osprey.start()
import osprey.prep

warnings.simplefilter('ignore')

def complex_renamer(complex_pdb: str, desired_target_id: str, desired_design_id: str):

    print("\nRenaming %s IDs. Target ID: %s, Design ID: %s.\n" % (complex_pdb, desired_target_id, desired_design_id))

    # get the old IDs
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", complex_pdb)
    model = structure[0]

    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)
    old_target_chain = chain_ids[0]
    old_peptide_chain = chain_ids[1]

    # change both IDs
    for chain in model:
        if chain.id == old_target_chain:
            chain.id = desired_target_id
        elif chain.id == old_peptide_chain:
            chain.id = desired_design_id

    # write out the new PDB with suffix -renamed.pdb
    io = PDBIO()
    io.set_structure(structure)
    filename = complex_pdb[:-4] + "-renamed.pdb"
    io.save(filename)
    print("Saved renamed PDB to %s" % filename)

    return filename

def reflect_design(complex_pdb: str, targetID: str, designID: str):
    print("Reflecting chain %s" % designID)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", complex_pdb)
    model = structure[0]

    for chain in model:
        if chain.id == designID:
            counter = 1
            for res in chain:
                for a in res:
                    a.coord[2] = a.coord[2] * -1
                res.id = (' ', counter, ' ')
                counter += 1

    model.detach_child(targetID)
    io = PDBIO()
    io.set_structure(structure)
    filename = complex_pdb[:-4] + "-inverted-design.pdb"
    io.save(filename)
    print("Saved inverted PDB to %s" % filename)

    return filename


def extract_design(complex_pdb: str, targetID: str, designID: str):

    print("Extracting chain %s" % designID)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", complex_pdb)
    model = structure[0]

    for chain in model:
        if chain.id == designID:
            counter = 1
            for res in chain:
                # start numbering from 1
                res.id = (' ', counter, ' ')
                counter += 1

    model.detach_child(targetID)
    io = PDBIO()
    io.set_structure(structure)
    filename = complex_pdb[:-4] + "-design.pdb"
    io.save(filename)
    print("Saved extracted PDB to %s" % filename)

    return filename

def chain_relabel(input_directory: str, chainID: str):

    for pdb_path in os.listdir(input_directory):
        parser = PDBParser(PERMISSIVE=1)
        full_pdb_path = os.path.join(input_directory, pdb_path)
        M_pep_structure = parser.get_structure('match', full_pdb_path)
        M_peptide_model = M_pep_structure[0]

        for chain in M_peptide_model:
            chain.id = chainID

        io = PDBIO()
        io.set_structure(M_peptide_model)
        io.save(full_pdb_path)

        print("Renamed %s chain to %s" % (full_pdb_path, chainID))


# note: atom labels have not yet been corrected, so some side chain H are not correct
# we only need bb H atoms, so this is ok
def protonate_chains(matches_directory: str, target_pdb: str):

    # protonate the target
    print("\n\n--Protonating L-target--\n\n")
    loaded_target_pdb = osprey.prep.loadPDB(open(target_pdb, 'r').read())
    targetChain = loaded_target_pdb[0]
    D_DesignChain = loaded_target_pdb[1]

    with osprey.prep.LocalService():
        osprey.prep.deprotonate(targetChain)
        protonated_atoms = osprey.prep.inferProtonation(targetChain)
        for protonated_atom in protonated_atoms:
            protonated_atom.add()
        print('Added %d hydrogens to %s' % (len(protonated_atoms), targetChain))

    # writing PDB via OSPREY is a pain, so we'll do some messy manual work to add chain sequentially
    open(target_pdb, 'w').write(osprey.prep.savePDB(targetChain))
    open(target_pdb, 'a').write(osprey.prep.savePDB(D_DesignChain))
    for line in fileinput.input(target_pdb, inplace=True):
        if "END" in line:
            print("TER")
        elif "REMARK" in line:
            continue
        else:
            print(line, end='')
    print("saved protonated PDB to %s" % target_pdb)

    # protonate each match (in-place)
    print("\n\n--Protonating L-matches--\n\n")
    for pdb_path in os.listdir(matches_directory):
        full_pdb_path = os.path.join(matches_directory, pdb_path)
        prep_pdb = osprey.prep.loadPDB(open(full_pdb_path, 'r').read())

        designChain = prep_pdb[0]

        with osprey.prep.LocalService():
            osprey.prep.deprotonate(designChain)
            protonated_atoms = osprey.prep.inferProtonation(designChain)
            for protonated_atom in protonated_atoms:
                protonated_atom.add()
            print('Added %d hydrogens to %s' % (len(protonated_atoms), designChain))

        open(full_pdb_path, 'w').write(osprey.prep.savePDB(designChain))
        print("saved protonated PDB to %s" % full_pdb_path)


def scaffold_generator(complex_pdb: str, designID: str, input_directory: str,
                       output_directory: str, desired_chirality: str):

    print("\n\n------ running scaffold generator ------\n\n")

    print("Storing outputs in %s" % output_directory)
    os.mkdir(output_directory)

    total_disjoint = 0

    for pdb_path in os.listdir(input_directory):

        # get the complex
        parser = PDBParser(PERMISSIVE=1)
        complex_structure = parser.get_structure('complex', complex_pdb)
        complex_model = complex_structure[0]

        # get the MASTER L-peptide
        full_pdb_path = os.path.join(input_directory, pdb_path)
        M_pep_structure = parser.get_structure('match', full_pdb_path)
        M_peptide_model = M_pep_structure[0]

        # skip the match if it has a disjoint segment
        has_disjoint = False
        last_res = 0
        for chain in M_peptide_model:
            for res in chain:
                if last_res == 0:
                    last_res = res.id[1]
                elif (res.id[1] - last_res) != 1:
                    has_disjoint = True
                    print("SKIP: %s has a disjoint segment and will not be prepared" % pdb_path)
                    total_disjoint += 1
                    break
                else:
                    last_res = res.id[1]
        if has_disjoint:
            continue

        # flip each MASTER L-match to D-space (if requested), and change numbering and ID if needed
        for chain in M_peptide_model:
            counter = 1
            for res in chain:
                for a in res:
                    if desired_chirality == 'D':
                        a.coord[2] = a.coord[2] * -1
                res.id = (' ', counter, ' ')
                counter += 1

        # MASTER already aligns, so just add the peptide
        for chain in M_peptide_model:
            complex_model.add(chain)

        # remove the old peptide from the complex
        complex_model.detach_child(designID)

        # print new PDB to output_directory
        io = PDBIO()
        io.set_structure(complex_model)
        pdb_name = pdb_path[:-4] + "-complex.pdb"
        savename = os.path.join(output_directory, pdb_name)
        io.save(savename)

        # update terminal
        print("------ completed %s ------" % pdb_path)

    print("\n\nscaffold generation completed")
    print("Deleted %s matches due to disjoint segments\n\n" % total_disjoint)


def submit_MASTER(pdb_filename: str, database_location: str, number_matches: str):

    # make the protein data structure (pds) file
    print("Creating PDS...")
    pds_request = subprocess.run(["./resources/createPDS", "--type", "query", "--pdb", pdb_filename, '--pds',
                                  "query.pds"], capture_output=True, text=True)

    # log + catch any errors
    print(pds_request.stdout)
    if pds_request.stderr:
        print("ERROR creating PDS for %s" % pds_request, file=sys.stderr)
        print(pds_request.stdout, file=sys.stderr)
        sys.exit(1)

    # run MASTER (usually takes ~10 minutes to complete)
    print("Running MASTER...")
    MASTER_request = subprocess.run(["./resources/master", "--query", "query.pds", "--targetList",
                                     database_location, "--rmsdCut", "10.0", "--topN", number_matches,
                                     "--outType", "match", "--seqOut", "./matches.txt", "--structOut", "./matches"],
                                    capture_output=True, text=True)

    # catch any errors (don't log stdout, it's thousands of lines)
    if MASTER_request.stderr:
        print("ERROR running MASTER for %s" % pds_request, file=sys.stderr)
        print(pds_request.stdout, file=sys.stderr)
        sys.exit(1)
    print("%s MASTER search completed!" % pdb_filename)


def run_MASTER(pdb_name: str, input_chirality: str, output_chirality: str, MASTER_database: str,
               number_matches: str):

    # use A/z naming to avoid issues with PDBs having identical chain names when using MASTER
    renamed_pdb = complex_renamer(pdb_name, 'A', 'z')

    if input_chirality == 'L' and output_chirality == 'D':
        print("Preparing to generate D-peptide scaffolds from L-peptide binder")
        query_pdb = reflect_design(renamed_pdb, 'A', 'z')

    elif input_chirality == 'D' and output_chirality == 'D':
        print("Preparing to generate D-peptide scaffolds from D-peptide binder")
        query_pdb = reflect_design(renamed_pdb, 'A', 'z')

    elif input_chirality == 'L' and output_chirality == 'L':
        print("Preparing to generate L-peptide scaffolds from L-peptide binder")
        query_pdb = extract_design(renamed_pdb, 'A', 'z')

    else:
        print("Preparing to generate L-peptide scaffolds from D-peptide binder")
        query_pdb = extract_design(renamed_pdb, 'A', 'z')

    # run the MASTER search
    submit_MASTER(query_pdb, MASTER_database, number_matches)

    # rename MASTER returns (in-place) to B (avoid name clash when generating scaffolds)
    chain_relabel('matches', 'B')

    # add protons to these L-space backbones and the L-target (for later hull placement)
    # We haven't corrected atom labels so some sc protons will be wrong, but this is ok. We only care about bb protons.
    protonate_chains('matches', renamed_pdb)

    # use MASTER returns to generate D or L design:L-protein PDBs
    scaffold_generator(renamed_pdb, 'z', "matches", "scaffolds", output_chirality)


def cleanup_MASTER(file: str):

    # remove the copied filed
    os.remove(file)

    # create directory for output
    out_directory = "%s-MASTER" % file.split('.')[0]
    os.mkdir(out_directory)

    # copy over renamed files
    pdb_files = file.split('.')[0] + '*'
    for renamed in glob.glob(pdb_files):
        shutil.move(renamed, out_directory)

    # copy over MASTER outputs and scaffolds
    shutil.move('./query.pds', out_directory)
    shutil.move('matches.txt', out_directory)
    shutil.move('matches', out_directory)
    shutil.move('scaffolds', out_directory)

    print("Moved all MASTER outputs to %s" % out_directory)


def target_flex_MONTAGE(pdb_name: str, design_chirality: str, max_flex: int):

    print("Finding target flexibility for MONTAGE")
    with tempfile.TemporaryDirectory() as tmpdir:

        intrachain_contacts, interchain_contacts = SCOPE(pdb_name, tmpdir, 'B', ['VAL'],
                                                         True, design_chirality, [])

        # ranks by decreasing overlap
        target_flex = rank_flex_overlap('B', 'A', intrachain_contacts,
                                        interchain_contacts, tmpdir)

        print("Interchain pair (design, target): overlap area = %s" % target_flex)

    ordered_interchain = [t for (d, t) in target_flex]
    selected_flex = []
    for res in ordered_interchain:
        if res not in selected_flex:
            selected_flex.append(res)
    selected_flex = selected_flex[:max_flex]

    print("Selected target residues for flex: %s" % selected_flex)

    return selected_flex


def get_MONTAGE_gmec(id):

    out_directory = "%s-MONTAGE_GMEC" % id
    os.mkdir(out_directory)
    print("Storing outputs in %s" % out_directory)

    ordered_scores = {}

    # check the logfile and find out which runs were successful (> 0 K* score)
    logfile_locations = "%s-MONTAGE/match*-MONTAGE/kstar*/submit.out" % id

    for logfile in glob.glob(logfile_locations):

        match_num = logfile.split('/')[1]

        print("Checking scores for %s" % logfile)

        f = open(logfile, 'r')
        for line in f:

            if line.startswith("1,"):
                score = line.split(',')[2]
                if score.lower() == "none" or float(score) <= 0:
                    print("SKIP: negative or 0 K* score")
                    ordered_scores[match_num] = 0
                else:
                    print("GOOD: the K* score is %s" % score)
                    ordered_scores[match_num] = float(score)
                    pdb_location = (logfile[:-11] + "/ensembles/*pdb").replace("[MONTAGE]", "*")
                    save_location = out_directory + '/' + (logfile.split('-')[1].split('/')[1] + "-complex.pdb")
                    for file in glob.glob(pdb_location):
                        shutil.copy(file, save_location)

    print("Saved all valid matches for %s to %s\n" % (id, out_directory))
    print("MONTAGE matches ordered by K* score: %s" % dict(sorted(ordered_scores.items(), key=lambda item: item[1],
                                                                  reverse=True)))


def check_MONTAGE_done():

    # count the number of matches for all PDBs
    total_matches = 0
    for matchfile in glob.glob("*-MONTAGE/match*-MONTAGE"):
        total_matches += 1

    # count finished matches
    finished_matches = 0
    for logfile in glob.glob("*-MONTAGE/match*-MONTAGE/kstar*/submit.out"):
        f = open(logfile, 'r')
        for line in f:
            if "completed" in line:
                finished_matches += 1
                break

    # return status
    if total_matches == finished_matches:
        return True
    else:
        print("%s out of %s matches are still running..." % ((total_matches - finished_matches), total_matches))
        return False


def cluster_runner_MONTAGE():

    # copy over the submit script
    shutil.copy("resources/montage_cluster_runner.sh", "./montage_cluster_runner.sh")

    # submit to the cluster
    MONTAGE_request = subprocess.run(["./montage_cluster_runner.sh"], capture_output=True, text=True)

    # log + catch any errors
    print(MONTAGE_request.stdout)
    if MONTAGE_request.stderr:
        print("ERROR submitting MONTAGE to cluster!", file=sys.stderr)
        print(MONTAGE_request.stdout, file=sys.stderr)
        sys.exit(1)

    # delete the copied file
    os.remove("./montage_cluster_runner.sh")

    # check K* runs status
    montage_done = False
    while not montage_done:
        montage_done = check_MONTAGE_done()
        time.sleep(300)
    print("All matches finished running MONTAGE!")

    # get the GMEC PDBs for ARISE
    for id in glob.glob("*-MONTAGE/"):
        pdb_id = id.split("-")[0]
        get_MONTAGE_gmec(pdb_id)


def run_MONTAGE(input_pdb_directory: str, input_chirality: str, output_chirality: str, MASTER_matches: int,
                max_flex: int):

    # generate MASTER scaffolds from a directory containing PDB files
    for file in os.listdir(input_pdb_directory):
        print("\n--- Now creating scaffolds for %s ---\n" % file)

        # copy over the file to home for easy access
        full_filepath = os.path.join(input_pdb_directory, file)
        shutil.copy(full_filepath, file)

        # create scaffolds using MASTER
        run_MASTER(file, input_chirality, output_chirality, "./resources/db.txt", str(MASTER_matches))

        # organize outputs
        cleanup_MASTER(file)

    # run K* to update chemical contacts of MASTER returns to lower energy
    for scaff_dir in glob.glob("*-MASTER/scaffolds"):

        print("\n--- Now prepping K* files for %s ---\n" % scaff_dir)

        pdb_name = scaff_dir.split('/')[0].split('-')[0]

        # setup output dir for this pdb
        out_pdb_foldername = ("%s-MONTAGE" % pdb_name)
        os.mkdir(out_pdb_foldername)

        # find contacts and prepare OSPREY K* files for each scaffold
        for scaff in os.listdir(scaff_dir):

            out_match_foldername = ("%s-MONTAGE" % scaff.split('-')[0])
            os.mkdir(out_match_foldername)

            print("\nNow preparing %s\n" % scaff.split('-')[0])

            # copy over the scaffold pdb for easy access
            full_filepath = os.path.join(scaff_dir, scaff)
            shutil.copy(full_filepath, scaff)

            # find flex on L-target
            target_flex = target_flex_MONTAGE(scaff, output_chirality, max_flex)

            # maintain GLY and PRO residues
            gly_pro = find_gly_pro_index(scaff, "B")

            # prep K* files
            prep_status = osprey_fileprep_kstar(scaff, out_match_foldername, [],
                                                [ConfSpaceSpecs(doublet=[], flexset=target_flex, mutations=[],
                                                                graph_path=[])],
                                                'B', output_chirality, False, True,
                                                True, True, 'resources/K_bash.sh',
                                                gly_pro)

            # delete everything if this match isn't viable
            if prep_status == 'DELETE':
                print("\n--Deleting all prep files for %s--\n" % scaff.split('-')[0])
                shutil.rmtree(out_match_foldername)
                for pdb in glob.glob("%s*.pdb" % scaff.split("-")[0]):
                    os.remove(pdb)

            # cleanup match PDBs
            else:
                os.remove(scaff)
                for pdb in glob.glob("%s*.pdb" % scaff.split("-")[0]):
                    shutil.move(pdb, out_match_foldername)
                    print("Moved %s to %s" % (pdb, out_match_foldername))

        # move all match prep files into a single directory for that PDB
        for matchdir in glob.glob("match*-MONTAGE"):
            shutil.move(matchdir, out_pdb_foldername)

    # run these file on the cluster. When finished, get the GMEC PDBs.
    cluster_runner_MONTAGE()
