import math
import sys
import time

from Bio.PDB import PDBParser, PDBIO
from scipy.spatial import ConvexHull as ScipyConvexHull
from scipy.spatial import Delaunay

import os
import shutil

from Find_Doublets import SCOPE, pdb_to_coords, find_volume_overlap
from Confspace_Combiner import combine_Confspaces

# startup OSPREY (if not already running)
import jpype
if not jpype.isJVMStarted():
    import osprey
    osprey.start()
import osprey.prep

# BioPython throws errors that aren't really errors, so we'll ignore
import warnings

warnings.simplefilter('ignore')

# class that holds confspace info specs
class ConfSpaceSpecs:
    def __init__(self, doublet, flexset, mutations, graph_path):
        self.doublet = doublet
        self.flexset = flexset
        self.mutations = mutations
        self.graph_path = graph_path


def doublet_confspace_info(pdb_filename, designID, design_muts: list, design_chirality, outHulls, FixedResidues: list):
    print("--Obtaining confspace info using hulls--")

    doublet, flex = SCOPE(pdb_filename, outHulls, designID, design_muts, True, design_chirality, FixedResidues)

    all_confspace = []

    for p in doublet:
        if len(p) == 2:
            doub = list(p)
            res1_flex = flex[doub[0] - 1]
            res2_flex = flex[doub[1] - 1]
            flex_set = list(set(res1_flex).union(res2_flex))
            newspace = ConfSpaceSpecs(doub, flex_set, design_muts, doublet)
            all_confspace.append(newspace)
        elif len(p) == 1:
            doub = list(p)
            res_flex = flex[doub[0] - 1]
            newspace = ConfSpaceSpecs(doub, res_flex, design_muts, doublet)
            all_confspace.append(newspace)
        else:
            print("ERROR! 0 or > 2 residues in doublet", file=sys.stderr)

    return all_confspace


def array_hull_coords(pdb_name: str):

    xyz_coords = []

    f = open(pdb_name, "r")
    for line in f:
        if "TER" in line:
            continue
        coords = line.split('     ')[2]
        coords_format = ""
        counter = -1
        for i in coords:
            if i == ' ' and counter == 0:
                coords_format += ","
                counter += 1
            if i in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ".", "-"]:
                counter = 0
                coords_format += i

        all_coords = coords_format.split(',')
        atom_coords = []
        for i in all_coords:
            atom_coords.append(float(i))

        xyz_coords.append(atom_coords)

    return xyz_coords


def reduce_doublet(designID: str, targetID: str, spec, max_flex: int, hull_folder: str):

    # exit early if we are under the flex threshold
    if len(spec.flexset) < max_flex:
        print("Doublet %s has %s flexible residues, so will not be reduced" % (spec.doublet, len(spec.flexset)))
        return [spec]

    # use volume overlap to rank flex
    print("Doublet %s has %s flexible residues which is > max flex of %s residues." % (spec.doublet,
                                                                                       len(spec.flexset),
                                                                                       max_flex))
    hulls_loc = hull_folder + "/Chain%sRes%s.pdb"
    design_hulls = {}
    target_hulls = {}
    ordered_flex = {}

    # get the CH coords for doublet residues
    for res in spec.doublet:
        hull_loc = hulls_loc % (designID, res)
        coords = pdb_to_coords(hull_loc)
        design_hulls[res] = coords

    # get the CH for all flex residues
    for res in spec.flexset:
        hull_loc = hulls_loc % (targetID, res)
        coords = pdb_to_coords(hull_loc)
        target_hulls[res] = coords

    # compute mesh overlap using vtk boolean
    for tar_res, tar_hull in target_hulls.items():
        for des_res, des_hull in design_hulls.items():
            volume_overlap = find_volume_overlap(tar_hull, des_hull, True)
            if volume_overlap > 0.0:
                overlap_pair = (des_res, tar_res)
                ordered_flex[overlap_pair] = volume_overlap

    # order
    ordered_flex = dict(sorted(ordered_flex.items(), key=lambda item: item[1], reverse=True))
    print("Interchain pair (design, target): overlap area = %s" % ordered_flex)

    # prune
    ordered_interchain = [t for (d, t) in ordered_flex]
    selected_flex = []
    for res in ordered_interchain:
        if res not in selected_flex:
            selected_flex.append(res)
    selected_flex = selected_flex[:max_flex]

    # update doublet
    red_doublet = ConfSpaceSpecs(doublet=spec.doublet, flexset=selected_flex, mutations=[], graph_path=[])
    print("Updated doublet: mutate %s, flex %s" % (red_doublet.doublet, red_doublet.flexset))

    return [red_doublet]


def euclidean_dist(coord1, coord2):

    x_diff = (coord1[0] - coord2[0]) ** 2
    y_diff = (coord1[1] - coord2[1]) ** 2
    z_diff = (coord1[2] - coord2[2]) ** 2

    return math.sqrt((x_diff + y_diff + z_diff))


def get_CA_coords(chainID: str, resID: int, pdb_filename: str):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_filename)
    model = structure[0]

    for chain in model:
        if chain.id == chainID:
            for residue in chain:
                if residue.id[1] == resID:
                    for atom in residue:
                        if atom.name == 'CA':
                            coords = atom.get_vector()
                            return list(coords)

    print("ERROR! Couldn't find CA for residue %s in chain %s" % (resID, chainID), file=sys.stderr)
    exit()

def error_check_pdb(pdb_name: str):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    # return error if N/C-term PRO present (OSPREY anchor chords can't handle these currently)
    for chain in model:
        chain_length = len(chain)
        counter = 0
        for res in chain:
            counter += 1
            if counter == 1 or counter == chain_length:
                atype = res.get_resname()
                if atype == "PRO":
                    print("WARNING: PDB %s contains an N or C-term proline, which can't be used in OSPREY" % pdb_name)
                    print("Halting process for this file")
                    return "ERROR_TERM_PRO"

    # change CD -> CD1 labelling, which is common mislabel in PDBs for these residues
    for chain in model:
        for res in chain:
            resname = res.get_resname()
            resnum = res.id[1]
            if resname in ("ILE", "LEU", "PHE", "TRP", "TYR"):
                for a in res:
                    if a.fullname == ' CD ':
                        a.fullname = ' CD1'
                        print("Mislabel found with carbon label. Changing residue %s%s to CD1." % (resname, resnum))

    new_pdb_name = pdb_name.split(".")[0] + "-corrected.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(new_pdb_name)
    print("Now using corrected PDB file %s for future fileprep" % new_pdb_name)

    return new_pdb_name


def combine_PDBs(pdb_list, outfile_name, delete_old: bool, designID: str, flip_design: bool):
    with open(outfile_name, 'w') as outfile:
        for fname in pdb_list:
            with open(fname) as infile:
                for line in infile:
                    if "REMARK" in line:
                        pass
                    elif "END\n" in line:
                        outfile.write("TER\n")
                    else:
                        outfile.write(line)

    if delete_old:
        for pdb_file in pdb_list:
            os.remove(pdb_file)

    if flip_design:
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure('complex', outfile_name)
        model = structure[0]
        d_pep = model[designID]
        for r in d_pep:
            for a in r:
                a.coord[2] = a.coord[2] * -1

        io = PDBIO()
        io.set_structure(structure)
        io.save(outfile_name)
        print("Flipped chain back to D space and saved to %s" % outfile_name)


def osprey_fileprep_preprocess(pdb_name, designID: str, chirality: str, mindesign: bool, add_protons: bool):

    changed_pdb = False

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)

    if chirality == 'D':
        print("Have D-space chain, so created PDB in L-space for fileprep")
        d_pep = model[designID]
        for r in d_pep:
            for a in r:
                a.coord[2] = a.coord[2] * -1

        io = PDBIO()
        io.set_structure(structure)
        L_filename = pdb_name.split(".")[0] + "-flipped-design.pdb"
        io.save(L_filename)
        print("Saved flipped chain to %s" % L_filename)
        pdb_name = L_filename
        changed_pdb = True

    prep_pdb = osprey.prep.loadPDB(open(pdb_name, 'r').read())
    if designID == chain_ids[1]:
        target = prep_pdb[0]
        design = prep_pdb[1]
    elif designID == chain_ids[0]:
        target = prep_pdb[1]
        design = prep_pdb[0]
    else:
        print("ERROR! Was unable to find chain %s in PDB %s" % (designID, pdb_name), file=sys.stderr)

    mols = [target, design]

    with osprey.prep.LocalService():
        for mol in mols:
            for group in osprey.prep.duplicateAtoms(mol):
                for atomi in range(1, len(group.getAtoms())):
                    group.remove(atomi)
                    changed_pdb = True
                    print('Removed duplicate atom %s' % group)

        for mol in mols:
            for missing_atom in osprey.prep.inferMissingAtoms(mol):
                missing_atom.add()
                changed_pdb = True
                print('Added missing atom: %s' % missing_atom)

        if add_protons:
            for mol in mols:
                osprey.prep.deprotonate(mol)
                protonated_atoms = osprey.prep.inferProtonation(mol)
                for protonated_atom in protonated_atoms:
                    protonated_atom.add()
                    changed_pdb = True
                print('Added %d hydrogens to %s' % (len(protonated_atoms), mol))

    if changed_pdb:
        new_pdb_design = pdb_name.split(".")[0] + "-design-temp.pdb"
        new_pdb_target = pdb_name.split(".")[0] + "-target-temp.pdb"
        new_pdb_complex = pdb_name.split("-")[0] + "-processed.pdb"
        open(new_pdb_design, 'w').write(osprey.prep.savePDB(design))
        open(new_pdb_target, 'w').write(osprey.prep.savePDB(target))
        if chirality == 'D':
            combine_PDBs([new_pdb_target, new_pdb_design], new_pdb_complex, True, designID, True)
        elif chirality != 'D':
            combine_PDBs([new_pdb_target, new_pdb_design], new_pdb_complex, True, designID, False)
        print("Note: corrected or added atoms in PDB %s. Created new PDB %s for future fileprep." % (pdb_name, new_pdb_complex))
        return new_pdb_complex

    # TODO minimize the target chain if requested

    return pdb_name


def make_omol_file(pdb_name, isDesign: bool, want_complex: bool):
    processed_pdb = osprey.prep.loadPDB(open(pdb_name, 'r').read())

    target = processed_pdb[0]
    design = processed_pdb[1]
    mols = [target, design]

    with osprey.prep.LocalService():
        for mol in mols:
            bonds = osprey.prep.inferBonds(mol)
            for bond in bonds:
                mol.getBonds().add(bond)
            print('added %d bonds to %s' % (len(bonds), mol))

    if want_complex:
        omol_name = pdb_name.split(".")[0] + "-complex.omol"
        open(omol_name, 'w').write(osprey.prep.saveOMOL([target, design]))
    elif isDesign:
        omol_name = pdb_name.split(".")[0] + "-design.omol"
        open(omol_name, 'w').write(osprey.prep.saveOMOL([design]))
    else:
        omol_name = pdb_name.split(".")[0] + "-target.omol"
        open(omol_name, 'w').write(osprey.prep.saveOMOL([target]))


    print("Saved to %s" % omol_name)
    return omol_name


def get_omol_res_chain(omol_name: str, resnum: str, location: int):
    resType = ""
    chainID = ""

    # for complex omol, we need to specify target index (0) or design index (1)
    search_name = ("[molecule.%s.polymer]" % location)

    f = open(omol_name, 'r')
    resnum = resnum.strip()
    found_start = False
    found_chain = False

    for line in f:
        if found_start:
            chainID = line[1]
            found_chain = True
            found_start = False
        elif found_chain:
            id_region = line.split(",")[0]
            curr_resnum = id_region.split("\"")[1].strip()
            if curr_resnum == resnum:
                type_region = line.split(",")[1]
                curr_type = type_region.split("\"")[1]
                resType = curr_type
                break
        elif search_name in line:
            found_start = True

    return resType, chainID


def make_target_confspace(confspec: ConfSpaceSpecs, omol_name: str):
    osprey_omol = osprey.prep.loadOMOL(open(omol_name, 'r').read())
    target_conf = osprey_omol[0]
    target_conf_space = osprey.prep.ConfSpace(osprey_omol)
    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
    target_conf_space.getConflibs().add(lovell2000)

    flex_res = confspec.flexset
    for resnum in flex_res:
        resType, chainID = get_omol_res_chain(omol_name, str(resnum), 0)
        new_flex = target_conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_conf, chainID, str(resnum)))
        target_conf_space.addMutations(new_flex, resType)

    print('Target confspace for doublet %s:' % confspec.doublet)
    for pos in target_conf_space.positions():
        print('\t%6s flexing: %s' % (pos.getName(), target_conf_space.getMutations(pos)))

    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            target_conf_space.addConformationsFromLibraries(pos, mutation)
        # if pos.getType() in target_conf_space.getMutations(pos):
        #     target_conf_space.addWildTypeConformation(pos)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            for conf_info in target_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    print(
        'Target conformation space describes %s conformations for doublet %s' % (target_conf_space.countConformations(),
                                                                                 confspec.doublet))

    if len(confspec.doublet) == 0:
        doublet_print = "[MONTAGE]"
    else:
        doublet_print = ("[%s_%s]" % (confspec.doublet[0], confspec.doublet[1]))

    target_conf_path = omol_name.split(".")[0] + "-" + doublet_print + ".confspace"
    open(target_conf_path, 'w').write(osprey.prep.saveConfSpace(target_conf_space))
    print("Saved target confspace to %s" % target_conf_path)

    return target_conf_path

def get_omol_length(omol_name: str, location: int):

    search_name = ("[molecule.%s.polymer]" % location)

    f = open(omol_name, 'r')
    at_start = False
    num_res = 0
    for line in f:
        if at_start:
            if "{" in line:
                num_res += 1
        if search_name in line:
            at_start = True

    if num_res == 0:
        print("ERROR! Unable to find chain length from OMOL file", file=sys.stderr)
        exit()

    return num_res


def make_design_confspace(confspec: ConfSpaceSpecs, omol_name: str, mutTypes: list, translate_rotate: bool,
                          otherResAla: bool, fixed_residues: list):
    osprey_omol = osprey.prep.loadOMOL(open(omol_name, 'r').read())
    design_conf = osprey_omol[0]
    design_conf_space = osprey.prep.ConfSpace(osprey_omol)

    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'D-lovell2000-osprey3').load()
    design_conf_space.getConflibs().add(lovell2000)

    preHis_muts = len(mutTypes)
    try:
        mutTypes.remove('HIP')
        mutTypes.remove('HID')
        mutTypes.remove('HIE')
    except:
        pass
    postHis_muts = len(mutTypes)
    if preHis_muts != postHis_muts:
        mutTypes.append("HIS")

    for resnum in confspec.doublet:
        design_length = get_omol_length(omol_name, 0)
        resType, chainID = get_omol_res_chain(omol_name, str(resnum), 0)
        new_flex = design_conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, chainID, str(resnum)))
        if resnum in fixed_residues:
            print("Design chain %s%s has a fixed identity, so only getting WT flexibility" % (resType, resnum))
            design_conf_space.addMutations(new_flex, resType)
        elif resnum != 1 and resnum != design_length:
            design_conf_space.addMutations(new_flex, mutTypes)
        else:
            print("N or C terminus in doublet, so excluding PRO from mutations (if requested)")
            no_PRO = []
            for a in mutTypes:
                if a != 'PRO':
                    no_PRO.append(a)
            design_conf_space.addMutations(new_flex, no_PRO)

    if otherResAla:
        print("Now changing non-mutant design chain residues to ALA")
        design_length = get_omol_length(omol_name, 0)
        non_muts = []
        for i in range(1, design_length+1):
            if len(confspec.doublet) == 0:  # if doublet empty (for MONTAGE), just add to list
                if i not in fixed_residues:
                    non_muts.append(i)
            elif i != confspec.doublet[0] and i != confspec.doublet[1]:
                non_muts.append(i)
        for n in non_muts:
            resType, chainID = get_omol_res_chain(omol_name, "1", 0)
            new_ala = design_conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, chainID, str(n)))
            design_conf_space.addMutations(new_ala, 'ALA')

    print('Design confspace for doublet %s:' % confspec.doublet)
    for pos in design_conf_space.positions():
        print('\t%6s mutating: %s' % (pos.getName(), design_conf_space.getMutations(pos)))

    for pos in design_conf_space.positions():
        for mutation in design_conf_space.getMutations(pos):
            design_conf_space.addConformationsFromLibraries(pos, mutation)
        # if pos.getType() in design_conf_space.getMutations(pos):
        #     design_conf_space.addWildTypeConformation(pos)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in design_conf_space.positions():
        for mutation in design_conf_space.getMutations(pos):
            for conf_info in design_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    if translate_rotate:
        design_conf_space.addMotion(osprey.prep.moleculeTranslationRotation(design_conf))

    print(
        'Design conformation space describes %s conformations for doublet %s' % (design_conf_space.countConformations(),
                                                                                 confspec.doublet))

    if len(confspec.doublet) == 0:
        doublet_print = "[MONTAGE]"
    else:
        doublet_print = ("[%s_%s]" % (confspec.doublet[0], confspec.doublet[1]))

    design_conf_path = omol_name.split(".")[0] + "-" + doublet_print + ".confspace"
    open(design_conf_path, 'w').write(osprey.prep.saveConfSpace(design_conf_space))
    print("Saved design confspace to %s" % design_conf_path)

    return design_conf_path


def compile_confspaces(spaces: list):
    for s in spaces:
        confspace = osprey.prep.loadConfSpace(open(s, 'r').read())
        save_path = s.split(".")[0] + ".ccsx"

        compiler = osprey.prep.ConfSpaceCompiler(confspace)

        compiler.getForcefields().add(osprey.prep.Forcefield.Amber96)
        compiler.getForcefields().add(osprey.prep.Forcefield.EEF1)

        print('Compiling %s' % s)
        progress = compiler.compile()
        progress.printUntilFinish(10000)
        report = progress.getReport()

        if report.getError() is not None:
            raise Exception('Compilation failed', report.getError())

        open(save_path, 'wb').write(osprey.prep.saveCompiledConfSpace(report.getCompiled()))
        print('Saved compiled confspace to %s' % save_path)


def organize_kstar_files(out_directory: str, doublet: ConfSpaceSpecs, target_name: str, design_name: str, complex_name: str, KStarBash: str):

    if len(doublet.doublet) == 0:
        doublet_name = "[MONTAGE]"
    else:
        doublet_name = ("[%s_%s]" % (doublet.doublet[0], doublet.doublet[1]))

    new_directory = ("%s/kstar-%s" % (out_directory, doublet_name))
    try:
        os.mkdir(new_directory)
    except Exception as e:
        print("ERROR! Unable to make directory %s" % new_directory, file=sys.stderr)
        print(e)
        exit()

    t_dest = ("%s/%s" % (new_directory, target_name))
    print("Moving %s to %s" % (target_name, t_dest))
    os.rename(target_name, t_dest)

    d_dest = ("%s/%s" % (new_directory, design_name))
    print("Moving %s to %s" % (design_name, d_dest))
    os.rename(design_name, d_dest)

    c_dest = ("%s/%s" % (new_directory, complex_name))
    print("Moving %s to %s" % (complex_name, c_dest))
    os.rename(complex_name, c_dest)

    target_ccsx = target_name.split(".")[0] + ".ccsx"
    target_ccsx_path = ("%s/%s" % (new_directory, "target.ccsx"))
    print("Moving %s to %s" % (target_ccsx, target_ccsx_path))
    os.rename(target_ccsx, target_ccsx_path)

    design_ccsx = design_name.split(".")[0] + ".ccsx"
    design_ccsx_path = ("%s/%s" % (new_directory, "design.ccsx"))
    print("Moving %s to %s" % (design_ccsx, design_ccsx_path))
    os.rename(design_ccsx, design_ccsx_path)

    complex_ccsx = complex_name.split(".")[0] + ".ccsx"
    complex_ccsx_path = ("%s/%s" % (new_directory, "complex.ccsx"))
    print("Moving %s to %s" % (complex_ccsx, complex_ccsx_path))
    os.rename(complex_ccsx, complex_ccsx_path)

    try:
        e_direc = ("%s/%s" % (new_directory, "ensembles"))
        os.mkdir(e_direc)
    except Exception as e:
        print("ERROR! Unable to make directory %s" % e_direc, file=sys.stderr)
        print(e)
        exit()

    shutil.copy(KStarBash, new_directory)


def make_complex_confspace(doublet, omol_file, design_muts, trans_rot, otherResAla, FixedRes):

    # load omol
    osprey_omol = osprey.prep.loadOMOL(open(omol_file, 'r').read())
    target_conf = osprey_omol[0]
    design_conf = osprey_omol[1]
    conf_space = osprey.prep.ConfSpace(osprey_omol)

    # load lovell library
    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
    conf_space.getConflibs().add(lovell2000)

    # do HIS stuff
    preHis_muts = len(design_muts)
    try:
        design_muts.remove('HIP')
        design_muts.remove('HID')
        design_muts.remove('HIE')
    except:
        pass
    postHis_muts = len(design_muts)
    if preHis_muts != postHis_muts:
        design_muts.append("HIS")

    # set target flexibility
    flex_res = doublet.flexset
    for resnum in flex_res:
        resType, chainID = get_omol_res_chain(omol_file, str(resnum), 0)
        new_flex = conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_conf, chainID, str(resnum)))
        conf_space.addMutations(new_flex, resType)

    # mutate design chain
    for resnum in doublet.doublet:
        design_length = get_omol_length(omol_file, 1)
        resType, chainID = get_omol_res_chain(omol_file, str(resnum), 1)
        new_flex = conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, chainID, str(resnum)))
        if resnum in FixedRes:
            print("Design chain %s%s has a fixed identity, so only getting WT flexibility" % (resType, resnum))
            conf_space.addMutations(new_flex, resType)
        elif resnum != 1 and resnum != design_length:
            conf_space.addMutations(new_flex, design_muts)
        else:
            print("N or C terminus in doublet, so excluding PRO from mutations (if requested)")
            no_PRO = []
            for a in design_muts:
                if a != 'PRO':
                    no_PRO.append(a)
            conf_space.addMutations(new_flex, no_PRO)

    # if requested, mutate other res to ALA
    if otherResAla:
        print("Now changing non-mutant design chain residues to ALA")
        design_length = get_omol_length(omol_file, 1)
        non_muts = []
        for i in range(1, design_length+1):
            if len(doublet.doublet) == 0:  # if doublet empty (for MONTAGE), just add to list
                if i not in FixedRes:
                    non_muts.append(i)
            elif i != doublet.doublet[0] and i != doublet.doublet[1]:
                non_muts.append(i)
        for n in non_muts:
            resType, chainID = get_omol_res_chain(omol_file, "1", 1)
            new_ala = conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, chainID, str(n)))
            conf_space.addMutations(new_ala, 'ALA')

    # add motion to confspace
    print('Design confspace for doublet %s:' % doublet.doublet)
    for pos in conf_space.positions():
        print('\t%6s mutating: %s' % (pos.getName(), conf_space.getMutations(pos)))

    for pos in conf_space.positions():
        for mutation in conf_space.getMutations(pos):
            conf_space.addConformationsFromLibraries(pos, mutation)
        # if pos.getType() in design_conf_space.getMutations(pos):
        #     design_conf_space.addWildTypeConformation(pos)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in conf_space.positions():
        for mutation in conf_space.getMutations(pos):
            for conf_info in conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    # add translation/rotation if requested
    if trans_rot:
        conf_space.addMotion(osprey.prep.moleculeTranslationRotation(design_conf))

    print(
        'Design conformation space describes %s conformations for doublet %s' % (conf_space.countConformations(),
                                                                                 doublet.doublet))

    # save each .confspace file
    if len(doublet.doublet) == 0:  # handle MONTAGE
        doublet_print = "[MONTAGE]"
    else:
        doublet_print = ("[%s_%s]" % (doublet.doublet[0], doublet.doublet[1]))

    complex_conf_path = omol_file.split(".")[0].split('-')[0] + "-complex-" + doublet_print + ".confspace"
    target_conf_path = omol_file.split(".")[0].split('-')[0] + "-target-" + doublet_print + ".confspace"
    design_conf_path = omol_file.split(".")[0].split('-')[0] + "-design-" + doublet_print + ".confspace"

    open(complex_conf_path, 'w').write(osprey.prep.saveConfSpace(conf_space))
    print("Saved complex confspace to %s" % complex_conf_path)

    target_conf_space = conf_space.copy(target_conf)
    open(target_conf_path, 'w').write(osprey.prep.saveConfSpace(target_conf_space))
    print("Saved target confspace to %s" % target_conf_path)

    design_conf_space = conf_space.copy(design_conf)
    open(design_conf_path, 'w').write(osprey.prep.saveConfSpace(design_conf_space))
    print("Saved design confspace to %s" % design_conf_path)

    # return names for later file organization
    return target_conf_path, design_conf_path, complex_conf_path


def find_gly_pro_index(pdb_name: str, chain_id: str):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    gly_pro_index = []

    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                resname = residue.resname
                resid = residue.id[1]
                if resname == 'GLY':
                    gly_pro_index.append(resid)
                elif resname == 'PRO':
                    gly_pro_index.append(resid)

    print("Not mutating residues %s to ALA because GLY or PRO" % gly_pro_index)

    return gly_pro_index


def save_doublets(specs: list, filename: str):

    all_doublets = ""

    for i in range(0, len(specs)):
        if i != (len(specs) - 1):
            all_doublets += str(specs[i].doublet)
            all_doublets += ','
        else:
            all_doublets += str(specs[i].doublet)

    with open(filename, "w") as f:
        f.write(all_doublets)
        f.write('\n')
        for s in specs:
            new = ""
            new += str(s.doublet)
            new += ' '
            new += str(s.flexset)
            f.write(new)
            f.write('\n')

    print("Saved doublet info to %s" % filename)

    f.close()


# function for file prepping designs (D:L or L:L) for MONTAGE and ARISE
def osprey_fileprep_kstar(pdb_input, out_directory: str, design_muts: list, confspecs: list, designID: str,
                          design_chirality: str, minimize_design: bool, add_protons: bool, translation_rotation: bool,
                          otherResAla: bool, KStarBash: str, FixedResidues: list):
    print("Creating OSPREY files for your design")

    print("--Checking input PDB file for any labelling or content issues--")
    correct_pdb = error_check_pdb(pdb_input)
    if correct_pdb == "ERROR_TERM_PRO":
        return "DELETE"

    print("--Preprocessing PDB file--")
    preprocess_pdb = osprey_fileprep_preprocess(correct_pdb, designID, design_chirality, minimize_design, add_protons)

    print("--Making OMOL files--")
    if design_chirality == 'D':
        design_omol = make_omol_file(preprocess_pdb, True, False)
        target_omol = make_omol_file(preprocess_pdb, False, False)
    elif design_chirality == 'L':
        complex_omol = make_omol_file(preprocess_pdb, False, True)

    print("--Making OSPREY confspace files--")
    for spec in confspecs:

        print("Now making confspace files for doublet %s with flexset %s" % (spec.doublet, spec.flexset))

        if design_chirality == 'L':

            target_conf_name, design_conf_name, complex_conf_name = make_complex_confspace(spec, complex_omol,
                                                                                           design_muts,
                                                                                           translation_rotation,
                                                                                           otherResAla, FixedResidues)

        elif design_chirality == 'D':

            target_conf_name = make_target_confspace(spec, target_omol)
            design_conf_name = make_design_confspace(spec, design_omol, design_muts, translation_rotation,
                                                 otherResAla, FixedResidues)
            complex_conf_name = combine_Confspaces(target_conf_name, design_conf_name)

        print("--Compiling OSPREY confspace files for doublet %s--" % spec.doublet)
        with osprey.prep.LocalService():
            compile_confspaces([target_conf_name, design_conf_name, complex_conf_name])

        print("--Organizing files for %s into directory %s--" % (spec.doublet, out_directory))
        organize_kstar_files(out_directory, spec, target_conf_name, design_conf_name, complex_conf_name, KStarBash)

    if design_chirality == 'L':
        complex_omol_path = ("%s/%s" % (out_directory, complex_omol))
        os.rename(complex_omol, complex_omol_path)
        print("Moved complex omol to %s" % complex_omol_path)

    elif design_chirality == 'D':
        design_omol_path = ("%s/%s" % (out_directory, design_omol))
        target_omol_path = ("%s/%s" % (out_directory, target_omol))
        os.rename(design_omol, design_omol_path)
        print("Moved design omol to %s" % design_omol_path)
        os.rename(target_omol, target_omol_path)
        print("Moved target omol to %s" % target_omol_path)

    return "READY"

