from Bio.PDB import PDBParser
from pathlib import Path
import confTools
import subprocess
import itertools
import osprey

import confTools
import constructMaskFromList


# Allen McBride
# June 5, 2024

# This script automates several steps in one of our common workflows for
# protein redesign:
#
# * Strip hydrogens from a given structure.
#
# * Re-add hydrogens according to the Reduce tool from the Richardson lab.
#
# * Add missing atoms according to the pdb4amber tool from AmberTools.
#
# * Restore orginial chain identifiers (which are stripped by pdb4amber).
#
# * Identify residues meeting certain distance criteria according to the
# confTools module. (Currently this requires the osprey module, although the
# code could be refactored in the future to avoid this requirement.)
#
# * Use the Probe tool from the Richardson lab to identify residues involved in
# clashes of the "Bad Overlap" category.
#
# * Use the tleap and sander tools from AmberTools to perform molecular
# dynamics minimization, restraining all atoms not in the above residues (the
# number of steps is currently hard-coded by this script, in the nSteps
# variable).
#
# * Restore original chain identifiers once again.
#
# As these steps are performed, all intermediate PDB files are saved to the
# current directory. If a file already exists with the filename that this
# script would have written, assume the existing file is the one desired and
# proceed with subsequent steps.
#
# Each function that (potentially) generates a new PDB file returns a Path
# object corresponding to the desired PDB file.
#
# Command-line parameters, all required, are as follows:
#
# --pdbname <filename> (<filename> names a PDB file describing a single model
# of a starting structure.)
#
# --cld <list> (<list> is a space-separated list of zero or more chain
# identifiers we wish to redesign.)
#
# --clb <list> (<list> is a space-separated list of zero or more chain
# identifiers we wish to bind.)
#
# --cutoffWM <X> (<X> is a floating-point distance in Angstroms. We wish to
# consider mutations in all residues in the chains indicated by the --cld flag
# that have side-atoms within X Angstroms of any atom in the chains indicated
# by the --clb flag.
#
# --cutoffN <X> (<X> is a floating-point distance in Angstroms. We wish to
# consider alternative rotamers in all residues that have side-chain atoms
# within X Angstroms of any side-chain atom in a residue for which mutation is
# to be considered.)
#
# Example invokation: python ./prepComplex.py --pdbname test.pdb --cld H L
# --clb B --cutoffWM 6 --cutoffN 4


# Find "Bad Overlap" clashes in a structure and return a "restrainmask" in the
# format required by AmberTools' sander. This restraintmask restrains during
# minimization all atoms not part of residues involved in these clashes.
# "pdbpath" is a Path to a PDB file, and "model" is the same structure in the
# form of an object of Biopython's Model class.
def findClashesAndMask(pdbpath, model, mutableSetInfos):

    allInvolved = set().union(*[i.allMut | i.allAla | i.allFlex for i in mutableSetInfos])
    probeoutput = subprocess.run(
            ["probe", "-4H", "-mc", "-SElf", "-Unformated", "-CONdense", "-ONLYBADOUT",
                confTools.setToProbe(allInvolved), str(pdbpath)],
            stdout = subprocess.PIPE, stderr = subprocess.DEVNULL, text = True)
    resInvolvedInClashes = set()
    for line in probeoutput.stdout.split('\n'):
        byColons = line.split(':')
        if len(byColons) > 4:
            resInvolvedInClashes.add(byColons[3][3:6])
            resInvolvedInClashes.add(byColons[4][3:6])
    maxResNum = max(res.get_id()[1] for res in itertools.chain(*model))
    with Path('involvedResidues.txt').open(mode='w') as involved:
        involved.write('All residues involved in any conformation space, mutable or flexible:\n')
        involved.write(confTools.setToPyMol(allInvolved) + '\n')
    return constructMaskFromList.findMask(resInvolvedInClashes, maxResNum)


# If an output PDB file already exists, assume the existing file is what is
# desired and proceed, using this existing file.
def checkfile(path):
    if path.exists():
        if path.is_file():
            if path.stat().st_size > 0:
                print(f'Path for output {path} exists and is nonempty. Assuming correct.')
                return True
            else:
                print(f'Path for output {path} exists but is empty. Deleting.')
                path.unlink()
                return False
        else:
            print(f'Path {path} exists but is not a file. Exiting.')
            sys.exit()
    else:
        return False


def his(inpath):
    outpath = inpath.with_stem(inpath.stem + '_his')
    count = 0
    if not checkfile(outpath):
        print(f'Renaming HID, HIE, HIP to HIS in {inpath}. Output path: {outpath}')
        with inpath.open(mode='r') as infile, outpath.open(mode='x') as hisFile:
            for line in infile:
                if len(line) >= 20 and line[0:4] == 'ATOM' and (line[17:20] == 'HID' or line[17:20] == 'HIE' or line[17:20] == 'HIP'):
                    hisFile.write(line[0:17] + 'HIS' + line[20:])
                    count += 1
                else:
                    hisFile.write(line)
    return outpath, count


# Invoke the Reduce tool to remove existing hydrogen atoms.
def strip(inpath, force = False):
    outpath = inpath.with_stem(inpath.stem + '_strip')
    if force or not checkfile(outpath):
        print(f'Stripping {inpath}. Output path: {outpath}.')
        modechar = 'w' if force else 'x'
        with outpath.open(mode=modechar) as stripFile:
            subprocess.run(['reduce', '-Trim', str(inpath)], stdout = stripFile)
    return outpath


# Invoke the Reduce tool to add hydrogen atoms.
def reduce(inpath, force = False):
    outpath = inpath.with_stem(inpath.stem + '_reduce')
    if force or not checkfile(outpath):
        print(f'Reducing {inpath}. Output path: {outpath}.')
        modechar = 'w' if force else 'x'
        with outpath.open(mode=modechar) as reduceFile:
            subprocess.run(['reduce', '-BUILD', str(inpath)], stdout = reduceFile)
    return outpath


# Invoke the pdb4amber tool to do most of its default work, but do not include
# --add-missing-atoms, which seems to interfere with choosing alternates.  This
# strips chain identifiers, renumbers all residues consecutively, and renames
# CYS to CYX and HIS to HID, HIE, or HIP as appropriate.
def pdb4amber(inpath):
    outpath = inpath.with_stem(inpath.stem + '_prep')
    if not checkfile(outpath):
        print(f'Running pdb4amber on {inpath}. Output path: {outpath}.')
        subprocess.run(['pdb4amber', '-i',str(inpath), '-o', str(outpath)])
    return outpath


# Invoke the pdb4amber tool with the --add-missing-atoms flag.
def addmissing(inpath):
    outpath = inpath.with_stem(inpath.stem + '_add')
    if not checkfile(outpath):
        print(f'Running pdb4amber to add missing atoms on {inpath}. Output path: {outpath}.')
        subprocess.run(['pdb4amber', '--add-missing-atoms', '-i',str(inpath), '-o', str(outpath)])
    return outpath


# Based on p4a-undo.py from Osprey. Use the '_renum.txt' file emitted by
# pdb4amber to restore stripped chain identifiers, but do NOT restore original
# residue numbering scheme. Also remove CONECT records.
def rechain(inpath, renumpathbase):
    outpath = inpath.with_stem(inpath.stem + '_rc')
    if not checkfile(outpath):
        renumpath = renumpathbase.with_stem(renumpathbase.stem + '_renum').with_suffix('.txt')
        print(f'Restoring chain IDs for {inpath} based on {renumpath}. Output path: {outpath}.')
        with renumpath.open(mode='r') as renum_lines:
            renum_dict = { line.split()[4] : line.split()[1] for line in renum_lines }
        with outpath.open(mode='x') as newfile, inpath.open(mode='r') as oldfile:
            for line in oldfile:
                if not line.strip(): # Remove blank lines
                    continue
                if line[:6] == 'CONECT': # Remove CONECT records
                    continue
                if line[:4] != 'ATOM' and line[:4] != 'HETA' and line[:3] != 'TER':
                    print(line, end='', file = newfile)
                    continue
                res_num = line[22:26].strip()
                chain = renum_dict[res_num]
                print(f'{line[:21]}{chain}{line[22:]}'.strip(), file = newfile)
    return outpath


# Use the tleap and sander molecular dynamics tools from AmberTools to minimize
# a structure, along with some preliminary steps:
# 
# * Identify residues meeting certain distance criteria according to the
# confTools module.
#
# * Use the Probe tool from the Richardson lab to identify residues involved in
# clashes of the "Bad Overlap" category.
#
# * Use the tleap and sander tools from AmberTools to perform molecular
# dynamics minimization, restraining all atoms not in the above residues (the
# number of steps is currently hard-coded by this script, in the nSteps
# variable).
def minimize(inpath, args):
    nSteps = 100 # Number of minimization steps to perform
    outpath = inpath.with_stem(inpath.stem + f'_steps{nSteps}')
    if not checkfile(outpath):
        print(f'Running sander on {inpath}. Output path: {outpath}.')

        model = PDBParser().get_structure("X", inpath)[0]
        chainListDesign = [model[chainid] for chainid in args.cld]
        chainListBinder = [model[chainid] for chainid in args.clb]

        osprey.start(heapSizeMiB=args.heapSize)
        mol = osprey.readPdb(inpath)
        ffparams = osprey.ForcefieldParams()
        templateLib = osprey.TemplateLibrary(ffparams.forcefld)

        # mutableSetInfos = confTools.findConfSpacesSet(inpath, mol, templateLib,
        #         model, chainListDesign, chainListBinder, [], args.cutoffWM,
        #         args.cutoffN, [], [], 1, {}, False, False)

        confSpaces = confTools.ConfSpaces(
            name = 'B781',
            csDesign = None,
            csBinder = None,
            csComplex = None
        )
        chain_id = "B"
        resseqMut = 781
        resseqFlex1 = 776
        resseqFlex2 = 735
        icode = " "  # insertion code, use " " if none
        hetfield = " "  # " " means it's a standard amino acid (not a heteroatom)

        residueMut = model[chain_id][(hetfield, resseqMut, icode)]
        residueFlex1 = model[chain_id][(hetfield, resseqFlex1, icode)]
        residueFlex2 = model[chain_id][(hetfield, resseqFlex2, icode)]


        mutableSetInfos = {confTools.MutableSetInfo(
            confSpaces=confSpaces,
            allMut={residueMut},
            allAla=set(),
            allFlex={residueFlex1,residueFlex2},
            bbFlex=set(),
            tr=[]
        )}


        for mutableSetInfo in mutableSetInfos:
            print("mutableSetInfo", mutableSetInfo.printinfo())


        fixedContents1 = \
"""Minimizing input structure to resolve clashes, generated by prepComplex.py
&cntrl
 imin=1,
 maxcyc="""
        fixedContents2 = """,
 ntpr=1,
 ntwx=10,
 ntxo=1,
 ntr=1,
 ioutfm=0,
 igb=1,
 restraint_wt=5.0,
 restraintmask="""

        restraintMask = findClashesAndMask(inpath, model, mutableSetInfos)
        allContents = fixedContents1 + str(nSteps) + fixedContents2 + f'"{restraintMask}"\n/\n'
        sanderControlPath = Path('sanderControl.mdin')
        with sanderControlPath.open(mode='w') as sanderControl:
            sanderControl.write(allContents)

        tleapInPath = Path('tleap.in')
        prmtop = inpath.with_suffix('.prmtop')
        inpcrd = inpath.with_suffix('.inpcrd')
        with tleapInPath.open(mode='w') as tleapIn:
            tleapIn.write('source leaprc.protein.ff14SB\n')
            tleapIn.write(f'mol = loadpdb {str(inpath)}\n')
            tleapIn.write('check mol\n')
            tleapIn.write(f'saveamberparm mol {str(prmtop)} {str(inpcrd)}\n')
            tleapIn.write(f'quit\n')

        tleapOutPath = Path('tleap.out')
        with tleapOutPath.open(mode='x') as tleapOut:
            subprocess.run(['tleap', '-f', str(tleapInPath)], stdout = tleapOut)
            
            
        temp_str = [
            'sander', '-O',
            '-i', str(sanderControlPath),
            '-o', 'sanderOutput.mdout',
            '-p', str(prmtop),
            '-c', str(inpcrd),
            '-ref', str(inpcrd),
            '-x', 'coords.mdcrd'
            ]
            

        subprocess.run(temp_str)

        with outpath.open(mode='x') as minimized:
            subprocess.run(['ambpdb', '-p', str(prmtop), '-c', 'restrt'], stdout = minimized)

    return outpath
