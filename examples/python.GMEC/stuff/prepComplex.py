from argparse import ArgumentParser
from pathlib import Path

import prepTools


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


# Parse command-line arguments and execute workflow.
if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("--pdbname", type=str, required=True)
    parser.add_argument("--cld", default=[], type=str, nargs='*', required=True)
    parser.add_argument("--clb", default=[], type=str, nargs='*', required=True)
    parser.add_argument("--cutoffWM", type=float)
    parser.add_argument("--cutoffN", type=float, required=True)
    parser.add_argument("--heapSize",type=float)
    args = parser.parse_args()

    orig = Path(args.pdbname)
    stripped = prepTools.strip(orig)
    reduced = prepTools.reduce(stripped)
    prepped = prepTools.pdb4amber(reduced)
    rechained = prepTools.rechain(prepped, prepped)
    addedmissing = prepTools.addmissing(rechained)
    rechained2 = prepTools.rechain(addedmissing, addedmissing)
    minimized = prepTools.minimize(rechained2, args)
    final = prepTools.rechain(minimized, addedmissing)
