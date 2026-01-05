from Bio.PDB import PDBParser
from collections import defaultdict
from pathlib import Path
import math
import sys
import osprey
import csv


# Allen McBride
# May 30, 2024
# 
# This is an inefficient library for two general types of tasks:
#
#  * Compute sets of flexible and/or mutable residues for typical
#    Osprey searches
# 
#  * Create conformation spaces from that information, for use with
#    Osprey's Python API
#
# This process could be made more efficient by using the KDTree class from the
# scipy.spatial module.


DEBUG = False  # Set true to print debugging information
CB_AS_BACKBONE = False  # Set true to treat beta carbon atoms as backbone atoms


allIDs = ["ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET", "SER", "THR", "LYS", "ARG", "HID", "HIE", "HIP", "ASP", "GLU", "ASN", "GLN"]


if CB_AS_BACKBONE:
    bbAtomsList = ['CA', 'HA', 'N', 'HN', 'H', 'C', 'O', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3', 'OXT', 'CB']
else:
    bbAtomsList = ['CA', 'HA', 'N', 'HN', 'H', 'C', 'O', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3', 'OXT']


# These five function convert among the residue formats for different software.
# In particular, Osprey names residues as <chain ID><residue number>. These
# three functions

def ospreyToPDB(model, ospreyRes):
    return model[ospreyRes[0]][int(ospreyRes[1:])]

def ospreyFormat(res):
    return res.get_parent().get_id() + str(res.get_id()[1])

def residuesToOsprey(ress):
    assert [ospreyFormat(res) for res in ress] == list(map(ospreyFormat, ress))
    return list(map(ospreyFormat, ress))
    
def setToProbe(resSet):
    return ','.join([str(res.get_id()[1]) for res in resSet])

def setToPyMol(resSet):
    return 'resi ' + '+'.join([str(res.get_id()[1]) for res in resSet])


# When working across multiple structures with different residue numberings,
# it's helpful to have a "canonical" numbering to work with; these three
# functions help with that.

class CanonicalTools:
    def __init__(self, renumfn):
        self.addForCanonical = {}
        with Path(renumfn).open(mode='r') as renumFile:
            renumReader = csv.DictReader(renumFile)
            for row in renumReader:
                self.addForCanonical[row['verWord']] = row

    def addForCanonLkup(self, verWord, chainID):
        byChainID = self.addForCanonical[verWord]
        string = self.addForCanonical[verWord].get(chainID, '0')
        if string == '':
            return 0
        else:
            return int(string)
            
    def resToCanonicalID(self, res, verWord):
        chainID = res.get_parent().get_id()
        canonicalNumber = res.get_id()[1] + self.addForCanonLkup(verWord, chainID)
        return chainID + str(canonicalNumber)

    def canonicalIDToRes(self, canonicalID, model, verWord):
        chainID = canonicalID[0]
        canonicalResNum = int(canonicalID[1:])
        offset = self.addForCanonLkup(verWord, chainID)
        origResNum = canonicalResNum - offset
        return model[chainID][origResNum]

    def canonicalNumToNum(self, canonicalResNum, chainID, verWord):
        return canonicalResNum - self.addForCanonLkup(verWord, chainID)

    def numToCanonicalNum(self, resNum, chainID, verWord):
        return resNum + self.addForCanonLkup(verWord, chainID)

    def resIDToCanonicalID(self, resID, verWord):
        chainID = resID[0]
        resNum = int(resID[1:])
        return chainID + str(self.numToCanonicalNum(resNum, chainID, verWord))

    def canonicalIDToResID(self, canonicalID, verWord):
        chainID = canonicalID[0]
        canonicalNum = int(canonicalID[1:])
        return chainID + str(self.canonicalNumToNum(canonicalNum, chainID, verWord))


# These two classes are just data structures for convenience

class ConfSpaces:
    def __init__(self, name, csDesign, csBinder, csComplex):
        self.name = name
        self.csDesign = csDesign
        self.csBinder = csBinder
        self.csComplex = csComplex

class MutableSetInfo:
    def __init__(self, confSpaces, allMut, allAla, allFlex, bbFlex, tr):
        self.confSpaces = confSpaces
        self.allMut = allMut
        self.allAla = allAla
        self.allFlex = allFlex
        self.bbFlex = bbFlex
        self.tr = tr
    def printinfo(self):
        print('Name: ', self.confSpaces.name)
        print('Mutable residues: ', ' '.join(residuesToOsprey(self.allMut)))
        if self.allAla:
            print('Residues mutated to alanine: ', ' '.join(residuesToOsprey(self.allAla)))
        print('Other flexible residues: ', ' '.join(residuesToOsprey(self.allFlex)), f'(count: {len(self.allFlex)})')
        print('Backbone flexibility: ', ' '.join([f'({ospreyFormat(first)}, {ospreyFormat(second)})' for first, second in self.bbFlex]))
        print('Overall residue count:', len(self.allMut) + len(self.allAla) + len(self.allFlex))
        print('Chains with translation and rotation:', ' '.join(self.tr))


def findResSetName(resSet):
    name = str()
    namelist = []
    for res in resSet:
        namelist.append(res.get_parent().get_id() + str(res.get_id()[1]))
    return '_'.join(sorted(namelist))


# More than one other module needs to put the raw CLI argument list for maxFlex
# into a form more useful for confTools. But note this is still a hack for
# point mutations; more work is needed to specify non-singleton sets.
def maxFlexFromArg(maxFlexArg):
    if len(maxFlexArg) == 1 and maxFlexArg[0].split(',')[0] == 'all':
        maxFlexAll = int(maxFlexArg[0].split(',')[1])
        maxFlexMap = None
        print(f'No more than {maxFlexAll} flexible residues will be allowed in findFlexible.')
    else:
        maxFlexAll = None
        maxFlexMap = { frozenset([ospreyToPDB(model, res)]):int(maxFlex) for (res, maxFlex) in [tuple(pairString.split(',')) for pairString in maxFlexArg] }
    return maxFlexMap, maxFlexAll


# Given residue X and residue Y, compute the minimum among any distance from a
# side chain atom in X to a side chain atom in Y.
def SCToSCMinDist(resMut, resFix):
    minDist = math.inf
    for a in (a for a in resFix if a.get_id() not in bbAtomsList):
        for b in (b for b in resMut if b.get_id() not in bbAtomsList):
            dist = a - b
            if dist < minDist:
                minDist = dist
    return minDist


# Given residue X and residue Y, return True iff there is any side chain atom
# in X that is within cutoff Angstroms of any side chain atom in Y. The
# exception is that if (X, Y) or (Y, X) is found in the notNeighbors list of
# pairs, return False regardless.
def SCToSC(resMut, resFix, cutoff, notNeighbors):
    if (resMut, resFix) in notNeighbors or (resFix, resMut) in notNeighbors:
        return False

    for a in (a for a in resFix if a.get_id() not in bbAtomsList):
        for b in (b for b in resMut if b.get_id() not in bbAtomsList):
            if a - b <= cutoff:
                if DEBUG:
                    print(f'{ospreyFormat(resMut)} has a side chain within {cutoff} of side chain of {ospreyFormat(resFix)} ({b.get_id()} to {a.get_id()}, {a - b})')
                return True
    return False


# Find the set of residues in chainListDesign that have side-chain atoms within
# cutoff Angstroms of any atom in chainListBinder, so long as it is not glycine
# or proline.
def findWantMutable(chainListDesign, chainListBinder, cutoff, onlyMutate):

    # Similar to SCToSC, but we are interested in distances from any side-chain
    # atom in resMut to ANY atom in chain.
    def SCToAnyInChain(resMut, chain, cutoff):
        for resFix in chain:
            for a in resFix:
                for b in (b for b in resMut if b.get_id() not in bbAtomsList):
                    if a - b <= cutoff:
                        if DEBUG:
                            print(f'{ospreyFormat(resMut)} has a side chain within {cutoff} of {ospreyFormat(resFix)} ({b.get_id()} to {a.get_id()}, {a - b})')
                            print(f'{ospreyFormat(resMut)} has a side chain within {cutoff} of {ospreyFormat(res)}')
                        return True
        return False

    wantMutable = set()
    for chainDesign in chainListDesign:
        for chainBinder in chainListBinder:
            for resDesign in chainDesign:
                if (
                        (not onlyMutate or resDesign in onlyMutate) and
                        resDesign.get_resname() not in {"GLY", "PRO"} and 
                        SCToAnyInChain(resDesign, chainBinder, cutoff)
                        ):
                    wantMutable.add(resDesign)
    return wantMutable


# Return a set of sets; each inner set contains positions in which we are
# interested in potentially simulataneous mutations, for purposes of a
# mutational scan.  First find the set of all residues for which we wish to
# consider mutations.  Then, for each residue in this set, constuct a set
# containing it and each other want-mutable residue with a side-chain atom
# within cutoffNeighbors Angstroms of a side-chain atom in the focal residue.
# If this construction results in more than maxSimultaneous elements, order the
# set according to each residue's minimum side-chain-atom-to-side-chain-atom
# distance with the focal residue, and remove all but the first maxSimultaneous
# of this list (including the focal residue).  Return the set of all such
# unique sets.
def findMutableSets(chainListDesign, chainListBinder, cutoffWantMutable, cutoffNeighbors, notNeighbors, onlyMutate, maxSimultaneous):
    mutableSets = set()
    wantMutable = findWantMutable(chainListDesign, chainListBinder, cutoffWantMutable, onlyMutate)
    print('All residues we want to mutate: ', setToPyMol(wantMutable))
    for res in wantMutable:
        mutableSetDists = {res: 0.0}
        for other in wantMutable - {res}:
            if SCToSC(res, other, cutoffNeighbors, notNeighbors):
                mutableSetDists[other] = SCToSCMinDist(res, other)
        mutableSets.add(frozenset(sorted(mutableSetDists.keys(), key = lambda e: mutableSetDists[e])[:maxSimultaneous]))
    return mutableSets


# For a given set of mutable residues, find the set of all other residues, on
# any chain, that have side-chain-atom-to-side-chain-atom distances within
# cutoff of some mutable residue in the set. Except that if a mutable set is
# found in the maxFlexMap dictionary, only accept the closest N such residues,
# where N is the corresponding key in the dictionary.
def findFlexible(mutableSet, chainList, cutoff, notNeighbors, maxFlexMap, maxFlexAll, mustBeFlex):

    def SCToSCSet(res, resSet, cutoff, notNeighbors):
        for resFromSet in resSet:
            if SCToSC(res, resFromSet, cutoff, notNeighbors):
                return True
        return False

    flexDists = {}
    for chain in chainList:
        for res in chain:
            if res not in mutableSet:
                if SCToSCSet(res, mutableSet, cutoff, notNeighbors):
                    minDist = min(SCToSCMinDist(res, resMut) for resMut in mutableSet)
                    flexDists[res] = minDist
    foundFlex = set(flexDists.keys())
    allFlex = foundFlex | mustBeFlex
    #print('allFlex: ', [ospreyFormat(r) for r in allFlex])
    ordered = sorted(foundFlex - mustBeFlex, key = lambda e: flexDists[e])
    if maxFlexAll or frozenset(mutableSet) in maxFlexMap:
        maxFlex = maxFlexAll if maxFlexAll else maxFlexMap[frozenset(mutableSet)]
        nToRemove = max(0, len(allFlex) - maxFlex)
    else:
        maxFlex = None
        nToRemove = 0

    nToKeep = len(ordered) - nToRemove
    if nToKeep < 0:
        print(f'Warning: Cannot remove more than {len(ordered)} flexible residues')
        nToKeep = 0

    if nToRemove > 0:
        print(f'For mutable set {{{" ".join(residuesToOsprey(mutableSet))}}}, {nToRemove} residues not being set flexible, despite being within {cutoff} angstroms:')
        for res in ordered[nToKeep:]:
            print(f'    Residue {ospreyFormat(res)} is within {flexDists[res]}.')
        print('Residues set flexible, and distances:\n    ', end="")
        for res in ordered[:nToKeep]:
            print(f'{ospreyFormat(res)}({flexDists[res]:.2f}); ', end="")
        print()
    return set(ordered[:nToKeep]), nToRemove


# Given a set of residues, return a dictionary whose keys are the chains
# involved and whose values are the subsets of the set that are part of a given
# chain.
def splitResiduesByChainID(residues):
    byChainID = defaultdict(set)
    for res in residues:
        byChainID[res.get_parent().get_id()].add(res)
    return byChainID


# Given a set of mutable residues, return a set of pairs where each pair
# represents the staring and ending residues of a segment for which we want
# backbone flexibility. These segments consist of any consecutive sequence of
# mutable residues, including their neighbors (by chain position) on either
# side, if such exist, are not N-termini, and are not proline. Do not include
# any segments shorter than three residues.
def findBBFlexPairsPerChain(mutableSet, chain, isCATS=False):

    def addBBFlexPair(start, end):
        # Assumes that start and end are potentially valid members of a flexible backbone segment
        
        def bbFlexCheck(index):
            #return index in chain and (not isCATS or chain[index].get_resname() != "PRO")
            return (index in chain 
                    and index - 1 in chain
                    and chain[index].get_resname() != "PRO")

        startIndex = start.get_id()[1]
        endIndex = end.get_id()[1]
        chain = start.get_parent()
        strandStartIndex = startIndex - 1 if bbFlexCheck(startIndex - 1) else startIndex
        strandEndIndex = endIndex + 1 if bbFlexCheck(endIndex + 1) else endIndex
        if strandEndIndex - strandStartIndex > 1:
            pairs.add((chain[strandStartIndex], chain[strandEndIndex]))
            #print('adding:')
            #print(f'    chain: {chain.get_id()}; res: {strandStartIndex, strandEndIndex}')
            #print(f'pairs is now:')
            #for pair in pairs:
            #    print(f'    {pair[0].get_id()[1]}, {pair[1].get_id()[1]}')

    byChainID = splitResiduesByChainID(mutableSet)
    #print(f'finding bb flex pairs, chain {chain.get_id()}')
    resSet = byChainID[chain.get_id()]
    pairs = set()
    start = None
    prev = None
    for res in sorted(resSet, key = lambda res : res.get_id()[1]):
        #print(f'considering res {res.get_id()[1]}')
        if start is None:
            start = res
        elif res.get_id()[1] != prev.get_id()[1] + 1:
            addBBFlexPair(start, prev)
            start = res
        prev = res
    if (start is not None and prev is not None):
        addBBFlexPair(start, prev)
    return pairs


# Write a .pert file in the format required by Osprey's DEEPer-related classes.
def writePertFile(model, mutableSet, chainList, chain, bbFlexPairs, flexSetByChainID, nRemoved):

    # Format for specifying residues in DEEPer is just numbers
    def resPairsToNumPairs(resPairs):
        return {(first.get_id()[1], last.get_id()[1]) for first, last in resPairs}

    def countBackrubs():
        return sum([last - first - 1 for first, last in bbFlexPairsNum])

    def countResidues():
        return sum([last - first + 1 for first, last in bbFlexPairsNum])

    if not bbFlexPairs:
        #print(f'No backrub pairs for chain list {", ".join(c.get_id() for c in chainList)}; chain {chain.get_id()}.')
        return None

    bbFlexPairsNum = resPairsToNumPairs(bbFlexPairs)

    residuesString = '_'.join([ospreyFormat(res) for res in mutableSet])
    nRemovedTag = '-rem' + str(nRemoved) if nRemoved > 0 else ''
    pertfilename = f'{residuesString}-{"".join(c.get_id() for c in chainList)}-{chain.get_id()}{nRemovedTag}-DEEPer.pert'
    pertfilepath = Path(pertfilename)
    if not pertfilepath.exists():
        with pertfilepath.open(mode='w') as pertfile:
            pertfile.write('PERTURBATIONS\n')
            pertfile.write(str(countBackrubs()) + '\n')
            pertIndex = 0
            pertSetByResByChain = defaultdict(dict)
            for first, last in bbFlexPairsNum:
                pertSet = set()
                for startIndex in range(first, last - 1):
                    pertfile.write('BACKRUB\n')
                    pertfile.write(' '.join([str(i) for i in range(startIndex, startIndex + 3)]) + '\n')
                    pertfile.write('1 states\n-2.5 2.5\n')
                    pertSet.add(pertIndex)
                    pertIndex += 1
                for resi in range(first + 1, last): # Not including endpoints; doesn't seem to be necessary
                    pertSetByResByChain[chain][resi] = pertSet

            # Turns out we need to go through other flexible residues and specify
            # which perturbations they belong to (that is, none of them).  Mutable
            # residues are by definition flexible, and we need to include those not
            # on present chain
            mutByChainID = splitResiduesByChainID(mutableSet)
            for c in chainList:
                flexSetForChain = flexSetByChainID[c.get_id()] | mutByChainID[c.get_id()]
                for res in flexSetForChain:
                    resi = res.get_id()[1]
                    if resi not in pertSetByResByChain[c]:
                        pertSetByResByChain[c][resi] = set()

            # It is crucial that the residues be given in the order of the internal
            # "positions" list, which means by chain according to the order those
            # chains are passed into the conformation space factory, and
            # numerically within a chain.
            for c in chainList:
                pertSetByRes = pertSetByResByChain[c]
                for res, pertSet in sorted(pertSetByRes.items()):
                    #print(f'Residue number {res} has perturbation set: {pertSet}')
                    pertfile.write(f'RES\n{1 if pertSet else 0} states\nPERTURBATIONS ')
                    pertfile.write(' '.join([str(i) for i in pertSet]) + '\n')
                    if pertSet:
                        pertfile.write(' '.join(['0' for i in pertSet]) + '\n')
            #print("------")

    return pertfilename


# Return a MutableSetInfo object for a given set of mutable residues. The most
# important part of this object is its ConfSpaces object, which contains a
# conformation space for the designed chains, for the non-designed chains, and
# for the entire complex.
def findMutableSetInfo(name, identitiesByResidue, pdbfile, mol, templateLib, model, chainListDesign, chainListBinder, chainListTR, cutoff, notNeighbors, maxFlexMap, maxFlexAll, alaIfNotMut = set(), backrubs = False):
    
    mutableSet = set(identitiesByResidue.keys())
    mutByChainID = splitResiduesByChainID(mutableSet)

    if backrubs:
        bbFlexPairsByChain = {chain: findBBFlexPairsPerChain(mutableSet, chain) for chain in model}
        allBBFlex = set().union(*bbFlexPairsByChain.values())
    else:
        allBBFlex = set()

    #mustBeFlex = mutableSet | {res for res in pair for pair in allBBFlex}
    mustBeFlex = mutableSet | {res for pair in allBBFlex for res in pair}
    flexSet, nRemoved = findFlexible(mutableSet, list(model), cutoff, notNeighbors, maxFlexMap, maxFlexAll, mustBeFlex)
    flexByChainID = splitResiduesByChainID(flexSet)

    alaIfFlex = alaIfNotMut - mutableSet

    chainListComplex = chainListDesign + chainListBinder
    chainListList = [chainListDesign, chainListBinder, chainListComplex]

    if backrubs:

        # This loop looks for backrub endpoint residues that are neither
        # flexible nor mutable, and adds desired flexibility to them.
        # Surprisingly, the backrub still works even if they're not made
        # flexible at all, but some flexibility is probably good to allow them
        # to respond to the CA-CB vector rotation incurred by the backrub.
        #
        # It's okay to ignore alanine stuff here; residues mutated to alanine
        # will be a subset of those in flexSet.
        #
        # This will need to be slightly more complex if we decide we want the
        # backrub endpoints to have a voxel around their wildtype discrete
        # rotamer but not flexibility across other discrete rotamers.
        #
        for chain in model:
            cID = chain.get_id()
            anyFlex = mutByChainID[cID] | flexByChainID[cID]
            for pair in bbFlexPairsByChain[chain]:
                for res in pair: # That is, for each of the two endpoint residues of the backrub
                    if res not in anyFlex:
                        flexByChainID[cID].add(res)

        pertfilenames = {}
        for chainList in chainListList:
            pertfilenames[tuple(chainList)] = {}
            for chain in chainList:
                pertfilename = writePertFile(model, mutableSet, chainList, chain, bbFlexPairsByChain[chain], flexByChainID, nRemoved)
                pertfilenames[tuple(chainList)][chain] = pertfilename

    allAla = set()
    allFlex = set()
    strandListByChainTupleByChainID = defaultdict(dict)
    tr = []
    for chain in model:
        firstRes = min(chain, key = lambda res : res.get_id()[1])
        lastRes = max(chain, key = lambda res : res.get_id()[1])
        strand = osprey.Strand(mol, templateLib=templateLib, residues=[ospreyFormat(firstRes), ospreyFormat(lastRes)])

        for res in mutByChainID[chain.get_id()]:
            strand.flexibility[ospreyFormat(res)].setLibraryRotamers(osprey.WILD_TYPE, *identitiesByResidue[res]).addWildTypeRotamers().setContinuous()

        flexForChain = flexByChainID[chain.get_id()]
        alaForChain = flexForChain & alaIfFlex

        for res in alaForChain:
            strand.flexibility[ospreyFormat(res)].setLibraryRotamers('ALA').setContinuous()
            allAla.add(res)

        for res in flexForChain - alaForChain:
            strand.flexibility[ospreyFormat(res)].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
            allFlex.add(res)

        strandList = [strand]

        if chain in chainListTR:
            strandList.append(osprey.c.confspace.StrandFlex.TranslateRotate())
            tr.append(chain.get_id())

        for chainList in chainListList:
            if chain in chainList:
                strandListForChainList = strandList.copy()
                if backrubs and pertfilenames[tuple(chainList)][chain]:
                    deeperThing = osprey.DEEPerStrandFlex(strand, pertfilenames[tuple(chainList)][chain], [], pdbfile)
                    strandListForChainList.append(deeperThing)
                strandListByChainTupleByChainID[tuple(chainList)][chain.get_id()] = strandListForChainList

        # If we want CATS instead of DEEPer, it's far simpler.
        #for pair in bbFlexPairs:
        #    strandList.append(osprey.c.confspace.CATSStrandFlex(strand, ospreyFormat(pair[0]), ospreyFormat(pair[1])))

    confSpaces = ConfSpaces(
        name = name,
        csDesign = chainListToConfSpace(strandListByChainTupleByChainID, chainListDesign),
        csBinder = chainListToConfSpace(strandListByChainTupleByChainID, chainListBinder),
        csComplex = chainListToConfSpace(strandListByChainTupleByChainID, chainListComplex)
        )

    return MutableSetInfo(confSpaces, mutableSet, allAla, allFlex, allBBFlex, tr)


# Actually create Osprey ConfSpace objects corresponding to each chain in a
# given list.
def chainListToConfSpace(strandListByChainTupleByChainID, chainList):
    strandListList = []
    for chain in chainList:
        strandListList.append(strandListByChainTupleByChainID[tuple(chainList)][chain.get_id()])
    return osprey.ConfSpace(strandListList)


# Given a file containing a sequence space in a simple space-separated format,
# construct an Osprey ConfSpace object and return it in the form of a
# MutableSetInfo object.
def findConfSpaceFromSeqSpace(verWord, canonTools, mutationsFile, pdbfile, mol, templateLib, model, chainListDesign, chainListBinder, chainListTR, cutoff, notNeighbors, maxFlexMap = None, maxFlexAll = None, backrubs = False):

    #addForCanonical = findAddForCanonical(renumfn)

    identitiesByResidue = dict()
    for residueLine in mutationsFile:
        chainID = residueLine.split()[0]
        canonicalResNum = int(residueLine.split()[1])
        #resNum = canonicalResNum - int(addForCanonical[verWord][chainID])
        resNum = canonTools.canonicalNumToNum(canonicalResNum, chainID, verWord)
        identitiesByResidue[model[chainID][resNum]] = residueLine.split()[2:]
    
    return findMutableSetInfo(findResSetName(set(identitiesByResidue.keys())), identitiesByResidue, pdbfile, mol, templateLib, model, chainListDesign, chainListBinder, chainListTR, cutoff, notNeighbors, maxFlexMap, maxFlexAll, backrubs=backrubs)


# Prepare for a mutational scan by returning a set of MutableSetInfo objects,
# each containing an Osprey ConfSpace object defining the desired conformation
# space surrounding a given subset of mutable residues. Each such ConfSpace
# object is the basis for a separate Osprey search.
def findConfSpacesSet(pdbfile, mol, templateLib, model, chainListDesign, chainListBinder, chainListTR, cutoffWantMutable, cutoffNeighbors, notNeighbors, onlyMutate, maxSimultaneous = None, maxFlexMap = None, maxFlexAll = None, invAla2 = False, backrubs = False):
    
    mutableSets = findMutableSets(chainListDesign, chainListBinder, cutoffWantMutable, cutoffNeighbors, notNeighbors, onlyMutate, maxSimultaneous)
    alaIfNotMut = set().union(*mutableSets) if invAla2 else set()
    mutableSetInfos = set()
    for mutableSet in mutableSets:
        allIDsByResidue = {res: allIDs for res in mutableSet}
        mutableSetInfos.add(findMutableSetInfo(findResSetName(mutableSet), allIDsByResidue, pdbfile, mol, templateLib, model, chainListDesign, chainListBinder, chainListTR, cutoffNeighbors, notNeighbors, maxFlexMap, maxFlexAll, alaIfNotMut=alaIfNotMut, backrubs=backrubs))
    return mutableSetInfos
