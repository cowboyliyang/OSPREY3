
import osprey

osprey.start()

# define a strand
# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol1 = osprey.readPdb('3btd_clean_no_dup_strip_reduce_prep_rc_add_rc.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol1, templateLib=templateLib, residues=['E168', 'E191'])
protein.flexibility['E169'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['E170'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['E172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['E175'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['E190'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make the conf space
confSpace = osprey.ConfSpace(protein)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confEcalc)
astar = osprey.AStarMPLP(emat, confSpace)

# find the best sequence and rotamers
gmec = osprey.GMECFinder(astar, confEcalc).find()

# write the rigid GMEC to a pdb
gmecStructure = confSpace.makeMolecule(gmec.getAssignments())
osprey.writePdb(gmecStructure,  "gmec_wild.pdb")