import osprey
osprey.start(heapSizeMiB=768000)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('6dv2_strip_reduce_prep_rc_add_rc.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['B734', 'B782'])
protein.flexibility['B781'].setLibraryRotamers(osprey.WILD_TYPE, 'LYS').addWildTypeRotamers().setContinuous()
protein.flexibility['B780'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['B782'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

protein.flexibility['B735'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['B776'].setLibraryRotamers(osprey.WILD_TYPE,).addWildTypeRotamers().setContinuous()

protein.flexibility['B734'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['B775'].setLibraryRotamers(osprey.WILD_TYPE,).addWildTypeRotamers().setContinuous()

protein.flexibility['B736'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['B777'].setLibraryRotamers(osprey.WILD_TYPE,).addWildTypeRotamers().setContinuous()

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# configure PAStE
paste = osprey.Paste(
    complexConfSpace,
    numPDBs = 15,
    epsilon=0.1, # you proabably want something more precise in your real designs
    maxNumPfConfs=50000,
    useWindowCriterion=True,
    writeSequencesToConsole=True,
    writeSequencesToFile='paste.results.tsv'
#     mutFile='mut.txt'
)

# configure PAStE inputs for the conf space

# how should we define energies of conformations?
eref = osprey.ReferenceEnergies(paste.protein.confSpace, ecalc)
paste.protein.confEcalc = osprey.ConfEnergyCalculator(paste.protein.confSpace, ecalc, referenceEnergies=eref)

# compute the energy matrix
emat = osprey.EnergyMatrix(paste.protein.confEcalc, cacheFile='emat.%s.dat' % paste.protein.id)

# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
def makeAStar(rcs, emat=emat):
    return osprey.AStarTraditional(emat, rcs, showProgress=False)
paste.protein.confSearchFactory = osprey.Paste.ConfSearchFactory(makeAStar)

# run PAStE
paste.run()
