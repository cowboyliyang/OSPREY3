import osprey
osprey.start(heapSizeMiB=768000)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# 创建包含acetyl-lysine (ALY/K03) 和 succinyl-lysine (SLL/K04) 的自定义模板库
# 必须在创建Strand之前定义模板库
customizedTemplateLib = osprey.TemplateLibrary(
    extraTemplates=['K03_template.in', 'K04_template.in'],
    extraTemplateCoords=['K03_coords.in', 'K04_coords.in'],
    extraRotamers=['K03_rotlib.dat', 'K04_rotlib.dat']
)

# define the protein strand - 直接用PDB文件名，传入自定义模板库
protein = osprey.Strand('6dv2_strip_reduce_prep_rc_add_rc.pdb', templateLib=customizedTemplateLib, residues=['B734', 'B782'])
# B781位置: 允许野生型(LYS)、acetyl-lysine(ALY)、succinyl-lysine(SLL)
protein.flexibility['B781'].setLibraryRotamers(osprey.WILD_TYPE, 'ALY', 'SLL').addWildTypeRotamers().setContinuous()
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
parallelism = osprey.Parallelism(cpuCores=20)
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
