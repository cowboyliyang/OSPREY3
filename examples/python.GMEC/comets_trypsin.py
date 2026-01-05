import osprey
osprey.start(heapSizeMiB=768000)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol1 = osprey.readPdb('3btd.df.min_strip_reduce_prep_rc_add_rc_steps100_rc.pdb')
mol2 = osprey.readPdb('3bte.df.min_strip_reduce_prep_rc_add_rc_steps100_rc.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol1, templateLib=templateLib, residues=['A226', 'A248'])
protein.flexibility['A227'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA').addWildTypeRotamers().setContinuous()
protein.flexibility['A228'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['A230'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['A233'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['A247'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# protein = osprey.Strand(mol1, templateLib=templateLib, residues=['A226', 'A248'])
# protein.flexibility['A227'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
#                                                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL').addWildTypeRotamers().setContinuous()
# protein.flexibility['A228'].setLibraryRotamers(osprey.WILD_TYPE,'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
# 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL').addWildTypeRotamers().setContinuous()
# protein.flexibility['A230'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
#                                                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL').addWildTypeRotamers().setContinuous()
# protein.flexibility['A233'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
#                                                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL').addWildTypeRotamers().setContinuous()
# protein.flexibility['A247'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
#                                                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL').addWildTypeRotamers().setContinuous()


protein2 = osprey.Strand(mol2, templateLib=templateLib, residues=['A170', 'A192'])
protein2.flexibility['A171'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA').addWildTypeRotamers().setContinuous()#asp
protein2.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein2.flexibility['A174'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein2.flexibility['A177'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein2.flexibility['A191'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous() #val



# define the ligand strand
ligand1 = osprey.Strand(mol1, templateLib=templateLib, residues=['I11', 'I14'])
ligand1.flexibility['I12'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand1.flexibility['I13'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

ligand2 = osprey.Strand(mol2, templateLib=templateLib, residues=['I234', 'I237'])
ligand2.flexibility['I235'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand2.flexibility['I236'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()


# define the COMETS states
bound1 = osprey.COMETS_State('Bound', osprey.ConfSpace([protein, ligand1]))
bound2 = osprey.COMETS_State('Bound', osprey.ConfSpace([protein2, ligand2]))
unbound1 = osprey.COMETS_State('Unbound', osprey.ConfSpace(protein))
unbound2 = osprey.COMETS_State('Unbound', osprey.ConfSpace(protein2))

# configure COMETS
comets = osprey.COMETS(
	objective = osprey.COMETS_LME({
		bound1: 1.0,
		bound2: -1.0,
		unbound1: -1.0,
        unbound2: 1.0
	}),
	logFile = 'comets.tsv'
)


# Configure energy calculations
conf_spaces = [state.confSpace for state in comets.states]
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(conf_spaces, ffparams, parallelism=parallelism)

# Configure each state
for state in comets.states:
    # Reference energies and energy calculator
    eref = osprey.ReferenceEnergies(state.confSpace, ecalc)
    conf_ecalc = osprey.ConfEnergyCalculator(
        state.confSpace,
        ecalc,
        referenceEnergies=eref
    )

    # Energy matrix configuration
    emat = osprey.EnergyMatrix(
        conf_ecalc,
        cacheFile=f'emat_{state.name}.dat'
    )

    # Configure search strategy
    def make_astar(rcs, emat=emat):
        return osprey.AStarTraditional(emat, rcs, showProgress=False)

    # Assign state properties
    state.confEcalc = conf_ecalc
    state.fragmentEnergies = emat
    state.confTreeFactory = osprey.COMETS_ConfSearchFactory(make_astar)

# Run COMETS optimization
comets.findBestSequences(5)  # Find top 5 sequences
print("COMETS design completed successfully!")