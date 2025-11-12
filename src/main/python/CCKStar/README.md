# Convex Closure K*

CCK* is a suite of three algorithms (SCOPE, MONTAGE, ARISE) for the de novo design of L and D peptide binders.
Examples for running these algorithms can be found in main.py.

# SCOPE

The purpose of SCOPE (Side Chain Orientation and Position Evaluation) is to determine side chain contacts by configuration space (C-space) intersections based on rotameric geometry. This is useful for
prioritization of mutant and conformation modeling, such as setting flexible residues in OSPREY.

Workflow:
1. Based on the chirality, generate convex hulls for the design chain based on all rotamers for all specified amino acids
2. For the target (no mutations), generate convex hulls for the wildtype identity rotamers
3. Compute intra-chain (design chain residue intersections) and inter-chain (design residues intersecting with
   target residues) intersections
4. Return the residue C-space intersections

Terminology:
- Doublets: intra-chain intersections, returned as pairs 
- Islands: design chain residues with no intra-chain intersections. These may have inter-chain intersections.

Prerequisites:
1. A protonated PDB file containing only 2 chains (target and design ligand)
2. The python PyVista library

Variables:
- pdb_name: the input PDB 
- outfolder: where to save the convex hulls, written out in PDB format (manually created dir)
- designID: chain ID of the design chain 
- design_AA_type: the amino acids that will be used to generate the design chain convex hulls. Use 'wt' to include the wildtype identity.
- savePDB: boolean that determines if hulls are written to outfolder 
- design_chirality: the chirality of the design chain 
- fixed_identity: design chain residues indices that will NOT have a hull created from the identities specified in design_AA_type. Any residues in this list
will have only a convex hull comprised of the wildtype rotamers

# MONTAGE

The purpose of the MONTAGE (Motif-Oriented Noncanonical Template Assembly and Generation Engine) algorithm is to generate L or D polyalanine scaffolds for L targets.

Workflow:
1. Perform the require reflection operations and run MASTER
2. Mutate all design residues to ALA (except for GLY and PRO residues)
3. A VAL hull (approximately an ALA residue + VdW radius) is placed on each design residue
4. A wildtype hull (hull of WT identity rotamers) is placed on each L-target residue
5. Inter-chain hull intersection = flexibility on participating target residue
6. Run the K* algorithm with flexible target residues (from inter-chain overlaps) and translation+rotation of design chain
7. Scaffolds are ranked and returned by designability

Prerequisites:
- A MASTER database (see https://grigoryanlab.org/master/ for instructions)
- Filepaths in resources/db.txt must point to your MASTER database pds files
- Each input PDB must contain only two chains (target and design). The target must be listed first in the PDB (target chain ID (e.g., A) < design chain ID (e.g., B))
- Targets in the PDB must be L-space. Design chains can be L or D.
- The input directory can only contain PDB files with same ligand chirality workflow (e.g., L-peptide -> D-peptide)
- A conda env with OSPREY installed (required for K* searches). Specify the name in main.sh (placeholder is AmberTools23)
- The python PyVista library

Variables:
- input_pdb_directory = the directory containing all the two-chain PDBs 
- input_chirality = the chirality of the design chains in input_pdb_directory 
- output_chirality = the desired output chirality of the design chains 
- MASTER_matches = how many matches from MASTER to use to generate scaffolds 
- max_flex = the maximum number of target residues set as flexible during the K* search. If the number of inter-chain intersects > max_flex, residue flexibility will be selected in order of hull volume overlap.

# ARISE

The purpose of ARISE (Affinity-driven Rational Iteration for Sequence Engineering) is to generate high-affinity L or D peptides for L-targets.

Workflow:
1. Construct an undirected graph of potential chemical contacts using SCOPE
2. Based on this graph, use K* to mutate 2 connected design chain residues (400 sequences) with flexibility on all neighboring target residues
3. Get the GMEC PDB of the highest-affinity pair
4. Recompute the edges of the graph using SCOPE, and repeat the mutant assignment search until a full design chain sequence is found

Prerequisites:
- An input directory named MONTAGE_GMEC containing two-chain PDBs (see MONTAGE above)
- A conda env with OSPREY installed (required for K* searches). Specify the name in main.sh (placeholder is AmberTools23)
- The python PyVista library

Variables:
- round_number: the round you are starting on (if 1 assumes input dir is MONTAGE_GMEC)
- length_chain: the length of the design chain 
- finished_matches: a set of matches that are fully assigned. Pass an empty set to start from round 1. 
- visited_doubs: a dictionary of assigned doublet (pairs of design chain residues). For maintaining gly/pro (starting from round 1), pass a dictionary with all the gly/pro reside indexes (can be computed using maintain_gly_pro). If starting from round > 1, paste the dict from the previous round here.
- designID: the chain ID of the design chain 
- targetID: the chain ID of the target chain 
- final_designs_outfolder: where to save fully assigned designs (will save design in complex with the target)
- apo_tolerance: delta of apo wildtype partition function to any apo mutant partition function. Set close to 0 to prevent improving predicted binding (K* ratio) driven by destabilizing the apo ligand. 
- design_chirality: chirality of the design chain (D or L)