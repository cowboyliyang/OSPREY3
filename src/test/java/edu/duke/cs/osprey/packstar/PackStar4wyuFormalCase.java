package edu.duke.cs.osprey.packstar;

import cern.colt.matrix.DoubleMatrix1D;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Experiment harness for the 4WYU A+D PACK* P_C case (not production feature code).
 *
 * <p>The runner intentionally owns the complete nine-sequence case instead of
 * delegating sequence generation to the general design command.  This keeps
 * every requested sequence explicit, records the three state partition
 * functions, and makes it possible to compare the sampled result with an
 * independent full enumeration when the RC count is within the frozen limit.</p>
 */
public final class PackStar4wyuFormalCase {

    private static final String EVENT_NAME = "pdz34_two_pocket_contacts_v1";
    private static final double FIRST_MAX_DISTANCE = 3.5;
    private static final double SECOND_MAX_DISTANCE = 4.0;
    private static final double[] THRESHOLD_DELTAS = {-0.25, 0.0, 0.25};
    private static final int MAX_REFERENCE_ASSIGNMENTS = 65_536;
    private static final double RT = BoltzmannCalculator.RClassic
            * BoltzmannCalculator.TClassic;

    private PackStar4wyuFormalCase() {
    }

    public static void main(String[] args) throws Exception {
        Path proteinPdb = requiredPath("pc4wyu.proteinPdb");
        Path ligandPdb = requiredPath("pc4wyu.ligandPdb");
        Path outputDir = requiredPath("pc4wyu.output");
        Files.createDirectories(outputDir);

        int threads = integerProperty("pc4wyu.threads", 1);
        int maxNumConfs = integerProperty("pc4wyu.maxNumConfs", Integer.MAX_VALUE);
        double targetEpsilon = doubleProperty("pc4wyu.epsilon", 0.1);
        validateProtocol(requiredPath("pc4wyu.protocol"), targetEpsilon, threads);
        if (threads < 1) {
            throw new IllegalArgumentException("pc4wyu.threads must be positive");
        }
        if (!(targetEpsilon > 0.0) || !Double.isFinite(targetEpsilon)) {
            throw new IllegalArgumentException("pc4wyu.epsilon must be finite and positive");
        }

        CaseSpaces spaces = makeSpaces(proteinPdb, ligandPdb);
        List<Sequence> sequences = spaces.complex.seqSpace().getSequences(2);
        if (sequences.size() != 9) {
            throw new IllegalStateException(
                    "expected nine explicit A25/A47 sequences, got " + sequences.size());
        }

        List<String> stateRows = new ArrayList<>();
        for (Sequence sequence : sequences) {
            stateRows.add(sequence.toString() + "\t"
                    + count(sequence, spaces.protein) + "\t"
                    + count(sequence, spaces.ligand) + "\t"
                    + count(sequence, spaces.complex));
        }
        writeLines(outputDir.resolve("sequence_rc_counts.tsv"),
                List.of("sequence\tprotein_rcs\tligand_rcs\tcomplex_rcs"),
                stateRows);
        enforceReferenceLimit(spaces, sequences);

        writeManifest(outputDir, proteinPdb, ligandPdb, targetEpsilon,
                threads, maxNumConfs, sequences, spaces);
        if (Boolean.getBoolean("pc4wyu.preflight")) {
            System.out.println("[4WYU-PC] PREFLIGHT PASS: input, protocol, nine sequences and RC limits");
            return;
        }

        Parallelism parallelism = Parallelism.makeCpu(threads);
        ForcefieldParams forcefieldParams = new ForcefieldParams();
        try (EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(
                spaces.complex, forcefieldParams)
                .setType(EnergyCalculator.Type.Cpu)
                .setIsMinimizing(true)
                .setParallelism(parallelism)
                .build();
             EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
                     .setIsMinimizing(false)
                     .build()) {

            StateModel protein = makeStateModel(
                    spaces.protein, minimizingEcalc, rigidEcalc, parallelism, "protein");
            StateModel ligand = makeStateModel(
                    spaces.ligand, minimizingEcalc, rigidEcalc, parallelism, "ligand");
            StateModel complex = makeStateModel(
                    spaces.complex, minimizingEcalc, rigidEcalc, parallelism, "complex");

            List<Row> rows = new ArrayList<>();
            Computed ligandComputed = null;
            ReferenceSummary ligandReference = null;
            try (var contexts = minimizingEcalc.tasks.contextGroup()) {
                for (Sequence sequence : sequences) {
                    System.out.println("[4WYU-PC] sequence=" + sequence);
                    Computed proteinComputed = computeState(
                            protein, sequence.filter(spaces.protein.seqSpace()),
                            null, targetEpsilon, maxNumConfs, contexts);
                    if (ligandComputed == null) {
                        ligandComputed = computeState(
                            ligand, sequence.filter(spaces.ligand.seqSpace()),
                            null, targetEpsilon, maxNumConfs, contexts);
                    }
                    PackStarFunctionalEvent event = makeEvent(0.0);
                    Computed complexComputed = computeState(
                            complex, sequence.filter(spaces.complex.seqSpace()),
                            event, targetEpsilon, maxNumConfs, contexts);

                    minimizingEcalc.tasks.waitForFinish();
                    ReferenceSummary proteinReference = enumerate(
                            protein, sequence.filter(spaces.protein.seqSpace()), false);
                    if (ligandReference == null) {
                        ligandReference = enumerate(
                            ligand, sequence.filter(spaces.ligand.seqSpace()), false);
                    }
                    ReferenceSummary complexReference = enumerate(
                            complex, sequence.filter(spaces.complex.seqSpace()), true);

                    KStarScore kstar = new KStarScore(
                            proteinComputed.result,
                            ligandComputed.result,
                            complexComputed.result);
                    rows.add(new Row(
                            sequence.toString(),
                            proteinComputed, ligandComputed, complexComputed,
                            proteinReference, ligandReference, complexReference,
                            kstar));
                    writeCheckpoint(outputDir.resolve("completed_sequences.tsv"), rows);
                }
            }

            writeResults(outputDir.resolve("results.tsv"), rows);
            writeSelections(outputDir.resolve("selections.tsv"), rows);
            writeValidation(outputDir.resolve("validation.tsv"), rows);
        }

        Files.writeString(outputDir.resolve("COMPLETE"),
                "Execution completed at " + Instant.now()
                        + "; scientific diagnostics are in validation.tsv (not a validation PASS)\n",
                StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
    }

    /** This fixed case fails closed if its recorded protocol drifts from the harness. */
    static void validateProtocol(Path path, double epsilon, int threads) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        JsonNode config = mapper.readTree(path.toFile());
        JsonNode expected = mapper.readTree("""
            {
              "chains":{"protein":"A","peptide":"D","excluded":["B","C"]},
              "protein_design_positions":{"A25":["SER","THR","ALA"],"A47":["SER","THR","ALA"]},
              "protein_flexible_wild_type_positions":["A79","A167"],
              "peptide_flexible_wild_type_positions":["D205","D206"],
              "sequence_count":9,"max_simultaneous_mutations":2,
              "rotamer_policy":"A167_D205_native_rotamer_only_continuous; other_flexible_positions_library_plus_wild_type_continuous",
              "interaction_graph":"complete_flexible_interactions_with_static_shell; shell_shell_constant_omitted",
              "event":{"name":"pdz34_two_pocket_contacts_v1",
                "semantics":"polar_atom_distance_surrogate; not a hydrogen_bond or activity predicate",
                "conditions":[
                  {"residue_1":"A79","atom_1":"NE2","residue_2":"D206","atom_2":"OG1","max_distance_A":3.5},
                  {"residue_1":"A167","atoms_1":["OE1","OE2"],"residue_2":"D205","atom_2":"NE2","max_distance_A":4.0}],
                "logical_operator":"AND","lower_distance_bound_A":0.0,"threshold_deltas_A":[-0.25,0.0,0.25]},
              "thermal_model":{"temperature_K":298.15,"gas_constant_kcal_per_mol_K":0.0019891},
              "reference":{"measure":"all_unpruned_RC_assignments","max_assignments_per_state_and_sequence":65536,
                "energy":"full_minimized_conf_energy",
                "geometry":"independent_coordinate_distance_evaluation_after_minimization",
                "zero_weight_assignments":"positive_infinity_included_as_zero_weight; NaN_and_negative_infinity_fatal"},
              "packstar":{"target_epsilon":0.1,"samples":1000,"confidence_delta":0.05,
                "random_seed":42,"sampling_gpu":false,"sampling_threads":32}
            }
            """);
        var fields = expected.fields();
        while (fields.hasNext()) {
            var field = fields.next();
            if (!field.getValue().equals(config.get(field.getKey()))) {
                throw new IllegalArgumentException("protocol mismatch: " + field.getKey());
            }
        }
        validateRuntimeSettings(epsilon, threads);
    }

    static void validateRuntimeSettings(double epsilon, int threads) {
        if (epsilon != 0.1 || threads != 32
                || integerProperty("packstar.pac.samples", 1000) != 1000
                || integerProperty("packstar.pac.randomSeed", 42) != 42
                || doubleProperty("packstar.pac.confidence", 0.05) != 0.05
                || Boolean.parseBoolean(System.getProperty("packstar.pac.sampling.gpu", "false"))
                || BoltzmannCalculator.TClassic != 298.15
                || BoltzmannCalculator.RClassic != 0.0019891) {
            throw new IllegalArgumentException("runtime settings differ from frozen case protocol");
        }
    }

    private static CaseSpaces makeSpaces(Path proteinPdb, Path ligandPdb) {
        ForcefieldParams forcefieldParams = new ForcefieldParams();
        Molecule proteinMolecule = PDBIO.readFile(proteinPdb.toString());
        Molecule ligandMolecule = PDBIO.readFile(ligandPdb.toString());
        if (proteinMolecule.residues.isEmpty() || ligandMolecule.residues.isEmpty()) {
            throw new IllegalArgumentException(
                    "4WYU input contains no residues: protein=" + proteinPdb
                            + ", ligand=" + ligandPdb);
        }

        Strand proteinStrand = makeProteinStrand(proteinMolecule, forcefieldParams);
        Strand ligandStrand = makeLigandStrand(ligandMolecule, forcefieldParams);

        SimpleConfSpace protein = new SimpleConfSpace.Builder()
                .addStrand(proteinStrand)
                .build();
        SimpleConfSpace ligand = new SimpleConfSpace.Builder()
                .addStrand(ligandStrand)
                .build();
        SimpleConfSpace complex = new SimpleConfSpace.Builder()
                .addStrand(proteinStrand)
                .addStrand(ligandStrand)
                .build();

        requireWildType(protein, "A25", "SER");
        requireWildType(protein, "A47", "SER");
        requireWildType(protein, "A79", "HID");
        requireWildType(protein, "A167", "GLU");
        requireWildType(ligand, "D205", "GLN");
        requireWildType(ligand, "D206", "THR");
        return new CaseSpaces(protein, ligand, complex);
    }

    private static Strand makeProteinStrand(Molecule molecule,
                                             ForcefieldParams forcefieldParams) {
        ResidueTemplateLibrary templates = new ResidueTemplateLibrary.Builder(
                forcefieldParams.forcefld).build();
        Strand.Builder builder = new Strand.Builder(molecule)
                .setTemplateLibrary(templates)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames);
        builder.setResidueMutability("A25", List.of("SER", "THR", "ALA"), true, true);
        builder.setResidueMutability("A47", List.of("SER", "THR", "ALA"), true, true);
        builder.setResidueMutability("A79", List.of(), true, true);
        builder.setResidueMutability("A167", List.of(), true, true);
        Strand strand = builder.build();
        strand.flexibility.get("A167").setNoRotamers().addWildTypeRotamers().setContinuous();
        return strand;
    }

    private static Strand makeLigandStrand(Molecule molecule,
                                            ForcefieldParams forcefieldParams) {
        ResidueTemplateLibrary templates = new ResidueTemplateLibrary.Builder(
                forcefieldParams.forcefld).build();
        Strand.Builder builder = new Strand.Builder(molecule)
                .setTemplateLibrary(templates)
                .setTemplateMatchingMethod(Residue.TemplateMatchingMethod.AtomNames);
        builder.setResidueMutability("D205", List.of(), true, true);
        builder.setResidueMutability("D206", List.of(), true, true);
        Strand strand = builder.build();
        strand.flexibility.get("D205").setNoRotamers().addWildTypeRotamers().setContinuous();
        return strand;
    }

    private static void requireWildType(SimpleConfSpace space,
                                         String residueNumber,
                                         String expected) {
        SimpleConfSpace.Position position = space.positions.stream()
                .filter(pos -> pos.resNum.equals(residueNumber))
                .findFirst()
                .orElseThrow(() -> new IllegalStateException(
                        "missing configured position " + residueNumber));
        String actual = position.resFlex.wildType;
        if (!expected.equalsIgnoreCase(actual)) {
            throw new IllegalStateException("wild type mismatch at " + residueNumber
                    + ": expected " + expected + ", got " + actual);
        }
    }

    private static StateModel makeStateModel(SimpleConfSpace space,
                                             EnergyCalculator minimizingEcalc,
                                             EnergyCalculator rigidEcalc,
                                             Parallelism parallelism,
                                             String name) {
        SimpleReferenceEnergies references = new SimpleReferenceEnergies.Builder(
                space, minimizingEcalc).build();
        ConfEnergyCalculator minimizingConfEcalc = new ConfEnergyCalculator.Builder(
                space, minimizingEcalc)
                .setReferenceEnergies(references)
                .build();
        ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator(
                minimizingConfEcalc, rigidEcalc);
        EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(
                rigidConfEcalc).build().calcEnergyMatrix();
        EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(
                minimizingConfEcalc).build().calcEnergyMatrix();
        UpdatingEnergyMatrix corrections = new UpdatingEnergyMatrix(
                space, minimizingEmat, minimizingConfEcalc);
        return new StateModel(space, minimizingConfEcalc, rigidEmat,
                minimizingEmat, corrections, parallelism, name);
    }

    private static Computed computeState(StateModel model,
                                         Sequence sequence,
                                         PackStarFunctionalEvent event,
                                         double targetEpsilon,
                                         int maxNumConfs,
                                         edu.duke.cs.osprey.parallelism.TaskExecutor.ContextGroup contexts) {
        RCs rcs = sequence.makeRCs(model.space);
        long start = System.nanoTime();
        PackStarFunctionalObservableResult observable;
        PartitionFunction.Result result;
        try (PackStarPartitionFunction pfunc = new PackStarPartitionFunction(
                model.space, model.rigidEmat, model.minimizingEmat,
                model.minimizingConfEcalc, rcs, model.parallelism,
                "4wyu-formal-" + model.name)) {
            pfunc.setCorrections(model.corrections);
            pfunc.setReduceMinimizations(true);
            pfunc.setCorrectionTighteningEnabled(true);
            pfunc.setInstanceId(instanceId(model.name));
            if (event != null) {
                pfunc.setFunctionalEvent(EVENT_NAME, event);
            }
            pfunc.init(targetEpsilon);
            pfunc.putTaskContexts(contexts);
            pfunc.compute(maxNumConfs);
            result = pfunc.makeResult();
            observable = ((PackStarResult) result).getFunctionalObservableResult();
        }
        double seconds = (System.nanoTime() - start) / 1.0e9;
        return new Computed(result, observable, seconds, rcs.getNumConformations().intValueExact());
    }

    private static int instanceId(String name) {
        return switch (name) {
            case "protein" -> 0;
            case "ligand" -> 1;
            case "complex" -> 2;
            default -> throw new IllegalArgumentException("unknown state " + name);
        };
    }

    private static ReferenceSummary enumerate(StateModel model,
                                              Sequence sequence,
                                              boolean withEvent) {
        RCs rcs = sequence.makeRCs(model.space);
        int count = rcs.getNumConformations().intValueExact();
        ReferenceAccumulator accumulator = new ReferenceAccumulator(
                model.minimizingConfEcalc, rcs, withEvent, count);
        accumulator.visit(0, new int[rcs.getNumPos()]);
        return accumulator.finish();
    }

    private static PackStarFunctionalEvent makeEvent(double delta) {
        return PackStarFunctionalEvents.allOf(
                PackStarFunctionalEvents.atomDistanceWithin(
                        "A79", "NE2", "D206", "OG1",
                        0.0, FIRST_MAX_DISTANCE + delta),
                PackStarFunctionalEvents.anyOf(
                        PackStarFunctionalEvents.atomDistanceWithin(
                                "A167", "OE1", "D205", "NE2",
                                0.0, SECOND_MAX_DISTANCE + delta),
                        PackStarFunctionalEvents.atomDistanceWithin(
                                "A167", "OE2", "D205", "NE2",
                                0.0, SECOND_MAX_DISTANCE + delta)));
    }

    private static boolean independentEvent(Molecule molecule, double delta) {
        double first = distance(molecule, "A79", "NE2", "D206", "OG1");
        double secondOne = distance(molecule, "A167", "OE1", "D205", "NE2");
        double secondTwo = distance(molecule, "A167", "OE2", "D205", "NE2");
        return first <= FIRST_MAX_DISTANCE + delta
                && Math.min(secondOne, secondTwo) <= SECOND_MAX_DISTANCE + delta;
    }

    private static double distance(Molecule molecule,
                                   String residueOne, String atomOne,
                                   String residueTwo, String atomTwo) {
        Residue firstResidue = molecule.getResByPDBResNumberOrNull(residueOne);
        Residue secondResidue = molecule.getResByPDBResNumberOrNull(residueTwo);
        if (firstResidue == null || secondResidue == null) {
            throw new IllegalStateException("reference event residue missing");
        }
        Atom first = firstResidue.getAtomByName(atomOne);
        Atom second = secondResidue.getAtomByName(atomTwo);
        if (first == null || second == null) {
            throw new IllegalStateException("reference event atom missing");
        }
        double[] a = first.getCoords();
        double[] b = second.getCoords();
        if (a == null || b == null || a.length != 3 || b.length != 3) {
            throw new IllegalStateException("reference event coordinates missing");
        }
        double squared = 0.0;
        for (int i = 0; i < 3; i++) {
            if (!Double.isFinite(a[i]) || !Double.isFinite(b[i])) {
                throw new IllegalStateException("reference event coordinates non-finite");
            }
            double difference = a[i] - b[i];
            squared += difference * difference;
        }
        double out = Math.sqrt(squared);
        if (!Double.isFinite(out)) {
            throw new IllegalStateException("reference event distance non-finite");
        }
        return out;
    }

    private static void restore(DoubleMatrix1D params,
                                edu.duke.cs.osprey.confspace.ParametricMolecule pmol) {
        if (params == null) {
            if (!pmol.dofs.isEmpty()) {
                throw new IllegalStateException("reference minimization returned no DOF vector");
            }
            return;
        }
        if (params.size() != pmol.dofs.size()) {
            throw new IllegalStateException("reference DOF vector length mismatch");
        }
        for (int i = 0; i < params.size(); i++) {
            double value = params.get(i);
            if (!Double.isFinite(value)) {
                throw new IllegalStateException("reference DOF vector contains non-finite value");
            }
        }
        for (int i = 0; i < params.size(); i++) {
            pmol.dofs.get(i).apply(params.get(i));
        }
    }

    private static void enforceReferenceLimit(CaseSpaces spaces,
                                              List<Sequence> sequences) {
        for (Sequence sequence : sequences) {
            for (SimpleConfSpace space : List.of(spaces.protein, spaces.ligand, spaces.complex)) {
                int count = count(sequence, space);
                if (count > MAX_REFERENCE_ASSIGNMENTS) {
                    throw new IllegalStateException("reference RC count exceeds frozen limit: "
                            + sequence + " " + count + " > " + MAX_REFERENCE_ASSIGNMENTS);
                }
            }
        }
    }

    private static int count(Sequence sequence, SimpleConfSpace space) {
        return sequence.makeRCs(space).getNumConformations().intValueExact();
    }

    private static void writeManifest(Path outputDir,
                                      Path proteinPdb,
                                      Path ligandPdb,
                                      double targetEpsilon,
                                      int threads,
                                      int maxNumConfs,
                                      List<Sequence> sequences,
                                      CaseSpaces spaces) throws IOException {
        List<String> lines = new ArrayList<>();
        lines.add("key\tvalue");
        lines.add("case\t4wyu_A_D_pc_formal_v1");
        lines.add("protein_pdb\t" + proteinPdb.toAbsolutePath());
        lines.add("ligand_pdb\t" + ligandPdb.toAbsolutePath());
        lines.add("event\t" + EVENT_NAME);
        lines.add("event_semantics\tpolar_distance_surrogate_not_hydrogen_bond_or_activity");
        lines.add("first_threshold_A\t" + FIRST_MAX_DISTANCE);
        lines.add("second_threshold_A\t" + SECOND_MAX_DISTANCE);
        lines.add("threshold_deltas_A\t-0.25,0.0,0.25");
        lines.add("temperature_K\t" + BoltzmannCalculator.TClassic);
        lines.add("gas_constant_kcal_per_mol_K\t" + BoltzmannCalculator.RClassic);
        lines.add("rt_kcal_per_mol\t" + RT);
        lines.add("target_epsilon\t" + targetEpsilon);
        lines.add("pac_samples_property\t" + System.getProperty("packstar.pac.samples", "1000"));
        lines.add("pac_confidence_delta\t" + System.getProperty("packstar.pac.confidence", "0.05"));
        lines.add("pac_random_seed\t" + System.getProperty("packstar.pac.randomSeed", "42"));
        lines.add("threads\t" + threads);
        lines.add("max_num_confs\t" + maxNumConfs);
        lines.add("reference_max_assignments\t" + MAX_REFERENCE_ASSIGNMENTS);
        lines.add("sequence_count\t" + sequences.size());
        lines.add("protein_positions\tA25,A47,A79,A167");
        lines.add("ligand_positions\tD205,D206");
        lines.add("protein_positions_count\t" + spaces.protein.positions.size());
        lines.add("ligand_positions_count\t" + spaces.ligand.positions.size());
        lines.add("complex_positions_count\t" + spaces.complex.positions.size());
        writeLines(outputDir.resolve("run_manifest.tsv"), lines, List.of());
    }

    private static void writeCheckpoint(Path path, List<Row> rows) throws IOException {
        List<String> lines = new ArrayList<>();
        lines.add("sequence\tpackstar_kstar_log10\tpc_status\tpc_probability\tref_kstar_log10");
        for (Row row : rows) {
            lines.add(row.sequence + "\t" + format(row.kstar.scoreLog10())
                    + "\t" + row.complexComputed.observable.getStatus()
                    + "\t" + format(row.complexComputed.observable.getProbability())
                    + "\t" + format(row.referenceKstarLog10()));
        }
        writeLines(path, lines, List.of());
    }

    private static void writeResults(Path path, List<Row> rows) throws IOException {
        List<Row> byKstar = new ArrayList<>(rows);
        byKstar.sort(Comparator.comparingDouble(Row::packstarKstarRankValue).reversed());
        List<Row> byPc = new ArrayList<>(rows);
        byPc.sort(Comparator.comparingDouble(Row::pcRankValue).reversed());
        List<Row> byReference = new ArrayList<>(rows);
        byReference.sort(Comparator.comparingDouble(Row::referenceKstarRankValue).reversed());
        for (int i = 0; i < rows.size(); i++) {
            rows.get(i).packstarKstarRank = rank(
                    byKstar, rows.get(i), Row::packstarKstarRankValue);
            rows.get(i).pcRank = rank(byPc, rows.get(i), Row::pcRankValue);
            rows.get(i).referenceKstarRank = rank(
                    byReference, rows.get(i), Row::referenceKstarRankValue);
        }

        List<String> lines = new ArrayList<>();
        lines.add(String.join("\t", List.of(
                "sequence", "packstar_kstar_log10", "packstar_kstar_lower_log10",
                "packstar_kstar_upper_log10", "packstar_protein_status",
                "packstar_ligand_status", "packstar_complex_status",
                "packstar_complex_p_c_status", "packstar_complex_p_c",
                "packstar_complex_g_restrict", "packstar_p_c_samples",
                "packstar_p_c_hits", "packstar_p_c_ess", "packstar_p_c_max_weight",
                "packstar_complex_seconds", "reference_log10_z_protein",
                "reference_log10_z_ligand", "reference_log10_z_complex",
                "reference_kstar_log10", "reference_gmec_protein",
                "reference_gmec_ligand", "reference_gmec_complex",
                "reference_p_c_minus", "reference_p_c_base", "reference_p_c_plus",
                "packstar_kstar_rank", "p_c_rank", "reference_kstar_rank")));
        for (Row row : rows) {
            PackStarFunctionalObservableResult observable = row.complexComputed.observable;
            lines.add(String.join("\t", List.of(
                    row.sequence,
                    format(row.kstar.scoreLog10()),
                    format(row.kstar.lowerBoundLog10()),
                    format(row.kstar.upperBoundLog10()),
                    row.proteinComputed.result.status.toString(),
                    row.ligandComputed.result.status.toString(),
                    row.complexComputed.result.status.toString(),
                    observable.getStatus().toString(),
                    format(observable.getProbability()),
                    format(observable.getRestrictionFreeEnergy()),
                    Integer.toString(observable.getSampleCount()),
                    Integer.toString(observable.getHitCount()),
                    format(observable.getEffectiveSampleSize()),
                    format(observable.getMaxNormalizedWeight()),
                    format(row.complexComputed.seconds),
                    format(row.proteinReference.log10Z),
                    format(row.ligandReference.log10Z),
                    format(row.complexReference.log10Z),
                    format(row.referenceKstarLog10()),
                    format(row.proteinReference.gmec),
                    format(row.ligandReference.gmec),
                    format(row.complexReference.gmec),
                    format(row.complexReference.probabilities[0]),
                    format(row.complexReference.probabilities[1]),
                    format(row.complexReference.probabilities[2]),
                    Integer.toString(row.packstarKstarRank),
                    Integer.toString(row.pcRank),
                    Integer.toString(row.referenceKstarRank))));
        }
        writeLines(path, lines, List.of());
    }

    private static void writeSelections(Path path, List<Row> rows) throws IOException {
        List<String> lines = new ArrayList<>();
        lines.add("criterion\tsequence\tvalue\tselection_scope");
        select(lines, rows, "packstar_kstar", Row::packstarKstarRankValue, true);
        select(lines, rows, "packstar_p_c", Row::pcRankValue, true);
        select(lines, rows, "reference_kstar", Row::referenceKstarLog10, true);
        select(lines, rows, "reference_p_c", r -> r.complexReference.probabilities[1], true);
        select(lines, rows, "reference_gmec_complex", r -> r.complexReference.gmec, false);
        writeLines(path, lines, List.of());
    }

    private static void select(List<String> lines, List<Row> rows, String criterion,
                               java.util.function.ToDoubleFunction<Row> score, boolean maximize) {
        double best = maximize ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        int valid = 0;
        for (Row row : rows) {
            double value = score.applyAsDouble(row);
            if (!Double.isFinite(value)) continue;
            valid++;
            best = maximize ? Math.max(best, value) : Math.min(best, value);
        }
        if (valid == 0) {
            lines.add(criterion + "\tNA\tNA\tUNRESOLVED");
            return;
        }
        for (Row row : rows) {
            double value = score.applyAsDouble(row);
            if (Double.compare(value, best) == 0) {
                lines.add(criterion + "\t" + row.sequence + "\t" + format(value)
                        + "\t" + (valid == rows.size() ? "ALL_CANDIDATES_POINT_ESTIMATE"
                        : "RESOLVED_SUBSET_ONLY"));
            }
        }
    }

    private static void writeValidation(Path path, List<Row> rows) throws IOException {
        List<String> lines = new ArrayList<>();
        lines.add("sequence\tpc_status\tpc_abs_error\tthreshold_monotone\tkstar_reference_in_interval"
                + "\tpc_diagnostic\tinterpretation");
        for (Row row : rows) {
            var obs = row.complexComputed.observable;
            double[] p = row.complexReference.probabilities;
            boolean monotone = p[0] <= p[1] && p[1] <= p[2];
            if (!monotone) throw new IllegalStateException("reference thresholds not monotone");
            Double lower = row.kstar.lowerBoundLog10();
            Double upper = row.kstar.upperBoundLog10();
            String covered = lower == null || upper == null ? "NA"
                    : Boolean.toString(lower <= row.referenceKstarLog10()
                            && row.referenceKstarLog10() <= upper);
            lines.add(row.sequence + "\t" + obs.getStatus() + "\t"
                    + (obs.isResolved() ? format(Math.abs(obs.getProbability() - p[1])) : "NA")
                    + "\t" + monotone + "\t" + covered + "\t"
                    + obs.getDiagnostic().replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
                    + "\tPOINT_ESTIMATE_ONLY_NO_PC_CONFIDENCE_INTERVAL");
        }
        writeLines(path, lines, List.of());
    }

    private static Row bestByRank(List<Row> rows, boolean kstar, String criterion) {
        Row best = rows.stream()
                .filter(row -> (kstar ? row.packstarKstarRank : row.pcRank) > 0)
                .min(Comparator.comparingInt(row -> kstar
                        ? row.packstarKstarRank : row.pcRank))
                .orElse(null);
        if (best == null) {
            System.out.println("[4WYU-PC] no resolved candidate for " + criterion);
        }
        return best;
    }

    private static int rank(List<Row> rows, Row target,
                            java.util.function.ToDoubleFunction<Row> value) {
        if (!Double.isFinite(value.applyAsDouble(target))) {
            return 0;
        }
        int finiteRank = 1;
        for (Row row : rows) {
            if (!Double.isFinite(value.applyAsDouble(row))) {
                continue;
            }
            if (value.applyAsDouble(row) > value.applyAsDouble(target)) finiteRank++;
        }
        return finiteRank;
    }

    private static double rankValue(Double value) {
        return value != null && Double.isFinite(value) ? value : Double.NEGATIVE_INFINITY;
    }

    private static String format(Double value) {
        return value == null || !Double.isFinite(value) ? "NA"
                : String.format(Locale.ROOT, "%.12g", value);
    }

    private static String format(double value) {
        return Double.isFinite(value) ? String.format(Locale.ROOT, "%.12g", value) : "NA";
    }

    private static void writeLines(Path path, List<String> headerAndRows,
                                   List<String> extraRows) throws IOException {
        List<String> lines = new ArrayList<>(headerAndRows);
        lines.addAll(extraRows);
        Files.write(path, lines, StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING);
    }

    private static Path requiredPath(String property) {
        String value = System.getProperty(property);
        if (value == null || value.trim().isEmpty()) {
            throw new IllegalArgumentException("missing -D" + property);
        }
        return Path.of(value).toAbsolutePath();
    }

    private static int integerProperty(String property, int defaultValue) {
        String value = System.getProperty(property);
        return value == null ? defaultValue : Integer.parseInt(value.trim());
    }

    private static double doubleProperty(String property, double defaultValue) {
        String value = System.getProperty(property);
        return value == null ? defaultValue : Double.parseDouble(value.trim());
    }

    private static final class ReferenceAccumulator {
        private final ConfEnergyCalculator confEcalc;
        private final RCs rcs;
        private final boolean withEvent;
        private final double[] energies;
        private final boolean[][] hits;
        private int next;
        private double minimum = Double.POSITIVE_INFINITY;

        private ReferenceAccumulator(ConfEnergyCalculator confEcalc, RCs rcs,
                                     boolean withEvent, int count) {
            this.confEcalc = confEcalc;
            this.rcs = rcs;
            this.withEvent = withEvent;
            this.energies = new double[count];
            this.hits = withEvent ? new boolean[THRESHOLD_DELTAS.length][count] : null;
        }

        private void visit(int position, int[] assignment) {
            if (position == rcs.getNumPos()) {
                int[] copy = assignment.clone();
                var energy = confEcalc.calcEnergy(new RCTuple(copy));
                if (Double.isNaN(energy.energy)
                        || energy.energy == Double.NEGATIVE_INFINITY) {
                    throw new IllegalStateException("invalid reference energy");
                }
                energies[next] = energy.energy;
                if (Double.isFinite(energy.energy)) {
                    minimum = Math.min(minimum, energy.energy);
                }
                if (withEvent && Double.isFinite(energy.energy)) {
                    if (energy.pmol == null) {
                        throw new IllegalStateException(
                                "finite reference energy returned no molecule");
                    }
                    restore(energy.params, energy.pmol);
                    for (int i = 0; i < THRESHOLD_DELTAS.length; i++) {
                        hits[i][next] = independentEvent(
                                energy.pmol.mol, THRESHOLD_DELTAS[i]);
                    }
                }
                next++;
                if (next % 1000 == 0 || next == energies.length) {
                    System.out.println("[4WYU-PC] reference " + next + "/" + energies.length);
                }
                return;
            }
            for (int rc : rcs.get(position)) {
                assignment[position] = rc;
                visit(position + 1, assignment);
            }
        }

        private ReferenceSummary finish() {
            if (next != energies.length || !Double.isFinite(minimum)) {
                throw new IllegalStateException("reference enumeration count mismatch");
            }
            double total = 0.0;
            double[] eventMass = withEvent ? new double[THRESHOLD_DELTAS.length] : null;
            for (int i = 0; i < energies.length; i++) {
                double weight = Math.exp((minimum - energies[i]) / RT);
                total += weight;
                if (withEvent) {
                    for (int j = 0; j < THRESHOLD_DELTAS.length; j++) {
                        if (hits[j][i]) {
                            eventMass[j] += weight;
                        }
                    }
                }
            }
            if (!(total > 0.0) || !Double.isFinite(total)) {
                throw new IllegalStateException("reference partition sum is invalid");
            }
            double logZ = -minimum / RT + Math.log(total);
            double[] probabilities = withEvent
                    ? new double[]{eventMass[0] / total, eventMass[1] / total, eventMass[2] / total}
                    : new double[]{Double.NaN, Double.NaN, Double.NaN};
            return new ReferenceSummary(energies.length, minimum, logZ,
                    logZ / Math.log(10.0), probabilities);
        }
    }

    private static final class CaseSpaces {
        private final SimpleConfSpace protein;
        private final SimpleConfSpace ligand;
        private final SimpleConfSpace complex;

        private CaseSpaces(SimpleConfSpace protein, SimpleConfSpace ligand,
                           SimpleConfSpace complex) {
            this.protein = protein;
            this.ligand = ligand;
            this.complex = complex;
        }
    }

    private static final class StateModel {
        private final SimpleConfSpace space;
        private final ConfEnergyCalculator minimizingConfEcalc;
        private final EnergyMatrix rigidEmat;
        private final EnergyMatrix minimizingEmat;
        private final UpdatingEnergyMatrix corrections;
        private final Parallelism parallelism;
        private final String name;

        private StateModel(SimpleConfSpace space,
                           ConfEnergyCalculator minimizingConfEcalc,
                           EnergyMatrix rigidEmat,
                           EnergyMatrix minimizingEmat,
                           UpdatingEnergyMatrix corrections,
                           Parallelism parallelism,
                           String name) {
            this.space = space;
            this.minimizingConfEcalc = minimizingConfEcalc;
            this.rigidEmat = rigidEmat;
            this.minimizingEmat = minimizingEmat;
            this.corrections = corrections;
            this.parallelism = parallelism;
            this.name = name;
        }
    }

    private static final class Computed {
        private final PartitionFunction.Result result;
        private final PackStarFunctionalObservableResult observable;
        private final double seconds;
        private final int rcCount;

        private Computed(PartitionFunction.Result result,
                         PackStarFunctionalObservableResult observable,
                         double seconds, int rcCount) {
            this.result = result;
            this.observable = observable;
            this.seconds = seconds;
            this.rcCount = rcCount;
        }
    }

    private static final class ReferenceSummary {
        private final int assignmentCount;
        private final double gmec;
        private final double logZ;
        private final double log10Z;
        private final double[] probabilities;

        private ReferenceSummary(int assignmentCount, double gmec,
                                 double logZ, double log10Z,
                                 double[] probabilities) {
            this.assignmentCount = assignmentCount;
            this.gmec = gmec;
            this.logZ = logZ;
            this.log10Z = log10Z;
            this.probabilities = probabilities;
        }
    }

    private static final class Row {
        private final String sequence;
        private final Computed proteinComputed;
        private final Computed ligandComputed;
        private final Computed complexComputed;
        private final ReferenceSummary proteinReference;
        private final ReferenceSummary ligandReference;
        private final ReferenceSummary complexReference;
        private final KStarScore kstar;
        private int packstarKstarRank;
        private int pcRank;
        private int referenceKstarRank;

        private Row(String sequence,
                    Computed proteinComputed,
                    Computed ligandComputed,
                    Computed complexComputed,
                    ReferenceSummary proteinReference,
                    ReferenceSummary ligandReference,
                    ReferenceSummary complexReference,
                    KStarScore kstar) {
            this.sequence = sequence;
            this.proteinComputed = proteinComputed;
            this.ligandComputed = ligandComputed;
            this.complexComputed = complexComputed;
            this.proteinReference = proteinReference;
            this.ligandReference = ligandReference;
            this.complexReference = complexReference;
            this.kstar = kstar;
        }

        private double referenceKstarLog10() {
            return complexReference.log10Z - proteinReference.log10Z - ligandReference.log10Z;
        }

        private double packstarKstarRankValue() {
            return rankValue(kstar.scoreLog10());
        }

        private double pcRankValue() {
            return complexComputed.observable.isResolved()
                    ? rankValue(complexComputed.observable.getProbability())
                    : Double.NEGATIVE_INFINITY;
        }

        private double referenceKstarRankValue() {
            return referenceKstarLog10();
        }
    }
}
