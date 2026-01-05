package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar.ConfSpaces;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.io.IOException;

/**
 * 实际运行MARK*来分析minimize顺序 vs 贡献
 * 使用更紧的epsilon来确保minimize更多conformations
 */
public class TestActualOrderVsContribution {

    @Test
    public void testWithTightEpsilon() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 分析: Minimize顺序 vs Partition Function贡献");
        System.out.println("使用tight epsilon来确保minimize更多conformations");
        System.out.println("=".repeat(80) + "\n");

        // 2个flexible残基，但是用很紧的epsilon
        ConfSpaces confSpaces = makeSmall1GUA(2);

        runWithTightEpsilon(confSpaces, 0.01); // 更紧的epsilon
    }

    private void runWithTightEpsilon(ConfSpaces confSpaces, double epsilon) throws IOException {
        Parallelism parallelism = Parallelism.makeCpu(4);

        System.out.println("设置:");
        System.out.println("  - Flexible positions: " + confSpaces.complex.positions.size());
        System.out.println("  - Target epsilon: " + epsilon);
        System.out.println();

        // 创建energy calculators
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
            .setParallelism(parallelism)
            .setIsMinimizing(true)
            .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
            .setIsMinimizing(false)
            .build();

        // 创建conf energy calculator
        ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpaces.complex, minimizingEcalc)
            .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaces.complex, minimizingEcalc)
                .build()
                .calcReferenceEnergies())
            .build();

        System.out.println("计算energy matrices...");

        // 计算energy matrices
        EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
            .build()
            .calcEnergyMatrix();

        EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(
            new ConfEnergyCalculator.Builder(confSpaces.complex, rigidEcalc)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaces.complex, rigidEcalc)
                    .build()
                    .calcReferenceEnergies())
                .build())
            .build()
            .calcEnergyMatrix();

        // 创建RCs
        RCs rcs = new RCs(confSpaces.complex);

        System.out.println("运行MARK* with tracking...\n");

        // 使用tracker来运行MARK*并捕获输出
        MARKStarTracker tracker = new MARKStarTracker();
        PartitionFunction.Result result = tracker.runAndTrack(
            confSpaces.complex, rigidEmat, minimizingEmat, confEcalc, rcs, parallelism, epsilon);

        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 计算完成!");
        System.out.println("=".repeat(80));
        System.out.println("Status: " + result.status);
        System.out.println("Lower bound (q*): " + result.values.qstar);
        System.out.println("Upper bound (p*): " + result.values.pstar);
        System.out.println("Effective epsilon: " + result.values.getEffectiveEpsilon());
        System.out.println("Conformations minimized: " + result.numConfs);

        // 分析minimize顺序 vs 贡献
        String outputFile = "markstar_order_vs_contribution_" +
            confSpaces.complex.positions.size() + "pos_eps" +
            String.format("%.3f", epsilon).replace(".", "p") + ".csv";

        tracker.analyzeCorrelation(outputFile);
    }

    private static ConfSpaces makeSmall1GUA(int numFlex) {
        ConfSpaces confSpaces = new ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .addMoleculeForWildTypeRotamers(mol)
            .build();

        // Protein strand
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("1", "180")
            .build();

        int start = 21;
        for (int i = start; i < start + numFlex; i++) {
            protein.flexibility.get(i + "")
                .setLibraryRotamers(Strand.WildType)
                .addWildTypeRotamers()
                .setContinuous();
        }

        // Ligand strand
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("181", "215")
            .build();
        ligand.flexibility.get("209")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }
}
