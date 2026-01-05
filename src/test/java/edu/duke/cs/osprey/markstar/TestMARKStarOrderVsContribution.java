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
 * 实际运行MARK*并分析：
 * 1. Minimize顺序
 * 2. 最终贡献大小
 * 3. 两者的相关性
 */
public class TestMARKStarOrderVsContribution {

    @Test
    public void testSmallProblem() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("测试: MARK* Minimize顺序 vs Partition Function贡献");
        System.out.println("=".repeat(80) + "\n");

        // 创建一个小问题（3个flexible残基）
        ConfSpaces confSpaces = makeSmall1GUA(3);

        runMARKStarWithTracking(confSpaces, 0.68);
    }

    @Test
    public void testTinyProblem() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("测试: MARK* Minimize顺序 vs Partition Function贡献 (Very Small)");
        System.out.println("=".repeat(80) + "\n");

        // 创建一个更小的问题（2个flexible残基）- 更快完成
        ConfSpaces confSpaces = makeSmall1GUA(2);

        runMARKStarWithTracking(confSpaces, 0.68);
    }

    @Test
    public void testMediumProblem() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("测试: MARK* Minimize顺序 vs Partition Function贡献 (Medium)");
        System.out.println("=".repeat(80) + "\n");

        // 中等大小问题（4个flexible残基）
        ConfSpaces confSpaces = makeSmall1GUA(4);

        runMARKStarWithTracking(confSpaces, 0.68);
    }

    private void runMARKStarWithTracking(ConfSpaces confSpaces, double epsilon) throws IOException {
        Parallelism parallelism = Parallelism.makeCpu(4);

        System.out.println("设置: " + confSpaces.complex.positions.size() + " flexible positions");
        System.out.println("Target epsilon: " + epsilon);
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
            confSpaces.complex.positions.size() + "pos.csv";

        tracker.analyzeCorrelation(outputFile);
    }

    private static ConfSpaces makeSmall1GUA(int numFlex) {
        ConfSpaces confSpaces = new ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .addMoleculeForWildTypeRotamers(mol)
            .build();

        // Protein strand - small number of flexible residues
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

    /**
     * 用于快速演示 - 只有1个flexible残基
     */
    @Test
    public void testMinimalProblem() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("测试: MARK* Minimize顺序 vs Partition Function贡献 (Minimal - 演示用)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = makeSmall1GUA(1);
        runMARKStarWithTracking(confSpaces, 0.9);
    }
}
