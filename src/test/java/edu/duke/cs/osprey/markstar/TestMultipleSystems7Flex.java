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
 * 在5个不同的系统配置上测试7个flexible residues
 * 每个配置使用不同的起始残基位置
 *
 * 目的：验证SubtreeUpper vs ErrorBound策略在多个系统上的普适性
 */
public class TestMultipleSystems7Flex {

    /**
     * System 1: Residues 21-27 (原始配置)
     */
    @Test
    public void testSystem1_Residues21to27() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("System 1: 7 flexible residues (21-27)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(21, 7);
        runAndAnalyze(confSpaces, 0.01, "system1_7flex");
    }

    /**
     * System 2: Residues 30-36 (不同区域)
     */
    @Test
    public void testSystem2_Residues30to36() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("System 2: 7 flexible residues (30-36)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(30, 7);
        runAndAnalyze(confSpaces, 0.01, "system2_7flex");
    }

    /**
     * System 3: Residues 40-46 (中间区域)
     */
    @Test
    public void testSystem3_Residues40to46() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("System 3: 7 flexible residues (40-46)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(40, 7);
        runAndAnalyze(confSpaces, 0.01, "system3_7flex");
    }

    /**
     * System 4: Residues 50-56 (另一个区域)
     */
    @Test
    public void testSystem4_Residues50to56() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("System 4: 7 flexible residues (50-56)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(50, 7);
        runAndAnalyze(confSpaces, 0.01, "system4_7flex");
    }

    /**
     * System 5: Residues 60-66 (更后面的区域)
     */
    @Test
    public void testSystem5_Residues60to66() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("System 5: 7 flexible residues (60-66)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(60, 7);
        runAndAnalyze(confSpaces, 0.01, "system5_7flex");
    }

    /**
     * 运行所有5个系统的测试（串行执行避免冲突）
     */
    @Test
    public void testAll5Systems() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("运行所有5个系统的7-flex测试");
        System.out.println("=".repeat(80) + "\n");

        int[][] configs = {
            {21, 1},  // System 1: Residues 21-27
            {30, 2},  // System 2: Residues 30-36
            {40, 3},  // System 3: Residues 40-46
            {50, 4},  // System 4: Residues 50-56
            {60, 5}   // System 5: Residues 60-66
        };

        for (int[] config : configs) {
            int startResidue = config[0];
            int systemNum = config[1];
            String tag = "system" + systemNum + "_7flex";

            System.out.println("\n" + "=".repeat(80));
            System.out.println("System " + systemNum + ": 7 flexible residues (" + startResidue + "-" + (startResidue+6) + ")");
            System.out.println("=".repeat(80) + "\n");

            ConfSpaces confSpaces = make1GUAWithFlexibleResidues(startResidue, 7);
            runAndAnalyze(confSpaces, 0.01, tag);

            System.out.println("\n✓ System " + systemNum + " 完成\n");
        }

        System.out.println("\n" + "=".repeat(80));
        System.out.println("所有5个系统测试完成！");
        System.out.println("=".repeat(80));
    }

    /**
     * 运行MARK*并分析结果
     */
    private void runAndAnalyze(ConfSpaces confSpaces, double epsilon, String tag) throws IOException {
        Parallelism parallelism = Parallelism.makeCpu(4);

        System.out.println("配置:");
        System.out.println("  - Flexible positions: " + confSpaces.complex.positions.size());
        System.out.println("  - Target epsilon: " + epsilon);

        // 估算conformation space大小
        int totalPositions = confSpaces.complex.positions.size();
        System.out.println("  - Total positions: " + totalPositions);
        System.out.println("  - Estimated conf space size: ~10^" + (totalPositions * 2) + " (假设每个位置~100个rotamers)");
        System.out.println();

        // 创建energy calculators
        System.out.println("创建energy calculators...");
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
        long startTime = System.currentTimeMillis();

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

        long ematTime = System.currentTimeMillis() - startTime;
        System.out.println("Energy matrix计算完成 (耗时: " + ematTime / 1000.0 + "s)\n");

        // 创建RCs
        RCs rcs = new RCs(confSpaces.complex);

        System.out.println("运行MARK* with tracking...");
        System.out.println("这可能需要几分钟到几十分钟，取决于系统大小...\n");

        startTime = System.currentTimeMillis();

        // 使用tracker来运行MARK*并捕获输出
        MARKStarTracker tracker = new MARKStarTracker();
        PartitionFunction.Result result = tracker.runAndTrack(
            confSpaces.complex, rigidEmat, minimizingEmat, confEcalc, rcs, parallelism, epsilon);

        long pfuncTime = System.currentTimeMillis() - startTime;

        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 计算完成!");
        System.out.println("=".repeat(80));
        System.out.println("Partition function计算耗时: " + pfuncTime / 1000.0 + "s");
        System.out.println("Status: " + result.status);
        System.out.println("Lower bound (q*): " + result.values.qstar);
        System.out.println("Upper bound (p*): " + result.values.pstar);
        System.out.println("Effective epsilon: " + result.values.getEffectiveEpsilon());
        System.out.println("Conformations minimized: " + result.numConfs);
        System.out.println();

        // 分析minimize顺序 vs 贡献
        String outputFile = "markstar_order_vs_contribution_" + tag + "_eps0p010.csv";

        System.out.println("开始分析minimize顺序 vs 贡献相关性...");
        tracker.analyzeCorrelation(outputFile);

        System.out.println("\n" + "=".repeat(80));
        System.out.println("完成! 详细结果已保存到: " + outputFile);
        System.out.println("=".repeat(80));
    }

    /**
     * 创建1GUA系统，指定起始残基和flexible残基数量
     */
    private static ConfSpaces make1GUAWithFlexibleResidues(int startResidue, int numFlexProtein) {
        ConfSpaces confSpaces = new ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .addMoleculeForWildTypeRotamers(mol)
            .build();

        // Protein strand - 从指定residue开始，添加numFlexProtein个flexible残基
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("1", "180")
            .build();

        System.out.println("设置protein flexible残基:");
        for (int i = startResidue; i < startResidue + numFlexProtein; i++) {
            protein.flexibility.get(i + "")
                .setLibraryRotamers(Strand.WildType)
                .addWildTypeRotamers()
                .setContinuous();
            System.out.println("  - Residue " + i + ": wild-type rotamers + continuous");
        }

        // Ligand strand - 固定1个flexible残基
        Strand ligand = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("181", "215")
            .build();

        ligand.flexibility.get("209")
            .setLibraryRotamers(Strand.WildType)
            .addWildTypeRotamers()
            .setContinuous();
        System.out.println("设置ligand flexible残基:");
        System.out.println("  - Residue 209: wild-type rotamers + continuous");
        System.out.println();

        confSpaces.protein = new SimpleConfSpace.Builder().addStrand(protein).build();
        confSpaces.ligand = new SimpleConfSpace.Builder().addStrand(ligand).build();
        confSpaces.complex = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        return confSpaces;
    }
}
