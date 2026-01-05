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
 * 测试更大的系统（更多flexible residues）
 * 以确保minimize足够多的conformations来分析相关性
 *
 * Epsilon保持0.01 (1%)
 */
public class TestLargerSystem {

    @Test
    public void test5FlexibleResidues() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 大规模测试: 5个flexible残基");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(5);
        runAndAnalyze(confSpaces, 0.01, "5flex");
    }

    @Test
    public void test6FlexibleResidues() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 大规模测试: 6个flexible残基");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(6);
        runAndAnalyze(confSpaces, 0.01, "6flex");
    }

    @Test
    public void test7FlexibleResidues() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 大规模测试: 7个flexible残基");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(7);
        runAndAnalyze(confSpaces, 0.01, "7flex");
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
        String outputFile = "markstar_order_vs_contribution_" + tag + "_eps0p01.csv";

        System.out.println("开始分析minimize顺序 vs 贡献相关性...");
        tracker.analyzeCorrelation(outputFile);

        System.out.println("\n" + "=".repeat(80));
        System.out.println("完成! 详细结果已保存到: " + outputFile);
        System.out.println("=".repeat(80));
    }

    /**
     * 创建1GUA系统，指定flexible残基数量
     */
    private static ConfSpaces make1GUAWithFlexibleResidues(int numFlexProtein) {
        ConfSpaces confSpaces = new ConfSpaces();
        confSpaces.ffparams = new ForcefieldParams();

        Molecule mol = PDBIO.read(FileTools.readResource("/1gua_adj.min.pdb"));

        ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(confSpaces.ffparams.forcefld)
            .addMoleculeForWildTypeRotamers(mol)
            .build();

        // Protein strand - 从residue 21开始，添加numFlexProtein个flexible残基
        Strand protein = new Strand.Builder(mol)
            .setTemplateLibrary(templateLib)
            .setResidues("1", "180")
            .build();

        int start = 21;
        System.out.println("设置protein flexible残基:");
        for (int i = start; i < start + numFlexProtein; i++) {
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

    /**
     * 快速测试 - 4个flexible残基（预计几分钟完成）
     */
    @Test
    public void test4FlexibleResiduesTouchstone() throws IOException {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("MARK* 快速测试: 4个flexible残基 (预计5-10分钟)");
        System.out.println("=".repeat(80) + "\n");

        ConfSpaces confSpaces = make1GUAWithFlexibleResidues(4);
        runAndAnalyze(confSpaces, 0.01, "4flex");
    }
}
