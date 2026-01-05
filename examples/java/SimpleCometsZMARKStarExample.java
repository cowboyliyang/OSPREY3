/*
** 简化的 COMETSZ + MARKStar 示例
**
** 这是一个简化版本，基于现有的测试代码 TestCometsZWithBBKStarAndMARKStar
** 可以直接运行，演示如何使用 COMETSZ + MARKStar 进行多状态蛋白质设计
**
** 使用场景：蛋白质-配体结合亲和力优化
*/

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.CometsZ;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SimpleCometsZMARKStarExample {

    public static void main(String[] args) {
        System.out.println("\n" + "=".repeat(80));
        System.out.println("简化的 COMETSZ + MARKStar 多状态蛋白质设计示例");
        System.out.println("=".repeat(80) + "\n");

        // 运行示例
        runExample();

        System.out.println("\n示例运行完成！");
    }

    public static void runExample() {

        // ===== 第一步：定义系统 =====
        System.out.println("【步骤 1】定义蛋白质-配体系统\n");

        // 读取 PDB 结构
        Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");
        ForcefieldParams ffparams = new ForcefieldParams();

        // 定义蛋白质链（设计位点：G649，其他位置保持野生型）
        Strand protein = new Strand.Builder(mol)
            .setResidues("G648", "G654")
            .build();
        protein.flexibility.get("G649")
            .setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU")
            .addWildTypeRotamers()
            .setContinuous();
        protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // 定义配体链（全部保持野生型，只考虑构象柔性）
        Strand ligand = new Strand.Builder(mol)
            .setResidues("A155", "A194")
            .build();
        ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
        ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

        // 创建构象空间
        SimpleConfSpace proteinSpace = new SimpleConfSpace.Builder().addStrand(protein).build();
        SimpleConfSpace ligandSpace = new SimpleConfSpace.Builder().addStrand(ligand).build();
        SimpleConfSpace complexSpace = new SimpleConfSpace.Builder().addStrands(protein, ligand).build();

        System.out.println("  ✓ 蛋白质构象空间: " + proteinSpace.positions.size() + " 个灵活位点");
        System.out.println("  ✓ 配体构象空间: " + ligandSpace.positions.size() + " 个灵活位点");
        System.out.println("  ✓ 复合物构象空间: " + complexSpace.positions.size() + " 个灵活位点\n");

        // ===== 第二步：配置 COMETSZ =====
        System.out.println("【步骤 2】配置 COMETSZ 多状态设计\n");

        // 创建三个状态
        CometsZ.State stateProtein = new CometsZ.State("Protein", proteinSpace);
        CometsZ.State stateLigand = new CometsZ.State("Ligand", ligandSpace);
        CometsZ.State stateComplex = new CometsZ.State("Complex", complexSpace);

        // 定义目标函数：最小化结合自由能
        // ΔG_binding = G_complex - G_protein - G_ligand
        CometsZ.LMFE objective = new CometsZ.LMFE.Builder()
            .addState(stateComplex, 1.0)
            .addState(stateProtein, -1.0)
            .addState(stateLigand, -1.0)
            .build();

        System.out.println("  目标函数: ΔG_binding = G(complex) - G(protein) - G(ligand)");
        System.out.println("  优化目标: 找到结合自由能最低的序列\n");

        // 创建 COMETSZ 实例
        CometsZ cometsZ = new CometsZ.Builder(objective)
            .setEpsilon(0.95)                        // 95% 置信度
            .setMaxSimultaneousMutations(1)          // 每次只考虑一个位点的突变
            .setObjectiveWindowSize(100.0)           // 能量窗口大小
            .setObjectiveWindowMax(100.0)            // 能量窗口最大值
            .setLogFile(new File("cometsz.results.tsv"))
            .build();

        System.out.println("  ✓ Epsilon (精度): 0.95");
        System.out.println("  ✓ 最大同时突变数: 1");
        System.out.println("  ✓ 结果保存至: cometsz.results.tsv\n");

        // ===== 第三步：初始化 MARKStar =====
        System.out.println("【步骤 3】初始化 MARKStar 配分函数计算器\n");

        Map<CometsZ.State, EnergyMatrix> rigidEmats = new HashMap<>();
        Map<CometsZ.State, EnergyMatrix> minimizingEmats = new HashMap<>();

        // 创建共享的能量计算器
        List<SimpleConfSpace> confSpaces = List.of(proteinSpace, ligandSpace, complexSpace);
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpaces, ffparams)
            .setParallelism(Parallelism.makeCpu(4))
            .setIsMinimizing(true)
            .build();

        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
            .setIsMinimizing(false)
            .build();

        System.out.println("  正在为每个状态计算能量矩阵...");

        try {
            // 为每个状态设置能量矩阵和 MARKStar
            for (CometsZ.State state : cometsZ.states) {

                System.out.println("    处理状态: " + state.name);

                // 设置柔性能量计算器
                state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, minimizingEcalc)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, minimizingEcalc)
                        .build()
                        .calcReferenceEnergies())
                    .build();

                // 计算柔性能量矩阵
                EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
                    .build()
                    .calcEnergyMatrix();
                minimizingEmats.put(state, minimizingEmat);

                // 计算刚性能量矩阵
                ConfEnergyCalculator rigidConfEcalc = new ConfEnergyCalculator.Builder(state.confSpace, rigidEcalc)
                    .setReferenceEnergies(state.confEcalc.eref)
                    .build();
                EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                    .build()
                    .calcEnergyMatrix();
                rigidEmats.put(state, rigidEmat);

                // 设置 MARKStar 配分函数工厂
                state.pfuncFactory = (rcs) -> {
                    MARKStarBound markstar = new MARKStarBound(
                        state.confSpace,
                        rigidEmats.get(state),
                        minimizingEmats.get(state),
                        state.confEcalc,
                        rcs,
                        Parallelism.makeCpu(4)
                    );
                    // 设置校正矩阵（必需）
                    markstar.setCorrections(new edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix(
                        state.confSpace,
                        minimizingEmats.get(state),
                        state.confEcalc
                    ));
                    return markstar;
                };

                // 设置构象树工厂
                state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(minimizingEmats.get(state), rcs)
                    .setTraditional()
                    .build();

                // 设置片段能量
                state.fragmentEnergies = minimizingEmats.get(state);
            }

            System.out.println("\n  ✓ 能量矩阵计算完成");
            System.out.println("  ✓ MARKStar 配置完成\n");

            // ===== 第四步：运行 COMETSZ =====
            System.out.println("【步骤 4】运行 COMETSZ 序列优化\n");
            System.out.println("  正在搜索最优序列...");
            System.out.println("  （这可能需要几分钟，取决于系统大小）\n");
            System.out.println("=".repeat(80) + "\n");

            // 创建新的能量计算器用于运行
            try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
                .setParallelism(Parallelism.makeCpu(4))
                .build()) {

                // 刷新状态的能量计算器
                for (CometsZ.State state : cometsZ.states) {
                    state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
                }

                // 运行 COMETSZ，寻找最优的 5 个序列
                List<CometsZ.SequenceInfo> sequences = cometsZ.findBestSequences(5);

                // ===== 第五步：输出结果 =====
                System.out.println("\n" + "=".repeat(80));
                System.out.println("【结果】找到 " + sequences.size() + " 个最优序列");
                System.out.println("=".repeat(80) + "\n");

                for (int i = 0; i < sequences.size(); i++) {
                    CometsZ.SequenceInfo seq = sequences.get(i);
                    System.out.println("序列 #" + (i + 1) + ":");
                    System.out.println("  " + seq.sequence);
                    System.out.println("  结合自由能 (ΔG): [" +
                        String.format("%.4f", seq.objective.lower) + ", " +
                        String.format("%.4f", seq.objective.upper) + "] kcal/mol");

                    // 输出各状态的自由能
                    for (CometsZ.State state : cometsZ.states) {
                        var pfunc = seq.pfuncResults.get(state);
                        var fe = pfunc.values.calcFreeEnergyBounds();
                        System.out.println("    " + state.name + ": [" +
                            String.format("%.4f", fe.lower) + ", " +
                            String.format("%.4f", fe.upper) + "] kcal/mol");
                    }
                    System.out.println();
                }

                System.out.println("=".repeat(80));
                System.out.println("所有结果已保存至: cometsz.results.tsv");
                System.out.println("=".repeat(80) + "\n");
            }

        } finally {
            minimizingEcalc.close();
            rigidEcalc.close();
        }

        // ===== 总结 =====
        System.out.println("【总结】");
        System.out.println("  本示例展示了如何使用:");
        System.out.println("  1. COMETSZ - 多状态设计，优化结合自由能");
        System.out.println("  2. MARKStar - 快速配分函数计算（比传统方法快 5-100 倍）");
        System.out.println("  3. 刚性/柔性能量矩阵 - 提供严格的数学界限");
        System.out.println("\n  优势:");
        System.out.println("  ✓ 同时考虑结合态和非结合态");
        System.out.println("  ✓ 快速准确的配分函数计算");
        System.out.println("  ✓ 自动处理构象柔性");
        System.out.println("  ✓ 提供严格的置信区间");
    }
}
