package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.io.ByteArrayOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * 追踪MARK*的执行，记录：
 * 1. 每个被minimize的conformation的顺序
 * 2. 最终的能量和Boltzmann weight
 * 3. 分析minimize顺序vs贡献大小的相关性
 *
 * 通过捕获stdout来解析MARK*的输出
 */
public class MARKStarTracker {

    public static class MinimizedConf {
        public final int order;              // Minimize顺序 (1st, 2nd, 3rd...)
        public final String confString;      // Conformation string
        public final double score;           // Initial score
        public final double lowerBound;      // Initial lower bound
        public final double corrected;       // Corrected energy
        public final double energy;          // Final minimized energy
        public final double oldConfLower;    // OLD conf lower bound (before minimize)
        public final double oldConfUpper;    // OLD conf upper bound (before minimize)
        public final BigDecimal subtreeLower;    // Subtree lower bound (before minimize)
        public final BigDecimal subtreeUpper;    // Subtree upper bound (before minimize)
        public final BigDecimal errorBound;      // Error bound (before minimize)
        public final BigDecimal globalLowerZ;    // Global lower Z bound at this point
        public final BigDecimal globalUpperZ;    // Global upper Z bound at this point
        public final double delta;           // Epsilon at this point

        public MinimizedConf(int order, String confString, double score, double lowerBound,
                           double corrected, double energy, double oldConfLower, double oldConfUpper,
                           BigDecimal subtreeLower, BigDecimal subtreeUpper, BigDecimal errorBound,
                           BigDecimal globalLowerZ, BigDecimal globalUpperZ, double delta) {
            this.order = order;
            this.confString = confString;
            this.score = score;
            this.lowerBound = lowerBound;
            this.corrected = corrected;
            this.energy = energy;
            this.oldConfLower = oldConfLower;
            this.oldConfUpper = oldConfUpper;
            this.subtreeLower = subtreeLower;
            this.subtreeUpper = subtreeUpper;
            this.errorBound = errorBound;
            this.globalLowerZ = globalLowerZ;
            this.globalUpperZ = globalUpperZ;
            this.delta = delta;
        }

        public BigDecimal getBoltzmannWeight() {
            BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
            return bc.calc(energy);
        }
    }

    private final List<MinimizedConf> minimizedConfs = new ArrayList<>();

    /**
     * 运行MARK*并捕获输出来追踪conformations
     */
    public PartitionFunction.Result runAndTrack(SimpleConfSpace confSpace,
                                                EnergyMatrix rigidEmat,
                                                EnergyMatrix minimizingEmat,
                                                ConfEnergyCalculator confEcalc,
                                                RCs rcs,
                                                Parallelism parallelism,
                                                double epsilon) {

        // 捕获stdout
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos);
        PrintStream old = System.out;
        System.setOut(ps);

        MARKStarBound pfunc = null;
        try {
            pfunc = new MARKStarBound(confSpace, rigidEmat, minimizingEmat, confEcalc, rcs, parallelism);
            pfunc.setReportProgress(true); // 打开进度报告

            // 创建correction matrix (CRITICAL!)
            UpdatingEnergyMatrix correctionEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat);
            pfunc.setCorrections(correctionEmat);

            pfunc.init(epsilon);
            pfunc.compute();

        } finally {
            // 恢复stdout
            System.out.flush();
            System.setOut(old);
        }

        // 解析输出
        String output = baos.toString();
        parseOutput(output);

        // 打印捕获的输出（用于调试）
        System.out.println(output);

        return pfunc != null ? pfunc.makeResult() : null;
    }

    /**
     * 解析MARK*的输出，提取minimize的conformations
     */
    private void parseOutput(String output) {
        /*
         * 寻找这样的行：
         * [4, 11, 7, 10]conf:  10, score:  -24.535469, lower:  -24.658562, corrected:  -24.658562 energy:  -27.246422,
         * confBounds:[-24.658562, -24.535469], subtreeBounds:[1.2e+03, 1.3e+03], errorBound:1.0e+02,
         * globalBounds:[1.239166e+03, 1.254382e+03], delta:    0.012163, time:     58.34s
         */

        Pattern pattern = Pattern.compile(
            "\\[([\\d,\\s]+)\\]conf:\\s*(\\d+),\\s*score:\\s*([\\d.e+-]+),\\s*lower:\\s*([\\d.e+-]+),\\s*corrected:\\s*([\\d.e+-]+)\\s*energy:\\s*([\\d.e+-]+),\\s*" +
            "confBounds:\\[([\\d.e+-]+),\\s*([\\d.e+-]+)\\],\\s*subtreeBounds:\\[([\\d.e+-]+),\\s*([\\d.e+-]+)\\],\\s*errorBound:([\\d.e+-]+),\\s*" +
            "globalBounds:\\[([\\d.e+-]+),\\s*([\\d.e+-]+)\\],\\s*delta:\\s*([\\d.]+)"
        );

        int confNum = 0;
        for (String line : output.split("\n")) {
            Matcher m = pattern.matcher(line);
            if (m.find()) {
                confNum++;
                String confString = "[" + m.group(1) + "]";
                int order = Integer.parseInt(m.group(2));
                double score = Double.parseDouble(m.group(3));
                double lower = Double.parseDouble(m.group(4));
                double corrected = Double.parseDouble(m.group(5));
                double energy = Double.parseDouble(m.group(6));
                double oldConfLower = Double.parseDouble(m.group(7));
                double oldConfUpper = Double.parseDouble(m.group(8));
                BigDecimal subtreeLower = new BigDecimal(m.group(9));
                BigDecimal subtreeUpper = new BigDecimal(m.group(10));
                BigDecimal errorBound = new BigDecimal(m.group(11));
                BigDecimal globalLowerZ = new BigDecimal(m.group(12));
                BigDecimal globalUpperZ = new BigDecimal(m.group(13));
                double delta = Double.parseDouble(m.group(14));

                MinimizedConf conf = new MinimizedConf(
                    order, confString, score, lower, corrected,
                    energy, oldConfLower, oldConfUpper,
                    subtreeLower, subtreeUpper, errorBound,
                    globalLowerZ, globalUpperZ, delta
                );
                minimizedConfs.add(conf);
            }
        }

        System.out.println("\n解析到 " + minimizedConfs.size() + " 个minimized conformations");
    }

    /**
     * 分析minimize顺序 vs 贡献大小
     */
    public void analyzeCorrelation(String outputFile) throws IOException {
        if (minimizedConfs.isEmpty()) {
            System.out.println("没有找到任何minimized conformations!");
            return;
        }

        System.out.println("\n" + "=".repeat(80));
        System.out.println("分析: Minimize顺序 vs Partition Function贡献");
        System.out.println("=".repeat(80));

        // 计算total Z
        BigDecimal totalZ = BigDecimal.ZERO;
        for (MinimizedConf conf : minimizedConfs) {
            totalZ = totalZ.add(conf.getBoltzmannWeight());
        }

        System.out.println("\n总共minimize了 " + minimizedConfs.size() + " 个conformations");
        System.out.println("Total Z (from minimized confs): " + totalZ);

        // 按minimize顺序排序
        List<MinimizedConf> byOrder = new ArrayList<>(minimizedConfs);
        byOrder.sort(Comparator.comparingInt(c -> c.order));

        // 按贡献大小排序
        List<MinimizedConf> byContribution = new ArrayList<>(minimizedConfs);
        byContribution.sort((a, b) -> b.getBoltzmannWeight().compareTo(a.getBoltzmannWeight()));

        // 创建贡献排名map
        Map<String, Integer> contributionRank = new HashMap<>();
        for (int i = 0; i < byContribution.size(); i++) {
            contributionRank.put(byContribution.get(i).confString, i + 1);
        }

        // 显示前20个被minimize的conformations
        System.out.println("\n" + "=".repeat(80));
        System.out.println("前20个被minimize的conformations:");
        System.out.println("-".repeat(80));
        System.out.printf("%-6s | %-20s | %-10s | %-12s | %-10s | %-8s\n",
            "Order", "Conformation", "Energy", "Boltzmann", "% of Z", "Rank");
        System.out.println("-".repeat(80));

        double cumulative = 0;
        for (int i = 0; i < Math.min(20, byOrder.size()); i++) {
            MinimizedConf conf = byOrder.get(i);
            BigDecimal weight = conf.getBoltzmannWeight();
            double pct = weight.divide(totalZ, 10, java.math.RoundingMode.HALF_UP).doubleValue() * 100;
            cumulative += pct;
            int rank = contributionRank.get(conf.confString);

            System.out.printf("%-6d | %-20s | %10.4f | %12.6e | %9.4f%% | #%-7d\n",
                conf.order,
                truncate(conf.confString, 20),
                conf.energy,
                weight.doubleValue(),
                pct,
                rank);
        }

        System.out.println("\n前20个conformations的累积贡献: " + String.format("%.2f%%", cumulative));

        // 分析不同percentiles的累积贡献
        System.out.println("\n" + "=".repeat(80));
        System.out.println("累积贡献分析:");
        System.out.println("-".repeat(80));

        int[] percentiles = {10, 25, 50, 75, 100};
        for (int pct : percentiles) {
            int nConfs = Math.max(1, (byOrder.size() * pct) / 100);
            BigDecimal cumulativeZ = BigDecimal.ZERO;
            for (int i = 0; i < nConfs && i < byOrder.size(); i++) {
                cumulativeZ = cumulativeZ.add(byOrder.get(i).getBoltzmannWeight());
            }

            double pctContribution = cumulativeZ.divide(totalZ, 10, java.math.RoundingMode.HALF_UP).doubleValue() * 100;
            System.out.printf("前 %3d%% (%4d confs) 贡献: %6.2f%% of total Z\n",
                pct, nConfs, pctContribution);
        }

        // 计算Spearman相关系数
        System.out.println("\n" + "=".repeat(80));
        System.out.println("相关性分析:");
        System.out.println("-".repeat(80));

        double spearman = calculateSpearman(byOrder, contributionRank);
        System.out.printf("Spearman秩相关系数: %.4f\n", spearman);

        // Top-K overlap分析
        int[] ks = {5, 10, 20, Math.min(50, byOrder.size())};
        for (int k : ks) {
            if (k > byOrder.size()) continue;

            Set<String> topKByOrder = new HashSet<>();
            for (int i = 0; i < k; i++) {
                topKByOrder.add(byOrder.get(i).confString);
            }

            Set<String> topKByContribution = new HashSet<>();
            for (int i = 0; i < k; i++) {
                topKByContribution.add(byContribution.get(i).confString);
            }

            topKByOrder.retainAll(topKByContribution);
            double overlap = (topKByOrder.size() * 100.0) / k;

            System.out.printf("Top-%d overlap: %d/%d (%.1f%%)\n", k, topKByOrder.size(), k, overlap);
        }

        // 写入CSV文件
        try (FileWriter writer = new FileWriter(outputFile)) {
            writer.write("MinimizeOrder,Conformation,Score,LowerBound,Corrected,FinalEnergy,OldConfLower,OldConfUpper," +
                "SubtreeLower,SubtreeUpper,ErrorBound,BoltzmannWeight,PercentContribution,ContributionRank\n");

            for (MinimizedConf conf : byOrder) {
                BigDecimal weight = conf.getBoltzmannWeight();
                double pct = weight.divide(totalZ, 10, java.math.RoundingMode.HALF_UP).doubleValue() * 100;
                int rank = contributionRank.get(conf.confString);

                writer.write(String.format("%d,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%e,%e,%e,%e,%.6f,%d\n",
                    conf.order,
                    conf.confString,
                    conf.score,
                    conf.lowerBound,
                    conf.corrected,
                    conf.energy,
                    conf.oldConfLower,
                    conf.oldConfUpper,
                    conf.subtreeLower.doubleValue(),
                    conf.subtreeUpper.doubleValue(),
                    conf.errorBound.doubleValue(),
                    weight.doubleValue(),
                    pct,
                    rank));
            }
        }

        System.out.println("\n详细数据已保存到: " + outputFile);

        // 最终结论
        System.out.println("\n" + "=".repeat(80));
        System.out.println("结论:");
        System.out.println("=".repeat(80));

        double top10PctContribution = 0;
        int top10Count = Math.max(1, byOrder.size() / 10);
        BigDecimal top10Z = BigDecimal.ZERO;
        for (int i = 0; i < top10Count; i++) {
            top10Z = top10Z.add(byOrder.get(i).getBoltzmannWeight());
        }
        top10PctContribution = top10Z.divide(totalZ, 10, java.math.RoundingMode.HALF_UP).doubleValue() * 100;

        if (spearman > 0.7 || top10PctContribution > 50) {
            System.out.println("✓ 优先minimize的conformations确实对Z贡献最大!");
            System.out.println("  - 前10%的confs贡献了 " + String.format("%.1f%%", top10PctContribution) + " 的Z");
            System.out.println("  - Spearman相关系数: " + String.format("%.3f", spearman));
        } else if (spearman > 0.3) {
            System.out.println("○ 有一定相关性，但不是非常强");
            System.out.println("  - 前10%的confs贡献了 " + String.format("%.1f%%", top10PctContribution) + " 的Z");
            System.out.println("  - Spearman相关系数: " + String.format("%.3f", spearman));
        } else {
            System.out.println("✗ minimize顺序和贡献大小关系不大");
            System.out.println("  - 前10%的confs只贡献了 " + String.format("%.1f%%", top10PctContribution) + " 的Z");
            System.out.println("  - Spearman相关系数: " + String.format("%.3f", spearman));
        }
        System.out.println("=".repeat(80));
    }

    private double calculateSpearman(List<MinimizedConf> byOrder, Map<String, Integer> contributionRank) {
        if (byOrder.size() < 2) return 0;

        double sumDSquared = 0;
        for (int i = 0; i < byOrder.size(); i++) {
            int orderRank = i + 1;
            int contribRank = contributionRank.get(byOrder.get(i).confString);
            int d = orderRank - contribRank;
            sumDSquared += d * d;
        }

        int n = byOrder.size();
        return 1.0 - (6.0 * sumDSquared) / (n * (n * n - 1));
    }

    private String truncate(String s, int maxLen) {
        if (s.length() <= maxLen) return s;
        return s.substring(0, maxLen - 3) + "...";
    }

    public List<MinimizedConf> getMinimizedConfs() {
        return new ArrayList<>(minimizedConfs);
    }
}
