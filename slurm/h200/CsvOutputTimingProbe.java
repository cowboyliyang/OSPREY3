import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.bench.GenericPDBBench;
import edu.duke.cs.osprey.tools.MathTools;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/** Timing fixtures reconstructed from rounded CSV bounds, not scientific reruns. */
public class CsvOutputTimingProbe {
    private static BigDecimal decimal(String text) {
        if (text.isEmpty()) return null;
        double log = Double.parseDouble(text);
        if (log == Double.NEGATIVE_INFINITY) return BigDecimal.ZERO;
        if (log == Double.POSITIVE_INFINITY) return MathTools.BigPositiveInfinity;
        if (Double.isNaN(log)) return MathTools.BigNaN;
        int exponent = (int)Math.floor(log);
        return BigDecimal.valueOf(Math.pow(10, log - exponent)).scaleByPowerOfTen(exponent);
    }

    private static PartitionFunction.Result result(String[] row, int offset) {
        if (row[offset + 2].equals("N/A")) return null;
        BigDecimal lower = decimal(row[offset]);
        BigDecimal upper = decimal(row[offset + 1]);
        var values = new PartitionFunction.Values();
        values.qstar = lower == null ? MathTools.BigNaN : lower;
        values.qprime = MathTools.bigSubtract(upper == null ? MathTools.BigNaN : upper,
                values.qstar, PartitionFunction.decimalPrecision);
        var result = new PartitionFunction.Result(PartitionFunction.Status.valueOf(row[offset + 2]),
                values, Integer.parseInt(row[offset + 4]));
        String[] stats = {"numConfsScored", "numPartialMinimizations"};
        for (int i = 0; i < stats.length; i++) result.setStat(stats[i], Long.parseLong(row[offset + 5 + i]));
        return result;
    }

    private static Sequence sequence(String text) throws Exception {
        var constructor = SeqSpace.class.getDeclaredConstructor();
        constructor.setAccessible(true);
        var space = constructor.newInstance();
        var makePos = SeqSpace.class.getDeclaredMethod("makePos", String.class, String.class, List.class);
        makePos.setAccessible(true);
        String[] types = text.trim().split("\\s+");
        for (int i = 0; i < types.length; i++) makePos.invoke(space, "A" + i, types[i], List.of(types[i]));
        return space.makeWildTypeSequence();
    }

    public static void main(String[] args) throws Exception {
        Path audit = Path.of(args[0]);
        Path output = Path.of(args[1]);
        var writer = GenericPDBBench.class.getDeclaredMethod("writeKStarResults", List.class,
                String.class, String.class, String.class, double.class, double.class);
        writer.setAccessible(true);
        List<String> records = Files.readAllLines(audit);
        var header = Arrays.asList(records.get(0).split("\t", -1));
        try (var report = new PrintWriter(Files.newBufferedWriter(output.resolve("csv_timings.tsv")))) {
            report.println("index\tdesign_id\trows\tfirst_export_ms\trepeat_export_ms\toriginal_gap_s\tselected_extreme");
            for (String line : records.subList(1, records.size())) {
                String[] item = line.split("\t", -1);
                if (!item[header.indexOf("stage")].equals("CSV_WRITTEN")) continue;
                String design = item[header.indexOf("design_id")];
                Path input = Path.of(item[header.indexOf("path")]).resolve(design + "_packstar.csv");
                var rows = Files.readAllLines(input);
                var csvHeader = Arrays.asList(rows.get(0).split(",", -1));
                var scores = new ArrayList<KStar.ScoredSequence>();
                for (String raw : rows.subList(1, rows.size())) {
                    String[] row = raw.split(",", -1);
                    var score = new KStarScore(decimal(row[4]), decimal(row[5]), decimal(row[6]));
                    // Fixtures only: attach original-state bounds to the immutable score holder.
                    for (int i = 0; i < 3; i++) {
                        var field = KStarScore.class.getField(new String[]{"protein", "ligand", "complex"}[i]);
                        field.setAccessible(true);
                        int offset = csvHeader.indexOf(new String[]{"prot", "lig", "comp"}[i]
                                + "_qstar_lb_log10");
                        if (offset < 0) throw new IllegalArgumentException("Missing state columns: " + input);
                        field.set(score, result(row, offset));
                    }
                    scores.add(new KStar.ScoredSequence(sequence(row[1]), score));
                }
                double[] times = new double[2];
                for (int repeat = 0; repeat < 2; repeat++) {
                    Collections.shuffle(scores, new Random(42));
                    long start = System.nanoTime();
                    writer.invoke(null, scores, design + "_probe" + repeat, "packstar", output.toString(), 0.683, 0.0);
                    times[repeat] = (System.nanoTime() - start) / 1e6;
                    Path csv = output.resolve(design + "_probe" + repeat + "_packstar.csv");
                    if (Files.readAllLines(csv).size() != rows.size()) throw new AssertionError("Missing CSV rows: " + design);
                }
                report.printf(Locale.ROOT, "%s\t%s\t%d\t%.3f\t%.3f\t%s\t%s%n",
                        item[0], design, scores.size(), times[0], times[1],
                        item[header.indexOf("output_gap_s")], item[header.indexOf("any_extreme")]);
            }
        }
    }
}
