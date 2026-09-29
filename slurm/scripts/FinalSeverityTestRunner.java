import edu.duke.cs.osprey.packstar.PackStarFrequencySeverityPAC;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import static org.junit.platform.engine.discovery.DiscoverySelectors.selectClass;

/** Focused tests and an independent numerical replay; execute through Slurm. */
public class FinalSeverityTestRunner {
    private static double number(String value) {
        if (value.equals("nan")) return Double.NaN;
        if (value.equals("inf")) return Double.POSITIVE_INFINITY;
        return Double.parseDouble(value);
    }

    public static void main(String[] args) throws Exception {
        var request = LauncherDiscoveryRequestBuilder.request()
                .selectors(selectClass("edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC"),
                           selectClass("edu.duke.cs.osprey.packstar.TestPackStarConfig"))
                .build();
        var listener = new SummaryGeneratingListener();
        var launcher = LauncherFactory.create();
        launcher.registerTestExecutionListeners(listener);
        launcher.execute(request);
        var summary = listener.getSummary();
        summary.printTo(new java.io.PrintWriter(System.out));
        summary.printFailuresTo(new java.io.PrintWriter(System.out));
        if (summary.getTestsFoundCount() < 18 || summary.getTestsFailedCount() != 0
                || summary.getTestsAbortedCount() != 0 || summary.getTestsSkippedCount() != 0) {
            throw new IllegalStateException("focused severity tests failed or did not run");
        }
        int count = 0;
        int rejected = 0;
        boolean saw2rl0 = false;
        boolean saw3cal = false;
        boolean saw4u3s = false;
        boolean saw4wyu = false;
        try (var reader = Files.newBufferedReader(Path.of(args[0]))) {
            String header = reader.readLine();
            if (!header.startsWith("case_id\tlog_clip\tcap\talpha\ttail_count\texpected_log_p\texpected_rejected\tlog_weights")) {
                throw new IllegalArgumentException("unexpected replay header");
            }
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t", -1);
                double[] weights = Arrays.stream(fields[7].split(","))
                        .mapToDouble(FinalSeverityTestRunner::number).toArray();
                var test = PackStarFrequencySeverityPAC.testConditionalMeanExcess(
                        weights, number(fields[1]), number(fields[2]), number(fields[3]));
                double expectedLogP = number(fields[5]);
                if (test.exceedanceCount != Integer.parseInt(fields[4])
                        || test.rejected != Boolean.parseBoolean(fields[6])
                        || Double.isNaN(test.logPValue) != Double.isNaN(expectedLogP)
                        || (!Double.isNaN(expectedLogP)
                            && Math.abs(test.logPValue - expectedLogP) > 1.0e-8)) {
                    throw new AssertionError("Python/Java replay mismatch: " + fields[0]);
                }
                count++;
                if (test.rejected) {
                    rejected++;
                    saw2rl0 |= fields[0].startsWith("rb2:2rl0:");
                    saw3cal |= fields[0].startsWith("rb1:3cal:");
                    saw4u3s |= fields[0].startsWith("rb1:4u3s:");
                    saw4wyu |= fields[0].startsWith("rb1:4wyu:");
                    System.out.println("REJECT " + fields[0] + " tail=" + test.exceedanceCount
                            + " p=" + test.pValueUpper);
                }
            }
        }
        if (count != 6065 || rejected != 4 || !saw2rl0 || !saw3cal || !saw4u3s || !saw4wyu) {
            throw new AssertionError("unexpected frozen replay coverage or decisions");
        }
        System.out.println("PASS: Java replay agrees with independent Python calculation for "
                + count + " archived final batches; rejected=" + rejected);
    }
}
