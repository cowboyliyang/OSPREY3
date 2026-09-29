import edu.duke.cs.osprey.kstar.KStarScore;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import static org.junit.platform.engine.discovery.DiscoverySelectors.selectClass;

/** Run through Slurm, using only the patched KStarScore over the frozen runtime. */
public class Log10OutputRegression {
    public static void main(String[] args) throws Exception {
        var request = LauncherDiscoveryRequestBuilder.request()
                .selectors(selectClass("edu.duke.cs.osprey.kstar.TestKStarScore"),
                           selectClass("edu.duke.cs.osprey.kstar.TestKStarScoreLog10"))
                .build();
        var listener = new SummaryGeneratingListener();
        var launcher = LauncherFactory.create();
        launcher.registerTestExecutionListeners(listener);
        launcher.execute(request);
        var summary = listener.getSummary();
        summary.printTo(new java.io.PrintWriter(System.out));
        summary.printFailuresTo(new java.io.PrintWriter(System.out));
        if (summary.getTestsSucceededCount() < 6 || summary.getTestsFailedCount() != 0
                || summary.getTestsAbortedCount() != 0 || summary.getTestsSkippedCount() != 0) {
            throw new AssertionError("KStarScore regression tests failed");
        }
        int count = 0;
        long start = System.nanoTime();
        for (String line : Files.readAllLines(Path.of(args[0]))) {
            if (line.startsWith("decimal\t")) continue;
            String[] fields = line.split("\t");
            double actual = KStarScore.scoreToLog10(new BigDecimal(fields[0]));
            double expected = Double.parseDouble(fields[1]);
            if (!Double.isFinite(actual) || Math.abs(actual - expected) > 1e-6) {
                throw new AssertionError("Archived log mismatch: " + line + " actual=" + actual);
            }
            count++;
        }
        if (count == 0) throw new AssertionError("No archived Z values replayed");
        System.out.printf("PASS archived Z replay: %d bounds in %.3f ms%n", count, (System.nanoTime() - start)/1e6);
        System.out.println("PRODUCTION_CLASS=" + KStarScore.class.getProtectionDomain().getCodeSource().getLocation());
    }
}
