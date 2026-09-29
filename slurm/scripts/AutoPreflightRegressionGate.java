import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import static org.junit.platform.engine.discovery.DiscoverySelectors.selectClass;

/** Run the focused regressions against an incremental, isolated class overlay. */
public class AutoPreflightRegressionGate {
    public static void main(String[] args) {
        var request = LauncherDiscoveryRequestBuilder.request();
        for (String type : new String[]{
                "packstar.TestPackStarAutoPreflight", "packstar.TestPackStarAdmissionDecision",
                "branchdp.TestBranchDpAdmission", "branchdp.TestBranchDecompositionStrategies",
                "packstar.TestPackStarTripleDecompositionCosts"})
            request.selectors(selectClass("edu.duke.cs.osprey." + type));
        var listener = new SummaryGeneratingListener();
        var launcher = LauncherFactory.create();
        launcher.registerTestExecutionListeners(listener);
        launcher.execute(request.build());
        var summary = listener.getSummary();
        summary.printTo(new java.io.PrintWriter(System.out));
        summary.printFailuresTo(new java.io.PrintWriter(System.out));
        if (summary.getTestsSucceededCount() < 42 || summary.getTestsFailedCount() != 0
                || summary.getTestsAbortedCount() != 0)
            throw new IllegalStateException("All automatic-preflight regressions must pass");
    }
}
