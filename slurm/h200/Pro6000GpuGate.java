import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import static org.junit.platform.engine.discovery.DiscoverySelectors.selectMethod;
public class Pro6000GpuGate {
    public static void main(String[] args) {
        var request = LauncherDiscoveryRequestBuilder.request();
        for (String method : new String[]{"gpuExactTriples_nonLeafRegularMatchesJava", "gpuExactTriples_childSlicedMatchesJava", "gpuExactTriples_hybridMatchesJava", "gpuExactTriples_outOfCoreMatchesJava", "gpuOutOfCore_multiGpuMatchesJava"})
            request.selectors(selectMethod("edu.duke.cs.osprey.branchdp.TestGpuFullDP", method));
        var listener = new SummaryGeneratingListener();
        var launcher = LauncherFactory.create();
        launcher.registerTestExecutionListeners(listener);
        launcher.execute(request.build());
        var summary = listener.getSummary();
        summary.printTo(new java.io.PrintWriter(System.out));
        summary.printFailuresTo(new java.io.PrintWriter(System.out));
        if (summary.getTestsSucceededCount() != 5 || summary.getTestsFailedCount() != 0 || summary.getTestsAbortedCount() != 0)
            throw new IllegalStateException("All five actual GPU checks must pass");
    }
}
