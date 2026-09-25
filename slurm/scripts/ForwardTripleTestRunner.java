import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import static org.junit.platform.engine.discovery.DiscoverySelectors.selectClass;

/** Slurm-only focused verification of the final compiled overlay. */
public class ForwardTripleTestRunner {
    public static void main(String[] args) {
        var request = LauncherDiscoveryRequestBuilder.request()
            .selectors(selectClass("edu.duke.cs.osprey.packstar.TestPackStarForwardTripleSelection"),
                       selectClass("edu.duke.cs.osprey.packstar.TestPackStarTripleDecompositionCosts"),
                       selectClass("edu.duke.cs.osprey.packstar.TestPackStarProposalLearning"),
                       selectClass("edu.duke.cs.osprey.packstar.TestPackStarTripleEtaCorrections"),
                       selectClass("edu.duke.cs.osprey.packstar.TestPackStarFrequencySeverityPAC"),
                       selectClass("edu.duke.cs.osprey.packstar.TestPackStarConfig"),
                       selectClass("edu.duke.cs.osprey.branchdp.TestBranchDecompositionStrategies"))
            .build();
        var listener = new SummaryGeneratingListener();
        var launcher = LauncherFactory.create();
        launcher.registerTestExecutionListeners(listener);
        launcher.execute(request);
        listener.getSummary().printTo(new java.io.PrintWriter(System.out));
        listener.getSummary().printFailuresTo(new java.io.PrintWriter(System.out));
        if (listener.getSummary().getTestsFoundCount() != 58
                || listener.getSummary().getTestsFailedCount() != 0)
            throw new IllegalStateException("focused forward-selector tests failed");
    }
}
