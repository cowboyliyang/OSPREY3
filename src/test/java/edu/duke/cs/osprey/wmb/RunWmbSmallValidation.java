package edu.duke.cs.osprey.wmb;

import org.junit.platform.engine.discovery.DiscoverySelectors;
import org.junit.platform.launcher.Launcher;
import org.junit.platform.launcher.LauncherDiscoveryRequest;
import org.junit.platform.launcher.core.LauncherDiscoveryRequestBuilder;
import org.junit.platform.launcher.core.LauncherFactory;
import org.junit.platform.launcher.listeners.SummaryGeneratingListener;
import org.junit.platform.launcher.listeners.TestExecutionSummary;

public class RunWmbSmallValidation {

	public static void main(String[] args) {
		LauncherDiscoveryRequest request = LauncherDiscoveryRequestBuilder.request()
				.selectors(DiscoverySelectors.selectMethod(
						TestWeightedMiniBucket.class,
						"looseProposalImportanceSamplingMatchesEnumeratedPartitionFunction"))
				.build();

		SummaryGeneratingListener listener = new SummaryGeneratingListener();
		Launcher launcher = LauncherFactory.create();
		launcher.registerTestExecutionListeners(listener);
		launcher.execute(request);

		TestExecutionSummary summary = listener.getSummary();
		System.out.println("WMB-IS small validation tests found: " + summary.getTestsFoundCount());
		System.out.println("WMB-IS small validation tests succeeded: " + summary.getTestsSucceededCount());
		System.out.println("WMB-IS small validation tests failed: " + summary.getTestsFailedCount());
		for (TestExecutionSummary.Failure failure : summary.getFailures()) {
			System.err.println(failure.getTestIdentifier().getDisplayName());
			failure.getException().printStackTrace(System.err);
		}

		if (summary.getTestsFailedCount() > 0 || summary.getTestsSucceededCount() == 0) {
			System.exit(1);
		}
	}
}
