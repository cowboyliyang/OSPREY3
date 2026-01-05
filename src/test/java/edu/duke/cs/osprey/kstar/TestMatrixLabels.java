package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.stream.Collectors;

public class TestMatrixLabels {

	@Test
	public void testMinimalMatrixLabels() {

		// Use a very small system - just the 2RL0 test case
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		// Create minimal states
		CometsZ.State protein = new CometsZ.State("Protein", confSpaces.protein);
		CometsZ.State ligand = new CometsZ.State("Ligand", confSpaces.ligand);
		CometsZ.State complex = new CometsZ.State("Complex", confSpaces.complex);

		List<CometsZ.State> states = List.of(complex, protein, ligand);

		List<SimpleConfSpace> confSpaceList = states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaceList, confSpaces.ffparams)
			.setParallelism(Parallelism.makeCpu(4))  // Use fewer threads
			.build()
		) {

			for (CometsZ.State state : states) {

				System.out.println("\n========================================");
				System.out.println("Processing state: " + state.name);
				System.out.println("========================================");

				// Set context label for reference energies
				SimplerEnergyMatrixCalculator.setContextLabel(state.name);

				// Calculate reference energies
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// Calculate minimizing energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				SimplerEnergyMatrixCalculator.clearContextLabel();

				System.out.println("Finished state: " + state.name);
			}

			System.out.println("\n========================================");
			System.out.println("All matrices computed successfully!");
			System.out.println("========================================");
		}
	}
}
