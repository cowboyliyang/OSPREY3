package edu.duke.cs.osprey.minimization;

/**
 * A lean line searcher optimized for GridDP warm-started CCD.
 *
 * When starting from a good initial point (e.g., GridDP output),
 * we're already in the correct basin and only need fine-tuning.
 * This searcher does pure quadratic fitting without:
 * - Wall jumping (unnecessary when already in basin)
 * - Surfing (minimal improvement near minimum)
 *
 * Result: 3 energy evaluations per DOF per iteration (vs 5+ for SurfingLineSearcher).
 */
public class QuadraticLineSearcher implements LineSearcher {

	private static final double ShapeEpsilon = 1e-12;

	private ObjectiveFunction.OneDof f;

	@Override
	public void init(ObjectiveFunction.OneDof f) {
		this.f = f;
	}

	@Override
	public double search(double xd) {

		double xdmin = f.getXMin();
		double xdmax = f.getXMax();

		// Eval 1: current position
		double fxd = f.getValue(xd);

		// Step size: just use the initial step size, no adaptive scaling needed
		double step = f.getInitialStepSize();

		// Shrink step if it would go out of bounds on both sides
		while (xd - step < xdmin && xd + step > xdmax) {
			step /= 2;
		}

		// Eval 2 & 3: positive and negative neighbors
		double xdp = xd + step;
		double xdm = xd - step;
		double fxdp = Double.POSITIVE_INFINITY;
		if (xdp <= xdmax) {
			fxdp = f.getValue(xdp);
		}
		double fxdm = Double.POSITIVE_INFINITY;
		if (xdm >= xdmin) {
			fxdm = f.getValue(xdm);
		}

		// Fit quadratic: q(x) = fx + a*(x-xd)^2 + b*(x-xd)
		double shape = fxdp + fxdm - 2 * fxd;
		double xdstar;

		if (shape < -ShapeEpsilon || Double.isNaN(shape) || Double.isInfinite(shape)) {
			// Concave down or invalid: pick the better endpoint
			xdstar = (fxdm < fxdp) ? xdm : xdp;
		} else if (shape <= ShapeEpsilon) {
			// Flat: stay put
			xdstar = xd;
		} else {
			// Concave up: step to quadratic minimum
			xdstar = xd + (fxdm - fxdp) * step / (2 * shape);
		}

		// Clamp to bounds
		if (xdstar < xdmin) xdstar = xdmin;
		if (xdstar > xdmax) xdstar = xdmax;

		// If xdstar is different from all 3 sampled points, evaluate it
		// Otherwise we already know the energy at one of the sampled points
		double fxdstar;
		if (xdstar == xd) {
			fxdstar = fxd;
		} else if (xdstar == xdm) {
			fxdstar = fxdm;
		} else if (xdstar == xdp) {
			fxdstar = fxdp;
		} else {
			fxdstar = f.getValue(xdstar); // Eval 4 (only when needed)
		}

		// Only accept if we improved
		if (fxdstar >= fxd) {
			// Quadratic step didn't help — check if either neighbor is better
			if (fxdm < fxd && fxdm <= fxdp) {
				xdstar = xdm;
			} else if (fxdp < fxd && fxdp < fxdm) {
				xdstar = xdp;
			} else {
				xdstar = xd;
			}
		}

		// Leave the molecule at the final position (LineSearcher contract)
		f.setX(xdstar);
		return xdstar;
	}
}
