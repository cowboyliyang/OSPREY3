package edu.duke.cs.osprey.minimization;

/**
 * A minimal line searcher that only does wall jumping.
 *
 * For each DOF, evaluates at ±1 degree (or unit) from current position
 * and jumps there if energy is lower. Only 2 energy evaluations per DOF
 * per iteration.
 *
 * Designed for multi warm start where GridDP already provides a good
 * starting point and we just need to hop over nearby energy barriers.
 */
public class WallJumpingLineSearcher implements LineSearcher {

	private ObjectiveFunction.OneDof f;

	@Override
	public void init(ObjectiveFunction.OneDof f) {
		this.f = f;
	}

	@Override
	public double search(double xd) {

		double xdmin = f.getXMin();
		double xdmax = f.getXMax();

		double fxd = f.getValue(xd);
		double xdstar = xd;
		double fxdstar = fxd;

		// Jump ±1 degree from current position
		double xdm = xd - 1;
		double xdp = xd + 1;

		if (xdm >= xdmin) {
			double fxdm = f.getValue(xdm);
			if (fxdm < fxdstar) {
				xdstar = xdm;
				fxdstar = fxdm;
			}
		}

		if (xdp <= xdmax) {
			double fxdp = f.getValue(xdp);
			if (fxdp < fxdstar) {
				xdstar = xdp;
				fxdstar = fxdp;
			}
		}

		f.setX(xdstar);
		return xdstar;
	}
}
