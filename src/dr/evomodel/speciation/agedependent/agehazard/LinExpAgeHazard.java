package dr.evomodel.speciation.agedependent.agehazard;

import dr.inference.model.Parameter;

/**
 * @author Frederik M. Andersen
 *
 * Linear-exponential age hazard:
 *     h(a, epoch) = (1 + r[epoch] * gamma[epoch] * a) * exp(-gamma[epoch] * a)
 *
 * r=0 or gamma=0 collapses to a constant-rate process.
 * r > 1 produces a hump with peak at a* = (r-1)/(r*gamma)
 * of amplitude r * exp(1/r - 1).
 *
 * Both r and gamma are Parameters with dimension 1 (shared across all epochs)
 * or dimension numEpochs (one value per epoch).
 */
public class LinExpAgeHazard extends AgeHazard {

    private final Parameter r;
    private final Parameter gamma;

    public LinExpAgeHazard(Parameter r, Parameter gamma) {
        super("linExpAgeHazard");
        this.r = r;
        this.gamma = gamma;
        addVariable(r);
        if (gamma != r) addVariable(gamma);
    }

    @Override
    public int getEpochCount() {
        return Math.max(r.getDimension(), gamma.getDimension());
    }

    private int idx(Parameter p, int epoch) {
        return p.getDimension() == 1 ? 0 : epoch;
    }

    @Override
    public double evaluate(double age, int epoch) {
        double rv = r.getParameterValue(idx(r, epoch));
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        return (1.0 + rv * gv * age) * Math.exp(-gv * age);
    }

    @Override
    public double maxHazard(double originTime, int epoch) {
        double rv = r.getParameterValue(idx(r, epoch));
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        double b  = rv * gv;
        double hMax = Math.max(1.0, evaluate(originTime, epoch));
        if (b > 0.0 && gv > 0.0) {
            double aStar = (b - gv) / (gv * b);
            if (aStar > 0.0 && aStar < originTime) {
                hMax = Math.max(hMax, evaluate(aStar, epoch));
            }
        }
        return hMax;
    }

    /** Overrides the default loop: reads parameters once rather than per-cell. */
    @Override
    public void evaluate(int Na, double da, double[] result, int epoch) {
        double rv  = r.getParameterValue(idx(r, epoch));
        double gv  = gamma.getParameterValue(idx(gamma, epoch));
        double lin = rv * gv;
        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            result[j] = (1.0 + lin * a) * Math.exp(-gv * a);
        }
    }
}
