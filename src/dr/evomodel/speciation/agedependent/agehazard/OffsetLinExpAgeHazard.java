package dr.evomodel.speciation.agedependent.agehazard;

import dr.inference.model.Parameter;

/**
 * @author Frederik M. Andersen
 *
 * Offset linear-exponential age hazard:
 *     h(a, epoch) = h0[epoch] + (1 + r[epoch] * gamma[epoch] * a) * exp(-gamma[epoch] * a)
 *
 * h0 >= 0 is a sustained baseline ("plateau"): as a -> infinity, h -> h0, so the hazard
 * retains a non-zero long-age tail (sustained transmission over the chronic phase) on top
 * of the early linExp hump/excess. h0 = 0 recovers {@link LinExpAgeHazard}.
 *
 * Note: a non-zero h0 means the hazard never falls below h0, so the rate-zero age-axis
 * truncation in the PDE model will not shrink the age grid (NaTrunc stays at Na). That is
 * the correct behaviour for a sustained tail, but it removes the truncation speed-up.
 *
 * h0, r and gamma are Parameters with dimension 1 (shared across all epochs) or dimension
 * numEpochs (one value per epoch); mixed dimensions are permitted.
 */
public class OffsetLinExpAgeHazard extends AgeHazard {

    private final Parameter h0;
    private final Parameter r;
    private final Parameter gamma;

    public OffsetLinExpAgeHazard(Parameter h0, Parameter r, Parameter gamma) {
        super("offsetLinExpAgeHazard");
        this.h0 = h0;
        this.r = r;
        this.gamma = gamma;
        addVariable(h0);
        if (r != h0) addVariable(r);
        if (gamma != h0 && gamma != r) addVariable(gamma);
    }

    @Override
    public int getEpochCount() {
        return Math.max(h0.getDimension(), Math.max(r.getDimension(), gamma.getDimension()));
    }

    @Override
    public double evaluate(double age, int epoch) {
        double h  = h0.getParameterValue(idx(h0, epoch));
        double rv = r.getParameterValue(idx(r, epoch));
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        return h + (1.0 + rv * gv * age) * Math.exp(-gv * age);
    }

    @Override
    public double maxHazard(double originTime, int epoch) {
        double h  = h0.getParameterValue(idx(h0, epoch));
        double rv = r.getParameterValue(idx(r, epoch));
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        double b  = rv * gv;
        double linMax = Math.max(1.0, (1.0 + b * originTime) * Math.exp(-gv * originTime));
        if (b > 0.0 && gv > 0.0) {
            double aStar = (b - gv) / (gv * b);
            if (aStar > 0.0 && aStar < originTime) {
                linMax = Math.max(linMax, (1.0 + b * aStar) * Math.exp(-gv * aStar));
            }
        }
        return h + linMax;
    }

    /** Overrides the default loop: reads parameters once rather than per-cell. */
    @Override
    public void evaluate(int Na, double da, double[] result, int epoch) {
        double h   = h0.getParameterValue(idx(h0, epoch));
        double rv  = r.getParameterValue(idx(r, epoch));
        double gv  = gamma.getParameterValue(idx(gamma, epoch));
        double lin = rv * gv;
        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            result[j] = h + (1.0 + lin * a) * Math.exp(-gv * a);
        }
    }
}
