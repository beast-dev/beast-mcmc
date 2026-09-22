package dr.evomodel.speciation.agedependent.agehazard;

import dr.inference.model.Parameter;

/**
 * @author Frederik M. Andersen
 *
 * Exponentially decaying age hazard:
 *     h(a, epoch) = exp(-gamma[epoch] * a)
 *
 * The hazard starts at 1 at age 0 and decays monotonically with rate gamma, so a
 * lineage is most active immediately after it arises. gamma = 0 collapses to a
 * constant-rate process; negative gamma gives a monotonically increasing hazard.
 * This is the r = 0 special case of {@link LinExpAgeHazard}, kept separate so the
 * one-parameter form can be used without an r Parameter to fix or sample.
 *
 * gamma is a Parameter with dimension 1 (shared across all epochs) or dimension
 * numEpochs (one value per epoch).
 */
public class ExpAgeHazard extends AgeHazard {

    private final Parameter gamma;

    public ExpAgeHazard(Parameter gamma) {
        super("expAgeHazard");
        this.gamma = gamma;
        addVariable(gamma);
    }

    @Override
    public int getEpochCount() {
        return gamma.getDimension();
    }

    @Override
    public double evaluate(double age, int epoch) {
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        return Math.exp(-gv * age);
    }

    @Override
    public double maxHazard(double originTime, int epoch) {
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        return Math.max(1.0, Math.exp(-gv * originTime));
    }

    @Override
    public void evaluate(int Na, double da, double[] result, int epoch) {
        double gv = gamma.getParameterValue(idx(gamma, epoch));
        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            result[j] = Math.exp(-gv * a);
        }
    }
}
