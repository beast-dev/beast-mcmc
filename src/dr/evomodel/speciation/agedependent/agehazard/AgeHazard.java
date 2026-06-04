package dr.evomodel.speciation.agedependent.agehazard;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Frederik M. Andersen
 *
 * Pluggable age-dependent hazard profile h(a, epoch).
 * The full rate is scale(t) * h(a, epoch).
 *
 * A single instance covers all epochs. Implementations hold their own Parameters
 * (with dimension 1 for a hazard shared across all epochs, or numEpochs for
 * per-epoch hazards) and fire model-changed events upward so the parent
 * likelihood model can call makeDirty().
 */
public abstract class AgeHazard extends AbstractModel {

    public AgeHazard(String name) {
        super(name);
    }

    /** Hazard at a single age for the given epoch. */
    public abstract double evaluate(double age, int epochIndex);

    /** Upper bound on h(a, epochIndex) over [0, originTime]. */
    public abstract double maxHazard(double originTime, int epochIndex);

    /**
     * Fill result[0..Na] with h(j*da, epochIndex) for j=0..Na.
     * Override in subclasses for efficiency (avoids per-call parameter reads).
     */
    public void evaluate(int Na, double da, double[] result, int epochIndex) {
        for (int j = 0; j <= Na; j++) {
            result[j] = evaluate(j * da, epochIndex);
        }
    }

    /**
     * Number of epochs this hazard covers: 1 means it is constant across all epochs,
     * N means it has distinct parameters for each of N epochs.
     */
    public abstract int getEpochCount();

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index,
                                               Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {}
    protected void restoreState() {}
    protected void acceptState() {}
}
