package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractAdaptableTreeOperator extends AbstractTreeOperator implements AdaptableMCMCOperator {

    public final AdaptationMode mode;
    private final double targetAcceptanceProbability;
    private long adaptationCount = 0;

    public AbstractAdaptableTreeOperator(AdaptationMode mode) {
        this(mode, DEFAULT_ADAPTATION_TARGET);
    }

    public AbstractAdaptableTreeOperator(AdaptationMode mode, double targetAcceptanceProbability) {
        this.mode = mode;
        if (System.getProperty("mcmc.adaptation_target") != null) {
            this.targetAcceptanceProbability = Double.parseDouble(System.getProperty("mcmc.adaptation_target"));
        } else {
            this.targetAcceptanceProbability = targetAcceptanceProbability;
        }
    }

    @Override
    public void setAdaptableParameter(double value) {
        adaptationCount ++;
        setAdaptableParameterValue(value);
    }

    @Override
    public double getAdaptableParameter() {
        return getAdaptableParameterValue();
    }

    @Override
    public long getAdaptationCount() {
        return adaptationCount;
    }

    @Override
    public void setAdaptationCount(long count) {
        adaptationCount = count;
    }

    /**
     * Sets the adaptable parameter value.
     *
     * @param value the value to set the adaptable parameter to
     */
    protected abstract void setAdaptableParameterValue(double value);

    /**
     * Gets the adaptable parameter value.
     *
     * @returns the value
     */
    protected abstract double getAdaptableParameterValue();

    @Override
    public AdaptationMode getMode() {
        return mode;
    }

    @Override
    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    @Override
    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    @Override
    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    @Override
    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    @Override
    public final String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public double getTargetAcceptanceProbability() {
        return targetAcceptanceProbability;
    }
}
