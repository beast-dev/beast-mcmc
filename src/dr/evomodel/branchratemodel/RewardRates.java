package dr.evomodel.branchratemodel;

import dr.inference.model.Parameter;

public final class RewardRates {

    private final Parameter values;
    private final Parameter fixedValues;
    private final Parameter varyingValues;
    private final Parameter stateIndices;

    public RewardRates(Parameter values, Parameter fixedValues, Parameter varyingValues, Parameter stateIndices) {
        if (values == null) {
            throw new IllegalArgumentException("values must be non-null");
        }
        if (varyingValues == null) {
            throw new IllegalArgumentException("varyingValues must be non-null");
        }
        if (stateIndices == null) {
            throw new IllegalArgumentException("stateIndices must be non-null");
        }

        this.values = values;
        this.fixedValues = fixedValues;
        this.varyingValues = varyingValues;
        this.stateIndices = stateIndices;
    }

    public Parameter getValues() {
        return values;
    }

    public Parameter getFixedValues() {
        return fixedValues;
    }

    public Parameter getVaryingValues() {
        return varyingValues;
    }

    public Parameter getStateIndices() {
        return stateIndices;
    }

    /** Raw reward value for atomic state stateIndex, via the stateIndices -> values mapping. */
    public double getRawReward(final int stateIndex) {
        if (stateIndex < 0 || stateIndex >= stateIndices.getDimension()) {
            throw new IllegalArgumentException("stateIndex out of range: " + stateIndex);
        }
        final int rewardRateIndex = (int) Math.round(stateIndices.getParameterValue(stateIndex));
        if (rewardRateIndex < 0 || rewardRateIndex >= values.getDimension()) {
            throw new IllegalArgumentException(
                    "Reward-rate mapping for state " + stateIndex + " points outside rewardRates: " +
                            rewardRateIndex);
        }
        return values.getParameterValue(rewardRateIndex);
    }
}
