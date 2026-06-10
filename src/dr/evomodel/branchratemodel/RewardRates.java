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
}
