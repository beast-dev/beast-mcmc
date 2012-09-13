/*
 * DeltaMixOperator.java
 *
 * Copyright (C) 2010 Alexander V. Alekseyenko
 *
 */

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * A generic operator for use with a sum-constrained (possibly weighted) vector parameter. This proposal is
 * non-symmetrical, performs better than DeltaExchangeOperator when the components are on vastly different scales,
 * does not propose states outside the sum-constraints and thus does not require bounds checking on the proposed state
 *
 * @author Alexander V. Alekseyenko
 * @version $Id: DeltaMixOperator.java,v 1.6 2010/09/23 $
 */
public class DeltaMixOperator extends AbstractCoercableOperator {

    public DeltaMixOperator(Parameter parameter, double delta) {

        super(CoercionMode.COERCION_ON);

        this.parameter = parameter;
        this.delta = delta;
        setWeight(1.0);

        parameterWeights = new int[parameter.getDimension()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
    }

    public DeltaMixOperator(Parameter parameter, int[] parameterWeights, double delta, double weight, CoercionMode mode) {

        super(mode);

        this.parameter = parameter;
        this.delta = delta;
        setWeight(weight);
        this.parameterWeights = parameterWeights;
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     * performs a delta exchange operation between two scalars in the vector
     * and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        // get two dimensions
        final int dim = parameter.getDimension();
        final int dim1 = MathUtils.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = MathUtils.nextInt(dim);
        }

        double scalar1 = parameter.getParameterValue(dim1);
        double scalar2 = parameter.getParameterValue(dim2);

        // exchange a random fraction of dim1 component with dim2
        final double d = MathUtils.nextDouble() * delta * scalar1;
        scalar1 -= d;
        if (parameterWeights[dim1] != parameterWeights[dim2]) {
            scalar2 += d * (double) parameterWeights[dim1] / (double) parameterWeights[dim2];
        } else {
            scalar2 += d;
        }

        parameter.setParameterValue(dim1, scalar1);
        parameter.setParameterValue(dim2, scalar2);

        return Math.log(scalar2 / (scalar1 + d));
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return parameter.getParameterName();
    }

    public double getCoercableParameter() {
        return Math.log(delta);
    }

    public void setCoercableParameter(double value) {
        delta = Math.exp(value);
    }

    public double getRawParameter() {
        return delta;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double d = OperatorUtils.optimizeWindowSize(delta, parameter.getParameterValue(0) * 2.0, prob, targetProb);


        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing delta to about " + d;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing delta to about " + d;
        } else return "";
    }

    public String toString() {
        return getOperatorName() + "(windowsize=" + delta + ")";
    }

    // Private instance variables

    private Parameter parameter = null;
    private final int[] parameterWeights;
    private double delta = 0.2;
}