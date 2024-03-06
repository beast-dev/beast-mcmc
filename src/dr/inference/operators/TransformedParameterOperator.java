package dr.inference.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.model.Parameter;
import dr.inference.model.TransformedParameter;

public class TransformedParameterOperator extends AbstractAdaptableOperator {
    private final boolean isAdaptable;
    private final SimpleMCMCOperator subOperator;
    private final TransformedParameter parameter;
    private final boolean checkValid;
    private final BoundedSpace generalBounds;
    public static final String TRANSFORMED_OPERATOR = "transformedParameterOperator";

    public TransformedParameterOperator(SimpleMCMCOperator operator, BoundedSpace generalBounds) {

        this.subOperator = operator;
        setWeight(operator.getWeight());
        this.isAdaptable = operator instanceof AbstractAdaptableOperator;
        this.parameter = (TransformedParameter) operator.getParameter();

        this.generalBounds = generalBounds;
        this.checkValid = generalBounds != null;
    }


    @Override
    protected void setAdaptableParameterValue(double value) {
        if (isAdaptable) {
            ((AbstractAdaptableOperator) subOperator).setAdaptableParameterValue(value);
        }
    }

    @Override
    protected double getAdaptableParameterValue() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterValue();
        }
        return 0;
    }

    @Override
    public double getRawParameter() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getRawParameter();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getAdaptableParameterName() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterName();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getOperatorName() {
        return TRANSFORMED_OPERATOR + "." + subOperator.getOperatorName();
    }

    @Override
    public double doOperation() {
        double[] oldValues = parameter.getParameterUntransformedValues();
        double ratio = subOperator.doOperation();
        double[] newValues = parameter.getParameterUntransformedValues();


        if (checkValid) { // GH: below is sloppy, but best I could do without refactoring how Parameter handles bounds
            if (generalBounds == null && !parameter.isWithinBounds()) {
                return Double.NEGATIVE_INFINITY;
            } else if (!generalBounds.isWithinBounds(parameter.getParameterValues())) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Compute Jacobians
        ratio += parameter.diffLogJacobian(oldValues, newValues);

        return ratio;
    }

    @Override
    public Parameter getParameter() {
        return subOperator.getParameter();
    }
}
